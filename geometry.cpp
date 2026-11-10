#include "geometry.h"

#include <stdexcept>

namespace {

const int GEOMETRY_TAG = 4711;

// Load one array of an .npz file; returns false if the file has no such array.
bool loadOptional(const std::string& fileName, const std::string& name, cnpy::NpyArray& array)
{
    try {
        array = cnpy::npz_load(fileName, name);
    }
    catch (std::runtime_error& e) {
        if (std::string(e.what()).find("not found") != std::string::npos)
            return false;
        throw;
    }
    if (array.fortran_order)
        throw std::runtime_error("The array '" + name + "' in " + fileName + " is stored in Fortran order; C order is expected.");
    return true;
}

cnpy::NpyArray loadRequired(const std::string& fileName, const std::string& name)
{
    cnpy::NpyArray array;
    if (!loadOptional(fileName, name, array))
        throw std::runtime_error("The geometry file " + fileName + " has no '" + name + "' array. Re-run the preprocessor.");
    return array;
}

// Fill the locally owned blocks (bulk and envelope) of a distributed field with
// valueAt(globalIndex), which is evaluated on the main process only.
template <typename V, typename F>
void scatterToField(MultiScalarField3D<V>& field, F valueAt)
{
    MultiBlockManagement3D const& management = field.getMultiBlockManagement();
    std::map<plint, Box3D> const& bulks = management.getSparseBlockStructure().getBulks();
    ThreadAttribution const& threads = management.getThreadAttribution();
    plint envelope = management.getEnvelopeWidth();
    int rank = global::mpi().getRank();

    std::vector<V> buffer;
    for (auto const& entry : bulks) {
        plint blockId = entry.first;
        Box3D box = entry.second.enlarge(envelope);
        int owner = threads.getMpiProcess(blockId);
        if (rank != 0 && rank != owner)
            continue;

        buffer.assign(box.nCells(), V());
        if (rank == 0) {
            size_t i = 0;
            for (plint x = box.x0; x <= box.x1; ++x)
                for (plint y = box.y0; y <= box.y1; ++y)
                    for (plint z = box.z0; z <= box.z1; ++z, ++i)
                        if (x >= 0 && y >= 0 && z >= 0 && x < Nx && y < Ny && z < Nz)
                            buffer[i] = valueAt(gT(x, y, z));
            if (owner != 0)
                global::mpi().send(buffer.data(), (int)buffer.size(), owner, GEOMETRY_TAG);
        }
        else {
            global::mpi().receive(buffer.data(), (int)buffer.size(), 0, GEOMETRY_TAG);
        }

        if (rank == owner) {
            ScalarField3D<V>& block = field.getComponent(blockId);
            Dot3D location = block.getLocation();
            size_t i = 0;
            for (plint x = box.x0; x <= box.x1; ++x)
                for (plint y = box.y0; y <= box.y1; ++y)
                    for (plint z = box.z0; z <= box.z1; ++z, ++i)
                        block.get(x - location.x, y - location.y, z - location.z) = buffer[i];
        }
    }
    global::mpi().barrier();
}

}  // namespace


GeometryFile::GeometryFile(const std::string& fileName_)
    : nx(0), ny(0), nz(0), dx(0.0), hasOrigin(false), origin(0.0, 0.0, 0.0), hasLengthUnit(false), lengthUnit(1.0),
      hasStent(false), fileName(fileName_)
{
    // Small arrays, read on the main process and broadcast: sizes, dx, origin, openings
    std::vector<double> header;      // nx, ny, nz, dx, hasOrigin, ox, oy, oz, hasStent, numOpenings, hasLengthUnit, lengthUnit
    std::vector<double> openingData; // per opening: label, radius, nx, ny, nz
    std::string error;

    if (global::mpi().isMainProcessor()) {
        try {
            labels = loadRequired(fileName, "geometryFlag");
            if (labels.shape.size() != 3 || labels.word_size != sizeof(unsigned short))
                throw std::runtime_error("'geometryFlag' must be a 3D int16 array.");

            cnpy::NpyArray dxArray = loadRequired(fileName, "dx");
            cnpy::NpyArray originArray, unitArray;
            bool originFound = loadOptional(fileName, "origin", originArray);
            bool unitFound = loadOptional(fileName, "lengthUnit", unitArray);
            bool stentFound = loadOptional(fileName, "stent", stent);
            if (stentFound && (stent.shape != labels.shape || stent.word_size != sizeof(unsigned short)))
                throw std::runtime_error("The 'stent' array must be an int16 array of the same shape as 'geometryFlag'.");

            cnpy::NpyArray index = loadRequired(fileName, "openingIndex");
            cnpy::NpyArray radius = loadRequired(fileName, "openingRadius");
            cnpy::NpyArray normal = loadRequired(fileName, "openingNormal");
            size_t numOpenings = index.num_vals;
            if (radius.num_vals != numOpenings || normal.shape.size() != 2 || normal.shape[0] != numOpenings || normal.shape[1] != 3)
                throw std::runtime_error("The opening arrays in the geometry file have inconsistent sizes.");

            const double* o = originFound ? originArray.data<double>() : nullptr;
            header = {(double)labels.shape[0], (double)labels.shape[1], (double)labels.shape[2],
                      dxArray.data<double>()[0], originFound ? 1.0 : 0.0,
                      o ? o[0] : 0.0, o ? o[1] : 0.0, o ? o[2] : 0.0,
                      stentFound ? 1.0 : 0.0, (double)numOpenings,
                      unitFound ? 1.0 : 0.0, unitFound ? unitArray.data<double>()[0] : 1.0};
            for (size_t i = 0; i < numOpenings; i++) {
                const double* n = normal.data<double>() + 3 * i;
                openingData.insert(openingData.end(), {(double)index.data<unsigned short>()[i], radius.data<double>()[i], n[0], n[1], n[2]});
            }
        }
        catch (std::exception& e) {
            error = e.what();
        }
    }

    // Propagate errors so that all processes stop consistently
    int failed = error.empty() ? 0 : 1;
    global::mpi().bCast(&failed, 1);
    if (failed) {
        global::mpi().bCast(error);
        throw std::runtime_error(error);
    }

    header.resize(12);
    global::mpi().bCast(header.data(), (int)header.size());
    nx = (int)header[0]; ny = (int)header[1]; nz = (int)header[2];
    dx = header[3];
    hasOrigin = header[4] != 0.0;
    origin = Array<double, 3>(header[5], header[6], header[7]);
    hasStent = header[8] != 0.0;
    int numOpenings = (int)header[9];
    hasLengthUnit = header[10] != 0.0;
    lengthUnit = header[11];

    openingData.resize(5 * numOpenings);
    if (numOpenings > 0)
        global::mpi().bCast(openingData.data(), (int)openingData.size());
    for (int i = 0; i < numOpenings; i++) {
        openingLabels.push_back((int)openingData[5 * i]);
        openingRadius.push_back(openingData[5 * i + 1]);
        openingNormal.push_back(vec3d(openingData[5 * i + 2], openingData[5 * i + 3], openingData[5 * i + 4]));
    }

    // Opening voxels: one scan of the label array on the main process, then broadcast
    Nx = nx; Ny = ny; Nz = nz;
    std::map<int, std::vector<int> > coordinates;
    if (global::mpi().isMainProcessor()) {
        std::vector<int> slot(65536, -1);
        for (int i = 0; i < numOpenings; i++)
            slot[openingLabels[i]] = i;
        const unsigned short* data = labels.data<unsigned short>();
        for (int x = 0; x < nx; x++)
            for (int y = 0; y < ny; y++)
                for (int z = 0; z < nz; z++) {
                    int s = slot[data[gT(x, y, z)]];
                    if (s >= 0) {
                        std::vector<int>& c = coordinates[openingLabels[s]];
                        c.push_back(x); c.push_back(y); c.push_back(z);
                    }
                }
    }
    for (int label : openingLabels) {
        std::vector<int>& c = coordinates[label];
        int count = (int)c.size();
        global::mpi().bCast(&count, 1);
        c.resize(count);
        if (count > 0)
            global::mpi().bCast(c.data(), count);
        std::vector<Index3D>& nodes = openingNodes[label];
        for (int k = 0; k + 2 < count; k += 3)
            nodes.push_back({c[k], c[k + 1], c[k + 2]});
    }
}

void GeometryFile::fillLabels(MultiScalarField3D<int>& field) const
{
    const unsigned short* data = global::mpi().isMainProcessor() ? labels.data<unsigned short>() : nullptr;
    scatterToField(field, [data](long long id) { return (int)data[id]; });
}

void GeometryFile::fillOccupancy(MultiScalarField3D<int>& field) const
{
    const unsigned short* data = global::mpi().isMainProcessor() ? labels.data<unsigned short>() : nullptr;
    scatterToField(field, [data](long long id) { return data[id] > UNUSED ? 1 : 0; });
}

void GeometryFile::fillPorosity(MultiScalarField3D<T>& field) const
{
    bool main = global::mpi().isMainProcessor();
    const unsigned short* data = main ? labels.data<unsigned short>() : nullptr;
    const unsigned short* device = (main && hasStent) ? stent.data<unsigned short>() : nullptr;
    // Porous forcing only in fluid voxels, never in walls or outside voxels (the force
    // itself also skips boundary cells such as openings)
    scatterToField(field, [data, device](long long id) {
        return (device && device[id] > 0 && data[id] >= FLUID) ? (T)1.0 : (T)0.0;
    });
}

void GeometryFile::release()
{
    labels = cnpy::NpyArray();
    stent = cnpy::NpyArray();
}
