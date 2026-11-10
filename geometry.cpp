#include "geometry.h"

#include <algorithm>
#include <cmath>
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

// Store the values of one cell in a scalar or tensor field
template <typename V>
void setCell(V& cell, const V* values)
{
    cell = values[0];
}

template <typename V, pluint N>
void setCell(Array<V, N>& cell, const V* values)
{
    for (pluint k = 0; k < N; ++k)
        cell[k] = values[k];
}

// Fill the locally owned blocks (bulk and envelope) of a distributed field with N values
// per cell: valuesAt(globalIndex, values) writes them, on the main process only.
template <typename V, int N, typename MultiField, typename F>
void scatterComponents(MultiField& field, F valuesAt)
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

        buffer.assign(N * box.nCells(), V());
        if (rank == 0) {
            size_t i = 0;
            for (plint x = box.x0; x <= box.x1; ++x)
                for (plint y = box.y0; y <= box.y1; ++y)
                    for (plint z = box.z0; z <= box.z1; ++z, ++i)
                        if (x >= 0 && y >= 0 && z >= 0 && x < Nx && y < Ny && z < Nz)
                            valuesAt(gT(x, y, z), &buffer[N * i]);
            if (owner != 0)
                global::mpi().send(buffer.data(), (int)buffer.size(), owner, GEOMETRY_TAG);
        }
        else {
            global::mpi().receive(buffer.data(), (int)buffer.size(), 0, GEOMETRY_TAG);
        }

        if (rank == owner) {
            auto& block = field.getComponent(blockId);
            Dot3D location = block.getLocation();
            size_t i = 0;
            for (plint x = box.x0; x <= box.x1; ++x)
                for (plint y = box.y0; y <= box.y1; ++y)
                    for (plint z = box.z0; z <= box.z1; ++z, ++i)
                        setCell(block.get(x - location.x, y - location.y, z - location.z), &buffer[N * i]);
        }
    }
    global::mpi().barrier();
}

// One value per cell: valueAt(globalIndex), evaluated on the main process only.
template <typename V, typename F>
void scatterToField(MultiScalarField3D<V>& field, F valueAt)
{
    scatterComponents<V, 1>(field, [&valueAt](long long id, V* values) { values[0] = valueAt(id); });
}

// Value of a float32 or float64 array
double valueOf(const cnpy::NpyArray& array, size_t id)
{
    return array.word_size == sizeof(float) ? (double)array.data<float>()[id] : array.data<double>()[id];
}

}  // namespace


GeometryFile::GeometryFile(const std::string& fileName_)
    : nx(0), ny(0), nz(0), dx(0.0), hasOrigin(false), origin(0.0, 0.0, 0.0), hasLengthUnit(false), lengthUnit(1.0),
      hasStent(false), hasMultipliers(false), numPorous(0), porousBox(0, -1, 0, -1, 0, -1),
      multiplierMin(0.0, 0.0), multiplierMax(0.0, 0.0), multiplierMean(0.0, 0.0), fileName(fileName_)
{
    // Small arrays, read on the main process and broadcast: sizes, dx, origin, openings
    std::vector<double> header;      // nx, ny, nz, dx, hasOrigin, ox, oy, oz, hasStent, numOpenings, hasLengthUnit, lengthUnit
    std::vector<double> openingData; // per opening: label, radius, nx, ny, nz
    std::vector<double> deviceData;  // hasMultipliers, numPorous, porousBox (6), min, max, mean (linear, quadratic)
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
            bool linearFound = loadOptional(fileName, "linearMultiplier", linearMultiplier);
            bool quadraticFound = loadOptional(fileName, "quadraticMultiplier", quadraticMultiplier);
            if (linearFound != quadraticFound || (linearFound && !stentFound))
                throw std::runtime_error("The geometry file must contain 'linearMultiplier' and 'quadraticMultiplier' together, and only with a 'stent' array.");
            for (cnpy::NpyArray const* m : {&linearMultiplier, &quadraticMultiplier})
                if (linearFound && (m->shape != labels.shape || (m->word_size != sizeof(float) && m->word_size != sizeof(double))))
                    throw std::runtime_error("The resistance multipliers must be float32 or float64 arrays of the same shape as 'geometryFlag'.");
            if (stentFound)
                deviceData = summarizeDevice(linearFound);

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

    if (hasStent) {
        deviceData.resize(14);
        global::mpi().bCast(deviceData.data(), (int)deviceData.size());
        hasMultipliers = deviceData[0] != 0.0;
        numPorous = (long long)deviceData[1];
        porousBox = Box3D((plint)deviceData[2], (plint)deviceData[3], (plint)deviceData[4],
                          (plint)deviceData[5], (plint)deviceData[6], (plint)deviceData[7]);
        multiplierMin = Array<double, 2>(deviceData[8], deviceData[9]);
        multiplierMax = Array<double, 2>(deviceData[10], deviceData[11]);
        multiplierMean = Array<double, 2>(deviceData[12], deviceData[13]);
    }

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

// Porous voxels: flow diverter voxels in the fluid, never walls or outside voxels (the
// force itself also skips boundary cells such as openings)
bool GeometryFile::isPorous(size_t id) const
{
    return stent.data<unsigned short>()[id] > 0 && labels.data<unsigned short>()[id] >= FLUID;
}

// Main process: number, bounding box and multiplier range of the porous voxels
std::vector<double> GeometryFile::summarizeDevice(bool withMultipliers)
{
    size_t n1 = labels.shape[1], n2 = labels.shape[2];
    double count = 0.0;
    double box[6] = {1e30, -1e30, 1e30, -1e30, 1e30, -1e30};
    double lo[2] = {1e300, 1e300}, hi[2] = {-1e300, -1e300}, sum[2] = {0.0, 0.0};
    for (size_t id = 0; id < labels.num_vals; ++id) {
        if (!isPorous(id))
            continue;
        double c[3] = {(double)(id / (n1 * n2)), (double)((id / n2) % n1), (double)(id % n2)};
        for (int d = 0; d < 3; ++d) {
            box[2 * d] = std::min(box[2 * d], c[d]);
            box[2 * d + 1] = std::max(box[2 * d + 1], c[d]);
        }
        for (int k = 0; k < 2; ++k) {
            double m = withMultipliers ? valueOf(k == 0 ? linearMultiplier : quadraticMultiplier, id) : 1.0;
            if (!(m >= 0.0) || !std::isfinite(m))
                throw std::runtime_error("The resistance multipliers must be finite and >= 0 in every flow diverter voxel.");
            lo[k] = std::min(lo[k], m);
            hi[k] = std::max(hi[k], m);
            sum[k] += m;
        }
        count += 1.0;
    }
    if (count == 0.0)
        return {withMultipliers ? 1.0 : 0.0, 0.0, 0, -1, 0, -1, 0, -1, 0, 0, 0, 0, 0, 0};
    return {withMultipliers ? 1.0 : 0.0, count, box[0], box[1], box[2], box[3], box[4], box[5],
            lo[0], lo[1], hi[0], hi[1], sum[0] / count, sum[1] / count};
}

void GeometryFile::fillPorosity(MultiTensorField3D<T, 2>& field) const
{
    bool main = global::mpi().isMainProcessor();
    PLB_ASSERT(hasStent);
    scatterComponents<T, 2>(field, [this, main](long long id, T* values) {
        bool porous = main && isPorous(id);
        values[0] = porous ? (hasMultipliers ? (T)valueOf(linearMultiplier, id) : (T)1.0) : (T)0.0;
        values[1] = porous ? (hasMultipliers ? (T)valueOf(quadraticMultiplier, id) : (T)1.0) : (T)0.0;
    });
}

void GeometryFile::release()
{
    labels = cnpy::NpyArray();
    stent = cnpy::NpyArray();
    linearMultiplier = cnpy::NpyArray();
    quadraticMultiplier = cnpy::NpyArray();
}
