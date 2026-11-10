#include "opening.h"

#include <algorithm>
#include <fstream>
#include <map>
#include <sstream>

namespace {

double component(const vec3d& v, int axis)
{
    return axis == 0 ? v.x : (axis == 1 ? v.y : v.z);
}

const char* openingTypeName(OpeningType type)
{
    switch (type) {
        case OPENING_VELOCITY: return "velocity";
        case OPENING_MURRAY:   return "Murray outlet";
        case OPENING_PRESSURE: return "pressure";
        case OUTLET_FREEFLOW:  return "free outflow";
    }
    return "unknown";
}

std::string labelStr(int label)
{
    return "Opening (label " + std::to_string(label) + ")";
}

}  // namespace


OpeningHandler::OpeningHandler(const std::vector<Index3D>& nodes_, int label_, OpeningType type_, T radius_lb, vec3d dirVec)
    : label(label_), type(type_), nodes(nodes_), R(radius_lb), faceAxis(-1), faceSign(0),
      profileFlux(0.0), bcParameter(0.0),
      hasScaleFunction(false), cScale(1.0), cTimePos(0), cTimeVal(0.0)
{
    if (type < OPENING_VELOCITY || type > OUTLET_FREEFLOW)
        throw std::runtime_error(labelStr(label) + ": unknown opening type " + std::to_string((int)type) + " (valid types: 1-4).");

    int lo[3] = {Nx, Ny, Nz};
    int hi[3] = {-1, -1, -1};
    double sum[3] = {0.0, 0.0, 0.0};
    for (auto const& node : nodes) {
        int c[3] = {node.x, node.y, node.z};
        for (int a = 0; a < 3; a++) {
            lo[a] = std::min(lo[a], c[a]);
            hi[a] = std::max(hi[a], c[a]);
            sum[a] += c[a];
        }
        nodeSet.insert(gT(node.x, node.y, node.z));
    }

    if (nodes.empty())
        throw std::runtime_error(labelStr(label) + " has no voxels in the geometry file.");

    boundingBox = Box3D(lo[0], hi[0], lo[1], hi[1], lo[2], hi[2]);
    center.set(sum[0] / nodes.size(), sum[1] / nodes.size(), sum[2] / nodes.size());

    // The opening must be a planar cut on one face of the domain bounding box.
    int N[3] = {Nx, Ny, Nz};
    for (int a = 0; a < 3 && faceAxis < 0; a++) {
        if (lo[a] == hi[a] && (lo[a] == 0 || lo[a] == N[a] - 1)) {
            faceAxis = a;
            faceSign = (lo[a] == 0) ? 1 : -1;
        }
    }
    if (faceAxis < 0) {
        std::ostringstream msg;
        msg << labelStr(label) << " does not lie on a face of the domain (bounding box x " << lo[0] << "-" << hi[0]
            << ", y " << lo[1] << "-" << hi[1] << ", z " << lo[2] << "-" << hi[2]
            << "). Openings must be planar cuts on a face of the bounding box.";
        throw std::runtime_error(msg.str());
    }

    // Rows of adjacent nodes along one in-face axis; together they cover exactly the opening.
    {
        int a1 = (faceAxis + 1) % 3, a2 = (faceAxis + 2) % 3;
        std::map<int, std::vector<int> > rows;
        for (auto const& node : nodes) {
            int c[3] = {node.x, node.y, node.z};
            rows[c[a2]].push_back(c[a1]);
        }
        for (auto& row : rows) {
            std::vector<int>& u = row.second;
            std::sort(u.begin(), u.end());
            size_t start = 0;
            for (size_t k = 1; k <= u.size(); k++) {
                if (k == u.size() || u[k] != u[k - 1] + 1) {
                    plint b0[3], b1[3];
                    b0[faceAxis] = b1[faceAxis] = lo[faceAxis];
                    b0[a1] = u[start];
                    b1[a1] = u[k - 1];
                    b0[a2] = b1[a2] = row.first;
                    runs.push_back(Box3D(b0[0], b1[0], b0[1], b1[1], b0[2], b1[2]));
                    start = k;
                }
            }
        }
    }

    // Direction of the opening (from the centerline). It must point into the domain.
    direction = dirVec.getNormal();
    if (direction.norm() == 0.0) {
        pcout << "-> WARNING: no direction given for " << labelStr(label) << ", using the face normal." << std::endl;
        direction.set(faceAxis == 0 ? faceSign : 0, faceAxis == 1 ? faceSign : 0, faceAxis == 2 ? faceSign : 0);
    }
    T cosAngle = component(direction, faceAxis) * faceSign;
    if (cosAngle <= 0.0)
        throw std::runtime_error(labelStr(label) + ": the centerline direction points out of the domain or along the face. "
                                 "Check that the centerline ends at this opening.");
    if (cosAngle < 0.5)
        pcout << "-> WARNING: the centerline direction of " << labelStr(label) << " is more than 60 degrees off the face normal." << std::endl;
}

void OpeningHandler::printOpeningDetails(const SimPar& s) const
{
    static const char* faceNames[3] = {"x", "y", "z"};
    T cosAngle = component(direction, faceAxis) * faceSign;

    pcout << "---- Opening parameters -----" << std::endl;
    pcout << "-> Name: " << getName() << std::endl;
    pcout << "-> Geometry label: " << label << std::endl;
    pcout << "-> Opening type: " << type << " (" << openingTypeName(type) << ")" << std::endl;
    pcout << "-> Face: " << faceNames[faceAxis] << (faceSign > 0 ? " min" : " max") << std::endl;
    pcout << "-> Radius (centerline) [m]: " << R * s.C_l << std::endl;
    pcout << "-> Center [lb]: " << center.x << " " << center.y << " " << center.z << std::endl;
    pcout << "-> Direction (inward): " << direction.x << " " << direction.y << " " << direction.z
          << " (" << std::acos(std::min((T)1.0, cosAngle)) * 180.0 / std::acos(-1.0) << " deg from the face normal)" << std::endl;
    pcout << "-> Nodes: " << nodes.size() << " in " << runs.size() << " rows (area " << nodes.size() * s.C_l * s.C_l << " m^2)" << std::endl;
    pcout << "-> Scale function: " << (hasScaleFunction ? "yes" : "no");
    if (hasScaleFunction)
        pcout << ", " << scaleTime.size() << " samples, period " << getPeriod() << " s";
    pcout << std::endl;

    if (type == OPENING_VELOCITY) {
        T area = nodes.size() * s.C_l * s.C_l;
        T flowRate = profileFlux * bcParameter * s.C_l * s.C_l * s.C_l / s.C_t;   // unscaled by the waveform
        pcout << "-> Peak velocity (before scaling) [m/s]: " << bcParameter * s.C_l / s.C_t
              << " [lb]: " << bcParameter << std::endl;
        pcout << "-> Flow rate (before scaling) [m^3/s]: " << flowRate
              << " (mean velocity over the opening " << flowRate / area << " m/s)" << std::endl;
        if (std::fabs(bcParameter) > s.U_max_lb)
            pcout << "-> *WARNING*: the peak velocity exceeds " << s.U_max_lb << " in lattice units ("
                  << s.U_max_lb * s.C_l / s.C_t << " m/s). Consider decreasing the time step." << std::endl;
    }
    else if (type == OPENING_MURRAY) {
        pcout << "-> Flow rate: set every time step by Murray's law" << std::endl;
    }
    else if (type == OPENING_PRESSURE) {
        pcout << "-> Gauge pressure (before scaling) [Pa]: " << bcParameter / 3.0 * s.C_p
              << " (lattice density " << 1.0 + bcParameter << ")" << std::endl;
        if (std::fabs(bcParameter) > 0.05)
            pcout << "-> *WARNING*: the lattice density deviates by more than 5% from 1 (" << 1.0 + bcParameter
                  << "); compressibility errors will be large. Consider decreasing the time step." << std::endl;
    }
}

// Convert the SI boundary value to lattice units.
void OpeningHandler::setBCParameter(T parameter, const SimPar& s)
{
    if (type == OPENING_VELOCITY) {
        // Peak velocity [m/s] -> lattice units
        bcParameter = parameter * s.C_t / s.C_l;
    }
    else if (type == OPENING_PRESSURE) {
        // Gauge pressure [Pa] -> lattice density deviation: p_lb = c_s^2 (rho - 1), c_s^2 = 1/3
        bcParameter = 3.0 * parameter / s.C_p;
    }
    else {
        // Murray outlets are recomputed every time step; free outflow has no parameter.
        bcParameter = 0.0;
    }
}

T OpeningHandler::getImposedFlowRate(const SimPar& s) const
{
    return getImposedFluxLB() * s.C_l * s.C_l * s.C_l / s.C_t;
}

T OpeningHandler::measureFlowRate(MultiBlockLattice3D<T, DESCRIPTOR>& lattice, const SimPar& s) const
{
    OpeningFluxFunctional3D<T, DESCRIPTOR> functional(&nodeSet, faceAxis, faceSign);
    applyProcessingFunctional(functional, boundingBox, lattice);
    return functional.getFlux() * s.C_l * s.C_l * s.C_l / s.C_t;
}

// Parabolic profile with peak 1 (lattice units) in the centerline direction. The opening is
// approximated by an ellipse: the farthest node from the centre defines the major axis and the
// centerline radius the minor axis. Nodes outside the ellipse get zero velocity.
void OpeningHandler::createPoiseuilleProfile()
{
    pcout << "-> Creating direction-corrected Poiseuille velocity profile on label: " << label << std::endl;

    if (!isVelocityType())
        pcout << "WARNING! Setting a velocity profile on a non-velocity opening! This has no effect." << std::endl;

    // Farthest node from the centre
    T l_max = 0.0;
    vec3d v_max(0, 0, 0);
    for (auto const& v : nodes) {
        vec3d dist(v.x - center.x, v.y - center.y, v.z - center.z);
        T length = dist.norm();
        if (length > l_max) {
            l_max = length;
            v_max.set(dist.x, dist.y, dist.z);
        }
    }

    profile.clear();
    if (l_max == 0.0 || R <= 0.0) {
        pcout << "-> WARNING: degenerate opening, using a uniform profile instead." << std::endl;
        createBluntVelocityProfile();
        return;
    }

    // Minor axis, perpendicular to the major axis, scaled to the centerline radius
    vec3d v_min = v_max.cross(direction);
    v_min.normalize();
    v_min = v_min * R;
    T l_min = v_min.norm();
    T r_min2 = l_min * l_min;

    pcout << "-> Ellipse semi-axes [lb]: " << l_min << ", " << l_max << std::endl;

    int numClamped = 0;
    for (auto const& v : nodes) {
        vec3d cR(v.x - center.x, v.y - center.y, v.z - center.z);

        T proj_r_min = cR.dot(v_min) / l_min;
        T proj_r_max = cR.dot(v_max) / l_max * (l_min / l_max);  // Scaled to the minor axis

        T proj_r_len2 = proj_r_max * proj_r_max + proj_r_min * proj_r_min;
        T magnitude = (r_min2 - proj_r_len2) / r_min2;
        if (magnitude < 0.0) {
            magnitude = 0.0;  // Outside the ellipse: no backflow at the rim
            numClamped++;
        }
        profile[gT(v.x, v.y, v.z)] = direction * magnitude;
    }
    if (numClamped > 0)
        pcout << "-> " << numClamped << " rim node(s) outside the fitted ellipse set to zero velocity." << std::endl;

    profileFlux = 0.0;
    for (auto const& p : profile)
        profileFlux += component(p.second, faceAxis) * faceSign;
    if (profileFlux <= 0.0)
        throw std::runtime_error(labelStr(label) + ": the velocity profile has no flux through the opening.");
}

void OpeningHandler::createBluntVelocityProfile()
{
    pcout << "-> Creating blunt velocity profile on label: " << label << std::endl;

    if (!isVelocityType())
        pcout << "WARNING! Setting a velocity profile on a non-velocity opening! This has no effect." << std::endl;

    profile.clear();
    for (auto const& v : nodes)
        profile[gT(v.x, v.y, v.z)] = direction;

    profileFlux = 0.0;
    for (auto const& p : profile)
        profileFlux += component(p.second, faceAxis) * faceSign;
    if (profileFlux <= 0.0)
        throw std::runtime_error(labelStr(label) + ": the velocity profile has no flux through the opening.");
}

void OpeningHandler::createConstantPressureProfile()
{
    pcout << "-> Creating constant pressure on label: " << label << std::endl;
    if (type != OPENING_PRESSURE)
        pcout << "WARNING! Setting a pressure profile on a non-pressure opening! This has no effect." << std::endl;
    // The pressure is uniform over the opening; its value is bcParameter (scaled by the waveform).
}

// Waveform file format: first line N, followed by N lines "time value".
// Times must start at 0 and increase; the last time is the period.
void OpeningHandler::loadScaleFunction(const std::string& fileName)
{
    if (type == OPENING_MURRAY || type == OUTLET_FREEFLOW) {
        pcout << "-> *WARNING*: a " << openingTypeName(type) << " opening can't have a scale function. "
              << fileName << " is ignored." << std::endl;
        hasScaleFunction = false;
        return;
    }

    pcout << "-> Loading scale function: " << fileName << std::endl;

    // Read on the main process, then broadcast, so that all ranks see the same data and errors.
    int status = 0;  // 0 ok, 1 cannot open, 2 parse error
    std::vector<double> data;
    if (global::mpi().isMainProcessor()) {
        std::ifstream fin(fileName.c_str());
        if (!fin.is_open()) {
            status = 1;
        }
        else {
            int Ns = 0;
            if (!(fin >> Ns) || Ns < 0)
                status = 2;
            for (int i = 0; i < Ns && status == 0; i++) {
                double time, value;
                if (fin >> time >> value) {
                    data.push_back(time);
                    data.push_back(value);
                }
                else {
                    status = 2;
                }
            }
        }
    }
    global::mpi().bCast(&status, 1);
    if (status == 1)
        throw std::runtime_error("Cannot open scale function file " + fileName);
    if (status == 2)
        throw std::runtime_error("Cannot parse scale function file " + fileName +
                                 " (expected the number of samples N on the first line, followed by N lines 'time value').");

    int numValues = (int)data.size();
    global::mpi().bCast(&numValues, 1);
    data.resize(numValues);
    if (numValues > 0)
        global::mpi().bCast(&data[0], numValues);

    scaleTime.clear();
    scaleSignal.clear();
    for (int i = 0; i + 1 < numValues; i += 2) {
        scaleTime.push_back(data[i]);
        scaleSignal.push_back(data[i + 1]);
    }

    if (scaleTime.size() < 2)
        throw std::runtime_error("Scale function " + fileName + " needs at least 2 samples.");
    if (scaleTime.front() != 0.0)
        throw std::runtime_error("Scale function " + fileName + " must start at time 0.");
    for (size_t i = 1; i < scaleTime.size(); i++)
        if (scaleTime[i] <= scaleTime[i - 1])
            throw std::runtime_error("Scale function " + fileName + ": the times must be strictly increasing.");
    if (std::fabs(scaleSignal.front() - scaleSignal.back()) > 1e-6 * std::max((T)1.0, std::fabs(scaleSignal.front())))
        pcout << "-> *WARNING*: the first and last value of " << fileName
              << " differ, so the boundary value jumps at the start of every cycle." << std::endl;

    hasScaleFunction = true;
    cTimePos = 0;
    cTimeVal = 0.0;
    pcout << fileName << " loaded with " << scaleTime.size() << " data points, period " << scaleTime.back() << " s." << std::endl;
}

// The boundary condition is instantiated on the opening nodes only (row by row), so that
// wall voxels next to the opening and other openings on the same face are not touched.
void OpeningHandler::setBCType(MultiBlockLattice3D<T, DESCRIPTOR>& lattice)
{
    Box3D domain = lattice.getBoundingBox();
    if (type == OPENING_VELOCITY || type == OPENING_MURRAY) {
        OnLatticeBoundaryCondition3D<T, DESCRIPTOR>* bc = createLocalBoundaryCondition3D<T, DESCRIPTOR>();
        for (Box3D const& run : runs)
            bc->setVelocityConditionOnBlockBoundaries(lattice, domain, run, boundary::dirichlet);
        delete bc;
    }
    else if (type == OPENING_PRESSURE) {
        OnLatticeBoundaryCondition3D<T, DESCRIPTOR>* bc = createInterpBoundaryCondition3D<T, DESCRIPTOR>();
        for (Box3D const& run : runs)
            bc->setPressureConditionOnBlockBoundaries(lattice, domain, run, boundary::dirichlet);
        delete bc;
    }
    else if (type == OUTLET_FREEFLOW) {
        // Virtual outlet. The auxiliary fields must live as long as the lattice.
        MultiScalarField3D<T>* rhoBar = generateMultiScalarField<T>((MultiBlock3D&)lattice, 2).release();
        rhoBar->toggleInternalStatistics(false);

        MultiTensorField3D<T, 3>* j = generateMultiTensorField<T, 3>((MultiBlock3D&)lattice, 2).release();
        j->toggleInternalStatistics(false);

        std::vector<MultiBlock3D*> bcargs;
        bcargs.push_back(&lattice);
        bcargs.push_back(rhoBar);
        bcargs.push_back(j);

        for (Box3D const& run : runs)
            integrateProcessingFunctional(new VirtualOutlet<T, DESCRIPTOR>(1.0, domain, 1), run, bcargs, 2);
    }
}

void OpeningHandler::progressTime(MultiBlockLattice3D<T, DESCRIPTOR>& lattice, T dt)
{
    if (hasScaleFunction) {
        int len = scaleTime.size();
        T period = scaleTime[len - 1];

        cTimeVal += dt;
        if (cTimeVal >= period) {
            cTimeVal = std::fmod(cTimeVal, period);
            cTimePos = 0;
        }
        while (cTimePos < len - 2 && scaleTime[cTimePos + 1] < cTimeVal)
            cTimePos++;

        cScale = interpolate(scaleTime[cTimePos], scaleTime[cTimePos + 1], cTimeVal, scaleSignal[cTimePos], scaleSignal[cTimePos + 1]);
    }

    // Impose the profile scaled with the boundary parameter and the waveform
    if (isVelocityType()) {
        applyProcessingFunctional(new ImposeOpeningValues3D<T, DESCRIPTOR>(&nodeSet, &profile, bcParameter * cScale),
                                  boundingBox, lattice);
    }
    else if (type == OPENING_PRESSURE) {
        applyProcessingFunctional(new ImposeOpeningValues3D<T, DESCRIPTOR>(&nodeSet, 1.0 + bcParameter * cScale),
                                  boundingBox, lattice);
    }
    // Free outflow: nothing to impose.
}
