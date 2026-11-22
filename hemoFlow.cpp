#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "globals.h"
#include "helper.h"
#include "opening.h"
#include "porous.h"
#include "nonNewtonian.h"
#include "io.h"
#include "geometry.h"
#ifdef HEMOFLOW_GPU
#include "gpuSolver.h"
#include <unistd.h>
#endif

#include <chrono>

using namespace std;


/* ********** GLOBAL VARIABLES ************/
bool DEBUG = false;
bool CHECKPOINTING = false;   // Restart from a checkpoint

// Domain size
int Nx=0;
int Ny=0;
int Nz=0;

// Flow diverter coefficients
T linCoeff = 0.0;
T quadCoeff = 0.0;
T linCoeff_lb = 0.0;
T quadCoeff_lb = 0.0;

// Simulation parameters structure
SimPar sim;

// Fluid properties
T fluidKinematicViscosity = DEFAULT_BLOOD_KINEMATIC_VISCOSITY;  // [m^2/s], Newtonian model
T fluidDensity = DEFAULT_BLOOD_DENSITY;                          // [kg/m^3]
string fluidModel = "newtonian";                                 // newtonian | carreau
// Carreau-Yasuda parameters (SI), defaults: blood after Cho and Kensey (Carreau, a = 2)
struct CarreauParametersSI {
    T mu0 = DEFAULT_CARREAU_MU0;        // [Pa s]
    T muInf = DEFAULT_CARREAU_MU_INF;   // [Pa s]
    T lambda = DEFAULT_CARREAU_LAMBDA;  // [s]
    T n = DEFAULT_CARREAU_N;
    T a = 2.0;
} carreauSI;

// Technical simulation parameters
bool useCheckpoint = false;
bool useSparse = true;    // Sparse block decomposition: only blocks that contain the vessel are allocated
int blockSize = 24;       // Side length of the blocks of the sparse decomposition [LB nodes]
int envelopeWidth = 1;    // Set before the lattice is created (2 with pressure openings)
string outputFolder;
string workingFolder;
string outputCoordinates = "stl";   // Units of the output grid: stl (those of the input STL) | m
bool outputCoordinatesGiven = false;
T simLength;
int numCycles;
T cyclePeriod;
T saveFreqTime;
T checkpointFreqTime = 0.0;
string mode;
OutputOptions outputOptions;

// Vector of openings
vector<OpeningHandler*> openings;

// Simulation data structures
MultiBlockLattice3D<T, DESCRIPTOR> *lattice = nullptr;
MultiScalarField3D<int> *labelField = nullptr;   // Voxel labels (distributed)
MultiTensorField3D<T, 2> *porosityField = nullptr;   // (linear, quadratic) resistance multipliers in porous voxels, 0 elsewhere (distributed)
Box3D porousBox;   // Bounding box of the porous voxels

// Where the time loop runs (doc/dev/TODO-GPU.md): <simulation><device> auto | gpu | cpu
string deviceOption = "auto";
bool gpuMode = false;   // The time loop runs on the GPU (one dense block per GPU rank)
#ifdef HEMOFLOW_GPU
bool gpuRank = false;                          // This rank runs a block on its GPU (chooseGpu)
GpuSolver *gpu = nullptr;                      // The GPU lattice, once the initial state is set
MPI_Comm worldCommunicator = MPI_COMM_NULL;    // GPU mode on several ranks: all ranks of the launch
#endif


namespace {

// ********** Configuration helpers **********

std::string trim(const std::string& s)
{
    size_t b = s.find_first_not_of(" \t\r\n");
    if (b == std::string::npos)
        return "";
    size_t e = s.find_last_not_of(" \t\r\n");
    return s.substr(b, e - b + 1);
}

template <typename V>
bool parseValue(const std::string& text, V& value)
{
    std::istringstream ss(text);
    V tmp;
    if (!(ss >> tmp))
        return false;
    ss >> std::ws;
    if (!ss.eof())
        return false;
    value = tmp;
    return true;
}

template <>
bool parseValue<std::string>(const std::string& text, std::string& value)
{
    value = text;
    return true;
}

template <>
bool parseValue<bool>(const std::string& text, bool& value)
{
    std::string word = text;
    std::transform(word.begin(), word.end(), word.begin(), ::tolower);
    if (word == "true" || word == "1" || word == "yes" || word == "on") { value = true; return true; }
    if (word == "false" || word == "0" || word == "no" || word == "off") { value = false; return true; }
    return false;
}

// Text of a child element. Returns false if the element does not exist.
bool readText(XMLreaderProxy const& parent, const std::string& tag, std::string& text)
{
    try {
        XMLreaderProxy child = parent[tag];
        child.read(text);
    }
    catch (PlbIOException&) {
        return false;
    }
    text = trim(text);
    return true;
}

// Optional element: returns false if it is missing or empty, throws if it cannot be parsed.
template <typename V>
bool readOptional(XMLreaderProxy const& parent, const std::string& where, const std::string& tag, V& value)
{
    std::string text;
    if (!readText(parent, tag, text) || text.empty())
        return false;
    if (!parseValue(text, value))
        throw std::runtime_error("Cannot parse the value '" + text + "' of <" + where + "><" + tag + ">.");
    return true;
}

template <typename V>
void readRequired(XMLreaderProxy const& parent, const std::string& where, const std::string& tag, V& value)
{
    if (!readOptional(parent, where, tag, value))
        throw std::runtime_error("Missing required element <" + where + "><" + tag + "> in the configuration file.");
}

bool hasSection(XMLreader const& xml, const std::string& name)
{
    try {
        xml[name];
    }
    catch (PlbIOException&) {
        return false;
    }
    return true;
}

bool hasChild(XMLreaderProxy const& parent, const std::string& name)
{
    try {
        parent[name];
    }
    catch (PlbIOException&) {
        return false;
    }
    return true;
}

// Paths in the configuration are relative to the configuration file, unless absolute.
std::string resolvePath(const std::string& path)
{
    if (!path.empty() && path[0] == '/')
        return path;
    return workingFolder + "/" + path;
}

// One boundary condition to be created
struct OpeningSpec {
    int label;
    OpeningType type;
    std::string name;
    double parameter;          // Velocity: mean velocity [m/s]; pressure: gauge pressure [Pa]
    std::string waveformFile;  // Empty: constant
};

const char* VALID_MODES = "'' (no mode), 'murray', 'pzero'";


// ********** Dynamics of the GPU mode **********
//
// In GPU mode the lattices have NoDynamics(1.0), the dynamics of the outside voxels, as their
// background dynamics, shared by all cells that use it. Palabos gives every other cell its own copy
// of a dynamics object, and in the GPU build a heap allocation costs about 1.2 kB of host memory
// (CUDA managed memory) against about 50 bytes in a CPU build: with the fluid dynamics as
// background and a NoDynamics copy per outside voxel, NAP180 at 20 M voxels needed more than
// 80 GB. The CPU build keeps the fluid dynamics as background (the results are the same).
// The functionals cover bulk and envelope of every block (the label field has its envelope filled,
// GeometryFile), so that Palabos does not regenerate the dynamics of all envelope cells afterwards.

// A copy of the dynamics in every fluid, opening and device voxel (label other than UNUSED and WALL)
class FluidDynamicsFunctional3D : public BoxProcessingFunctional3D_LS<T, DESCRIPTOR, int> {
public:
    explicit FluidDynamicsFunctional3D(Dynamics<T, DESCRIPTOR>* dynamics_) : dynamics(dynamics_) { }
    FluidDynamicsFunctional3D(FluidDynamicsFunctional3D const& rhs) : dynamics(rhs.dynamics->clone()) { }
    FluidDynamicsFunctional3D& operator=(FluidDynamicsFunctional3D const&) = delete;
    ~FluidDynamicsFunctional3D() { delete dynamics; }

    virtual void process(Box3D domain, BlockLattice3D<T, DESCRIPTOR>& lattice, ScalarField3D<int>& labels)
    {
        Dot3D offset = computeRelativeDisplacement(lattice, labels);
        for (plint iX = domain.x0; iX <= domain.x1; ++iX)
            for (plint iY = domain.y0; iY <= domain.y1; ++iY)
                for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                    int label = labels.get(iX + offset.x, iY + offset.y, iZ + offset.z);
                    if (label != UNUSED && label != WALL)
                        lattice.attributeDynamics(iX, iY, iZ, dynamics->clone());
                }
    }
    virtual FluidDynamicsFunctional3D* clone() const { return new FluidDynamicsFunctional3D(*this); }
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const
    {
        modified[0] = modif::dataStructure;
        modified[1] = modif::nothing;
    }
    virtual BlockDomain::DomainT appliesTo() const { return BlockDomain::bulkAndEnvelope; }

private:
    Dynamics<T, DESCRIPTOR>* dynamics;
};

// The outside voxels (label UNUSED) back to the shared background dynamics
class OutsideToBackgroundFunctional3D : public BoxProcessingFunctional3D_LS<T, DESCRIPTOR, int> {
public:
    virtual void process(Box3D domain, BlockLattice3D<T, DESCRIPTOR>& lattice, ScalarField3D<int>& labels)
    {
        Dot3D offset = computeRelativeDisplacement(lattice, labels);
        Dynamics<T, DESCRIPTOR>* background = &lattice.getBackgroundDynamics();
        for (plint iX = domain.x0; iX <= domain.x1; ++iX)
            for (plint iY = domain.y0; iY <= domain.y1; ++iY)
                for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ)
                    if (labels.get(iX + offset.x, iY + offset.y, iZ + offset.z) == UNUSED &&
                        &lattice.get(iX, iY, iZ).getDynamics() != background)
                        lattice.attributeDynamics(iX, iY, iZ, background);
    }
    virtual OutsideToBackgroundFunctional3D* clone() const { return new OutsideToBackgroundFunctional3D(*this); }
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const
    {
        modified[0] = modif::dataStructure;
        modified[1] = modif::nothing;
    }
    virtual BlockDomain::DomainT appliesTo() const { return BlockDomain::bulkAndEnvelope; }
};

// The bulk dynamics of the fluid
Dynamics<T, DESCRIPTOR>* newBulkDynamics()
{
#if LES
    return new SmagorinskyRegularizedDynamics<T, DESCRIPTOR>(sim.omega, cSmago);
#else
    if (sim.nonNewtonian)
        return new CarreauBackgroundDynamics(sim.omega, sim.carreau);
    return new BackgroundDynamics(sim.omega);
#endif
}

}  // namespace


// *** LB parameters from the voxel size, time step and viscosity
void calcSimulationParameters(SimPar &sim, T dx, T dt, T nu, T rho, T U_max_LB_ = 0.1)
{
    sim.C_l = dx;
    sim.U_max_lb = U_max_LB_;

    T nu_lb;
    if(dt > 0.0){
        sim.C_t = dt;
        nu_lb = nu * sim.C_t / sim.C_l / sim.C_l;
        T tau = 3.0*nu_lb+0.5;
        sim.omega = 1.0 / tau;
    }
    else {
        sim.omega = 1.0;
        nu_lb =  0.5 / 3.0;
        sim.C_t = nu_lb * sim.C_l * sim.C_l / nu;
    }

    // Derived quantities using density. The density only scales the pressure and stress output.
    sim.C_r = rho;
    sim.C_p = sim.C_r * sim.C_l * sim.C_l / (sim.C_t * sim.C_t);
    sim.C_m = sim.C_r * sim.C_l * sim.C_l * sim.C_l;

    // Porous coefficients. The force density in a porous voxel is -(a u + b u|u|) / dx, so
    // a [Pa s/m] and b [Pa s^2/m^2] describe the pressure drop across ONE voxel-thick layer.
    // https://youtu.be/sOQMXxoKFQM?t=1168
    linCoeff_lb = linCoeff * sim.C_l*sim.C_l * sim.C_t / sim.C_m;
    quadCoeff_lb = quadCoeff * sim.C_l*sim.C_l * sim.C_l / sim.C_m;

    if(sim.U_max_lb > 0.1)
        pcout << "** Parameter sanity WARNING ** LBM velocity seems high:" << sim.U_max_lb << endl;
    if(sim.omega > 1.9)
        pcout << "** Parameter sanity WARNING ** LBM omega seems high: " << sim.omega
              << " (tau = " << 1.0 / sim.omega << "). The simulation may become unstable; consider a smaller dt." << endl;

}

// Carreau-Yasuda parameters in lattice units. sim.omega belongs to the infinite-shear
// viscosity, which calcSimulationParameters() was given for the non-Newtonian model.
//  Y.I. Cho and K.R. Kensey, Biorheology 28, 241 (1991); B.M. Johnston et al., J. Biomech. 37, 709 (2004)
void setCarreauModel(SimPar &sim, CarreauParametersSI const& p, T rho)
{
    sim.nonNewtonian = true;
    sim.carreau.nu0 = p.mu0 / rho * sim.C_t / sim.C_l / sim.C_l;
    sim.carreau.nuInf = p.muInf / rho * sim.C_t / sim.C_l / sim.C_l;
    sim.carreau.lambda = p.lambda / sim.C_t;
    sim.carreau.n = p.n;
    sim.carreau.a = p.a;
}

// Advance an opening by dt and impose its boundary value on the CPU lattice. With the GPU
// lattice, imposeOpenings() writes all values to the GPU at the end instead.
void advanceOpening(OpeningHandler& o, T dt)
{
    o.advanceTime(dt);
#ifdef HEMOFLOW_GPU
    if (gpu)
        return;
#endif
    o.imposeValues(*lattice);
}

// Advance the openings by dt and impose their boundary values.
void imposeOpenings(T dt)
{
    const T murrayExponent = 3.0;
    T sumInflow = 0.0;          // Net flow rate imposed on velocity openings [lattice units], inflow positive
    T sumDiameterPower = 0.0;   // Sum of d^3 over all other openings

    // Phase I: openings with a prescribed velocity
    for (auto &o : openings) {
        if (o->getOpeningType() == OPENING_VELOCITY) {
            advanceOpening(*o, dt);
            sumInflow += o->getImposedFluxLB();
        }
        else {
            sumDiameterPower += pow(o->getRadius() * 2, murrayExponent);
        }
    }

    // Phase II: distribute the outflow over the Murray outlets by d^3; pressure outlets take the rest.
    // C. Chnafa, O. Brina, V. M. Pereira, and D. A. Steinman, "Better Than Nothing: A Rational Approach for
    // Minimizing the Impact of Outflow Strategy on Cerebrovascular Simulations," AJNR 39(2), 337-343, 2018.
    for (auto &o : openings) {
        if (o->getOpeningType() == OPENING_MURRAY) {
            T flowRate = -sumInflow * pow(o->getRadius() * 2, murrayExponent) / sumDiameterPower;  // Outflow is negative
            o->setBCParameterLB(flowRate / o->getProfileFlux());  // Peak velocity of the discrete profile
            advanceOpening(*o, dt);
        }
        else if (o->getOpeningType() != OPENING_VELOCITY) {
            advanceOpening(*o, dt);   // Pressure (and free outflow) openings
        }
    }
#ifdef HEMOFLOW_GPU
    if (gpu)
        gpu->imposeOpenings(openings);
#endif
}

// GPU mode: the GPU state into the CPU lattice, with the current boundary values, for output,
// checkpoints and diagnostics (io.cpp works on the CPU lattice)
void syncCpuLattice()
{
#ifdef HEMOFLOW_GPU
    if (!gpu)
        return;
    auto start = std::chrono::steady_clock::now();
    gpu->copyTo(*lattice);
    for (auto &o : openings)
        o->imposeValues(*lattice);
    pcout << "GPU state copied to the CPU lattice in "
          << std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count() << " s." << endl;
#endif
}

T averageEnergy()
{
#ifdef HEMOFLOW_GPU
    if (gpu)
        return gpu->averageEnergy();
#endif
    return computeAverageEnergy(*lattice);
}

void collideAndStream()
{
#ifdef HEMOFLOW_GPU
    if (gpu) {
        gpu->collideAndStream();
        return;
    }
#endif
    lattice->collideAndStream();
}

void print_help()
{
    pcout << "Usage: mpirun -np N ./hemoFlow config.xml [OPTIONS]" << endl
          << "Options:" << endl
          << "  -d, --debug         Every 200 steps, print the imposed and measured flow rate of every opening" << endl
          << "  -c, --checkpoint    Restart from the checkpoint in the output directory" << endl
          << "  -h, --help          Show this help message" << endl;
}

// Build the list of boundary conditions from the configuration and the openings in the geometry file.
vector<OpeningSpec> buildOpeningSpecs(XMLreader const& xml, const vector<int>& labels)
{
    int numOpenings = (int)labels.size();
    auto has = [&labels](int label) { return std::find(labels.begin(), labels.end(), label) != labels.end(); };
    vector<OpeningSpec> specs;
    XMLreaderProxy geometry = xml["geometry"];

    if (mode == "murray" || mode == "pzero") {
        // Inlet: label 10 (the centerline inlet), parameters from <opening_0>.
        if (!has(FIRST_OPENING))
            throw std::runtime_error("Mode '" + mode + "' needs the inlet label 10, which is not in the geometry file.");
        if (!hasChild(geometry, "opening_0"))
            throw std::runtime_error("Mode '" + mode + "' needs <geometry><opening_0> with the inlet <parameter> and <timeScaleFunction>.");
        XMLreaderProxy inletXml = geometry["opening_0"];

        OpeningSpec inlet;
        inlet.label = FIRST_OPENING;
        inlet.type = OPENING_VELOCITY;
        inlet.name = "Velocity inlet";
        readRequired(inletXml, "geometry><opening_0", "parameter", inlet.parameter);
        inlet.waveformFile = "";
        readOptional(inletXml, "geometry><opening_0", "timeScaleFunction", inlet.waveformFile);
        specs.push_back(inlet);

        // Smallest outlet (label 11) is the 0 Pa outlet in Murray mode.
        const int pressureLabel = FIRST_OPENING + 1;
        if (mode == "murray" && numOpenings > 1 && !has(pressureLabel))
            throw std::runtime_error("Mode 'murray' needs the pressure outlet label 11, which is not in the geometry file.");

        for (int label : labels) {
            if (label == FIRST_OPENING)
                continue;
            OpeningSpec outlet;
            outlet.label = label;
            outlet.parameter = 0.0;
            outlet.waveformFile = "";
            if (mode == "pzero" || label == pressureLabel) {
                outlet.type = OPENING_PRESSURE;
                outlet.name = "Pressure outlet " + std::to_string(label);
            }
            else {
                outlet.type = OPENING_MURRAY;
                outlet.name = "Murray outlet " + std::to_string(label);
            }
            specs.push_back(outlet);
        }
        return specs;
    }

    // Default: every opening is described by <opening_i>, i = 0 .. N-1
    std::set<int> seenLabels;
    for (int o = 0; o < numOpenings; o++) {
        string tag = "opening_" + std::to_string(o);
        string where = "geometry><" + tag;
        if (!hasChild(geometry, tag))
            throw std::runtime_error("The geometry file has " + std::to_string(numOpenings) + " openings, but <geometry><" + tag +
                                     "> is missing. Describe every opening as <opening_0> ... <opening_" +
                                     std::to_string(numOpenings - 1) + "> or set a <mode>.");
        XMLreaderProxy openingXml = geometry[tag];

        OpeningSpec spec;
        int type = 0;
        readRequired(openingXml, where, "label", spec.label);
        readRequired(openingXml, where, "type", type);
        if (type < OPENING_VELOCITY || type > OUTLET_FREEFLOW)
            throw std::runtime_error("<" + where + "><type> is " + std::to_string(type) + "; valid types are 1 (velocity), 2 (Murray), 3 (pressure), 4 (free outflow).");
        spec.type = static_cast<OpeningType>(type);
        spec.name = "Opening " + std::to_string(spec.label);
        readOptional(openingXml, where, "name", spec.name);
        spec.parameter = 0.0;
        if (spec.type == OPENING_VELOCITY)
            readRequired(openingXml, where, "parameter", spec.parameter);
        else if (spec.type == OPENING_PRESSURE)
            readOptional(openingXml, where, "parameter", spec.parameter);
        spec.waveformFile = "";
        readOptional(openingXml, where, "timeScaleFunction", spec.waveformFile);

        if (!has(spec.label))
            throw std::runtime_error("<" + where + "> has label " + std::to_string(spec.label) + ", which is not in the geometry file.");
        if (!seenLabels.insert(spec.label).second)
            throw std::runtime_error("Label " + std::to_string(spec.label) + " is used by more than one <opening_i>.");
        specs.push_back(spec);
    }
    if (hasChild(geometry, "opening_" + std::to_string(numOpenings)))
        pcout << "WARNING: <opening_" << numOpenings << "> and higher are ignored; the geometry file has only "
              << numOpenings << " openings." << endl;
    return specs;
}

// *** Steps of main(), in this order

void printBanner()
{
    pcout   << "********************************* " << endl
            << "*        hemoFlow  v0.31        * " << endl
            << "********************************* " << endl;
    pcout << "Floating-point precision: " << (sizeof(T) == sizeof(float) ? "single" : "double") << endl;
#ifdef HEMOFLOW_GPU
    pcout << "GPU build (in development, doc/dev/TODO-GPU.md)." << endl;
#endif
}

// Command-line options. Returns false if the program has to stop, with its exit code.
bool parseCommandLine(int argc, char *argv[], int& exitCode)
{
    // *** Reading in command line arguments
    if (argc < 2)
    {
        pcout << "Not enough arguments" << endl;
        print_help();
        exitCode = -1;
        return false;
    }

    for (int i = 1; i < argc; ++i)
    {
        std::string arg = argv[i];

        if (arg == "-h" || arg == "--help")
        {
            print_help();
            exitCode = 0;
            return false;
        }
        else if (i == 1)
        {
            continue;   // Configuration file
        }
        else if (arg == "-d" || arg == "--debug")
        {
            DEBUG = true;
            pcout << "Debug mode enabled" << endl;
        }
        else if (arg == "-c" || arg == "--checkpoint")
        {
            CHECKPOINTING = true;
            pcout << "Restarting from checkpoint" << endl;
        }
        else
        {
            pcout << "Unknown argument: " << arg << endl;
            print_help();
            exitCode = 1;
            return false;
        }
    }
    return true;
}

// <simulation>: run control, output fields; <fluid>: fluid model and properties
void readSimulationSection(XMLreader const& xml)
{
    XMLreaderProxy simXml = xml["simulation"];

    // Main simulation parameters
    readRequired(simXml, "simulation", "outputDir", outputFolder);
    readRequired(simXml, "simulation", "cycle", numCycles);
    readRequired(simXml, "simulation", "saveFrequency", saveFreqTime);
    readOptional(simXml, "simulation", "blockSize", blockSize);
    readOptional(simXml, "simulation", "sparse", useSparse);
    useCheckpoint = readOptional(simXml, "simulation", "checkpointFrequency", checkpointFreqTime);
    if (!useCheckpoint)
        pcout << "Note: no <checkpointFrequency> in the config, checkpointing is disabled." << std::endl;
    mode = "";
    readOptional(simXml, "simulation", "mode", mode);
    if (mode != "" && mode != "murray" && mode != "pzero")
        throw std::runtime_error("Unknown boundary condition <mode> '" + mode + "'. Valid modes: " + VALID_MODES + ".");
    outputCoordinatesGiven = readOptional(simXml, "simulation", "outputCoordinates", outputCoordinates);
    if (outputCoordinates != "stl" && outputCoordinates != "m")
        throw std::runtime_error("Unknown <simulation><outputCoordinates> '" + outputCoordinates + "' (valid: stl, m).");
    readOptional(simXml, "simulation", "device", deviceOption);
    if (deviceOption != "auto" && deviceOption != "gpu" && deviceOption != "cpu")
        throw std::runtime_error("Unknown <simulation><device> '" + deviceOption + "' (valid: auto, gpu, cpu).");
#ifndef HEMOFLOW_GPU
    if (deviceOption == "gpu")
        throw std::runtime_error("<simulation><device>gpu</device> needs the GPU build (cmake -DHEMOFLOW_GPU=ON); this build runs on the CPU.");
#endif

    if (numCycles < 1)
        throw std::runtime_error("<simulation><cycle> must be at least 1.");
    if (saveFreqTime <= 0.0)
        throw std::runtime_error("<simulation><saveFrequency> must be positive.");
    if (blockSize < 4)
        throw std::runtime_error("<simulation><blockSize> must be at least 4.");

    // Output fields
    string fieldList;
    if (readOptional(simXml, "simulation", "outputFields", fieldList)) {
        outputOptions.velocity = outputOptions.pressure = outputOptions.stress = outputOptions.viscosity = false;
        std::istringstream fields(fieldList);
        string field;
        while (fields >> field) {
            if (field == "velocity") outputOptions.velocity = true;
            else if (field == "pressure") outputOptions.pressure = true;
            else if (field == "stress") outputOptions.stress = true;
            else if (field == "viscosity") outputOptions.viscosity = true;
            else throw std::runtime_error("Unknown output field '" + field + "' in <simulation><outputFields> (valid: velocity pressure stress viscosity).");
        }
    }

    // Fluid properties (optional)
    if (hasSection(xml, "fluid")) {
        XMLreaderProxy fluidXml = xml["fluid"];
        readOptional(fluidXml, "fluid", "model", fluidModel);
        readOptional(fluidXml, "fluid", "density", fluidDensity);
        bool hasNu = readOptional(fluidXml, "fluid", "kinematicViscosity", fluidKinematicViscosity);
        bool hasCarreau = false;
        hasCarreau |= readOptional(fluidXml, "fluid", "mu0", carreauSI.mu0);
        hasCarreau |= readOptional(fluidXml, "fluid", "muInf", carreauSI.muInf);
        hasCarreau |= readOptional(fluidXml, "fluid", "lambda", carreauSI.lambda);
        hasCarreau |= readOptional(fluidXml, "fluid", "n", carreauSI.n);
        hasCarreau |= readOptional(fluidXml, "fluid", "a", carreauSI.a);
        if (fluidModel == "carreau" && hasNu)
            throw std::runtime_error("<fluid><kinematicViscosity> is not used by the carreau model; set mu0 and muInf [Pa s] instead.");
        if (fluidModel == "newtonian" && hasCarreau)
            throw std::runtime_error("<fluid> has Carreau parameters (mu0, muInf, lambda, n, a), but <model> is not 'carreau'.");
    }
    if (fluidModel != "newtonian" && fluidModel != "carreau")
        throw std::runtime_error("Unknown <fluid><model> '" + fluidModel + "' (valid: newtonian, carreau).");
    if (fluidKinematicViscosity <= 0.0 || fluidDensity <= 0.0)
        throw std::runtime_error("The fluid viscosity and density must be positive.");
    if (fluidModel == "carreau") {
        if (!(carreauSI.muInf > 0.0 && carreauSI.mu0 >= carreauSI.muInf))
            throw std::runtime_error("Carreau model: 0 < muInf <= mu0 is required.");
        if (!(carreauSI.lambda >= 0.0))
            throw std::runtime_error("Carreau model: lambda must not be negative.");
        if (!(carreauSI.n > 0.0 && carreauSI.n <= 1.0))
            throw std::runtime_error("Carreau model: 0 < n <= 1 is required (shear-thinning).");
        if (!(carreauSI.a > 0.0))
            throw std::runtime_error("Carreau model: the Yasuda exponent a must be positive.");
    #if LES
        throw std::runtime_error("The carreau fluid model cannot be combined with LES (LES 1 in globals.h).");
    #endif
        // The time step and the stability checks use the lowest viscosity
        fluidKinematicViscosity = carreauSI.muInf / fluidDensity;
    }
}

// Output folder, created by the main process. Returns its path.
std::string createOutputFolder()
{
    // Output directory, created by the main process
    std::string outDir = resolvePath(outputFolder);
    int dirStatus = 1;
    if (global::mpi().isMainProcessor()) {
        if (dirExists(outDir) <= 0) {
            pcout << "Output folder " << outDir << " does not exist! Creating it...." << std::endl;
            dirStatus = (mkpath(outDir.c_str(), 0777) == 0) ? 1 : 0;
        }
    }
    global::mpi().bCast(&dirStatus, 1);
    global::mpi().barrier();
    if (!dirStatus)
        throw std::runtime_error("Cannot create the output folder " + outDir);
    pcout << "Output folder: " << outDir+"/" << std::endl;
    global::directories().setOutputDir(outDir+"/");
    return outDir;
}

// Domain size, output frame, flow diverter coefficients and lattice Boltzmann parameters
void setUpFromGeometry(XMLreader const& xml, GeometryFile const& geometry)
{
    XMLreaderProxy simXml = xml["simulation"];

    Nx = geometry.nx;
    Ny = geometry.ny;
    Nz = geometry.nz;
    pcout << "Domain size: " << Nx << " x " << Ny << " x " << Nz << std::endl;
    double sim_dx = geometry.dx;
    pcout << "Resolution [m]: " << sim_dx << std::endl;

    // Physical position of voxel (0,0,0), used for the output (older geometry files have none)
    if (geometry.hasOrigin) {
        outputOptions.origin = geometry.origin;
        pcout << "Origin [m]: " << geometry.origin[0] << " " << geometry.origin[1] << " " << geometry.origin[2] << std::endl;
    }
    else {
        pcout << "Note: the geometry file has no 'origin'; the output is written in the voxel frame (origin 0,0,0)." << std::endl;
    }

    // Units of the output grid: by default those of the input STL (e.g. mm), so that the
    // output overlays the STL. Only the grid (spacing, origin) is scaled; the fields stay SI.
    if (outputCoordinates == "stl" && geometry.hasLengthUnit) {
        outputOptions.lengthUnit = geometry.lengthUnit;
        pcout << "Output grid in the coordinates of the input STL (1 unit = " << geometry.lengthUnit
              << " m); the fields are in SI units." << std::endl;
    }
    else if (outputCoordinates == "stl" && outputCoordinatesGiven) {
        throw std::runtime_error("The geometry file has no 'lengthUnit' (older preprocessor), so <outputCoordinates>stl</outputCoordinates> "
                                 "is not possible. Re-run the preprocessor or use <outputCoordinates>m</outputCoordinates>.");
    }
    else {
        if (!geometry.hasLengthUnit)
            pcout << "Note: the geometry file has no 'lengthUnit' (older preprocessor); the output grid is in m." << std::endl;
        else
            pcout << "Output grid in m." << std::endl;
    }

    // Flow diverter (if there is one)
    if (geometry.hasStent) {
        pcout << "Found flow diverter information to load." << std::endl;
        if (!hasSection(xml, "flowdiverter"))
            throw std::runtime_error("The geometry contains a flow diverter, but the config has no <flowdiverter> section.");
        XMLreaderProxy fdXml = xml["flowdiverter"];
        readRequired(fdXml, "flowdiverter", "linCoeff", linCoeff);
        readRequired(fdXml, "flowdiverter", "quadCoeff", quadCoeff);
    }

    // Time step (optional)
    double sim_dt = -1.0;
    if (readOptional(simXml, "simulation", "dt", sim_dt) && sim_dt <= 0.0)
        throw std::runtime_error("<simulation><dt> must be positive.");

    pcout << "Setting LBM parameters..." << std::endl;
    calcSimulationParameters(sim, sim_dx, sim_dt, fluidKinematicViscosity, fluidDensity);
    if (fluidModel == "carreau")
        setCarreauModel(sim, carreauSI, fluidDensity);
}

// Boundary conditions to create, from the configuration and the openings of the geometry file
vector<OpeningSpec> readOpenings(XMLreader const& xml, GeometryFile const& geometry)
{
    // **** Openings in the geometry file ****
    pcout << "Processing openings..." << std::endl;
    const vector<int>& openingLabels = geometry.openingLabels;
    int numOpenings = (int)openingLabels.size();
    if (numOpenings < 2)
        throw std::runtime_error("The geometry has fewer than 2 openings.");

    pcout << "Number of openings in geometry: " << numOpenings << std::endl;
    pcout << "Opening labels: ";
    for (int label : openingLabels)
        pcout << label << " ";
    pcout << endl;
    pcout << "Boundary condition mode: " << (mode.empty() ? "xml (default)" : mode) << std::endl;

    return buildOpeningSpecs(xml, openingLabels);
}

// The lattice with the sparse or dense block decomposition and the bulk dynamics
void createLattice(GeometryFile& geometry, vector<OpeningSpec> const& specs)
{
    // **** Lattice ****
    // The finite-difference ("interp") pressure boundary reads data two cells deep across
    // block boundaries. With a one-cell envelope the result would depend on how the domain
    // is split into blocks wherever a block boundary cuts a pressure opening.
    bool hasPressureOpening = std::any_of(specs.begin(), specs.end(),
                                          [](OpeningSpec const& spec) { return spec.type == OPENING_PRESSURE; });
    envelopeWidth = hasPressureOpening ? 2 : 1;
    pcout << "Envelope width: " << envelopeWidth << endl;

    if (gpuMode) {
        // One dense block per GPU rank. The GPU-side packing of MPI messages needs the blocking communicator.
        defaultMultiBlockPolicy3D().toggleBlockingCommunication(true);
        // Development only (doc/dev/TODO-GPU.md, Phases 5-6): HEMOFLOW_GPU_BLOCK_SIZE=n puts the
        // sparse decomposition with blocks of n on the GPU, to test and time several blocks
        const char* gpuBlockSize = std::getenv("HEMOFLOW_GPU_BLOCK_SIZE");
        if (gpuBlockSize != nullptr && std::atoi(gpuBlockSize) > 0) {
            blockSize = std::atoi(gpuBlockSize);
            useSparse = true;
            pcout << "GPU mode: sparse decomposition with blocks of " << blockSize
                  << " (HEMOFLOW_GPU_BLOCK_SIZE, for development)." << endl;
        }
        else {
            if (useSparse)
                pcout << "GPU mode: one dense block per GPU (<sparse> and <blockSize> apply to CPU runs)." << endl;
            useSparse = false;
        }
    }

    // Background dynamics: the fluid's on the CPU, the outside voxels' in GPU mode (see FluidDynamicsFunctional3D)
    Dynamics<T, DESCRIPTOR>* bulkDynamics = gpuMode ? new NoDynamics<T, DESCRIPTOR>(1.0) : newBulkDynamics();

    if (useSparse) {
        pcout << "Setting simulation domain mask for sparse decomposition (block size " << blockSize << ")..." << endl;
        // Keep every block that contains a wall, fluid or opening voxel (label > UNUSED).
        MultiScalarField3D<int> occupancy(Nx, Ny, Nz);
        geometry.fillOccupancy(occupancy);

        pcout << "Creating sparse representation ..." << endl;
        MultiBlockManagement3D sparseBlockManagement =
            computeSparseManagement(*plb::reparallelize(occupancy, blockSize, blockSize, blockSize), envelopeWidth);

        lattice = new MultiBlockLattice3D<T, DESCRIPTOR>(sparseBlockManagement,
                                                    defaultMultiBlockPolicy3D().getBlockCommunicator(),
                                                    defaultMultiBlockPolicy3D().getCombinedStatistics(),
                                                    defaultMultiBlockPolicy3D().getMultiCellAccess<T,DESCRIPTOR>(),
                                                    bulkDynamics);
    }
    else {
        lattice = new MultiBlockLattice3D<T, DESCRIPTOR>(
            defaultMultiBlockPolicy3D().getMultiBlockManagement(Box3D(0, Nx - 1, 0, Ny - 1, 0, Nz - 1), envelopeWidth),
            defaultMultiBlockPolicy3D().getBlockCommunicator(),
            defaultMultiBlockPolicy3D().getCombinedStatistics(),
            defaultMultiBlockPolicy3D().getMultiCellAccess<T,DESCRIPTOR>(),
            bulkDynamics);
    }

    #if LES
        instantiateStaticSmagorinsky(*lattice, lattice->getBoundingBox(), cSmago);
    #endif

    pcout << getMultiBlockInfo(*lattice) << endl;
    plint numBlocks = lattice->getMultiBlockManagement().getSparseBlockStructure().getNumBlocks();
    if (useSparse && numBlocks < global::mpi().getSize())
        pcout << "Note: " << numBlocks << " blocks for " << global::mpi().getSize()
              << " processes, so some processes have no work. Use a smaller <blockSize>." << endl;
}

// Voxel labels, porous multipliers and the boundary conditions of the openings
void createFieldsAndOpenings(GeometryFile& geometry, vector<OpeningSpec> const& specs)
{
    const vector<int>& openingLabels = geometry.openingLabels;

    // Voxel labels and porosity as distributed fields with the block structure of the lattice
    labelField = generateMultiScalarField<int>(*lattice, envelopeWidth).release();
    geometry.fillLabels(*labelField);
    if (geometry.hasStent && geometry.numPorous == 0) {
        pcout << "Note: the flow diverter has no voxels in the lumen, so it has no effect." << std::endl;
    } else if (geometry.hasStent) {
        porosityField = generateMultiTensorField<T, 2>(*lattice, envelopeWidth).release();
        geometry.fillPorosity(*porosityField);
        porousBox = geometry.porousBox;
        pcout << "Flow diverter: " << geometry.numPorous << " porous voxels. Resistance multipliers ";
        if (geometry.hasMultipliers)
            pcout << "(min / mean / max): linear " << geometry.multiplierMin[0] << " / " << geometry.multiplierMean[0]
                  << " / " << geometry.multiplierMax[0] << ", quadratic " << geometry.multiplierMin[1] << " / "
                  << geometry.multiplierMean[1] << " / " << geometry.multiplierMax[1] << std::endl;
        else
            pcout << "1 (the geometry file has none: older preprocessor)." << std::endl;
    }
    geometry.release();

    // GPU mode: the fluid dynamics in the fluid voxels, before the openings wrap them
    if (gpuMode)
        applyProcessingFunctional(new FluidDynamicsFunctional3D(newBulkDynamics()), lattice->getBoundingBox(), *lattice, *labelField);

    // **** Boundary conditions ****
    for (auto const& spec : specs) {
        int idx = (int)(std::find(openingLabels.begin(), openingLabels.end(), spec.label) - openingLabels.begin());
        pcout << "Processing opening with label " << spec.label << std::endl;

        auto *opening = new OpeningHandler(geometry.openingNodes.at(spec.label), spec.label, spec.type,
                                           geometry.openingRadius[idx] / sim.C_l, geometry.openingNormal[idx]);
        opening->setName(spec.name);
        opening->setBCType(*lattice);

        if (spec.type == OPENING_VELOCITY)
            opening->setBCParameter(spec.parameter * 2, sim);   // Mean velocity -> peak velocity of a paraboloid
        else
            opening->setBCParameter(spec.parameter, sim);

        if (!spec.waveformFile.empty())
            opening->loadScaleFunction(resolvePath(spec.waveformFile));

        if (opening->isVelocityType())
            opening->createPoiseuilleProfile();   // Peak velocity 1 (lattice units)
        else if (spec.type == OPENING_PRESSURE)
            opening->createConstantPressureProfile();

        opening->printOpeningDetails(sim);
        vec3d c = opening->getCenter();
        double u = outputOptions.lengthUnit;
        pcout << "-> Center [output grid coordinates]: " << (outputOptions.origin[0] + c.x * sim.C_l) / u << " "
              << (outputOptions.origin[1] + c.y * sim.C_l) / u << " " << (outputOptions.origin[2] + c.z * sim.C_l) / u << std::endl;
        openings.push_back(opening);
    }
}

// Cardiac cycle period and simulated time
void setCyclePeriod()
{
    // Derive the cardiac cycle period from the first opening with a waveform
    cyclePeriod = -1.0;
    for (auto &o : openings) {
        if (o->getHasScaleFunction()) {
            if (cyclePeriod < 0.0)
                cyclePeriod = o->getPeriod();
            else if (std::fabs(o->getPeriod() - cyclePeriod) > 1e-12)
                pcout << "WARNING: the waveform of opening " << o->getLabel() << " has a period of " << o->getPeriod()
                      << " s; the cycle length is taken from the first waveform (" << cyclePeriod << " s)." << std::endl;
        }
    }
    if (cyclePeriod <= 0.0)
        throw std::runtime_error("No opening has a <timeScaleFunction>, so the cardiac cycle period for <cycle> is undefined. "
                                 "Give the inlet a waveform file (a two-point constant file for steady inflow).");
    simLength = numCycles * cyclePeriod;
    pcout << "Cardiac cycle period: " << cyclePeriod << " s. Running " << numCycles << " cycles -> simLength = " << simLength << " s." << std::endl;
}

void printSimulationParameters()
{
    pcout   << "*********** Simulation parameters *********** " << endl
            << "size [LU]:   " << Nx << "x" << Ny << "x" << Nz << endl
            << "dx [m]:  " << sim.C_l << endl
            << "dt [s]:  " << sim.C_t << endl
            << "dm [kg]:  " << sim.C_m << endl
            << "fluid model: " << fluidModel << endl;
    if (sim.nonNewtonian) {
        // omega lies between its values at the zero-shear and the infinite-shear viscosity
        T omega0 = 1.0 / (3.0 * sim.carreau.nu0 + 0.5);
        pcout   << "Carreau-Yasuda: mu0 = " << carreauSI.mu0 << " Pa s, muInf = " << carreauSI.muInf << " Pa s, lambda = "
                << carreauSI.lambda << " s, n = " << carreauSI.n << ", a = " << carreauSI.a << endl
                << "omega range: " << omega0 << " (zero shear, tau = " << 1.0 / omega0 << ") to "
                << sim.omega << " (infinite shear, tau = " << 1.0 / sim.omega << ")" << endl;
    }
    else {
        pcout   << "nu [m^2/s]: " << fluidKinematicViscosity << endl
                << "omega:  " << sim.omega << endl
                << "nu [lb]: " << 1./3. * (1./sim.omega - 0.5) << endl;
    }
    pcout   << "rho [kg/m^3]: " << fluidDensity << endl
            << "U_max [lbm]: " << sim.U_max_lb << endl
            << "U_max [m/s]: " << sim.U_max_lb * sim.C_l / sim.C_t << endl << endl;
}

#ifdef HEMOFLOW_GPU
// Why the case cannot run on the GPU yet (empty: it can). doc/dev/TODO-GPU.md.
std::string gpuUnsupportedReason(GeometryFile const& geometry, vector<OpeningSpec> const& specs)
{
    if (LES)
        return "LES is not implemented on the GPU";
    for (auto const& spec : specs) {
        if (spec.type == OUTLET_FREEFLOW)
            return "free-outflow openings are not implemented on the GPU";
    }
    // Opening nodes on an edge of the bounding box (geometry files of older preprocessors)
    int N[3] = {geometry.nx, geometry.ny, geometry.nz};
    for (auto const& entry : geometry.openingNodes)
        for (Index3D const& node : entry.second) {
            int c[3] = {node.x, node.y, node.z};
            int faces = 0;
            for (int a = 0; a < 3; a++)
                faces += (c[a] == 0 || c[a] == N[a] - 1) ? 1 : 0;
            if (faces > 1)
                return "opening " + std::to_string(entry.first) + " touches an edge of the bounding box (older geometry file)";
        }
    return "";
}

// Where the time loop runs (collective): on the GPU, one dense block per GPU rank, if the case
// is supported there, world rank 0 has a GPU and every block fits its GPU's free memory;
// otherwise on the CPU with the usual decomposition. <device>gpu</device> makes a missing GPU
// or an unsupported case an error, <device>cpu</device> skips the GPU.
//
// GPU ranks: on each node the ranks, in rank order, take the first GPU they can see that no
// lower rank of the node has taken; ranks without one wait. GPUs are told apart by their UUID,
// so this works whether every rank sees all GPUs of the node or the launcher gives each rank its
// own (CUDA_VISIBLE_DEVICES, Slurm --gpus-per-task=1; recommended, see the README). Development
// only: HEMOFLOW_GPU_RANKS=n makes world ranks 0..n-1 the GPU ranks even if they share a GPU.
bool chooseGpu(GeometryFile const& geometry, vector<OpeningSpec> const& specs)
{
    // The same on all ranks (the geometry and the configuration are)
    std::string reason = (deviceOption == "cpu") ? "<device>cpu</device>" : gpuUnsupportedReason(geometry, specs);
    if (!reason.empty()) {
        if (deviceOption == "gpu")
            throw std::runtime_error("<simulation><device>gpu</device>, but " + reason + ".");
        pcout << "Execution: CPU (" << reason << ")." << endl;
        return false;
    }
    MPI_Comm world = global::mpi().getGlobalCommunicator();
    int rank = 0, size = 1;
    MPI_Comm_rank(world, &rank);
    MPI_Comm_size(world, &size);

    // The GPUs of each rank, gathered on the first rank of every node
    struct Visible {
        int count;
        char uuid[16][36];
    };
    Visible mine {};
    std::vector<GpuDevice> devices = visibleGpus();
    mine.count = (int)std::min<std::size_t>(devices.size(), 16);
    for (int d = 0; d < mine.count; ++d)
        std::snprintf(mine.uuid[d], sizeof(mine.uuid[d]), "%s", devices[d].uuid.c_str());
    MPI_Comm node;
    MPI_Comm_split_type(world, MPI_COMM_TYPE_SHARED, rank, MPI_INFO_NULL, &node);
    int nodeRank = 0, nodeSize = 1;
    MPI_Comm_rank(node, &nodeRank);
    MPI_Comm_size(node, &nodeSize);
    std::vector<Visible> all(nodeRank == 0 ? nodeSize : 0);
    MPI_Gather(&mine, sizeof(Visible), MPI_BYTE, all.data(), sizeof(Visible), MPI_BYTE, 0, node);
    // Per rank of the node: its GPU (index in its own list of visible GPUs, -1: none)
    std::vector<int> assigned(nodeSize, -1);
    if (nodeRank == 0) {
        std::set<std::string> taken;
        for (int r = 0; r < nodeSize; ++r)
            for (int d = 0; d < all[r].count; ++d)
                if (taken.insert(all[r].uuid[d]).second) {
                    assigned[r] = d;
                    break;
                }
    }
    int device = -1;
    MPI_Scatter(assigned.data(), 1, MPI_INT, &device, 1, MPI_INT, 0, node);
    const char* forcedRanks = std::getenv("HEMOFLOW_GPU_RANKS");
    if (forcedRanks != nullptr && std::atoi(forcedRanks) > 0) {
        if (rank < std::atoi(forcedRanks) && device < 0 && mine.count > 0)
            device = 0;   // Shares a GPU
        else if (rank >= std::atoi(forcedRanks))
            device = -1;
        if (rank == 0)
            pcout << "GPU ranks: the first " << std::atoi(forcedRanks) << " (HEMOFLOW_GPU_RANKS, for development)." << endl;
    }
    int isGpu = device >= 0 ? 1 : 0;
    std::vector<int> gpuFlags(size, 0);
    MPI_Allgather(&isGpu, 1, MPI_INT, gpuFlags.data(), 1, MPI_INT, world);
    int numGpuRanks = std::accumulate(gpuFlags.begin(), gpuFlags.end(), 0);
    int gpuIndex = std::accumulate(gpuFlags.begin(), gpuFlags.begin() + rank, 0);   // Block of this rank
    // Start-up contexts: every rank has one on the first GPU it sees; waiting ranks release theirs
    std::vector<int> waitingFlags(nodeRank == 0 ? nodeSize : 0);
    int waiting = isGpu ? 0 : 1;
    MPI_Gather(&waiting, 1, MPI_INT, waitingFlags.data(), 1, MPI_INT, 0, node);

    // One block per GPU rank (the dense decomposition of createLattice), with its envelope (two
    // cells with pressure openings) and the list of porous cells (all of them, an upper bound),
    // plus 5 % and 256 MB for the CUDA context and the runtime
    bool hasPressure = std::any_of(specs.begin(), specs.end(), [](OpeningSpec const& spec) { return spec.type == OPENING_PRESSURE; });
    int e = hasPressure ? 2 : 1;
    struct RankGpu {
        int gpuIndex;
        char uuid[36], bus[16], name[64], host[64];
        double needed, freeBytes, totalBytes, released;
        long long nx, ny, nz;
    };
    RankGpu info {};
    info.gpuIndex = isGpu ? gpuIndex : -1;
    if (isGpu && numGpuRanks > 0) {
        SparseBlockStructure3D blocks = createRegularDistribution3D(Box3D(0, geometry.nx - 1, 0, geometry.ny - 1, 0, geometry.nz - 1), numGpuRanks);
        Box3D bulk = blocks.getBulks().at(gpuIndex);
        info.nx = bulk.getNx() + 2 * e;
        info.ny = bulk.getNy() + 2 * e;
        info.nz = bulk.getNz() + 2 * e;
        double porousCells = geometry.hasStent ? (double)geometry.numPorous : 0.0;
        info.needed = ((double)info.nx * info.ny * info.nz * gpuBytesPerCell() + porousCells * gpuBytesPerPorousCell()) * 1.05 + 256e6;
        std::string name;
        std::size_t freeBytes = 0, totalBytes = 0;
        if (selectGpu(device) && queryGpu(name, freeBytes, totalBytes)) {
            info.freeBytes = (double)freeBytes;
            info.totalBytes = (double)totalBytes;
            std::snprintf(info.name, sizeof(info.name), "%s", name.c_str());
        }
        std::snprintf(info.uuid, sizeof(info.uuid), "%s", devices[device].uuid.c_str());
        std::snprintf(info.bus, sizeof(info.bus), "%s", devices[device].pciBus.c_str());
    }
    // Released contexts per GPU of the node, credited to the GPU ranks there (100 MB each, measured about 120 MB)
    std::vector<std::string> startUpUuid(nodeRank == 0 ? nodeSize : 0);
    if (nodeRank == 0)
        for (int r = 0; r < nodeSize; ++r)
            startUpUuid[r] = all[r].count > 0 ? all[r].uuid[0] : "";
    {
        std::vector<double> credit(nodeSize, 0.0);
        if (nodeRank == 0)
            for (int r = 0; r < nodeSize; ++r)
                if (waitingFlags[r] && !startUpUuid[r].empty())
                    for (int g = 0; g < nodeSize; ++g)
                        if (assigned[g] >= 0 && all[g].uuid[assigned[g]] == startUpUuid[r]) {
                            credit[g] += 100e6;   // Once per GPU: the first GPU rank on it
                            break;
                        }
        MPI_Scatter(credit.data(), 1, MPI_DOUBLE, &info.released, 1, MPI_DOUBLE, 0, node);
    }
    gethostname(info.host, sizeof(info.host) - 1);
    MPI_Comm_free(&node);
    std::vector<RankGpu> ranks(rank == 0 ? size : 0);
    MPI_Gather(&info, sizeof(RankGpu), MPI_BYTE, ranks.data(), sizeof(RankGpu), MPI_BYTE, 0, world);

    // Decision on world rank 0: every GPU must hold the blocks of its ranks
    int use = 0;
    std::ostringstream message, details;
    if (rank == 0) {
        message.precision(2);
        details.precision(2);
        message << std::fixed;
        details << std::fixed;
        struct Usage {
            double needed = 0.0, available = 0.0, total = 0.0;
        };
        std::map<std::string, Usage> perGpu;   // host/uuid -> the ranks on that GPU
        for (RankGpu const& r : ranks)
            if (r.gpuIndex >= 0) {
                Usage& u = perGpu[std::string(r.host) + "/" + r.uuid];
                u.needed += r.needed;
                u.available = std::max(u.available, r.freeBytes + r.released);
                u.total = r.totalBytes;
            }
        std::string fullest;
        double worst = -1.0;
        for (auto const& p : perGpu)
            if (p.second.needed / std::max(p.second.available, 1.0) > worst) {
                worst = p.second.needed / std::max(p.second.available, 1.0);
                fullest = p.first;
            }
        if (reason.empty() && (ranks[0].gpuIndex < 0 || ranks[0].totalBytes == 0.0))
            reason = "no GPU found";
        if (reason.empty()) {
            Usage const& u = perGpu[fullest];
            std::ostringstream memory;
            memory.precision(2);
            memory << std::fixed << u.needed / 1e9 << " GB needed, " << u.available / 1e9 << " GB of " << u.total / 1e9
                   << " GB available";
            if (worst > 1.0)
                reason = (numGpuRanks == 1 ? "the domain does not fit the GPU memory in one block ("
                                           : "the blocks do not fit the GPU memory (fullest GPU: ") + memory.str() + ")";
            else if (numGpuRanks == 1)
                message << "Execution: GPU (" << ranks[0].name << "), one block of " << ranks[0].nx - 2 * e << " x "
                        << ranks[0].ny - 2 * e << " x " << ranks[0].nz - 2 * e << " cells, " << memory.str() << ".";
            else {
                message << "Execution: " << perGpu.size() << (perGpu.size() == 1 ? " GPU, " : " GPUs, ") << numGpuRanks << " ranks with one block each ("
                        << "fullest GPU: " << memory.str() << ").";
                if ((int)perGpu.size() < numGpuRanks)
                    message << " Several ranks share a GPU.";
                for (int r = 0; r < size; ++r)
                    if (ranks[r].gpuIndex >= 0)
                        details << std::endl << "  rank " << r << " on " << ranks[r].host << ": " << ranks[r].name << " (PCI "
                                << ranks[r].bus << "), block of " << ranks[r].nx - 2 * e << " x " << ranks[r].ny - 2 * e << " x "
                                << ranks[r].nz - 2 * e << " cells, " << ranks[r].needed / 1e9 << " GB";
            }
        }
        use = reason.empty() ? 1 : 0;
    }
    MPI_Bcast(&use, 1, MPI_INT, 0, world);
    gpuRank = use && isGpu;
    if (use) {
        pcout << message.str() << details.str() << endl;
        return true;
    }
    if (deviceOption == "gpu") {
        int length = (int)reason.size();
        MPI_Bcast(&length, 1, MPI_INT, 0, world);
        reason.resize(length);
        MPI_Bcast(&reason[0], length, MPI_CHAR, 0, world);
        throw std::runtime_error("<simulation><device>gpu</device>, but " + reason + ".");
    }
    pcout << "Execution: CPU (" << reason << ")." << endl;
    return false;
}

// GPU mode: Palabos continues on the GPU ranks, with their own communicator. The other ranks
// release their GPU and wait, without spinning, for the end of the run (finishRun); they do not
// return.
void enterGpuMode()
{
    gpuMode = true;
    MPI_Comm world = global::mpi().getGlobalCommunicator();
    int rank = 0, size = 1;
    MPI_Comm_rank(world, &rank);
    MPI_Comm_size(world, &size);
    if (size > 1) {
        MPI_Comm gpuCommunicator;
        MPI_Comm_split(world, gpuRank ? 0 : MPI_UNDEFINED, rank, &gpuCommunicator);
        if (!gpuRank) {
            releaseGpu();
            MPI_Request request;
            MPI_Ibarrier(world, &request);
            int done = 0;
            while (!done) {
                MPI_Test(&request, &done, MPI_STATUS_IGNORE);
                if (!done)
                    usleep(100000);
            }
            int code = 0;
            MPI_Bcast(&code, 1, MPI_INT, 0, world);
            MPI_Finalize();
            std::_Exit(code);
        }
        worldCommunicator = world;
        global::mpi().init(gpuCommunicator);
        int waiting = size - global::mpi().getSize();
        if (waiting > 0)
            pcout << "The other " << waiting << " MPI ranks wait for the end of the GPU run." << endl;
    }
    // One block per GPU rank in the dense decomposition (the policy may have counted the world)
    defaultMultiBlockPolicy3D().setNumProcesses(global::mpi().getSize());
}

// GPU mode: a CPU lattice with placeholder dynamics (gpuSolver.h) and the current state of the
// CPU lattice is copied to the GPU. The CPU lattice keeps the real dynamics for output and
// checkpoints. Collective; false (on all ranks) if the GPU lattice fails on any rank.
bool startGpu()
{
    pcout << "Copying the lattice to the GPU..." << endl;
    Box3D box = lattice->getBoundingBox();
    auto *placeholder = new MultiBlockLattice3D<T, DESCRIPTOR>(
        lattice->getMultiBlockManagement(), defaultMultiBlockPolicy3D().getBlockCommunicator(),
        defaultMultiBlockPolicy3D().getCombinedStatistics(), defaultMultiBlockPolicy3D().getMultiCellAccess<T, DESCRIPTOR>(),
        new NoDynamics<T, DESCRIPTOR>(1.0));   // Outside voxels (see FluidDynamicsFunctional3D)
    applyProcessingFunctional(new FluidDynamicsFunctional3D(new BGKdynamics<T, DESCRIPTOR>(sim.omega)), box, *placeholder, *labelField);
    for (auto &o : openings)
        o->setPlaceholderBCType(*placeholder, sim.omega);
    defineDynamics(*placeholder, *labelField, box, new BounceBack<T, DESCRIPTOR>(1.0), (int)WALL);
    copy(*lattice, box, *placeholder, box, modif::staticVariables);
    for (auto &o : openings)
        o->imposeValues(*placeholder);
    GpuPorousMedium porous;
    if (porosityField != nullptr) {
        porous.multipliers = porosityField;
        porous.box = porousBox;
        porous.linCoeff = linCoeff_lb;
        porous.quadCoeff = quadCoeff_lb;
    }
    int failed = 0;
    std::string error;
    try {
        gpu = new GpuSolver(*placeholder, openings, sim, porous);
    }
    catch (std::exception const& exception) {
        failed = 1;
        error = exception.what();
    }
    delete placeholder;
    int anyFailed = failed;
    global::mpi().reduceAndBcast(anyFailed, MPI_MAX);
    if (anyFailed) {
        if (failed)
            std::cerr << "ERROR on rank " << global::mpi().getRank() << ": " << error << std::endl;
        pcout << "ERROR: the GPU lattice could not be set up." << endl;
        return false;
    }

    std::string name;
    std::size_t freeBytes = 0, totalBytes = 0;
    if (queryGpu(name, freeBytes, totalBytes))
        pcout << "GPU memory in use: " << (totalBytes - freeBytes) / 1e9 << " GB of " << totalBytes / 1e9 << " GB." << endl;
    return true;
}
#endif

// GPU mode with idle ranks: release them, with the exit code of the run
int finishRun(int code)
{
#ifdef HEMOFLOW_GPU
    if (worldCommunicator != MPI_COMM_NULL) {
        MPI_Request request;
        MPI_Ibarrier(worldCommunicator, &request);
        MPI_Wait(&request, MPI_STATUS_IGNORE);
        MPI_Bcast(&code, 1, MPI_INT, 0, worldCommunicator);
        worldCommunicator = MPI_COMM_NULL;
    }
#endif
    return code;
}

// Initial state (or checkpoint), time loop, output and checkpoints. Returns the exit code.
int runSimulation(std::string const& outDir)
{
    const int totalSteps = (int)std::llround(simLength / sim.C_t);
    const int saveFrequency = (int)std::llround(saveFreqTime / sim.C_t);
    if (saveFrequency < 1) {
        pcout << "ERROR: <saveFrequency> (" << saveFreqTime << " s) is smaller than the time step (" << sim.C_t << " s)." << endl;
        return -1;
    }
    pcout << "Saving frequency set to every " << saveFreqTime << " s (" << saveFrequency << " steps)." << endl;

    int checkpointFrequency = 0;
    if (useCheckpoint) {
        checkpointFrequency = (int)std::llround(checkpointFreqTime / sim.C_t);
        if (checkpointFrequency < 1) {
            pcout << "ERROR: <checkpointFrequency> (" << checkpointFreqTime << " s) is smaller than the time step (" << sim.C_t << " s)." << endl;
            return -1;
        }
        pcout << "Checkpointing will happen every " << checkpointFreqTime << " s (" << checkpointFrequency << " steps)." << endl;
    }

    const int lastCycleStartIter = (int)std::llround((numCycles - 1) * cyclePeriod / sim.C_t);
    pcout << "Writing output only from iteration " << lastCycleStartIter << " onward (last cycle)." << endl;
    pcout << "Output fields:" << (outputOptions.velocity ? " velocity" : "") << (outputOptions.pressure ? " pressure" : "")
          << (outputOptions.stress ? " stress" : "") << (outputOptions.viscosity ? " viscosity" : "") << endl;

    // Checkpoint file names relative to the output folder
    string chkParamFile = outDir+"/checkpoint_parameters.dat";
    string chkDataFile = outDir+"/checkpoint_lattice.dat";
    string chkParamFileOld = outDir+"/checkpoint_parameters_old.dat";
    string chkDataFileOld = outDir+"/checkpoint_lattice_old.dat";

    // If there is a flow diverter, set up the porous forcing
    if(porosityField != nullptr) {
        pcout << "Setting up porous layer for flow diverter..." << std::endl;
        // Only over the bounding box of the porous voxels: elsewhere all multipliers are 0
        integrateProcessingFunctional( new PorousForceFunctional<T, DESCRIPTOR>(linCoeff_lb, quadCoeff_lb), porousBox, *lattice, *porosityField);
    }

    // Walls and outside voxels. This also overrides the boundary dynamics that the openings
    // instantiated on wall/outside voxels of their bounding rectangles.
    pcout << "Defining walls..." << std::endl;
    // NoDynamics with the reference density, so outside cells report 0 Pa (in GPU mode the background dynamics)
    if (gpuMode)
        applyProcessingFunctional(new OutsideToBackgroundFunctional3D(), lattice->getBoundingBox(), *lattice, *labelField);
    else
        defineDynamics(*lattice, *labelField, lattice->getBoundingBox(), new NoDynamics<T, DESCRIPTOR>(1.0), (int)UNUSED);
    defineDynamics(*lattice, *labelField, lattice->getBoundingBox(), new BounceBack<T, DESCRIPTOR>(1.0), (int)WALL);

    pcout << "Initializing lattice in equilibrium..." << std::endl;
    initializeAtEquilibrium (*lattice, lattice->getBoundingBox(), (T)1.0, Array<T,3>((T)0.,(T)0.,(T)0.) );

    pcout << "Finalizing lattice..." << std::endl;
    lattice->initialize();

    // Set all the boundaries in the initialized lattice
    imposeOpenings(0.0);

    // iteration counter
    int stat_cycle = 0;

    // Restart from a checkpoint
    if(CHECKPOINTING) {
        pcout << endl << "*********** Restoring checkpoint ***********" << endl;

        plb_ifstream ifile(chkParamFile.c_str());
        if(ifile.is_open()) {
            ifile >> stat_cycle;
            global::mpi().bCast(&stat_cycle, 1); // Broadcast the iteration counter
        }
        else {
            pcout << "ERROR reading from the checkpoint parameter file: " << chkParamFile << "! Exiting..." << std::endl;
            return -1;
        }

        loadLattice(*lattice, chkDataFile, gpuMode);

        // Bring the boundary waveforms to the time of the checkpoint
        imposeOpenings(stat_cycle * sim.C_t);

        pcout << "Checkpoint at iteration " << stat_cycle << " (" << stat_cycle * sim.C_t << " s) loaded successfully." << std::endl;
    }

#ifdef HEMOFLOW_GPU
    if (gpuMode && !startGpu())
        return -1;
#endif
    delete labelField;   // Not needed any more
    labelField = nullptr;

    pcout << "Starting computation..." << endl;
    const int firstStep = stat_cycle;
    auto loopStart = std::chrono::steady_clock::now();

    while(stat_cycle < totalSteps)
    {
        if(stat_cycle % 200 == 0) {
            T cE = averageEnergy();
            pcout << "\rTime: " << stat_cycle*sim.C_t << "s / " << simLength << "s" << " [" << stat_cycle << " / " << totalSteps << "] " << " - Energy: " << cE << endl;

            // Capture numerical divergence if appears
            if (std::isnan(cE)){
                pcout << "ERROR: NaN average energy! Saving state and stopping simulation" << std::endl;
                syncCpuLattice();
                writeVTK(*lattice, sim, stat_cycle, outputOptions);
                return 1;
            }
            if (DEBUG) {
                syncCpuLattice();
                T sumMeasured = 0.0;
                for (auto &o : openings) {
                    T measured = o->measureFlowRate(*lattice, sim);
                    sumMeasured += measured;
                    pcout << "  label " << o->getLabel() << " (" << o->getName() << "): measured flow rate " << measured << " m^3/s";
                    if (o->isVelocityType())
                        pcout << ", imposed " << o->getImposedFlowRate(sim) << " m^3/s";
                    pcout << ", d^3 " << pow(o->getRadius() * 2 * sim.C_l, 3) << " m^3" << std::endl;
                }
                pcout << "  net inflow (mass balance): " << sumMeasured << " m^3/s" << std::endl;
            }
        }

        // Impose boundary conditions with dt progress in time
        imposeOpenings(sim.C_t);

        // Calculate next step
        collideAndStream();

        // Advance time
        stat_cycle++;

        // Save output
        if(stat_cycle >= lastCycleStartIter && stat_cycle % saveFrequency == 0) {
            pcout << "Writing output at: " << stat_cycle << " (" << stat_cycle*sim.C_t << " s)." << endl;
            syncCpuLattice();
            writeVTK(*lattice, sim, stat_cycle, outputOptions);
        }

        if(useCheckpoint && (stat_cycle % checkpointFrequency == 0)) {
            // Keep the previous checkpoint as *_old in case writing the new one fails
            if(global::mpi().isMainProcessor()) {
                if(fileExists(chkDataFile)){
                    if(fileExists(chkDataFileOld)){
                        if( ( std::remove( chkDataFileOld.c_str() ) + std::remove( chkParamFileOld.c_str() ) ) != 0 )
                            pcout << "WARNING: cannot remove old checkpoint file!" << std::endl;
                    }
                    if (std::rename(chkParamFile.c_str(), chkParamFileOld.c_str()) || std::rename(chkDataFile.c_str(), chkDataFileOld.c_str() ))
                        pcout << "WARNING: cannot rename old checkpoint file!" << std::endl;
                }
            }

            global::mpi().barrier();

            plb_ofstream ofile(chkParamFile.c_str()); ofile << stat_cycle << endl;
            syncCpuLattice();
            saveLattice(*lattice, chkDataFile, gpuMode);
        }
    } // End of main loop

    double loopSeconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - loopStart).count();
    if (stat_cycle > firstStep)
        pcout << "Time loop: " << loopSeconds << " s for " << stat_cycle - firstStep << " steps ("
              << (double)Nx * Ny * Nz * (stat_cycle - firstStep) / loopSeconds * 1e-6
              << " M updates/s of the domain box, " << (gpuMode ? "GPU" : "CPU") << ")." << endl;

    pcout << endl << "Simulation done successfully :)" << endl;

    return 0;
}

// *** Main simulation entry point
int main(int argc, char *argv[])
{
#ifdef HEMOFLOW_GPU
    // HPC-X: without its CUDA transports UCX does not stage the solver's (managed) host memory
    // through the GPU, which made CPU runs of the GPU build several times slower. The single-GPU
    // mode sends no MPI messages. (doc/dev/TODO-GPU.md, Phase 2)
    setenv("UCX_TLS", "^cuda", 0);
    // HPC-X: UCC collectives keep CUDA resources that MPI_Finalize destroys; on the ranks that
    // released their CUDA context during a GPU run that crashed. Palabos uses few collectives.
    setenv("OMPI_MCA_coll_ucc_enable", "0", 0);
#endif
    plbInit(&argc, &argv);
    printBanner();

    int exitCode = 0;
    if (!parseCommandLine(argc, argv, exitCode))
        return exitCode;

    string paramXmlFileName = argv[1];
    string outDir; // Output directory

    // This also checks if there are folders in the path or not!
    size_t folderIdx = paramXmlFileName.find_last_of("/\\");
    if (std::string::npos == folderIdx)
        workingFolder = ".";
    else
        workingFolder = paramXmlFileName.substr(0, folderIdx);

    // *** Read the configuration and the geometry, set up the lattice and the boundary conditions
    try {
        XMLreader xml(paramXmlFileName);

        if (!hasSection(xml, "simulation"))
            throw std::runtime_error("The configuration file has no <simulation> section.");
        if (!hasSection(xml, "geometry"))
            throw std::runtime_error("The configuration file has no <geometry> section.");
        XMLreaderProxy geomXml = xml["geometry"];

        readSimulationSection(xml);
        outDir = createOutputFolder();

        // Geometry file: read on the main process; the other processes get only the
        // parts of the voxel fields they own (see geometry.h)
        pcout << "Loading in data file..." << std::endl;
        string npzFileName;
        readRequired(geomXml, "geometry", "file", npzFileName);
        npzFileName = resolvePath(npzFileName);
        GeometryFile geometry(npzFileName);
        setUpFromGeometry(xml, geometry);
        vector<OpeningSpec> specs = readOpenings(xml, geometry);
#ifdef HEMOFLOW_GPU
        if (chooseGpu(geometry, specs))
            enterGpuMode();   // The other ranks wait there until the end of the run
#endif
        createLattice(geometry, specs);
        createFieldsAndOpenings(geometry, specs);
        setCyclePeriod();
    }
    catch (PlbIOException& exception) {
        pcout << "ERROR while processing input file " << paramXmlFileName << ": " << exception.what() << std::endl;
        return finishRun(-1);
    }
    catch (std::exception& exception) {
        pcout << "ERROR: " << exception.what() << std::endl;
        return finishRun(-1);
    }

    printSimulationParameters();

    return finishRun(runSimulation(outDir));
}
