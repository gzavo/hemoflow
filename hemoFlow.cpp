#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace std;

#include "globals.h"
#include "helper.h"
#include "opening.h"
#include "porous.h"
#include "io.h"
#include "geometry.h"


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

// Fluid properties (Newtonian)
T fluidKinematicViscosity = DEFAULT_BLOOD_KINEMATIC_VISCOSITY;  // [m^2/s]
T fluidDensity = DEFAULT_BLOOD_DENSITY;                          // [kg/m^3]

// Technical simulation parameters
bool useCheckpoint = false;
bool useSparse = true;    // Sparse block decomposition: only blocks that contain the vessel are allocated
int blockSize = 24;       // Side length of the blocks of the sparse decomposition [LB nodes]
int envelopeWidth = 1;    // Set before the lattice is created (2 with pressure openings)
string outputFolder;
string workingFolder;
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
MultiScalarField3D<T> *porosityField = nullptr;   // 1 in porous (flow diverter) voxels (distributed)


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

    // Note: the fluid is Newtonian. Non-Newtonian (e.g. Carreau) rheology is not implemented.
    //  B.M. Johnston et al., J. Biomech. 37, 709-720 (2004); Y.I. Cho and K.R. Kensey, Biorheology 28, 241 (1991)
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
            o->progressTime(*lattice, dt);
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
            o->progressTime(*lattice, dt);
        }
        else if (o->getOpeningType() != OPENING_VELOCITY) {
            o->progressTime(*lattice, dt);   // Pressure (and free outflow) openings
        }
    }
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

// *** Main simulation entry point
int main(int argc, char *argv[])
{
    plbInit(&argc, &argv);

    pcout   << "********************************* " << endl
            << "*        hemoFlow  v0.31        * " << endl
            << "********************************* " << endl;

    // *** Reading in command line arguments
    if (argc < 2)
    {
        pcout << "Not enough arguments" << endl;
        print_help();
        return -1;
    }

    for (int i = 1; i < argc; ++i)
    {
        std::string arg = argv[i];

        if (arg == "-h" || arg == "--help")
        {
            print_help();
            return 0;
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
            return 1;
        }
    }

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
        XMLreaderProxy simXml = xml["simulation"];
        XMLreaderProxy geomXml = xml["geometry"];

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

        if (numCycles < 1)
            throw std::runtime_error("<simulation><cycle> must be at least 1.");
        if (saveFreqTime <= 0.0)
            throw std::runtime_error("<simulation><saveFrequency> must be positive.");
        if (blockSize < 4)
            throw std::runtime_error("<simulation><blockSize> must be at least 4.");

        // Output fields
        string fieldList;
        if (readOptional(simXml, "simulation", "outputFields", fieldList)) {
            outputOptions.velocity = outputOptions.pressure = outputOptions.stress = false;
            std::istringstream fields(fieldList);
            string field;
            while (fields >> field) {
                if (field == "velocity") outputOptions.velocity = true;
                else if (field == "pressure") outputOptions.pressure = true;
                else if (field == "stress") outputOptions.stress = true;
                else throw std::runtime_error("Unknown output field '" + field + "' in <simulation><outputFields> (valid: velocity pressure stress).");
            }
        }

        // Fluid properties (optional)
        if (hasSection(xml, "fluid")) {
            XMLreaderProxy fluidXml = xml["fluid"];
            readOptional(fluidXml, "fluid", "kinematicViscosity", fluidKinematicViscosity);
            readOptional(fluidXml, "fluid", "density", fluidDensity);
        }
        if (fluidKinematicViscosity <= 0.0 || fluidDensity <= 0.0)
            throw std::runtime_error("The fluid viscosity and density must be positive.");

        // Output directory, created by the main process
        outDir = resolvePath(outputFolder);
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

        // Geometry file: read on the main process; the other processes get only the
        // parts of the voxel fields they own (see geometry.h)
        pcout << "Loading in data file..." << std::endl;
        string npzFileName;
        readRequired(geomXml, "geometry", "file", npzFileName);
        npzFileName = resolvePath(npzFileName);
        GeometryFile geometry(npzFileName);
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

        vector<OpeningSpec> specs = buildOpeningSpecs(xml, openingLabels);

        // **** Lattice ****
        // The finite-difference ("interp") pressure boundary reads data two cells deep across
        // block boundaries. With a one-cell envelope the result would depend on how the domain
        // is split into blocks wherever a block boundary cuts a pressure opening.
        bool hasPressureOpening = std::any_of(specs.begin(), specs.end(),
                                              [](OpeningSpec const& spec) { return spec.type == OPENING_PRESSURE; });
        envelopeWidth = hasPressureOpening ? 2 : 1;
        pcout << "Envelope width: " << envelopeWidth << endl;

        Dynamics<T, DESCRIPTOR>* bulkDynamics =
        #if LES
            new SmagorinskyRegularizedDynamics<T, DESCRIPTOR>(sim.omega, cSmago);
        #else
            new BackgroundDynamics(sim.omega);
        #endif

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

        // Voxel labels and porosity as distributed fields with the block structure of the lattice
        labelField = generateMultiScalarField<int>(*lattice, envelopeWidth).release();
        geometry.fillLabels(*labelField);
        if (geometry.hasStent) {
            porosityField = generateMultiScalarField<T>(*lattice, envelopeWidth).release();
            geometry.fillPorosity(*porosityField);
        }
        geometry.release();

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
            openings.push_back(opening);
        }

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
    catch (PlbIOException& exception) {
        pcout << "ERROR while processing input file " << paramXmlFileName << ": " << exception.what() << std::endl;
        return -1;
    }
    catch (std::exception& exception) {
        pcout << "ERROR: " << exception.what() << std::endl;
        return -1;
    }

    pcout   << "*********** Simulation parameters *********** " << endl
            << "size [LU]:   " << Nx << "x" << Ny << "x" << Nz << endl
            << "dx [m]:  " << sim.C_l << endl
            << "dt [s]:  " << sim.C_t << endl
            << "dm [kg]:  " << sim.C_m << endl
            << "nu [m^2/s]: " << fluidKinematicViscosity << endl
            << "rho [kg/m^3]: " << fluidDensity << endl
            << "omega:  " << sim.omega << endl
            << "nu [lb]: " << 1./3. * (1./sim.omega - 0.5) << endl
            << "U_max [lbm]: " << sim.U_max_lb << endl
            << "U_max [m/s]: " << sim.U_max_lb * sim.C_l / sim.C_t << endl << endl;

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
          << (outputOptions.stress ? " stress" : "") << endl;

    // Checkpoint file names relative to the output folder
    string chkParamFile = outDir+"/checkpoint_parameters.dat";
    string chkDataFile = outDir+"/checkpoint_lattice.dat";
    string chkParamFileOld = outDir+"/checkpoint_parameters_old.dat";
    string chkDataFileOld = outDir+"/checkpoint_lattice_old.dat";

    // If there is a flow diverter, set up the porous forcing
    if(porosityField != nullptr) {
        pcout << "Setting up porous layer for flow diverter..." << std::endl;
        integrateProcessingFunctional( new PorousForceFunctional<T, DESCRIPTOR>(linCoeff_lb, quadCoeff_lb), lattice->getBoundingBox(), *lattice, *porosityField);
    }

    // Walls and outside voxels. This also overrides the boundary dynamics that the openings
    // instantiated on wall/outside voxels of their bounding rectangles.
    pcout << "Defining walls..." << std::endl;
    defineDynamics(*lattice, *labelField, lattice->getBoundingBox(), new NoDynamics<T, DESCRIPTOR>(1.0), (int)UNUSED);  // Reference density, so outside cells report 0 Pa
    defineDynamics(*lattice, *labelField, lattice->getBoundingBox(), new BounceBack<T, DESCRIPTOR>(1.0), (int)WALL);
    delete labelField;   // Not needed any more
    labelField = nullptr;

    pcout << "Initializing lattice in equilibrium..." << std::endl;
    initializeAtEquilibrium (*lattice, lattice->getBoundingBox(), 1.0, Array<T,3>((T)0.,(T)0.,(T)0.) );

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

        loadBinaryBlock(*lattice, chkDataFile);

        // Bring the boundary waveforms to the time of the checkpoint
        imposeOpenings(stat_cycle * sim.C_t);

        pcout << "Checkpoint at iteration " << stat_cycle << " (" << stat_cycle * sim.C_t << " s) loaded successfully." << std::endl;
    }

    pcout << "Starting computation..." << endl;

    while(stat_cycle < totalSteps)
    {
        if(stat_cycle % 200 == 0) {
            T cE = computeAverageEnergy(*lattice);
            pcout << "\rTime: " << stat_cycle*sim.C_t << "s / " << simLength << "s" << " [" << stat_cycle << " / " << totalSteps << "] " << " - Energy: " << cE << endl;

            // Capture numerical divergence if appears
            if (std::isnan(cE)){
                pcout << "ERROR: NaN average energy! Saving state and stopping simulation" << std::endl;
                writeVTK(*lattice, sim, stat_cycle, outputOptions);
                return 1;
            }
            if (DEBUG) {
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
        lattice->collideAndStream();

        // Advance time
        stat_cycle++;

        // Save output
        if(stat_cycle >= lastCycleStartIter && stat_cycle % saveFrequency == 0) {
            pcout << "Writing output at: " << stat_cycle << " (" << stat_cycle*sim.C_t << " s)." << endl;
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
            saveBinaryBlock(*lattice, chkDataFile);
        }
    } // End of main loop

    pcout << endl << "Simulation done successfully :)" << endl;

    return 0;
}
