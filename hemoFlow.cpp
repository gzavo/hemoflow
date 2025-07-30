#include <map>
#include <string>
#include <sstream>
//#include <algorithm>
#include <cstdlib>
//#include <iomanip>
#include <vector>
#include <cmath>

using namespace std;

#include "globals.h"
#include "helper.h"
#include "opening.h"
#include "porous.h"
#include "io.h"


/* ********** GLOBAL VARIABLES ************/

// Domain size
int Nx=0;
int Ny=0;
int Nz=0;

// Domain data
cnpy::NpyArray geometryFlag;
unsigned short* gfData = nullptr;

// Flow diverter (stent) data
cnpy::NpyArray stentFlag;
unsigned short* sfData = nullptr;
T linCoeff = 0.0;
T quadCoeff = 0.0;
T linCoeff_lb = 0.0;
T quadCoeff_lb = 0.0;

// Simulation parameters structure
SimPar sim;

// Technical simulation parameters
bool useCheckpoint = true;
bool saveInitState = true;
int blockSize;
int envelopeWidth = 1;
string outputFolder;
string workingFolder;
T simLength;
T saveFreqTime;
T checkpointFreqTime;
string mode;

// Vector of openings
vector<OpeningHandler*> openings;

// Simulation data structures
MultiBlockLattice3D<T, DESCRIPTOR> *lattice = nullptr;
MultiNTensorField3D<T> *porosityField = nullptr;

// Carreau parameters for human blood
//  B.M.  Johnston,  P.R.  Johnson,  S.  Corney,  and  D. Kilpatrick, “Non-Newtonian blood flow in human  right  coronary  arteries:  steady  state  simulations,” Journal  of  Biomechanics, 37, 709 – 720 (2004)
//  Y.I.  Cho  and  K.R.  Kensey,  “Effects  of  the  non-Newtonian viscosity of blood on flows in a   diseased   arterial   vessel.   Part   1:   steady   flows,” Biorheology28, 241 (1991)
// nu0 = 5.6e-5; nuInf = 3.5e-6;
T nu0 = 5.6e-5;     // [m^2/s]
T nuInf = 3.22e-6;   // [m^2/s]
T lambda = 3.331;
T n = 0.3568;

// Simple find index of a value in an array (where indices are unique)
int findIndex(const unsigned short *array, int arraySize, unsigned short itemToFind) {
    for(int i = 0; i < arraySize; i++)
        if(array[i] == itemToFind)
            return i;
    return -1;  // Not found
}

// *** Calculating LB parameters using Re on the inlet: Re = U_avg * D / nu
void calcSimulationParameters(SimPar &sim, T dx, T dt = -1, T U_max_LB_ = 0.1)
{   
    sim.C_l = dx;
    sim.U_max_lb = U_max_LB_;

    T nuInf_lb;
    if(dt > 0.0){
        sim.C_t = dt;
        nuInf_lb = nuInf * sim.C_t / sim.C_l / sim.C_l;
        T tau = 3.0*nuInf_lb+0.5;
        sim.omega = 1.0 / tau;
    } 
    else {
        sim.omega = 1.0;
        nuInf_lb =  0.5 / 3.0;
        sim.C_t = nuInf_lb * sim.C_l * sim.C_l / nuInf;
    }

    T nu0_lb = nu0 * sim.C_t / sim.C_l / sim.C_l;

    // Derived quantities using density
    sim.C_r = BLOOD_DENSITY;    // TODO IF we are simulating blood.... Note: only changes pressure output values, the simulation results are independent!
    sim.C_p = sim.C_r * sim.C_l * sim.C_l / (sim.C_t * sim.C_t);
    sim.C_m = sim.C_r * sim.C_l * sim.C_l * sim.C_l;

    // TODO: convert linCoeff and quadCoeff
    //!!!!NOTE it is multiplied with one lattice lenght unit
    //https://youtu.be/sOQMXxoKFQM?t=1168
    
    linCoeff_lb = linCoeff * sim.C_l*sim.C_l * sim.C_t / sim.C_m;       // [ kg / (m2 s) ] 
    quadCoeff_lb = quadCoeff * sim.C_l*sim.C_l * sim.C_l / sim.C_m;     // [ kg / m3 ]

    // TODO: add sanity check on parameters here
    if(sim.U_max_lb > 0.1)
        pcout << "** Parameter sanity WARNING ** LBM velocity seems high:" << sim.U_max_lb << endl;
    if(sim.omega > 1.9)
        pcout << "** Parameter sanity WARNING ** LBM omega seems high:" << sim.omega << endl;

    // TODO: using the Smagorinsky dynamics as a base add dynamic viscosity (Carreau and rheoModel)
    // (ps: solve the Fokker-Plank in rheoModel with finite difference?)
    global::CarreauParameters().setNu0(nu0_lb);
    global::CarreauParameters().setNuInf(nuInf_lb);
    global::CarreauParameters().setLambda(3.313);  //1.
    global::CarreauParameters().setExponent(0.357);   //0.3
}

// Progress timer for the openings and impose BC values
void imposeOpenings(T dt)
{
    T murrayExponent = 3.0;
    T murrayOutletTotalDiamt = 0.0;
    T sumInflowRate = 0.0;

    // Phase I - Defined BCs
    // Get the sum defined inflowrate and the sum undefined outflow surface

    for (auto &o : openings)
    {

        if (o->getOpeningType() == OPENING_VELOCITY)
        { // All non-Murray velocity BCs

            o->progressTime(lattice, dt);
            sumInflowRate += o->getScaledFlowRate();
        }
        else
        {
            murrayOutletTotalDiamt += pow(o->getRadius()*2, murrayExponent); // LBM units with Murray exponent
        }
    }

    // Phase II - Automatic BCs
    // Set the undefined outflow rates according to Murray's law
    // C. Chnafa, O. Brina, V. M. Pereira, and D. A. Steinman, “Better Than Nothing: A Rational Approach for Minimizing the Impact of Outflow Strategy on Cerebrovascular Simulations,” American Journal of Neuroradiology, vol. 39, no. 2, pp. 337–343, 2018, doi: 10.3174/ajnr.A5484.

    for (auto &o : openings)
    {
        if (o->getOpeningType() == OPENING_MURRAY)
        {
            T flowRate = -1 * ((pow(o->getRadius()*2, murrayExponent)* sumInflowRate) / murrayOutletTotalDiamt); //-1 cause it is an outlet
            T murrayVelocity = flowRate / (pow(o->getRadius(), 2)*3.14); // v=Q/A
            // The profile flow-rate is 0.5 (we give maximum velocity as a parameter) only if we have a parabolic profile. Let's assume it for performance reasons.
            // T profileFlowRate = o->getProfileFlowRate();     // Use this if not parabolic!
            T profileFlowRate = 0.5;
            o->setBCParameter(murrayVelocity / profileFlowRate);
            o->progressTime(lattice, dt);
        }
    }
}

// *** Main simulation entry point
int main(int argc, char *argv[])
{
    plbInit(&argc, &argv);

    pcout   << "********************************* " << endl
            << "*        hemoFlow  v0.31        * " << endl
            << "********************************* " << endl;

    // *** Reading in command line arguments
    if(global::argc() < 2) {
        pcout << "Not enough arguments; the syntax is: "
              << (std::string)global::argv(0) << " parameter-input-file.xml [-r]" << std::endl;
        return -1;
    }

    // Reading in the config file name
    string paramXmlFileName;
    try {
        global::argv(1).read(paramXmlFileName);
    }
    catch (PlbIOException& exception) {
        pcout << "Wrong input XML; the syntax is: "
              << (std::string)global::argv(0) << " parameter-input-file.xml [-r]" << std::endl;
        return -1;
    }
    
    // Check if we have the checkpoint flag
    string checkpointFlag;
    bool isCheckpointed = false;
    if(global::argc() > 2) {
        try {
            global::argv(2).read(checkpointFlag);
            if(checkpointFlag.compare("-r")==0) {
                isCheckpointed = true;
                pcout << std::endl << "Restart from checkpoint is requested! The checkpoint data will be loaded after the geometry setup." << std::endl << std::endl;
            }
            else
                pcout << "Unknown command line argument: " << checkpointFlag << std::endl;
        }
        catch (PlbIOException& exception) {
            // No flag, nothing to do
        }
    }

    string outDir; // Output directory

    XMLreader xml(paramXmlFileName);

    // This also checks if there are folders in the path or not!
    size_t folderIdx = paramXmlFileName.find_last_of("/\\");
    if (std::string::npos == folderIdx)
        workingFolder = ".";
    else
        workingFolder = paramXmlFileName.substr(0, folderIdx);
    
    // *** Load in data files
    try {
        pcout << "Loading in data file..." << std::endl;

        xml["simulation"]["outputDir"].read(outputFolder);
        
        outDir = workingFolder + "/" + outputFolder;

        // Check if output dir exists and accessible on the main processor
        if(global::mpi().isMainProcessor()) 
            if (dirExists(outDir) <= 0) {
                pcout << "Output folder " << outDir << " does not exist! Creating it...." << std::endl;
                mkpath(outDir.c_str(), 0777);
            }

        // Sync up after creating the directory by the master process.
        //global::mpi().barrier();

        pcout << "Output folder: " << outDir+"/" << std::endl;

        global::directories().setOutputDir(outDir+"/");

        // Reading main simulation parameters
        xml["simulation"]["blockSize"].read(blockSize);
        xml["simulation"]["simLength"].read(simLength);
        xml["simulation"]["saveFrequency"].read(saveFreqTime);
        
        // Check for optional checkpoint argument
        try {
            xml["simulation"]["checkpointFrequency"].read(checkpointFreqTime);
        }
        catch (PlbIOException& exception) {
            pcout << "Warning: checkpointing tag was not found in config, checkpointing will be disabled!" << std::endl;
            useCheckpoint = false;
        }

        // Loading the input file
        string npzFileName;
        xml["geometry"]["file"].read(npzFileName);
        cnpy::npz_t geom_npz = cnpy::npz_load(workingFolder + "/" + npzFileName);
        pcout << "Input data elements: " << geom_npz.size() << std::endl;
        
        for (auto const& array : geom_npz) 
            pcout << "->" << array.first << std::endl;

        // TODO - Consider hdf5-based input files (see the Turbulence branch)
        // Loading geometry
        geometryFlag = geom_npz["geometryFlag"];
        gfData = geometryFlag.data<unsigned short>();
        Nx = geometryFlag.shape[0];
        Ny = geometryFlag.shape[1];
        Nz = geometryFlag.shape[2];
        pcout << "Domain size: " << Nx << " x " << Ny << " x " << Nz << std::endl;

        // Reading dx = C_l from the geometry file
        cnpy::NpyArray dxA = geom_npz["dx"];
        double sim_dx = (dxA.data<double>())[0];

        pcout << "Resolution [m]: " << sim.C_l << std::endl;

        // Loading the geometry of the stent (if there is one)
        stentFlag = geom_npz["stent"];
        
        if(stentFlag.shape.size() > 1) {    // Check if there is data on the FD
            pcout << "Found flow diverter information to load." << std::endl;
            sfData = stentFlag.data<unsigned short>();
            // Also look for corresponding data in xml
            xml["flowdiverter"]["linCoeff"].read(linCoeff);
            xml["flowdiverter"]["quadCoeff"].read(quadCoeff);
        }

        // Check if time-step is specified
        double sim_dt = -1.0;
        try {
            xml["simulation"]["dt"].read(sim_dt);
        }
        catch (PlbIOException& exception) {}

        // Calcualte simulation parameters here!
        pcout << "Setting LBM parameters..." << std::endl;
        calcSimulationParameters(sim, sim_dx, sim_dt);  

        // **** Processing openings ****
        pcout << "Processing openings..." << std::endl;
        
        // Loading information on openings 
        cnpy::NpyArray openingIndex = geom_npz["openingIndex"];
        auto* oiData = openingIndex.data<unsigned short>();
        
        cnpy::NpyArray openingRadius = geom_npz["openingRadius"];
        auto* orData = openingRadius.data<double>();
        
        /* -- Not needed atm.
        cnpy::NpyArray openingQRatio = geom_npz["openingNormalizedQRatio"];
        double* oqData = openingQRatio.data<double>();
        
        cnpy::NpyArray openingCenter = geom_npz["openingCenter"];
        double* ocData = openingCenter.data<double>();
        */

        cnpy::NpyArray openingNormal = geom_npz["openingNormal"];
        auto* onData = openingNormal.data<double>();

        if(SPARSE) {
            pcout << "Setting simulation domain mask for sparse decomposition..." << endl;
            auto *flagMatrix = new MultiScalarField3D<int>(Nx,Ny,Nz);
            setToFunction(*flagMatrix, flagMatrix->getBoundingBox(), FlagMaskDomain3D<unsigned short>(gfData, 1));

            pcout << "Creating sparse representation ..." << endl;
            
            //Create sparse representation
            MultiBlockManagement3D sparseBlockManagement =
                        computeSparseManagement(*plb::reparallelize(*flagMatrix, blockSize, blockSize, blockSize), envelopeWidth);
                                                            
            #if LES
                lattice = new MultiBlockLattice3D<T, DESCRIPTOR> (sparseBlockManagement,
                                                            defaultMultiBlockPolicy3D().getBlockCommunicator(),
                                                            defaultMultiBlockPolicy3D().getCombinedStatistics(),
                                                            defaultMultiBlockPolicy3D().getMultiCellAccess<T,DESCRIPTOR>(),
                                                            new SmagorinskyRegularizedDynamics<T,DESCRIPTOR>(sim.omega, cSmago)); 
            #else
                lattice = new MultiBlockLattice3D<T, DESCRIPTOR> (sparseBlockManagement,
                                                            defaultMultiBlockPolicy3D().getBlockCommunicator(),
                                                            defaultMultiBlockPolicy3D().getCombinedStatistics(),
                                                            defaultMultiBlockPolicy3D().getMultiCellAccess<T,DESCRIPTOR>(),
                                                            new BackgroundDynamics(sim.omega));
            #endif
        }
        else {
            #if LES
                lattice = new MultiBlockLattice3D<T, DESCRIPTOR>(Nx, Ny, Nz, new SmagorinskyRegularizedDynamics<T,DESCRIPTOR>(sim.omega, cSmago));
            #else
                lattice = new MultiBlockLattice3D<T, DESCRIPTOR>(Nx, Ny, Nz, new BackgroundDynamics(sim.omega));
            #endif
        }

        #if LES
            instantiateStaticSmagorinsky(*lattice, lattice->getBoundingBox(), cSmago);
        #endif

        pcout << getMultiBlockInfo(*lattice) << endl;

        // Loop through the openings in the datafile 
        unsigned int numOpenings = openingRadius.shape[0];
        pcout << "Number of openings in geometry: " << numOpenings << std::endl;
        pcout << "Opening flags: ";

        for(unsigned int o=0; o < numOpenings; o++)
            pcout << oiData[o] << " ";
        pcout << endl;

        try
        {
            xml["simulation"]["mode"].read(mode);
        }
        catch (PlbIOException &exception)
        {
            pcout << "No mode is defined, all BCs must be declared in XML." << std::endl;
            mode = "";
        }
        // Loop through the openings following the xml config
        if (mode == "aneurysm")
        {
            pcout << "Aneurysm boundary condition mode activated";
            for (unsigned int o = 0; o < numOpenings; o++)
            {
                if (o == 0)
                {
                    int label = oiData[o];
                    int type_velocity = 1;
                    pcout << "Processing opening: " << o << std::endl;

                    string xmlTagOpening = "opening_" + std::to_string(o);

                    // Get the direction of the opening
                    int s = openingNormal.shape[1];
                    vec3d dir(onData[gT2D(s, o, 0)], onData[gT2D(s, o, 1)], onData[gT2D(s, o, 2)]);

                    // Create the opening
                    auto *opening = new OpeningHandler(gfData, static_cast<GeometryLabel>(label), static_cast<OpeningType>(type_velocity), orData[o] / sim.C_l, dir);

                    opening->setName("Velocity inlet");
                    opening->setBCType(lattice);

                    string parameterStr;
                    double parameter;
                    xml["geometry"][xmlTagOpening]["parameter"].read(parameterStr);
                    if (!parameterStr.empty())
                        parameter = std::stod(parameterStr);
                    
                    //*2 because we are putting in avg velo not parabolic
                    opening->setBCParameter(parameter*2, sim);

                    // Load scale function (fileName from XML)
                    string flowrateFunc;
                    xml["geometry"][xmlTagOpening]["timeScaleFunction"].read(flowrateFunc);
                    if (!flowrateFunc.empty())
                        opening->loadScaleFunction(workingFolder + "/" + flowrateFunc);

                    // Set profile
                    opening->createPoiseauilleProfile();
                    opening->printOpeningDetails(sim);
                    openings.push_back(opening);
                }
                //TODO: be consistent in opening labeling (area descending?->smallest is pressure)
                else if (o == numOpenings - 1)
                {
                    int label = oiData[o];
                    int type_pressure = 3;
                    int s = openingNormal.shape[1];
                    vec3d dir(onData[gT2D(s, o, 0)], onData[gT2D(s, o, 1)], onData[gT2D(s, o, 2)]);
                    // Create the opening
                    auto *opening = new OpeningHandler(gfData, static_cast<GeometryLabel>(label), static_cast<OpeningType>(type_pressure), orData[o] / sim.C_l, dir);

                    opening->setName("Pressure outlet");
                    opening->setBCType(lattice);
                    opening->setBCParameter(0, sim);
                    opening->createConstantPressureProfile();
                    opening->printOpeningDetails(sim);
                    openings.push_back(opening);
                }
                else
                {
                    int label = oiData[o];
                    int type_murray = 2;
                    string xmlTagOpening = "opening_" + std::to_string(o);

                    pcout << "Processing opening: " << o << std::endl;
                    // Get the direction of the opening
                    int s = openingNormal.shape[1];
                    vec3d dir(onData[gT2D(s, o, 0)], onData[gT2D(s, o, 1)], onData[gT2D(s, o, 2)]);

                    // Create the opening
                    auto *opening = new OpeningHandler(gfData, static_cast<GeometryLabel>(label), static_cast<OpeningType>(type_murray), orData[o] / sim.C_l, dir);

                    opening->setName("Murray outlet");
                    opening->setBCType(lattice);
                    opening->setBCParameter(0, sim);
                    opening->createPoiseauilleProfile();
                    opening->printOpeningDetails(sim);
                    openings.push_back(opening);
                }
            }
        }
        else
        {
            for (unsigned int o = 0; o < numOpenings; o++)
            {
                pcout << "Processing opening: " << o << std::endl;

                string xmlTagOpening = "opening_" + std::to_string(o);

                string name;
                xml["geometry"][xmlTagOpening]["name"].read(name);
                int type;
                xml["geometry"][xmlTagOpening]["type"].read(type);
                int label;
                xml["geometry"][xmlTagOpening]["label"].read(label);

                int openingIdx = findIndex(oiData, numOpenings, label);
                if (openingIdx == -1)
                    pcout << "ERROR: Opening label " << label << " was found in the config xml, but not in the geometry file!" << endl;

                // Get the direction of the opening
                int s = openingNormal.shape[1];
                vec3d dir(onData[gT2D(s, openingIdx, 0)], onData[gT2D(s, openingIdx, 1)], onData[gT2D(s, openingIdx, 2)]);

                // Create the opening
                auto *opening = new OpeningHandler(gfData, static_cast<GeometryLabel>(label), static_cast<OpeningType>(type), orData[openingIdx] / sim.C_l, dir);

                opening->setName(name);
                opening->setBCType(lattice);

                string parameterStr;
                double parameter;
                xml["geometry"][xmlTagOpening]["parameter"].read(parameterStr);
                if (!parameterStr.empty())
                    parameter = std::stod(parameterStr);

                opening->setBCParameter(parameter, sim);

                // Load scale function (fileName from XML)
                string flowrateFunc;
                xml["geometry"][xmlTagOpening]["timeScaleFunction"].read(flowrateFunc);
                if (!flowrateFunc.empty())
                    opening->loadScaleFunction(workingFolder + "/" + flowrateFunc);

                // Set profile
                if (type == OPENING_VELOCITY || type == OPENING_MURRAY || type == OUTLET_FREEFLOW)
                {
                    opening->setBCParameter(parameter*2, sim); //*2 because we are putting in avg velo not parabolic
                    opening->createPoiseauilleProfile(); // Normalized to max_vel = 1.0 (LBM units)
                    // opening->normalizeFlowRate();           // Normalize to Q=1 (LBM units)
                }
                else if (type == OPENING_PRESSURE)
                {
                    opening->createConstantPressureProfile(); // Contant pressure = 1.0 (LBM density)
                }

                opening->printOpeningDetails(sim);

                openings.push_back(opening);
            }
    }

        // Sanity check
        if(numOpenings != openings.size()) {
            pcout << "**WARNING** The number of opening definitions don't match between the config and the geometry file! So how many openings do we actually have?" << endl; 
        }

    }
    catch (PlbIOException& exception) {
        pcout << "Error while processing input file " << paramXmlFileName
              << ": " << exception.what() << std::endl;
        return -1;
    }
    
    pcout   << "*********** Simulation parameters *********** " << endl
            << "size [LU]:   " << Nx << "x" << Ny << "x" << Nz << endl
            << "dx [m]:  " << sim.C_l << endl
            << "dt [s]:  " << sim.C_t << endl
            << "dm [kg]:  " << sim.C_m << endl
            << "omega:  " << sim.omega << endl
            << "nu:     " << 1./3. * (1./sim.omega - 0.5) << endl
            << "U_max [lbm]: " << sim.U_max_lb << endl 
            << "U_max [m/s]: " << sim.U_max_lb * sim.C_l / sim.C_t << endl << endl;

    int saveFrequency;
    saveFrequency = (int)round(saveFreqTime/sim.C_t);
    pcout << "Saving frequency set to every " << saveFreqTime << " s (" << saveFrequency << " steps)." << endl;

    int checkpointFrequency = (int)round(checkpointFreqTime/sim.C_t);
    
    if(useCheckpoint)     
        pcout << "Chekpointing will happen every " << checkpointFreqTime << " s (" << checkpointFrequency << " steps)." << endl;
    
    // Checkpoint file names relative to the output folder
    string chkParamFile = outDir+"/checkpoint_parameters.dat";
    string chkDataFile = outDir+"/checkpoint_lattice.dat";
    string chkParamFileOld = outDir+"/checkpoint_parameters_old.dat";
    string chkDataFileOld = outDir+"/checkpoint_lattice_old.dat";

    // If there is data on porosity, set up porous layer in the simulation   
    if(sfData != nullptr) {
        pcout << "Setting up porous layer for flow diverter..." << std::endl;
        porosityField = defaultGenerateMultiNTensorField3D<T>(lattice->getMultiBlockManagement(), 1).release();
        applyProcessingFunctional(new InitializePorousField<T, unsigned short>(sfData), porosityField->getBoundingBox(), *porosityField);
        integrateProcessingFunctional( new PorousForceFunctional<T, DESCRIPTOR>(linCoeff_lb, quadCoeff_lb), lattice->getBoundingBox(), *lattice, *porosityField);
    }

    pcout << "Defining walls..." << std::endl;
    defineDynamics(*lattice, lattice->getBoundingBox(), new FlagMaskSingleDomain3D<unsigned short>(gfData, 0), new NoDynamics<T, DESCRIPTOR>);
    defineDynamics(*lattice, lattice->getBoundingBox(), new FlagMaskSingleDomain3D<unsigned short>(gfData, 1), new BounceBack<T, DESCRIPTOR>(1.0));

    // TODO: add some reparallelize here, check if it plays nice with checkpointing
    
    pcout << "Initializing lattice in equilibrium..." << std::endl;
    initializeAtEquilibrium (*lattice, lattice->getBoundingBox(), 1.0, Array<T,3>((T)0.,(T)0.,(T)0.) );

    pcout << "Finalizing lattice..." << std::endl;
    lattice->initialize();

    // Set all the boundaries in the initialized lattice
    imposeOpenings(0.0);

    // iteration counter
    int stat_cycle = 0;
    
    // Check if the simulation was checkpointed
    if(isCheckpointed) {
        pcout << endl << "*********** Restoring checkpoint ***********" << endl;
        
        // Load in the iteration counter
        string chkParamFile = outDir+"/checkpoint_parameters.dat";
        plb_ifstream ifile(chkParamFile.c_str());
        if(ifile.is_open()) {
            ifile >> stat_cycle;
            global::mpi().bCast(&stat_cycle, 1); // Broadcast the iteration counter
        }
        else {
            pcout << "ERROR reading from the checkpoint parameter file: checkpoint_parameters.dat! Exiting..." << std::endl;
            return -1;
        }
        
        // Load the lattice
        loadBinaryBlock(*lattice, outDir+"/checkpoint_lattice.dat");
        
        pcout << "Checkpoint at iteration " << stat_cycle << " loaded succesfully." << std::endl;
    }
    else { // If not, then let's chek the initial state and do a warm up.
        pcout << endl << "*********** Entering stationary warmup phase ***********" << endl;
           
        int convergenceSteps = 10*max(max(Nx, Ny), Nz);
        int rampupInterval = convergenceSteps/2;
        
        T minDE = 1e-12; T dE = 100; T prevE = 0;
    

        T cE = computeAverageEnergy(*lattice);
        if(isnan(cE)) {
            pcout << "WARNING: Energy (velocity) is NaN! Please check the simulation setup! Exiting..." << endl;
            return -1;
        }

        if(saveInitState) {
            pcout << "Energy at the initial state: "<< cE << endl;
            pcout << "Saving initial state with flow diverter..." << endl;
            writeVTK(*lattice, sim, -1, porosityField);
            //writeHDF5(*lattice, sim, -1, outDir, porosityField);
        }

        while(abs(dE) > minDE && stat_cycle < convergenceSteps )
        {
            for (auto &o : openings)
            {
                if (o->getOpeningType() == OPENING_VELOCITY || o->getOpeningType() == OPENING_MURRAY)
                {

                    o->progressWarmup(lattice, stat_cycle, rampupInterval);
                }
            }

            lattice->collideAndStream();

            T cE = computeAverageEnergy(*lattice);
            dE = cE - prevE; prevE = cE;
    
            if(stat_cycle % 500 == 0) {
                pcout << "Delta energy: " << abs(dE) << "/" << minDE << "  Cycle: [" << stat_cycle << "/" << convergenceSteps <<"]" << std::endl;
                writeVTK(*lattice, sim,stat_cycle);
            }
            
            stat_cycle++;
        }
        pcout << "Delta energy: " << abs(dE) << "/" << minDE << "  Cycle: [" << stat_cycle << "/" << convergenceSteps <<"]" << std::endl;
    
        pcout << endl << "*********** Entering transient simulation phase ***********" << endl;
    
        pcout << "Saving time step 0..." << endl;
        writeVTK(*lattice, sim, 0);
        // writeNPZ(*lattice, 0);
        //writeHDF5(*lattice, sim, 0, outDir);
        
        // Set the counter back
        stat_cycle = 0;
    }
    
    pcout << "Starting computation..." << endl;

    while(stat_cycle*sim.C_t <= simLength + sim.C_t)
    {
        
        if(stat_cycle % 200 == 0) {
            T cE = computeAverageEnergy(*lattice);
            pcout << "\rTime: " << stat_cycle*sim.C_t << "s / " << simLength << "s" << " [" << stat_cycle << " / " << std::round(simLength/sim.C_t) << "] " << " - Energy: " << cE << endl;
            
            // Capture numerical divergence if appears
            if (std::isnan(cE)){
                pcout << "ERROR: NaN average energy! Saving state and stopping simulation" << std::endl;
                writeVTK(*lattice, sim, stat_cycle);
                //writeHDF5(*lattice, sim, stat_cycle, outDir, porosityField);
                return 0;
            }
            for (auto &o : openings)
            {
                pcout << o->getName() << " flow rate SI: " << o->getFlowRate(sim) << " scaledVFR:" << o->getScaledFlowRate() << " with D^3:" << pow(o->getRadius() * 2, 3) << std::endl;
            }
        }

        // Impose boundary conditions with dt progress in time
        imposeOpenings(sim.C_t);

        // Calculate next step
        lattice->collideAndStream();

        // Advance time
        stat_cycle++;

        // Save output 
        if(stat_cycle % saveFrequency == 0) {
            pcout << "Writing output at: " << stat_cycle << " (" << stat_cycle*sim.C_t << " s)." << endl;
            writeVTK(*lattice, sim, stat_cycle);
            // writeNPZ(*lattice, sim, stat_cycle);
            //writeHDF5(*lattice, sim, stat_cycle, outDir, porosityField);   
        }
        
        if(useCheckpoint && (stat_cycle % checkpointFrequency == 0)) {
            // Overwriting previous checkpoint. Note: if failure happens during saving the checkpoint we cannot recover: TODO two step checkpoint
            if(global::mpi().isMainProcessor()) {
                if(fileExists(chkDataFile)){
                    // Remove prev-previous checkpoint            
                    if(fileExists(chkDataFileOld)){
                        if( ( std::remove( chkDataFileOld.c_str() ) + std::remove( chkParamFileOld.c_str() ) ) != 0 )
                            pcout << "WARNING: cannot remove old chekpoint file!" << std::endl;
                    }
                    // Rename previous checkpoint
                    if (std::rename(chkParamFile.c_str(), chkParamFileOld.c_str()) || std::rename(chkDataFile.c_str(), chkDataFileOld.c_str() )) 
                        pcout << "WARNING: cannot rename old chekpoint file!" << std::endl;
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
