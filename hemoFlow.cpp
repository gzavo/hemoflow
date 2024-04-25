#include <map>
#include <string>
#include <sstream>
#include <algorithm>
#include <cstdlib>
#include <iomanip>
#include <vector>
#include <math.h>

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
unsigned short* gfData = NULL;

// Flow diverter (stent) data
cnpy::NpyArray stentFlag;
unsigned short* sfData = NULL;
T linCoeff = 0.0;
T quadCoeff = 0.0;
T linCoeff_lb = 0.0;
T quadCoeff_lb = 0.0;

// Info on openings
cnpy::NpyArray openingIndex;
unsigned short* oiData = NULL;
cnpy::NpyArray openingRadius;
double* orData = NULL;
cnpy::NpyArray openingQRatio;
double* oqData = NULL;
cnpy::NpyArray openingCenter;
double* ocData = NULL;
cnpy::NpyArray openingTangent;
double* otData = NULL;

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

// Vector of openings
vector<OpeningHandler*> openings;

// Simulation data structures
MultiBlockLattice3D<T, DESCRIPTOR> *lattice = NULL;
MultiNTensorField3D<T> *porosityField = NULL;

// Carreau parameters
//  B.M.  Johnston,  P.R.  Johnson,  S.  Corney,  and  D. Kilpatrick, “Non-Newtonian blood flow in human  right  coronary  arteries:  steady  state  simulations,” Journal  of  Biomechanics, 37, 709 – 720 (2004)
//  Y.I.  Cho  and  K.R.  Kensey,  “Effects  of  the  non-Newtonian viscosity of blood on flows in a   diseased   arterial   vessel.   Part   1:   steady   flows,” Biorheology28, 241 (1991)
// nu0 = 5.6e-5; nuInf = 3.5e-6;
T nu0 = 5.6e-5;     // [m^2/s]
T nuInf = 3.22e-6;   // [m^2/s]
T lambda = 3.331;
T n = 0.3568;

// *** Calculating LB parameters using Re on the inlet: Re = U_avg * D / nu
void calcSimulationParameters(T D_m, SimPar &sim)
{   
    T D_lb = D_m  / sim.C_l;
    T U_avg = sim.Re * nuInf / D_m;
    
    sim.C_t =  U_AVG_LB / U_avg * sim.C_l;

    T nuInf_lb = U_AVG_LB * D_lb / sim.Re;
    T nu_ratio = nuInf_lb / nuInf;
    T nu0_lb = nu0 * nu_ratio;

    T tau = 3.0*nuInf_lb+0.5;

    sim.omega = 1.0 / tau;

    sim.C_r = BLOOD_DENSITY;    // TODO IF we are simulating blood.... Note: only changes pressure output values, the simulation results are independent!

    sim.C_p = sim.C_r * sim.C_l * sim.C_l / (sim.C_t * sim.C_t);
    sim.C_m = sim.C_r * sim.C_l * sim.C_l * sim.C_l;

    // TODO: convert linCoeff and quadCoeff
    linCoeff_lb = linCoeff * sim.C_l*sim.C_l * sim.C_t / sim.C_m;       // [ kg / (m2 s) ]
    quadCoeff_lb = quadCoeff * sim.C_l*sim.C_l * sim.C_l / sim.C_m;     // [ kg / m3 ]

    // TODO: add sanity check on parameters here

    // TODO: using the Smagorinsky dynamics as a base add dynamic viscosity (Carreau and rheoModel)
    // (ps: solve the Fokker-Plank in rheoModel with finite difference?)
    global::CarreauParameters().setNu0(nu0_lb);
    global::CarreauParameters().setNuInf(nuInf_lb);
    global::CarreauParameters().setLambda(3.313);  //1.
    global::CarreauParameters().setExponent(0.357);   //0.3
}

void processOpenings(string inletFlowrateFunc, T inletA)
{
    int numOpenings = openingRadius.shape[0];
    pcout << "-> Number of openings to process: " << numOpenings << std::endl;

    T Qin = U_AVG_LB * inletA;

    for(int i=0; i<numOpenings; i++){
        int flag = oiData[i];
        pcout << "Processing flag: " << flag << std::endl;

        if(flag < INLET){
            pcout << "WARNING! Wrong flag for an opening: " << flag << std::endl;
            continue;
        }

        int s = openingTangent.shape[1];
        vec3d dir(otData[gT2D(s,i,0)], otData[gT2D(s,i,1)], otData[gT2D(s,i,2)]);

        OpeningHandler *opening = new OpeningHandler(gfData, flag, orData[i] / sim.C_l, dir);
        
        if(flag == FIRST_OUTLET)
            opening->createConstantPressureProfile();
        else {
            if(flag == INLET)
                opening->createPoiseauilleProfile(U_AVG_LB);
                // opening->createBluntVelocityProfile(U_AVG_LB);
            else {
                // Calculate outflow velocity based on Murray's law
                T u_out = Qin * oqData[i] / (pow(orData[i] / sim.C_l, 2) * M_PI);
                opening->createPoiseauilleProfile(u_out);
                // opening->createBluntVelocityProfile(u_out);
            }
            
            if(!inletFlowrateFunc.empty())
                opening->loadScaleFunction(inletFlowrateFunc);
        }
        
        opening->printOpeningDetails();
        openings.push_back(opening);
    }

}


// *** Main simulation entry point
int main(int argc, char *argv[])
{
    plbInit(&argc, &argv);

    pcout   << "********************************* " << endl
            << "*       hemoFlowCFD  v0.3       * " << endl
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

        // Check if output dir exists and accessible
        if (dirExists(outDir) <= 0) {
            pcout << "Output folder " << outDir << " does not exist! Creating it...." << std::endl;

            if (global::mpi().isMainProcessor())
                mkpath(outDir.c_str(), 0777);
        }

        // Sync up after creating the directory by the master process.
        //global::mpi().barrier();

        pcout << "Output folder: " << outDir+"/" << std::endl;

        global::directories().setOutputDir(outDir+"/");

        xml["simulation"]["Re"].read(sim.Re);
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

        // Loading geometry
        geometryFlag = geom_npz["geometryFlag"];
        gfData = geometryFlag.data<unsigned short>();
        Nx = geometryFlag.shape[0];
        Ny = geometryFlag.shape[1];
        Nz = geometryFlag.shape[2];
        pcout << "Domain size: " << Nx << " x " << Ny << " x " << Nz << std::endl;

        // Reading dx = C_l
        cnpy::NpyArray dxA = geom_npz["dx"];
        sim.C_l = (dxA.data<double>())[0];

        pcout << "Resolution [m]: " << sim.C_l << std::endl;

        // Loading stent geometry
        stentFlag = geom_npz["stent"];
        
        if(stentFlag.shape.size() > 1) {    // Check if there is data on FD
            pcout << "Found flow diverter information to load." << std::endl;
            sfData = stentFlag.data<unsigned short>();
            // Also look for corresponding data in xml
            xml["flowdiverter"]["linCoeff"].read(linCoeff);
            xml["flowdiverter"]["quadCoeff"].read(quadCoeff);
        }

        // Loading information on openings
        openingIndex = geom_npz["openingIndex"];
        oiData = openingIndex.data<unsigned short>();
        
        openingRadius = geom_npz["openingRadius"];
        orData = openingRadius.data<double>();
        
        openingQRatio = geom_npz["openingNormalizedQRatio"];
        oqData = openingQRatio.data<double>();
        
        openingCenter = geom_npz["openingCenter"];
        ocData = openingCenter.data<double>();
        
        openingTangent = geom_npz["openingTangent"];
        otData = openingTangent.data<double>();


        pcout << "Number of openings: " << openingRadius.shape[0] << std::endl;

        pcout << "Processing openings..." << std::endl;

        T inletD = 2.0 * orData[0]; // [m]
        pcout << "-> Inlet radius [m]: " << orData[0] << std::endl;

        // Inlet area in lattice units
        T inletA = pow(orData[0] / sim.C_l, 2) * M_PI;


        string inletFlowrateFunc;
        xml["geometry"]["inletFlowrateFunc"].read(inletFlowrateFunc);
        processOpenings(workingFolder + "/" + inletFlowrateFunc, inletA);

        pcout << "Setting LBM parameters..." << std::endl;
        calcSimulationParameters(inletD, sim);
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
            << "omega:  " << sim.omega << endl
            << "nu:     " << 1./3. * (1./sim.omega - 0.5) << endl
            << "Re_inlet: " << sim.Re << endl
            << "U_avg(inlet) [lbm]: " << U_AVG_LB << endl 
            << "U_avg [m/s]: " << U_AVG_LB * sim.C_l / sim.C_t << endl << endl;

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

    // TODO: IMPORTANT! - Work out proper sparse mode, we waste up to 90% numerical cells. Take a hint from HemoCell.
    bool sparse = false;
    if(sparse) {
        pcout << "Setting simulation domain mask for sparse decomposition..." << endl;
        MultiScalarField3D<int> *flagMatrix = new MultiScalarField3D<int>(Nx,Ny,Nz);
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

    // If there is data on porosity, set up porous layer in the simulation   
    if(sfData != NULL) {
        pcout << "Setting up porous layer for flow diverter..." << std::endl;
        porosityField = defaultGenerateMultiNTensorField3D<T>(lattice->getMultiBlockManagement(), 1).release();
        applyProcessingFunctional(new InitializePorousField<T, unsigned short>(sfData), porosityField->getBoundingBox(), *porosityField);
        integrateProcessingFunctional( new PorousForceFunctional<T, DESCRIPTOR>(linCoeff_lb, quadCoeff_lb), lattice->getBoundingBox(), *lattice, *porosityField);
    }

    pcout << "Defining walls..." << std::endl;
    defineDynamics(*lattice, lattice->getBoundingBox(), new FlagMaskSingleDomain3D<unsigned short>(gfData, 0), new NoDynamics<T, DESCRIPTOR>);
    defineDynamics(*lattice, lattice->getBoundingBox(), new FlagMaskSingleDomain3D<unsigned short>(gfData, 1), new BounceBack<T, DESCRIPTOR>(1.0));

    // TODO: add some reparallelize here, check if it plays nice with checkpointing

    pcout << "Setting values on openings..." << std::endl;
    for(auto &o: openings){
        o->setBC(lattice);
    }

    pcout << "Initializing lattice in equilibrium..." << std::endl;
    initializeAtEquilibrium (*lattice, lattice->getBoundingBox(), 1.0, Array<T,3>((T)0.,(T)0.,(T)0.) );

    pcout << "Finalizing lattice..." << std::endl;
    lattice->initialize();

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
    else { // If not, then let's do a warm up.
        pcout << endl << "*********** Entering stationary warmup phase ***********" << endl;
           
        int convergenceSteps = 10*max(max(Nx, Ny), Nz);
        T minDE = 1e-11; T dE = 100; T prevE = 0;
    
        for(auto &o: openings)
        	o->imposeBC(lattice, 0.0);

        if(saveInitState) {
            pcout << "Saving initial state with flow diverter..." << endl;
            // writeVTK(*lattice, sim, -1, porosityField);
            writeHDF5(*lattice, sim, -1, outDir, porosityField);
        }

        while(abs(dE) > minDE && stat_cycle < convergenceSteps )
        {
            lattice->collideAndStream();
    
            T cE = computeAverageEnergy(*lattice);
            dE = cE - prevE; prevE = cE;
    
            if(stat_cycle % 500 == 0) {
                pcout << "Delta energy: " << abs(dE) << "/" << minDE << "  Cycle: [" << stat_cycle << "/" << convergenceSteps <<"]" << std::endl;
            }
            
            stat_cycle++;
        }
        pcout << "Delta energy: " << abs(dE) << "/" << minDE << "  Cycle: [" << stat_cycle << "/" << convergenceSteps <<"]" << std::endl;
    
        pcout << endl << "*********** Entering transient simulation phase ***********" << endl;
    
        pcout << "Saving time step 0..." << endl;
        // writeVTK(*lattice, 0);
        // writeNPZ(*lattice, 0);
        writeHDF5(*lattice, sim, 0, outDir);
        
        // Set the counter back
        stat_cycle = 0;
    }
    
    pcout << "Starting computation..." << endl;

    while(stat_cycle*sim.C_t <= simLength + sim.C_t)
    {
        
        if(stat_cycle % 200 == 0) {
            T cE = computeAverageEnergy(*lattice);
            pcout << "\rTime: " << stat_cycle*sim.C_t << "s / " << simLength << "s" << " [" << stat_cycle << " / " << std::round(simLength/sim.C_t) << "] " << " - Energy: " << cE <<"         ";
            
            // Capture numerical divergence if appears
            if (std::isnan(cE)){
                pcout << "ERROR: NaN average energy! Saving state and stopping simulation" << std::endl;
                writeHDF5(*lattice, sim, stat_cycle, outDir, porosityField);
                return 0;
            } 
        }

        // Impose boundary conditions
        for(auto &o: openings)
            o->imposeBC(lattice, sim.C_t);

        // Calculate next step
        lattice->collideAndStream();

        // Advance time
        stat_cycle++;

        // Save output 
        if(stat_cycle % saveFrequency == 0) {
            pcout << "Writing output at: " << stat_cycle << " (" << stat_cycle*sim.C_t << " s)." << endl;
            // writeVTK(*lattice, sim, stat_cycle);
            // writeNPZ(*lattice, sim, stat_cycle);
            writeHDF5(*lattice, sim, stat_cycle, outDir, porosityField);   
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
