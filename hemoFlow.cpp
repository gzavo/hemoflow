#include <map>
#include <string>
#include <sstream>
#include <algorithm>
#include <cstdlib>
#include <iomanip>
#include <vector>
#include <math.h>

#include <sys/types.h>
#include <sys/stat.h>

using namespace std;

#include "globals.h"
#include "helper.h"
#include "opening.h"
#include "porous.h"

// I/O
#include "cnpy.h"
#include "io/xdmfDataOutput.h"

/* ********** GLOBAL VARIABLES ************/

// Domain size
int Nx=0;
int Ny=0;
int Nz=0;

// Domain data
cnpy::NpyArray geometryFlag;
unsigned short* gfData = NULL;

// Stent geometry data
cnpy::NpyArray stentFlag;
unsigned short* sfData = NULL;

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

// Simulation parameters
T omega;
T dt;
T dx;
T Re;   

// Technical simulation parameters
int blockSize;
int envelopeWidth = 1;
string outputFolder;
string workingFolder;
T simLength;
T saveFreqTime;

// Openings
vector<OpeningHandler*> openings;

// Simulation data structures
MultiBlockLattice3D<T, DESCRIPTOR> *lattice = NULL;
OnLatticeBoundaryCondition3D<T, DESCRIPTOR> *boundaryCondition = NULL;
MultiNTensorField3D<T> *porosityField = NULL;

// Carreau parameters
//  B.M.  Johnston,  P.R.  Johnson,  S.  Corney,  and  D. Kilpatrick, “Non-Newtonian blood flow in human  right  coronary  arteries:  steady  state  simulations,” Journal  of  Biomechanics, 37, 709 – 720 (2004)
//  Y.I.  Cho  and  K.R.  Kensey,  “Effects  of  the  non-Newtonian viscosity of blood on flows in a   diseased   arterial   vessel.   Part   1:   steady   flows,” Biorheology28, 241 (1991)
T nu0 = 5.6e-5;     // [m^2/s]
T nuInf = 3.5e-6;   // [m^2/s]
T lambda = 3.331;
T n = 0.3568;

// Checks for a directory. Hopefuly a portable way. I want C++17...
int dirExists(string pathName)
{
    struct stat info;

    if( stat( pathName.c_str(), &info ) != 0 )
        return -1; // Cannot acces path
    else if( info.st_mode & S_IFDIR )  // S_ISDIR() doesn't exist on my windows 
        return 1;  // Path exists
    else
        return 0;  // Path does not exist
}

// *** Calculating LB parameters using Re on the inlet: Re = U_avg * D / nu
void calcSimulationParameters(T D_m)
{   
    T D_lb = D_m  / dx;
    T U_avg = Re * nuInf / D_m;
    
    dt =  U_AVG_LB / U_avg * dx;

    T nuInf_lb = U_AVG_LB * D_lb / Re;
    T nu_ratio = nuInf_lb / nuInf;
    T nu0_lb = nu0 * nu_ratio;

    T tau = 3.0*nuInf_lb+0.5;

    omega = 1.0 / tau;

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
        pcout << "-> Processing flag: " << flag << std::endl;

        if(flag < INLET){
            pcout << "WARNING! Wrong flag for an opening: " << flag << std::endl;
            continue;
        }

        int s = openingTangent.shape[1];
        vec3d dir(otData[gT2D(s,i,0)], otData[gT2D(s,i,1)], otData[gT2D(s,i,2)]);

        OpeningHandler *opening = new OpeningHandler(gfData, flag, orData[i] / dx, dir);
        
        if(flag == FIRST_OUTLET)
            opening->createConstantPressureProfile();
        else {
            if(flag == INLET)
                opening->createBluntVelocityProfile(U_AVG_LB);
            else {
                // Calculate outflow velocity based on Murray's law
                T u_out = Qin * oqData[i] / (pow(orData[i] / dx, 2) * M_PI);
                opening->createBluntVelocityProfile(u_out);
            }
            
            if(!inletFlowrateFunc.empty())
                opening->loadScaleFunction(inletFlowrateFunc);
        }
        
        opening->printOpeningDetails();
        openings.push_back(opening);
    }

}

// Write out data in vtk format
void writeVTK(MultiBlockLattice3D<T,DESCRIPTOR>& lattice, plint iter, MultiNTensorField3D<T> *field1 = NULL)
{
    VtkImageOutput3D<T> vtkOut(createFileName("vtk", iter, 6), dx);
    vtkOut.writeData<float>(*computeDensity(lattice), "density [Pa]", 1./3.* dx*dx/dt/dt);
    vtkOut.writeData<3,float>(*computeVelocity(lattice), "velocity [m/s]", dx/dt);
    vtkOut.writeData<6,float>(*computeShearStress(lattice), "sigma [1/m2s]", 1./(dx*dt*dt));
    vtkOut.writeData<float>(*computeSymmetricTensorNorm(*computeStrainRateFromStress(lattice)), "S_norm [1/s]", 1./dt );
    // TODO - output viscosity?
    
    vtkOut.writeData<float>(*computeDensity(lattice), "density [LBM]", 1.0);
    vtkOut.writeData<3,float>(*computeVelocity(lattice), "velocity [LBM]", 1.0);
    
    if (field1 != NULL)
       vtkOut.writeData<float>(*field1, "field1");
}

void writeHDF5(MultiBlockLattice3D<T,DESCRIPTOR>& lattice, plint iter, MultiNTensorField3D<T> *field1 = NULL)
{
    ParallelXdmfDataWriter3D xdmfOut(createFileName("hemoFlow_out", iter, 6));
    xdmfOut.writeDataField<float>(*computeDensity(lattice), "density");
    xdmfOut.writeDataField<float>(*computeVelocity(lattice), "velocity");
}

void writeNPZ(MultiBlockLattice3D<T,DESCRIPTOR>& lattice, plint iter)
{
    Box3D bb = lattice.getBoundingBox();
    long unsigned int nx = bb.getNx();
    long unsigned int ny = bb.getNy();
    long unsigned int nz = bb.getNz();

    TensorField3D<T,3> localVelocity(nx, ny, nz);
    copySerializedBlock(*computeVelocity(lattice), localVelocity);

    if(global::mpi().isMainProcessor()) {
        double *data = new double[3*nx*ny*nz];

        for(unsigned int i = 0; i < nx; i++) 
            for(unsigned int j = 0; j < ny; j++)
                for(unsigned int k = 0; k < nz; k++) {
                int idx = (i*nx*nz+j*nz+k)*3;

                data[idx]   = localVelocity.get(i, j, k)[0];
                data[idx+1] = localVelocity.get(i, j, k)[1];
                data[idx+2] = localVelocity.get(i, j, k)[2];
            }

        cnpy::npz_save(createFileName("output_", iter, 6) + ".npz", "velocity",&data[0],{3,nz,ny,nx},"w"); 
    }

    global::mpi().barrier();
}


// *** Main simulation entry point
int main(int argc, char *argv[])
{
    plbInit(&argc, &argv);

    // *** Reading in command line arguments
    string paramXmlFileName;
    try {
        global::argv(1).read(paramXmlFileName);
    }
    catch (PlbIOException& exception) {
        pcout << "Wrong parameters; the syntax is: "
              << (std::string)global::argv(0) << " parameter-input-file.xml" << std::endl;
        return -1;
    }

    XMLreader xml(paramXmlFileName);

    size_t folderIdx = paramXmlFileName.find_last_of("/\\");
    workingFolder = paramXmlFileName.substr(0, folderIdx);

    // *** Load in data files
    try {
        pcout << "Loading in data file..." << std::endl;

        xml["simulation"]["outputDir"].read(outputFolder);
        
        if (dirExists(workingFolder + "/" + outputFolder) == -1) {
        	pcout << "Output folder: " << workingFolder + "/" + outputFolder << " is not accessible! Exiting..." << std::endl;
        	return -1;
        }
        else if (dirExists(workingFolder + "/" + outputFolder) == 0) {
        	pcout << "Output folder needs to be created first! " << workingFolder + "/" + outputFolder << std::endl;
        }

        pcout << "Output folder: " << workingFolder + "/" + outputFolder+"/" << std::endl;

        global::directories().setOutputDir(workingFolder + "/" + outputFolder+"/");

        xml["simulation"]["Re"].read(Re);
        xml["simulation"]["blockSize"].read(blockSize);
        xml["simulation"]["simLength"].read(simLength);
        xml["simulation"]["saveFrequency"].read(saveFreqTime);

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

        // Reading dx
        cnpy::NpyArray dxA = geom_npz["dx"];
        dx = (dxA.data<double>())[0];

        pcout << "Resolution [m]: " << dx << std::endl;

        // Loading stent geometry
        stentFlag = geom_npz["stent"];
        
        if(stentFlag.shape.size() > 1) {    // Check if there is data on FD
            pcout << "Found flow diverter information to load." << std::endl;
            sfData = stentFlag.data<unsigned short>();
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
        T inletA = pow(orData[0] / dx, 2) * M_PI;


        string inletFlowrateFunc;
        xml["geometry"]["inletFlowrateFunc"].read(inletFlowrateFunc);
        processOpenings(workingFolder + "/" + inletFlowrateFunc, inletA);

        pcout << "Setting LBM parameters..." << std::endl;
        calcSimulationParameters(inletD);
    }
    catch (PlbIOException& exception) {
        pcout << "Error while processing input file " << paramXmlFileName
              << ": " << exception.what() << std::endl;
        return -1;
    }
    
    pcout   << "********************************* " << endl
            << "*     Simulation parameters     * " << endl
            << "********************************* " << endl
            << "size [LU]:   " << Nx << "x" << Ny << "x" << Nz << endl
            << "dx [m]:  " << dx << endl
            << "dt [s]:  " << dt << endl
            << "omega:  " << omega << endl
            << "nu:     " << 1./3. * (1./omega - 0.5) << endl
            << "Re_inlet: " << Re << endl
            << "U_avg(inlet) [lbm]: " << U_AVG_LB << endl 
            << "U_avg [m/s]: " << U_AVG_LB * dx / dt << endl << endl;

    int saveFrequency;
    saveFrequency = (int)round(saveFreqTime/dt);
    pcout << "Saving frequency set to every " << saveFreqTime << " s (" << saveFrequency << " steps)." << endl;

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

        lattice = new MultiBlockLattice3D<T, DESCRIPTOR> (sparseBlockManagement,
                                                          defaultMultiBlockPolicy3D().getBlockCommunicator(),
                                                          defaultMultiBlockPolicy3D().getCombinedStatistics(),
                                                          defaultMultiBlockPolicy3D().getMultiCellAccess<T,DESCRIPTOR>(),
                                                          // new GuoExternalForceBGKdynamics<T, DESCRIPTOR>(omega));
                                                          new ForcedCarreauDynamics<T, DESCRIPTOR>(omega));

    }
    else {
        lattice = new MultiBlockLattice3D<T, DESCRIPTOR>(Nx, Ny, Nz, new GuoExternalForceBGKdynamics<T, DESCRIPTOR>(omega) );
        // lattice = new MultiBlockLattice3D<T, DESCRIPTOR>(Nx, Ny, Nz, new ForcedCarreauDynamics<T, DESCRIPTOR>(omega) );
    }


    pcout << getMultiBlockInfo(*lattice) << endl;

    // If there is data on porosity, set up porous layer in the simulation   
    if(sfData != NULL) {
        pcout << "Setting up porous layer for flow diverter..." << std::endl;
        // TODO: calculate LB Darcy coefficient from config file values
        T DarcyCoeff = 0.5;

        porosityField = defaultGenerateMultiNTensorField3D<T>(lattice->getMultiBlockManagement(), 1).release();
        applyProcessingFunctional(new InitializePorousField<T, unsigned short>(sfData, DarcyCoeff), porosityField->getBoundingBox(), *porosityField);
        integrateProcessingFunctional( new DarcyPorousForceFunctional<T, DESCRIPTOR>, lattice->getBoundingBox(), *lattice, *porosityField);
    }

    boundaryCondition = createLocalBoundaryCondition3D<T,DESCRIPTOR>();

    pcout << "Defining walls..." << std::endl;
    defineDynamics(*lattice, lattice->getBoundingBox(), new FlagMaskSingleDomain3D<unsigned short>(gfData, 0), new NoDynamics<T, DESCRIPTOR>);
    defineDynamics(*lattice, lattice->getBoundingBox(), new FlagMaskSingleDomain3D<unsigned short>(gfData, 1), new BounceBack<T, DESCRIPTOR>(1.0));

    // TODO: add some reparallelize here

    pcout << "Setting values on openings..." << std::endl;
    for(auto &o: openings){
        o->setBC(lattice, boundaryCondition);
    }

    pcout << "Initializing lattice in equilibrium..." << std::endl;
    initializeAtEquilibrium (*lattice, lattice->getBoundingBox(), 1.0, Array<T,3>((T)0.,(T)0.,(T)0.) );

    pcout << "Finalizing lattice..." << std::endl;
    lattice->initialize();


    pcout << endl << "*********** Entering stationary warmup phase ***********" << endl;
    
    int stat_cycle = 0;
    int convergenceSteps = 10*max(max(Nx, Ny), Nz);
    T minDE = 1e-11; T dE = 100; T prevE = 0;

    for(auto &o: openings)
    	o->imposeBC(lattice, 0.0);

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

    pcout << "Saving initial state..." << endl;
    writeVTK(*lattice, 0, porosityField);
    // writeNPZ(*lattice, 0);
    //writeHDF5(*lattice, 0, porosityField);

    pcout << "Starting computation..." << endl;

    stat_cycle = 0;
    T currentTime = 0;
    while(currentTime <= simLength + dt)
    {
        
        if(stat_cycle % 200 == 0) {
            T cE = computeAverageEnergy(*lattice);
            pcout << "\rTime: " << currentTime << "s / " << simLength << "s" << " [" << stat_cycle << " / " << (int)(simLength/dt) << "] " << " - Energy: " << cE <<"         ";
        }

        // Impose boundary conditions
        for(auto &o: openings)
    		o->imposeBC(lattice, dt);

        // Calculate next step
        lattice->collideAndStream();

        // Advance time
        currentTime += dt; 
        stat_cycle++;

        // Save VTK output
        if(stat_cycle % saveFrequency == 0) {
            pcout << "Writing output at: " << stat_cycle << " (" << currentTime << " s)." << endl;
            writeVTK(*lattice, stat_cycle);
            // writeNPZ(*lattice, stat_cycle);
            //writeHDF5(*lattice, stat_cycle);
        }
    }

    pcout << endl << "Simulation done successfully :)" << endl;

    return 0;
}
