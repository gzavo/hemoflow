#include <map>
#include <string>
#include <sstream>
#include <algorithm>
#include <cstdlib>
#include <cstdint>
#include <iomanip>
#include <vector>
#include <filesystem>
#include <math.h>
#include <time.h>
#include <mpi.h>
#include <sys/mman.h>
#include <fcntl.h>

// For directory manipulations
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

int node_rank;

// Domain size
int Nx=0;
int Ny=0;
int Nz=0;

// Domain data
cnpy::NpyArray geometryFlag;
int8_t* gfData = NULL;

// Flow diverter (stent) data
cnpy::NpyArray stentFlag;
int8_t* sfData = NULL;
T linCoeff = 0.0;
T quadCoeff = 0.0;
T linCoeff_lb = 0.0;
T quadCoeff_lb = 0.0;


// Info on openings
cnpy::NpyArray openingIndex;
int8_t* oiData = NULL;
cnpy::NpyArray openingRadius;
double* orData = NULL;
cnpy::NpyArray openingQRatio;
double* oqData = NULL;
cnpy::NpyArray openingCenter;
double* ocData = NULL;
cnpy::NpyArray openingTangent;
double* otData = NULL;

string* FlowrateFuncs = NULL;

// Simulation parameters
T omega;
T C_l;  // Length conversion factor
T C_t;  // Time conversion factor
T C_r;  // Density conversion factor
T C_p;  // Pressure conversion factor (derived)
T C_m;  // Mass conversion factor (derived)
T Re;   

T tau_LB;

T U_AVG_LB = 0.05;

// Technical simulation parameters
bool useCheckpoint = true;
bool saveInitState = false;
int saveZeroState = 0;
int write2hdf5 = 0;
int blockSize;
int envelopeWidth = 1;
string outputFolder;
string workingFolder;
T simLength;
T saveFreqTime;
T checkpointFreqTime;

// Openings
vector<OpeningHandler*> openings;

// Simulation data structures
MultiBlockLattice3D<T, DESCRIPTOR> *lattice = NULL;
MultiNTensorField3D<T> *porosityField = NULL;

// Carreau parameters
//  B.M.  Johnston,  P.R.  Johnson,  S.  Corney,  and  D. Kilpatrick, “Non-Newtonian blood flow in human  right  coronary  arteries:  steady  state  simulations,” Journal  of  Biomechanics, 37, 709 – 720 (2004)
//  Y.I.  Cho  and  K.R.  Kensey,  “Effects  of  the  non-Newtonian viscosity of blood on flows in a   diseased   arterial   vessel.   Part   1:   steady   flows,” Biorheology28, 241 (1991)
// nu0 = 5.6e-5; nuInf = 3.5e-6;
T fluidDensity = 1055; // [kg/m3] Default blood density
T nu0 = 5.6e-5;     // [m^2/s]
T nuInf = 3.27e-6;   // [m^2/s]
T lambda = 3.331;
T n = 0.3568;
T nu_const;

//
// Directory handling routines here
//
// WARNING, really ugly! We need portable I/O

bool fileExists (const std::string& name) {
    ifstream f(name.c_str());
    return f.good();
}

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

/**
** mkpath - ensure all directories in path exist
** Algorithm takes the pessimistic view and works top-down to ensure
** each directory in path exists, rather than optimistically creating
** the last element and working backwards. It uses the custom makedir function below.
*/

// TODO: Unix specific, look for portable solution!
int do_mkdir(const char *path, mode_t mode)
{
    //Stat            st;
    struct stat st = {0};
    int    status = 0;

    if (stat(path, &st) != 0)
    {
        /* Directory does not exist. EEXIST for race condition */
        if (mkdir(path, mode) != 0 && errno != EEXIST)
            status = -1;
    }
    else if (!S_ISDIR(st.st_mode))
    {
        errno = ENOTDIR;
        status = -1;
    }

    return(status);
}

// mkpath(argv[i], 0777);
int mkpath(const char *path, mode_t mode)
{
    char           *pp;
    char           *sp;
    int             status;
    char           *copypath = strdup(path);

    status = 0;
    pp = copypath;
    while (status == 0 && (sp = strchr(pp, '/')) != 0)
    {
        if (sp != pp)
        {
            /* Neither root nor double slash in path */
            *sp = '\0';
            status = do_mkdir(copypath, mode);
            *sp = '/';
        }
        pp = sp + 1;
    }
    if (status == 0)
        status = do_mkdir(path, mode);
    free(copypath);
    return (status);
}

// *** Calculating LB parameters using Re on the inlet: Re = U_avg * D / nu
void calcSimulationParameters(T D_m)
{   

    T D_lb = D_m  / C_l;
    T U_avg = Re * nuInf / D_m;
    
    C_t =  U_AVG_LB / U_avg * C_l;

    T nuInf_lb = U_AVG_LB * D_lb / Re;
    T nu_ratio = nuInf_lb / nuInf;
    T nu0_lb = nu0 * nu_ratio;

    T tau = 3.0*nuInf_lb+0.5;

    omega = 1.0 / tau;

    C_r = fluidDensity;    // TODO IF we are simulating blood.... Note: only changes pressure output values, the simulation results are independent!

    C_p = C_r * C_l * C_l / (C_t * C_t);
    C_m = C_r * C_l * C_l * C_l;

    // TODO: convert linCoeff and quadCoeff
    linCoeff_lb = linCoeff * C_l*C_l * C_t / C_m;       // [ kg / (m2 s) ]
    quadCoeff_lb = quadCoeff * C_l*C_l * C_l / C_m;     // [ kg / m3 ]

    // TODO: add sanity check on parameters here

    // TODO: using the Smagorinsky dynamics as a base add dynamic viscosity (Carreau and rheoModel)
    // (ps: solve the Fokker-Plank in rheoModel with finite difference?)
    global::CarreauParameters().setNu0(nu0_lb);
    global::CarreauParameters().setNuInf(nuInf_lb);
    global::CarreauParameters().setLambda(3.313);  //1.
    global::CarreauParameters().setExponent(0.357);   //0.3

}

void processOpenings(string * FlowrateFuncs)
{
    int numOpenings = openingRadius.shape[0];
    pcout << "-> Number of openings to process: " << numOpenings << std::endl;

    for(int i=0; i<numOpenings; i++){
        int flag = int(oiData[i]);
        pcout << "Processing flag: " << flag << std::endl;

        if(flag < INLET_SVC){
            pcout << "WARNING! Wrong flag for an opening: " << flag << std::endl;
            continue;
        }

        int s = 3;
        vec3d dir(otData[gT2D(s,i,0)], otData[gT2D(s,i,1)], otData[gT2D(s,i,2)]);

        OpeningHandler *opening = new OpeningHandler(gfData, flag, orData[i] / C_l, dir);

        // LPA, RPA, and RPA side branches are all outlets
        // IVC, SVC, RHV, MHV, LHV have own flowrate funcs
        
        if (flag >= INLET_SVC && flag <= OUTLET_SIDE_V) {   // Inlets
            opening->createPoiseauilleProfile(U_AVG_LB);
            // opening->createBluntVelocityProfile(U_AVG_LB);
            if(flag==INLET_SVC)
                opening->loadScaleFunction(FlowrateFuncs[0], C_l);
            else if(flag==INLET_IVC)
                opening->loadScaleFunction(FlowrateFuncs[1], C_l);
            else if(flag==INLET_RHV)
                opening->loadScaleFunction(FlowrateFuncs[2], C_l);
            else if(flag==INLET_MHV)
                opening->loadScaleFunction(FlowrateFuncs[3], C_l);
            else if(flag==INLET_LHV)
                opening->loadScaleFunction(FlowrateFuncs[4], C_l);
            else if(flag==OUTLET_RPA)
                opening->loadScaleFunction(FlowrateFuncs[5], C_l);
            else if(flag==OUTLET_LPA)
                opening->loadScaleFunction(FlowrateFuncs[6], C_l);
            else if(flag==OUTLET_SIDE_V)
                opening->loadScaleFunction(FlowrateFuncs[7], C_l);
        }
        else
            // // Set pressure outlets with 1.0 pressure
            // opening->createConstantPressureProfile(1.0);

            // Set virtual outlet (copy velocity from last time step)
            opening->createVirtualOutletProfile(lattice);
        
        opening->printOpeningDetails();
        pcout << " " << std::endl;
        openings.push_back(opening);
    }

}


void writeHDF5(MultiBlockLattice3D<T,DESCRIPTOR>& lattice, plint iter, string outDir)
{

    T SaveTime = T();
    global::timer("SaveTime").restart();

    // Compute velocity in 3 dims, shear stress in 6 dims
    // Note the velovities are distributed on every processor
    MultiTensorField3D<double,3> DistributedVelocity = *computeVelocity(lattice);
    MultiScalarField3D<double> DistributedDensity = *computeDensity(lattice);
    MultiTensorField3D<double,6> DistributedShearStress = *computeShearStress(lattice);
    MultiScalarField3D<double> DistributedS_Norm = *computeSymmetricTensorNorm(*computeStrainRateFromStress(lattice));

    // Density/Velocity/... shared the same atomic block distribution!
    MultiBlockManagement3D VelocityBlockManagement = DistributedVelocity.getMultiBlockManagement();

    vector<plint> LocalBlockIDs = VelocityBlockManagement.getLocalInfo().getBlocks();

    // Start to count the writing time
    T FindAttributesTime = T();
    global::timer("FindAttributes").restart();

    // Again, Density/Velocity/Shear stress/S_norm... share the same distribution, so one ID vector is enough
    vector<vector<long unsigned int>> GlobalID;
    vector<float> VelocityX; vector<float> VelocityY; vector<float> VelocityZ; vector<float> Density;
    vector<float> SS1; vector<float> SS2; vector<float> SS3; vector<float> SS4; vector<float> SS5; vector<float> SS6;
    vector<float> SNorm;
    vector<int> this_rank;

    int RankID = global::mpi().getRank();

    // Now we loop through all local blocks on current MPI thread
    for(pluint iBlock=0; iBlock < LocalBlockIDs.size(); ++iBlock) {
        plint blockId = LocalBlockIDs[iBlock];

        // The "SmartBulk3D" object represents local atomic block in a global view, i.e. its bounding box coordinates are in global scale.
        // If you do not understand, go check the source codes of "MultiBlockManagement3D::findAllLocalRepresentations()"
        // Why we use it? Because we need to know which atomic blocks are stored on current MPI thread!
        SmartBulk3D LocalBulk(VelocityBlockManagement.getSparseBlockStructure(), VelocityBlockManagement.getEnvelopeWidth(), blockId);

        for(unsigned int i = LocalBulk.getBulk().x0; i <= LocalBulk.getBulk().x1; i++)
            for(unsigned int j = LocalBulk.getBulk().y0; j <= LocalBulk.getBulk().y1; j++)
                for(unsigned int k = LocalBulk.getBulk().z0; k <= LocalBulk.getBulk().z1; k++){
                    
                    // Now we convert the global scale coordinates to block local coordinates
                    unsigned int LocalX = LocalBulk.toLocalX(i);
                    unsigned int LocalY = LocalBulk.toLocalY(j);
                    unsigned int LocalZ = LocalBulk.toLocalZ(k);

                    GlobalID.push_back({k,j,i});

                    // Velocity
                    Array<double,3> const& foundVelocity = DistributedVelocity.getComponent(blockId).get(LocalX, LocalY, LocalZ);
                    // Note: Scale to physical unit before saving
                    float vel_scale = float(C_l/C_t);
                    VelocityX.push_back(float(foundVelocity[0])*vel_scale); VelocityY.push_back(float(foundVelocity[1])*vel_scale); VelocityZ.push_back(float(foundVelocity[2])*vel_scale);
                    // Density
                    double foundDensity = DistributedDensity.getComponent(blockId).get(LocalX, LocalY, LocalZ);
                    // Density.push_back(float(foundDensity)*1./3.*float(C_p));
                    Density.push_back(float(foundDensity));
                    // Shear Stress
                    Array<double,6> const& foundSS = DistributedShearStress.getComponent(blockId).get(LocalX, LocalY, LocalZ);
                    float SS_scale = C_m / (C_l*C_t*C_t);
                    SS1.push_back(float(foundSS[0])*SS_scale); SS2.push_back(float(foundSS[1])*SS_scale); SS3.push_back(float(foundSS[2])*SS_scale);
                    SS4.push_back(float(foundSS[3])*SS_scale); SS5.push_back(float(foundSS[4])*SS_scale); SS6.push_back(float(foundSS[5])*SS_scale);
                    // S_Norm
                    double foundS_Norm = DistributedS_Norm.getComponent(blockId).get(LocalX, LocalY, LocalZ);
                    SNorm.push_back(foundS_Norm*float(1./C_t));
                    // Rank of current mpi thread
                    this_rank.push_back(RankID);

                }
    }

    FindAttributesTime = global::timer("FindAttributes").stop();
    pcout << "Finding attributes time: " << FindAttributesTime << " sec" << endl;

    assert(GlobalID.size() > 0);

    ///////////////////////////// Saving HDF5 /////////////////////////////

    // Now save the partial local data to hdf5, if you dont understand, 
    // check (https://github.com/BlueBrain/HighFive/blob/master/src/examples/parallel_hdf5_collective_io.cpp)
    using namespace HighFive;
    
    FileAccessProps fapl;
    // Tell HDF5 to use MPI-IO
    fapl.add(MPIOFileAccess{MPI_COMM_WORLD, MPI_INFO_NULL});
    // Specify that we want all meta-data related operations to use MPI collective operations,
    // that is, all MPI ranks must participate in any HDF5 operations.
    fapl.add(MPIOCollectiveMetadata{});

    // Create the file as usual.
    std::string file_name = createFileName(outDir + "/output_", iter, 6);
    File file(file_name + ".h5", File::Truncate, fapl);

    // For compression
    DataSetCreateProps props;
    // Use chunking
    props.add(Chunking(std::vector<hsize_t>{100, 100, 100}));
    // Enable shuffle
    props.add(Shuffle());
    // Enable deflate
    props.add(Deflate(9));

    // Create the dataset as usual
    std::vector<size_t> Dims{(long unsigned int)Nz, (long unsigned int)Ny, (long unsigned int)Nx};
    DataSet velocity_x = file.createDataSet<float>("velocity_x", DataSpace(Dims), props);
    DataSet velocity_y = file.createDataSet<float>("velocity_y", DataSpace(Dims), props);
    DataSet velocity_z = file.createDataSet<float>("velocity_z", DataSpace(Dims), props);
    // Shear Stress
    DataSet SS_1 = file.createDataSet<float>("sigma_1", DataSpace(Dims), props);
    DataSet SS_2 = file.createDataSet<float>("sigma_2", DataSpace(Dims), props);
    DataSet SS_3 = file.createDataSet<float>("sigma_3", DataSpace(Dims), props);
    DataSet SS_4 = file.createDataSet<float>("sigma_4", DataSpace(Dims), props);
    DataSet SS_5 = file.createDataSet<float>("sigma_5", DataSpace(Dims), props);
    DataSet SS_6 = file.createDataSet<float>("sigma_6", DataSpace(Dims), props);
    // Density
    DataSet density = file.createDataSet<float>("density", DataSpace(Dims), props);
    // S_Norm
    DataSet S_Norm = file.createDataSet<float>("S_norm", DataSpace(Dims), props);
    // MPI rank
    DataSet Rank = file.createDataSet<int>("MPI_rank", DataSpace(Dims), props);

    auto xfer_props = DataTransferProps{};
    xfer_props.add(UseCollectiveIO{});

    // Each process writes the local attributes to the file
    velocity_x.select(ElementSet(GlobalID)).write(VelocityX, xfer_props);
    velocity_y.select(ElementSet(GlobalID)).write(VelocityY, xfer_props);
    velocity_z.select(ElementSet(GlobalID)).write(VelocityZ, xfer_props);
    // Shear Stress
    SS_1.select(ElementSet(GlobalID)).write(SS1, xfer_props);
    SS_2.select(ElementSet(GlobalID)).write(SS2, xfer_props);
    SS_3.select(ElementSet(GlobalID)).write(SS3, xfer_props);
    SS_4.select(ElementSet(GlobalID)).write(SS4, xfer_props);
    SS_5.select(ElementSet(GlobalID)).write(SS5, xfer_props);
    SS_6.select(ElementSet(GlobalID)).write(SS6, xfer_props);
    // Density
    density.select(ElementSet(GlobalID)).write(Density, xfer_props);
    // S_Norm
    S_Norm.select(ElementSet(GlobalID)).write(SNorm, xfer_props);
    // MPI Rank
    Rank.select(ElementSet(GlobalID)).write(this_rank, xfer_props);

    global::mpi().barrier();

    SaveTime = global::timer("SaveTime").stop();
    pcout << "Saving HDF5 time: " << SaveTime << " sec" << endl;

    T XDMFtime = T();
    global::timer("XDMFTime").restart();

    ///////////////////////////// Writing Xdmf /////////////////////////////
    if (global::mpi().isMainProcessor())
    {
        FILE *xmf = 0;

        /*
        * Open the file and write the header.
        */
        std::string xmf_name = createFileName(outDir + "/output_", iter, 6) + ".xmf";
        xmf = fopen(xmf_name.c_str(), "w");

        // HDF5 name
        // Find the last occurrence of the directory separator '/'
        size_t lastSlash = file_name.find_last_of('/');
        // Return the substring after the last '/'
        std::string h5_name = file_name.substr(lastSlash + 1);

        fprintf(xmf, "<?xml version=\"1.0\" ?>\n");
        fprintf(xmf, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
        fprintf(xmf, "<Xdmf Version=\"2.0\">\n");

        /*
        * Write the mesh description and the variables defined on the mesh.
        */
        fprintf(xmf, " <Domain>\n");

        fprintf(xmf, "   <Grid Name=\"mesh\" GridType=\"Uniform\">\n");
        // Regular mesh
        fprintf(xmf, "     <Topology TopologyType=\"3DCoRectMesh\" NumberOfElements=\"%d %d %d\"/>\n", Nz, Ny, Nx);
        fprintf(xmf, "     <Geometry GeometryType=\"Origin_DxDyDz\">\n");
        fprintf(xmf, "       <DataItem Name=\"Origin\" Dimensions=\"%d\" NumberType=\"Float\" Precision=\"4\" Format=\"XML\">\n", 3);
        fprintf(xmf, "          0 0 0\n");
        fprintf(xmf, "       </DataItem>\n");
        // Discretization step size
        fprintf(xmf, "       <DataItem Name=\"Spacing\" Dimensions=\"%d\" NumberType=\"Float\" Precision=\"4\" Format=\"XML\">\n", 3);
        fprintf(xmf, "          %f %f %f\n", C_l, C_l, C_l);
        fprintf(xmf, "       </DataItem>\n");
        fprintf(xmf, "     </Geometry>\n");
        fprintf(xmf, "     \n");
        // Density
        fprintf(xmf, "     <Attribute Name=\"Density [Pa]\" AttributeType=\"Scalar\" Center=\"Cell\">\n");
        fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", Nz, Ny, Nx);
        fprintf(xmf, "          %s.h5:/density\n", h5_name.c_str());
        fprintf(xmf, "       </DataItem>\n");
        fprintf(xmf, "     </Attribute>\n");
        fprintf(xmf, "     \n");
        // Velocities
        fprintf(xmf, "     <Attribute Name=\"Velocity-X [m/s]\" AttributeType=\"Scalar\" Center=\"Cell\">\n");
        fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", Nz, Ny, Nx);
        fprintf(xmf, "          %s.h5:/velocity_x\n", h5_name.c_str());
        fprintf(xmf, "       </DataItem>\n");
        fprintf(xmf, "     </Attribute>\n");
        fprintf(xmf, "     <Attribute Name=\"Velocity-Y [m/s]\" AttributeType=\"Scalar\" Center=\"Cell\">\n");
        fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", Nz, Ny, Nx);
        fprintf(xmf, "          %s.h5:/velocity_y\n", h5_name.c_str());
        fprintf(xmf, "       </DataItem>\n");
        fprintf(xmf, "     </Attribute>\n");
        fprintf(xmf, "     <Attribute Name=\"Velocity-Z [m/s]\" AttributeType=\"Scalar\" Center=\"Cell\">\n");
        fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", Nz, Ny, Nx);
        fprintf(xmf, "          %s.h5:/velocity_z\n", h5_name.c_str());
        fprintf(xmf, "       </DataItem>\n");
        fprintf(xmf, "     </Attribute>\n");
        fprintf(xmf, "     \n");
        // Shear Stress
        fprintf(xmf, "     <Attribute Name=\"Shear Stress 1 [1/m2s]\" AttributeType=\"Scalar\" Center=\"Cell\">\n");
        fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", Nz, Ny, Nx);
        fprintf(xmf, "          %s.h5:/sigma_1\n", h5_name.c_str());
        fprintf(xmf, "       </DataItem>\n");
        fprintf(xmf, "     </Attribute>\n");
        fprintf(xmf, "     <Attribute Name=\"Shear Stress 2 [1/m2s]\" AttributeType=\"Scalar\" Center=\"Cell\">\n");
        fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", Nz, Ny, Nx);
        fprintf(xmf, "          %s.h5:/sigma_2\n", h5_name.c_str());
        fprintf(xmf, "       </DataItem>\n");
        fprintf(xmf, "     </Attribute>\n");
        fprintf(xmf, "     <Attribute Name=\"Shear Stress 3 [1/m2s]\" AttributeType=\"Scalar\" Center=\"Cell\">\n");
        fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", Nz, Ny, Nx);
        fprintf(xmf, "          %s.h5:/sigma_3\n", h5_name.c_str());
        fprintf(xmf, "       </DataItem>\n");
        fprintf(xmf, "     </Attribute>\n");
        fprintf(xmf, "     \n");
        fprintf(xmf, "     <Attribute Name=\"Shear Stress 4 [1/m2s]\" AttributeType=\"Scalar\" Center=\"Cell\">\n");
        fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", Nz, Ny, Nx);
        fprintf(xmf, "          %s.h5:/sigma_4\n", h5_name.c_str());
        fprintf(xmf, "       </DataItem>\n");
        fprintf(xmf, "     </Attribute>\n");
        fprintf(xmf, "     <Attribute Name=\"Shear Stress 5 [1/m2s]\" AttributeType=\"Scalar\" Center=\"Cell\">\n");
        fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", Nz, Ny, Nx);
        fprintf(xmf, "          %s.h5:/sigma_5\n", h5_name.c_str());
        fprintf(xmf, "       </DataItem>\n");
        fprintf(xmf, "     </Attribute>\n");
        fprintf(xmf, "     <Attribute Name=\"Shear Stress 6 [1/m2s]\" AttributeType=\"Scalar\" Center=\"Cell\">\n");
        fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", Nz, Ny, Nx);
        fprintf(xmf, "          %s.h5:/sigma_6\n", h5_name.c_str());
        fprintf(xmf, "       </DataItem>\n");
        fprintf(xmf, "     </Attribute>\n");
        fprintf(xmf, "     \n");
        // S_Norm
        fprintf(xmf, "     <Attribute Name=\"S_Norm [1/s]\" AttributeType=\"Scalar\" Center=\"Cell\">\n");
        fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", Nz, Ny, Nx);
        fprintf(xmf, "          %s.h5:/S_norm\n", h5_name.c_str());
        fprintf(xmf, "       </DataItem>\n");
        fprintf(xmf, "     </Attribute>\n");
        fprintf(xmf, "     \n");
        // MPI Rank
        fprintf(xmf, "     <Attribute Name=\"MPI Rank\" AttributeType=\"Scalar\" Center=\"Cell\">\n");
        fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Int\" Precision=\"4\" Format=\"HDF\">\n", Nz, Ny, Nx);
        fprintf(xmf, "          %s.h5:/MPI_rank\n", h5_name.c_str());
        fprintf(xmf, "       </DataItem>\n");
        fprintf(xmf, "     </Attribute>\n");
        fprintf(xmf, "     \n");

        fprintf(xmf, "   </Grid>\n");
        fprintf(xmf, " </Domain>\n");

        /*
        * Write the footer and close the file.
        */
        fprintf(xmf, "</Xdmf>\n");
        fclose(xmf);
    }

    global::mpi().barrier();

    XDMFtime = global::timer("XDMFTime").stop();
    pcout << "Saving XDMF time: " << XDMFtime << " sec" << endl;

}


// *** Main simulation entry point
int main(int argc, char *argv[])
{

    // global::mpi().init(&argc, &argv);

    // Count the computation time
    T TotalTime = T();
    global::timer("TotalTime").restart();

    plbInit(&argc, &argv);

    pcout   << "********************************* " << endl
            << "*       hemoFlowCFD  v0.2       * " << endl
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

    MPI_Win win;
    
    // **************** Load in xml file
    try{

        // **************** Load in opening data
        pcout << "Loading in data file..." << std::endl;

        xml["simulation"]["outputDir"].read(outputFolder);
        
        outDir = outputFolder;

        // Check if output dir exists and accessible
        if (dirExists(outDir) <= 0) {
            pcout << "Output folder " << outDir << " does not exist! Creating it...." << std::endl;

            if (global::mpi().isMainProcessor())
                mkpath(outDir.c_str(), 0777);
        }

        // Sync up after creating the directory by the master process.
        global::mpi().barrier();

        pcout << "Output folder: " << outDir+"/" << std::endl;

        global::directories().setOutputDir(outDir+"/");

        // ***************************** XML informtion ***************************** //

        xml["simulation"]["blockSize"].read(blockSize);
        xml["simulation"]["simLength"].read(simLength);
        xml["simulation"]["saveFrequency"].read(saveFreqTime);
        int initflag;
        xml["simulation"]["saveInitState"].read(initflag);
        saveInitState = (bool)initflag;
        xml["simulation"]["saveZeroState"].read(saveZeroState);
        xml["simulation"]["write2hdf5"].read(write2hdf5);

        xml["fluid"]["Re"].read(Re);
        xml["fluid"]["U_AVG_LB"].read(U_AVG_LB);
        xml["fluid"]["tau_LB"].read(tau_LB);
        xml["fluid"]["kinematicViscosity"].read(nu_const);
        xml["fluid"]["fluidDensity"].read(fluidDensity);

        xml["flowdiverter"]["linCoeff"].read(linCoeff);
        xml["flowdiverter"]["quadCoeff"].read(quadCoeff);
        
        
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

        cnpy::npz_t geom_npz = cnpy::npz_load(npzFileName);
        pcout << "Input data elements: " << geom_npz.size() << std::endl;

        for (auto const& array : geom_npz) 
            pcout << "->" << array.first << std::endl;

        // ***************************** Computation domain information ***************************** //

        // Reading dx = C_l
        cnpy::NpyArray dxA = geom_npz["dx"];
        C_l = (dxA.data<double>())[0];

        // Reading Nx, Ny, Nz
        cnpy::NpyArray nx = geom_npz["Nx"];
        Nx = (nx.data<int>())[0];
        cnpy::NpyArray ny = geom_npz["Ny"];
        Ny = (ny.data<int>())[0];
        cnpy::NpyArray nz = geom_npz["Nz"];
        Nz = (nz.data<int>())[0];
        pcout << "Domain size: " << Nx << " x " << Ny << " x " << Nz << std::endl;

        pcout << "Resolution [m]: " << C_l << std::endl;

        // ***************************** Opening information ***************************** //

        // Loading information on openings
        openingIndex = geom_npz["openingIndex"];
        oiData = openingIndex.data<int8_t>();
        
        openingRadius = geom_npz["openingRadius"];
        orData = openingRadius.data<double>();
        
        openingQRatio = geom_npz["openingNormalizedQRatio"];
        oqData = openingQRatio.data<double>();
        
        openingCenter = geom_npz["openingCenter"];
        ocData = openingCenter.data<double>();
        
        openingTangent = geom_npz["openingTangent"];
        otData = openingTangent.data<double>();

        pcout << "Number of openings: " << openingRadius.shape[0] << std::endl;

        // ***************************** Flowrate functions (hardcoded) ***************************** //

        int num_funcs;
        xml["geometry"]["NumberofFlowrateFuncs"].read(num_funcs);

        FlowrateFuncs = new string[num_funcs];
        xml["geometry"]["inletFlowrateFunc_SVC"].read(FlowrateFuncs[0]);
        xml["geometry"]["inletFlowrateFunc_IVC"].read(FlowrateFuncs[1]);
        xml["geometry"]["inletFlowrateFunc_RHV"].read(FlowrateFuncs[2]);
        xml["geometry"]["inletFlowrateFunc_MHV"].read(FlowrateFuncs[3]);
        xml["geometry"]["inletFlowrateFunc_LHV"].read(FlowrateFuncs[4]);
        xml["geometry"]["outletFlowrateFunc_RPA"].read(FlowrateFuncs[5]);
        xml["geometry"]["outletFlowrateFunc_LPA"].read(FlowrateFuncs[6]);
        xml["geometry"]["outletFlowrateFunc_SIDE_V"].read(FlowrateFuncs[7]);
        xml["geometry"]["outletFlowrateFunc_SIDE_P"].read(FlowrateFuncs[8]); // pressure outlet
        xml["geometry"]["outletFlowrateFunc_OutFlow"].read(FlowrateFuncs[9]);


        // **************** Load in geometry flag data
        int world_rank, world_size;
        MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);

        int node_size;
        MPI_Comm node_comm;
        MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, world_rank, MPI_INFO_NULL, &node_comm);
        MPI_Comm_rank(node_comm, &node_rank);
        MPI_Comm_size(node_comm, &node_size);

        pcout << "Shared region node size: " << node_size << std::endl;

        MPI_Aint size_of_data;

        // ***************************** Geometry flag data ***************************** //

        // Only load the geometry flag data to the first thread on the same node
        if (node_rank == 0) {

            // Loading the geometry flag file (hdf5)
            string GfFileName;
            xml["geometry"]["GeometryFlag"].read(GfFileName);

            using namespace HighFive;
            // Create a new file using the default property lists.
            File file(GfFileName, File::ReadOnly);
            // let's create a dataset of this size
            auto dataset = file.getDataSet("geometryFlag");

            auto space = dataset.getStorageSize();
            pcout << "Geometry flag storage space: " << space << " bytes." << std::endl;

            // gfData = (int8_t*) malloc(sizeof(int8_t) * Nx*Ny*Nz);
            MPI_Win_allocate_shared(sizeof(int8_t) * Nx*Ny*Nz, sizeof(int8_t), MPI_INFO_NULL, node_comm, &gfData, &win);
            dataset.read(gfData);
            
        }
        else {
            int u_sh = sizeof(int8_t);
            MPI_Win_allocate_shared(0, sizeof(int8_t), MPI_INFO_NULL, node_comm, &gfData, &win);
            MPI_Win_shared_query(win, 0, &size_of_data, &u_sh, &gfData);
        }

        MPI_Win_fence(0, win);

        global::mpi().barrier();

        T inletD = 2.0 * orData[6]; // SVC [m]

        pcout << "Setting LBM parameters..." << std::endl;
        calcSimulationParameters(inletD);
        

    }
    catch (PlbIOException& exception) {
        pcout << "Error while processing input opening file " << paramXmlFileName
              << ": " << exception.what() << std::endl;
        return -1;
    }

    global::mpi().barrier();
    
    pcout   << "*********** Simulation parameters *********** " << endl
            << "size [LU]:   " << Nx << "x" << Ny << "x" << Nz << endl
            << "dx [m]:  " << C_l << endl
            << "dt [s]:  " << C_t << endl
            << "omega:  " << omega << endl
            << "nu:     " << 1./3. * (1./omega - 0.5) << endl
            << "Re_inlet: " << Re << endl
            << "U_avg(inlet) [lbm]: " << U_AVG_LB << endl 
            << "U_avg [m/s]: " << U_AVG_LB * C_l / C_t << endl << endl;

    int saveFrequency;
    saveFrequency = (int)round(saveFreqTime/C_t);
    pcout << "Saving frequency set to every " << saveFreqTime << " s (" << saveFrequency << " steps)." << endl;

    int checkpointFrequency = (int)round(checkpointFreqTime/C_t);
    
    if(useCheckpoint)     
        pcout << "Chekpointing will happen every " << checkpointFreqTime << " s (" << checkpointFrequency << " steps)." << endl;
    
    // Checkpoint file names relative to the output folder
    string chkParamFile = outDir+"/checkpoint_parameters.dat";
    string chkDataFile = outDir+"/checkpoint_lattice.dat";
    string chkParamFileOld = outDir+"/checkpoint_parameters_old.dat";
    string chkDataFileOld = outDir+"/checkpoint_lattice_old.dat";

    // TODO: IMPORTANT! - Work out proper sparse mode, we waste up to 90% numerical cells. Take a hint from HemoCell.
    bool sparse = true;
    T cSmago = 0.14;
    if(sparse) {
        pcout << "Setting simulation domain mask for sparse decomposition..." << endl;
        MultiScalarField3D<int> *flagMatrix = new MultiScalarField3D<int>(Nx,Ny,Nz);

        setToFunction(*flagMatrix, flagMatrix->getBoundingBox(), FlagMaskDomain3D<int8_t>(gfData, 1));

        pcout << "Creating sparse representation ..." << endl;
        
        // Create sparse representation
        MultiBlockManagement3D sparseBlockManagement =
                    computeSparseManagement(*plb::reparallelize(*flagMatrix, blockSize, blockSize, blockSize), envelopeWidth);

        lattice = new MultiBlockLattice3D<T, DESCRIPTOR> (sparseBlockManagement,
                                                          defaultMultiBlockPolicy3D().getBlockCommunicator(),
                                                          defaultMultiBlockPolicy3D().getCombinedStatistics(),
                                                          defaultMultiBlockPolicy3D().getMultiCellAccess<T,DESCRIPTOR>(),
                                                          new BackgroundDynamics(omega));
                                                        //   new SmagorinskyRegularizedDynamics<T,DESCRIPTOR>(omega, cSmago)); // Uncomment if simulation LES
    }
    else {
        lattice = new MultiBlockLattice3D<T, DESCRIPTOR>(Nx, Ny, Nz, 
                                                         new BackgroundDynamics(omega));
                                                        //  new SmagorinskyRegularizedDynamics<T,DESCRIPTOR>(omega, cSmago));
    }

    // Uncomment the following line to instantiate the Smagorinsky LES model,
    //   if the background-dynamics (i.e. the dynamics given to "lattice" in
    //   the constructor) is BGKdynamics instead of SmagorinskyBGKdynamics.
    
    // instantiateStaticSmagorinsky(*lattice, lattice->getBoundingBox(), cSmago);

    pcout << getMultiBlockInfo(*lattice) << endl;

    // // If there is data on porosity, set up porous layer in the simulation   
    // if(sfData != NULL) {
    //     pcout << "Setting up porous layer for flow diverter..." << std::endl;
    //     porosityField = defaultGenerateMultiNTensorField3D<T>(lattice->getMultiBlockManagement(), 1).release();
    //     applyProcessingFunctional(new InitializePorousField<T, int8_t>(sfData), porosityField->getBoundingBox(), *porosityField);
    //     integrateProcessingFunctional( new PorousForceFunctional<T, DESCRIPTOR>(linCoeff_lb, quadCoeff_lb), lattice->getBoundingBox(), *lattice, *porosityField);
    // }

    pcout << "Defining walls..." << std::endl;
    defineDynamics(*lattice, lattice->getBoundingBox(), new FlagMaskSingleDomain3D<int8_t>(gfData, 0), new NoDynamics<T, DESCRIPTOR>);
    defineDynamics(*lattice, lattice->getBoundingBox(), new FlagMaskSingleDomain3D<int8_t>(gfData, 1), new BounceBack<T, DESCRIPTOR>(1.0));

    pcout << "Wall defined" << std::endl;

    // gfData is very large whehn resolution is high, deallocate it to save more RAM space
    MPI_Win_free(&win);

    // TODO: add some reparallelize here, check if it plays nice with checkpointing

    // ********************************************** Openings **********************************************

    pcout << "Processing openings..." << std::endl;
    processOpenings(FlowrateFuncs);

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
        T LoadCheckpointTime = T();
        global::timer("LoadCheckpointTime").restart();

        loadBinaryBlock(*lattice, outDir+"/checkpoint_lattice.dat");

        LoadCheckpointTime = global::timer("LoadCheckpointTime").stop();
        
        pcout << "Checkpoint at iteration " << stat_cycle << " loaded succesfully." << "   Execution time: " << LoadCheckpointTime / 60 << " mins" << std::endl;
    }
    else { // If not, then let's do a warm up.
        pcout << endl << "*********** Entering stationary warmup phase ***********" << endl;
           
        int convergenceSteps = 10*max(max(Nx, Ny), Nz);
        T minDE = 1e-11; T dE = 100; T prevE = 0;
    
        for(auto &o: openings)
        	o->imposeBC(lattice, 0.0);

        if(saveInitState){
            pcout << "Saving initial state with flow diverter..." << endl;

            if(write2hdf5)  writeHDF5(*lattice, -1, outDir);

        }

        // Record the computing time of 500 iteration
        T WarmUpTime = T();
        global::timer("WarmUpTime").restart();
    
        while(abs(dE) > minDE && stat_cycle < convergenceSteps )
        {
            lattice->collideAndStream();

            T cE = computeAverageEnergy(*lattice);
            dE = cE - prevE; prevE = cE;
    
            if(stat_cycle % 500 == 0) {
                WarmUpTime = global::timer("WarmUpTime").stop();

                if (std::isnan(cE)){
                    pcout << "NaN average energy!" << std::endl;
                    if(write2hdf5)  writeHDF5(*lattice, -666, outDir);
                    return 0;
                } 

                pcout << "Delta energy: " << abs(dE) << "/" << minDE << "  Cycle: [" << stat_cycle << "/" \
                      << convergenceSteps <<"]      " << "Execution time: " << WarmUpTime / 60.0 << " mins" << std::endl;

                global::timer("WarmUpTime").restart();
            }
            
            stat_cycle++;
        }

        global::timer("WarmUpTime").stop();
    
        pcout << endl << "*********** Entering transient simulation phase ***********" << endl;

        if (saveZeroState){
            pcout << "Saving time step 0..." << endl;

                if(write2hdf5)  writeHDF5(*lattice, 0, outDir);

        }
        
        // Set the counter back
        stat_cycle = 0;
    }
    
    pcout << "Starting computation..." << endl;

    // Count the actual computing time
    T ComputeTime = T();
    global::timer("ComputeTime").restart();

    while(stat_cycle*C_t <= simLength + C_t)
    {
        
        // very inefficient! The idea is to capture the divergence point
        T cE = computeAverageEnergy(*lattice);

        if (std::isnan(cE)){
            pcout << "NaN average energy!" << std::endl;
            if(write2hdf5)  writeHDF5(*lattice, -666, outDir);
            return 0;
        } 

        
        if(stat_cycle % 1000 == 0) {
            pcout << "\rTime: " << stat_cycle*C_t << "s / " << simLength << "s" << " [" << stat_cycle << " / " << std::round(simLength/C_t) << "] " << " - Energy: " << cE << std::endl;
        }

        // Impose boundary conditions
        for(auto &o: openings)
            o->imposeBC(lattice, C_t);

        // Calculate next step
        lattice->collideAndStream();

        // Advance time
        stat_cycle++;

        // Save VTK output 
        if(stat_cycle % saveFrequency == 0) {

            global::timer("ComputeTime").stop();

            pcout << "Writing output at: " << stat_cycle << " (" << stat_cycle*C_t << " s)." << endl;

            if(write2hdf5)  writeHDF5(*lattice, stat_cycle, outDir);

            global::timer("ComputeTime").start();
            
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

            // Record the saving time
            T SaveCheckpointTime = T();
            global::timer("SaveCheckpointTime").restart();

            saveBinaryBlock(*lattice, chkDataFile);

            SaveCheckpointTime = global::timer("SaveCheckpointTime").stop();
        
            pcout << "Checkpoint at iteration " << stat_cycle << " saved succesfully." << "   Execution time: " << SaveCheckpointTime / 60 << " mins" << std::endl;
        }
    } // End of main loop

    ComputeTime = global::timer("ComputeTime").stop();

    TotalTime = global::timer("TotalTime").stop();

    pcout << endl << "Simulation done successfully :)" << endl;
    pcout << endl << "Execution time: " << TotalTime / 60 << " mins." << endl;
    pcout << endl << "Computation time: " << ComputeTime / 60 << " mins." << endl;

    return 0;
}
