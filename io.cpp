#include "io.h"

// ***********************************
// *** Directory handling routines ***
// ***********************************

// WARNING, not portable! We need portable I/O code in the future.

bool fileExists (const std::string& name) {
    ifstream f(name.c_str());
    return f.good();
}

// Checks for a directory. Hopefuly a portable way. TODO: replace with C++17 method.
int dirExists(const string& pathName)
{
    struct stat info{};

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


// ****************************
// *** Data saving routines ***
// ****************************

// Write out data in VTK image format (.vti) in SI units
void writeVTK(MultiBlockLattice3D<T,DESCRIPTOR>& lattice, const SimPar &sim, plint iter, const OutputOptions &out)
{
    VtkImageOutput3D<T> vtkOut(createFileName("vtk", iter, 6), sim.C_l, out.origin);

    if (out.velocity)
        vtkOut.writeData<3,float>(*computeVelocity(lattice), "velocity [m/s]", sim.C_l/sim.C_t);

    if (out.pressure) {
        // Gauge pressure p = c_s^2 (rho - 1) C_p. Subtract the reference density in double
        // precision before converting to float to avoid cancellation.
        std::unique_ptr<MultiScalarField3D<T> > rhoDeviation = subtract(*computeDensity(lattice), (T)1.0);
        vtkOut.writeData<float>(*rhoDeviation, "pressure [Pa]", (float)(sim.C_p / 3.0));
    }

    if (out.stress) {
        // Viscous stress sigma = -(1 - omega/2) Pi^neq (lattice units), times C_p for [Pa].
        vtkOut.writeData<6,float>(*computePiNeq(lattice), "stress [Pa]", (float)(-(1.0 - sim.omega / 2.0) * sim.C_p));
    }
}

#ifdef HEMOFLOW_WITH_HDF5
// TODO - too slow, optimize the arrays (MPI rank is now saved in every lattice?).
// TODO - Optimize chunk size.
// TODO - save as vectors and matrices instead of 3D scalar arrays (also modify xdmf) - https://github.com/BlueBrain/HighFive/blob/master/src/examples/create_dataset_double.cpp
void writeHDF5(MultiBlockLattice3D<T,DESCRIPTOR>& lattice, const SimPar &sim, plint iter, string outDir,
               const Array<double, 3> &origin, MultiNTensorField3D<T> *field1)
{

    T SaveTime = T();
    global::timer("SaveTime").restart();

    // Compute velocity in 3 dims, shear stress in 6 dims
    // Note the velovities are distributed on every processor
    MultiTensorField3D<double,3> DistributedVelocity = *computeVelocity(lattice);
    MultiScalarField3D<double> DistributedDensity = *computeDensity(lattice);
    // Viscous stress from the non-equilibrium momentum flux (computeShearStress is zero for the bulk dynamics used here)
    MultiTensorField3D<double,6> DistributedShearStress = *computePiNeq(lattice);
    MultiScalarField3D<double> DistributedS_Norm = *computeSymmetricTensorNorm(*computeStrainRateFromStress(lattice));
    // MultiScalarField3D<double> DistributedField1 = *field1; // Used for additional fields, e.g. porosity, do any necessary calculations here. 

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
    vector<float> Field1;
    vector<int> this_rank;

    int RankID = global::mpi().getRank();

    // Now we loop through all local blocks on current MPI thread
    for(long blockId : LocalBlockIDs) {
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
                    float vel_scale = float(sim.C_l/sim.C_t);
                    VelocityX.push_back(float(foundVelocity[0])*vel_scale); VelocityY.push_back(float(foundVelocity[1])*vel_scale); VelocityZ.push_back(float(foundVelocity[2])*vel_scale);
                    // Density
                    double foundDensity = DistributedDensity.getComponent(blockId).get(LocalX, LocalY, LocalZ);
                    // Gauge pressure [Pa]
                    Density.push_back(float((foundDensity - 1.0) * sim.C_p / 3.0));
                    // Shear Stress
                    Array<double,6> const& foundSS = DistributedShearStress.getComponent(blockId).get(LocalX, LocalY, LocalZ);
                    float SS_scale = float(-(1.0 - sim.omega / 2.0) * sim.C_p);  // Pi^neq -> viscous stress [Pa]
                    SS1.push_back(float(foundSS[0])*SS_scale); SS2.push_back(float(foundSS[1])*SS_scale); SS3.push_back(float(foundSS[2])*SS_scale);
                    SS4.push_back(float(foundSS[3])*SS_scale); SS5.push_back(float(foundSS[4])*SS_scale); SS6.push_back(float(foundSS[5])*SS_scale);
                    // S_Norm
                    double foundS_Norm = DistributedS_Norm.getComponent(blockId).get(LocalX, LocalY, LocalZ);
                    SNorm.push_back(foundS_Norm*float(1./sim.C_t));
                    // Additional field - No unit conversion!
                    if (field1 != nullptr) {
                        double foundField1 = *field1->getComponent(blockId).get(LocalX, LocalY, LocalZ);
                        Field1.push_back(foundField1); 
                    }
                    else {
                        Field1.push_back(0.0);
                    }    
                    // Rank of current mpi thread
                    this_rank.push_back(RankID);

                }
    }

    FindAttributesTime = global::timer("FindAttributes").stop();
    pcout << "Finding attributes time: " << FindAttributesTime << " sec" << endl;

    assert(!GlobalID.empty());

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
    props.add(Deflate(7));

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
    // Field1
    DataSet Field1_data = file.createDataSet<float>("Field1", DataSpace(Dims), props);
    
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
    // Field1
    Field1_data.select(ElementSet(GlobalID)).write(Field1, xfer_props);
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
        FILE *xmf = nullptr;

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
        fprintf(xmf, "          %g %g %g\n", origin[2], origin[1], origin[0]);
        fprintf(xmf, "       </DataItem>\n");
        // Discretization step size
        fprintf(xmf, "       <DataItem Name=\"Spacing\" Dimensions=\"%d\" NumberType=\"Float\" Precision=\"4\" Format=\"XML\">\n", 3);
        fprintf(xmf, "          %f %f %f\n", sim.C_l, sim.C_l, sim.C_l);
        fprintf(xmf, "       </DataItem>\n");
        fprintf(xmf, "     </Geometry>\n");
        fprintf(xmf, "     \n");
        // Density
        fprintf(xmf, "     <Attribute Name=\"Pressure [Pa]\" AttributeType=\"Scalar\" Center=\"Cell\">\n");
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
        fprintf(xmf, "     <Attribute Name=\"Stress 1 [Pa]\" AttributeType=\"Scalar\" Center=\"Cell\">\n");
        fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", Nz, Ny, Nx);
        fprintf(xmf, "          %s.h5:/sigma_1\n", h5_name.c_str());
        fprintf(xmf, "       </DataItem>\n");
        fprintf(xmf, "     </Attribute>\n");
        fprintf(xmf, "     <Attribute Name=\"Stress 2 [Pa]\" AttributeType=\"Scalar\" Center=\"Cell\">\n");
        fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", Nz, Ny, Nx);
        fprintf(xmf, "          %s.h5:/sigma_2\n", h5_name.c_str());
        fprintf(xmf, "       </DataItem>\n");
        fprintf(xmf, "     </Attribute>\n");
        fprintf(xmf, "     <Attribute Name=\"Stress 3 [Pa]\" AttributeType=\"Scalar\" Center=\"Cell\">\n");
        fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", Nz, Ny, Nx);
        fprintf(xmf, "          %s.h5:/sigma_3\n", h5_name.c_str());
        fprintf(xmf, "       </DataItem>\n");
        fprintf(xmf, "     </Attribute>\n");
        fprintf(xmf, "     \n");
        fprintf(xmf, "     <Attribute Name=\"Stress 4 [Pa]\" AttributeType=\"Scalar\" Center=\"Cell\">\n");
        fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", Nz, Ny, Nx);
        fprintf(xmf, "          %s.h5:/sigma_4\n", h5_name.c_str());
        fprintf(xmf, "       </DataItem>\n");
        fprintf(xmf, "     </Attribute>\n");
        fprintf(xmf, "     <Attribute Name=\"Stress 5 [Pa]\" AttributeType=\"Scalar\" Center=\"Cell\">\n");
        fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", Nz, Ny, Nx);
        fprintf(xmf, "          %s.h5:/sigma_5\n", h5_name.c_str());
        fprintf(xmf, "       </DataItem>\n");
        fprintf(xmf, "     </Attribute>\n");
        fprintf(xmf, "     <Attribute Name=\"Stress 6 [Pa]\" AttributeType=\"Scalar\" Center=\"Cell\">\n");
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
        // Field1
        fprintf(xmf, "     <Attribute Name=\"Additional Field [-]\" AttributeType=\"Scalar\" Center=\"Cell\">\n");
        fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", Nz, Ny, Nx);
        fprintf(xmf, "          %s.h5:/Field1\n", h5_name.c_str());
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

#endif  // HEMOFLOW_WITH_HDF5

void writeNPZ(MultiBlockLattice3D<T,DESCRIPTOR>& lattice, const SimPar &sim, plint iter)
{
    Box3D bb = lattice.getBoundingBox();
    size_t nx = bb.getNx();
    size_t ny = bb.getNy();
    size_t nz = bb.getNz();

    TensorField3D<T,3> localVelocity(nx, ny, nz);
    copySerializedBlock(*computeVelocity(lattice), localVelocity);

    if(global::mpi().isMainProcessor()) {
        const double velScale = sim.C_l / sim.C_t;
        std::vector<double> data(3 * nx * ny * nz);

        // C-ordered (Nx, Ny, Nz, 3) array, like the geometry arrays
        for(size_t i = 0; i < nx; i++)
            for(size_t j = 0; j < ny; j++)
                for(size_t k = 0; k < nz; k++) {
                    size_t idx = ((i * ny + j) * nz + k) * 3;
                    for (int d = 0; d < 3; d++)
                        data[idx + d] = localVelocity.get(i, j, k)[d] * velScale;
                }

        std::string fileName = global::directories().getOutputDir() + createFileName("output_", iter, 6) + ".npz";
        cnpy::npz_save(fileName, "velocity", data.data(), {nx, ny, nz, 3}, "w");
    }

    global::mpi().barrier();
}
