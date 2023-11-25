#include "opening.h"


OpeningHandler::OpeningHandler(int8_t* flagAray, int8_t flag_, T radius_lb, vec3d dirVec)
{
    flag = flag_;
    hasScaleFunction = false;
    cTimePos = 0;
    cTimeVal = 0;

    vector<int> xCoord; vector<int> yCoord; vector<int> zCoord;
    
    pcout << "-> Looping through flag array... " << std::endl;

    for(int i=0; i<Nx; i++)
        for(int j=0; j<Ny; j++)
            for(int k=0; k<Nz; k++) 
                if(flagAray[gT(i,j,k)]==flag) {
                    nodes.push_back({i,j,k});
                    xCoord.push_back(i); yCoord.push_back(j); zCoord.push_back(k);                
                }
            
    if(nodes.size() == 0) {
        pcout << "WARNING! Non-existing opening (zero size)!" << std::endl;
        return;
    }

    pcout << "-> Calculating bounding box and center point... " << std::endl;
    boundingBox = new Box3D(*min_element(xCoord.begin(), xCoord.end()),
                            *max_element(xCoord.begin(), xCoord.end()),
                            *min_element(yCoord.begin(), yCoord.end()),
                            *max_element(yCoord.begin(), yCoord.end()),
                            *min_element(zCoord.begin(), zCoord.end()),
                            *max_element(zCoord.begin(), zCoord.end()) );
    

    R = radius_lb; // Radius data from the voxelizer code (comes from the centerline caluclation).

    // Geometric centerpoint of the opening
    center.x = accumulate( xCoord.begin(), xCoord.end(), 0.0) / xCoord.size();
    center.y = accumulate( yCoord.begin(), yCoord.end(), 0.0) / yCoord.size();
    center.z = accumulate( zCoord.begin(), zCoord.end(), 0.0) / zCoord.size();

    // Direction in which the opening faces
    direction.x = dirVec.x; direction.y = dirVec.y; direction.z = dirVec.z;
}

void OpeningHandler::printOpeningDetails()
{
    pcout << "Opening parameters:" << std::endl;
    
    if(flag > FLUID && flag < OUTLET_RPA)
        pcout << "=> This is the inlet." << std::endl;
    if(flag == OUTLET_RPA)
        pcout << "=> This is the smallest pressure outlet." << std::endl;

    pcout << "-> flag: " << int(flag) << std::endl;
    pcout << "-> radius [lb]: " << R << std::endl;
    //pcout << "-> Q rate: " << Qr << std::endl;
    pcout << "-> Center [lb]: " << center.x << " " << center.y << " " << center.z << std::endl;
    pcout << "-> Normal: " << direction.x << " " << direction.y << " " << direction.z << std::endl;
    pcout << "-> Area [lb]: " << nodes.size() << std::endl;
    pcout << "-> Scale function: " << hasScaleFunction << " length: " << scaleSignal.size() << std::endl;

}

void OpeningHandler::createPoiseauilleProfile(T u_avg)
{
    pcout << "-> Creating direction-corrected Pouseuille velocity profile on flag: " << int(flag) << std::endl;

    // Sanity check
    if ( !(boundingBox->x0 == boundingBox->x1 || boundingBox->y0 == boundingBox->y1 || boundingBox->z0 == boundingBox->z1) ) {
        pcout << "!!! ERROR: The opening is not parallel to any major plane! This functionality is not implemented, the opening will not funxtion!" << std::endl;
    }

    // Paraboloid height from average
    T u_max = 2. * u_avg;

    // Search for the farthest point from the centerpoint of the opening.
    T l_max = 0.0;
    vec3d v_max(0,0,0);

    for(auto const& v: nodes) {
        vec3d dist(v.x-center.x, v.y-center.y, v.z-center.z); 
        T length = dist.norm();

        if(length > l_max) {
            l_max = length;
            v_max.set(dist.x, dist.y, dist.z);
        }
    }

    // The normal radial direction should be perpendicular to it
    vec3d v_min = v_max.cross(direction);
    v_min.normalize();
    v_min = v_min * this->R;

    pcout << "The scaling coordinate system (R=" << this->R << "): " << v_min.norm() << ", " << l_max << std::endl;

    // Now use the closes and farthest points as coordinate system to scale the paraboloid.
    // Note: these two direction vectors are supposed to be perpendicular!
    for(auto const& v: nodes) {
        vec3d cR(v.x-center.x, v.y-center.y, v.z-center.z); 

        T l_min = v_min.norm();

        T proj_r_min = cR.dot(v_min) / l_min;
        T proj_r_max = cR.dot(v_max) / l_max * (l_min /l_max); // Scale it to r_min

        // At this point both projections are in (0,l_min == R), so the projected vector length will also be in (0, R)

        T proj_r_len2 = proj_r_max * proj_r_max + proj_r_min * proj_r_min;
        T r_min2 = l_min*l_min;

        vec3d vel = direction * (u_max / r_min2 * (r_min2 - proj_r_len2));
        
        velArr[v.x][v.y][v.z].set(vel.x, vel.y, vel.z);

    }

}

void OpeningHandler::createBluntVelocityProfile(T u_avg)
{
    pcout << "-> Creating blunt velocity profile on flag: " << int(flag) << std::endl;

    // Sanity check
    if ( !(boundingBox->x0 == boundingBox->x1 || boundingBox->y0 == boundingBox->y1 || boundingBox->z0 == boundingBox->z1) ) {
        pcout << "!!! ERROR: The opening is not parallel to any major plane! This functionality is not implemented, the opening will not funxtion!" << std::endl;
    }

    vec3d vel = direction * u_avg;

    for(auto const& v: nodes) {
        velArr[v.x][v.y][v.z].set(vel.x, vel.y, vel.z);
    }
}

void OpeningHandler::createConstantPressureProfile(T density)
{
    pcout << "-> Creating constant pressure on flag: " << int(flag) << std::endl;

    // if(flag < FIRST_OUTLET) {
    //     pcout << "WARNING! Calculating pressure profile for a non-pressure opening!" << std::endl;
    // }

    for(auto const& v: nodes) {
        presArr[v.x][v.y][v.z] = density;
    }
}

void OpeningHandler::loadScaleFunction(string fileName, T C_l)
{
    pcout << "-> Loading scale function: " << fileName << std::endl;

    plb_ifstream finSign(fileName.c_str());
    //istream &finSign = pfinSign.getOriginalStream();

    if(!finSign.good()) {
        pcout << "WARNING!!! Flow rate scale file " << fileName << " is not readable!" << std::endl;
        // hasScaleFunction = false;
    }


    T time, Q;
    int Ns;

    finSign >> Ns;
    global::mpi().bCast(&Ns, 1);
    global::mpi().barrier();

    for(int i=0; i<Ns; i++) {
        finSign >> time; finSign >> Q;
        // A_p = dx^2 * A_lb (Unit casting)
        // Note we do not use A = pi * r^2 here because opening may not be a circle
        // T area = nodes.size() * C_l*C_l;
        T area = pow(R*C_l, 2) * M_PI;
        if (i == Ns-1) pcout << "cross-section circle area: " << area << std::endl;
        // The average velocity on opening
        T value = Q / area;
        global::mpi().bCast(&time, 1);
        global::mpi().bCast(&value, 1);
        global::mpi().barrier();
        scaleTime.push_back(time); scaleSignal.push_back(value);
    }

    finSign.close();

    hasScaleFunction = true;

    pcout << fileName << " loaded with " << scaleTime.size() << " data points." << std::endl;
}

void OpeningHandler::setBC(MultiBlockLattice3D<T, DESCRIPTOR> *lattice)
{
    // SET SMALLEST OUTLET AS PRESSURE OUTLET
    // Set all outlets as pressure outlets
    if (flag < OUTLET_RPA) {
        OnLatticeBoundaryCondition3D<T, DESCRIPTOR> *bc = createLocalBoundaryCondition3D<T,DESCRIPTOR>();
        // OnLatticeBoundaryCondition3D<T, DESCRIPTOR> *bc = createZouHeBoundaryCondition3D<T,DESCRIPTOR>();
        bc->setVelocityConditionOnBlockBoundaries(*lattice, *boundingBox, boundary::dirichlet);
    }
    else {
        // OnLatticeBoundaryCondition3D<T, DESCRIPTOR> *bc = createLocalBoundaryCondition3D<T,DESCRIPTOR>();
        // OnLatticeBoundaryCondition3D<T, DESCRIPTOR> *bc = createZouHeBoundaryCondition3D<T,DESCRIPTOR>();
        // OnLatticeBoundaryCondition3D<T, DESCRIPTOR> *bc = createInterpBoundaryCondition3D<T,DESCRIPTOR>();
        // OnLatticeBoundaryCondition3D<T, DESCRIPTOR> *bc = createEquilibriumBoundaryCondition3D<T,DESCRIPTOR>();
        // OnLatticeBoundaryCondition3D<T, DESCRIPTOR> *bc = createInterpBoundaryCondition3D<T,DESCRIPTOR>();
        // bc->setPressureConditionOnBlockBoundaries(*lattice, *boundingBox, boundary::dirichlet);

        // Virtual outlet
        MultiScalarField3D<T> *rhoBar = generateMultiScalarField<T>((MultiBlock3D&) *lattice, 2).release();
        rhoBar->toggleInternalStatistics(false);

        MultiTensorField3D<T,3> *j = generateMultiTensorField<T,3>((MultiBlock3D&) *lattice, 2).release();
        j->toggleInternalStatistics(false);

        std::vector<MultiBlock3D*> bcargs;
        bcargs.push_back(lattice);
        bcargs.push_back(rhoBar);
        bcargs.push_back(j);

        integrateProcessingFunctional(new VirtualOutlet<T,DESCRIPTOR>(1.0, lattice->getBoundingBox(), 1),
                *boundingBox, bcargs, 2);
        setBoundaryVelocity(*lattice, *boundingBox, Array<T,3>((T) 0, (T) 0, (T) 0));
    }
}

void OpeningHandler::setVelocityProfile(MultiBlockLattice3D<T, DESCRIPTOR> *lattice, field3D &velocityArr)
{
    // For flowrate, additional measurements required after lattice initialization
    if(flag < OUTLET_RPA) {
        // Set given profile
        setBoundaryVelocity(*lattice, *boundingBox, VelocityProfile3D<T,DESCRIPTOR>(&velocityArr, 1.0));
    }
    else {
        pcout << "WARNING, opening " << flag << " incorrectly addressed as velocity opening (instead of pressure)!" << std::endl;
    }
}

void OpeningHandler::setPressureProfile(MultiBlockLattice3D<T, DESCRIPTOR> *lattice, scalar3D &pressureArr)
{
    // For flowrate, additional measurements required after lattice initialization
    if(flag >= OUTLET_RPA) {
        // Set given profile
        setBoundaryDensity(*lattice, *boundingBox, PressureProfile3D<T,DESCRIPTOR>(&pressureArr, 1.0));
    }
    else {
        pcout << "WARNING, opening " << flag << " incorrectly addressed as pressure opening (instead of velocity)!" << std::endl;
    }
}

void OpeningHandler::imposeBC(MultiBlockLattice3D<T, DESCRIPTOR> *lattice, T dt)
{
    T scale = 1.0;

    if(hasScaleFunction)
    {
        int len = scaleTime.size();

        cTimeVal += dt;
        if(scaleTime[cTimePos] < cTimeVal)  // TODO : We might need to skip some positions if simulation dt is too large. (With LBM, heck no....)
            cTimePos++;

        if (cTimePos > len-1) {
            cTimePos = 0;
            cTimeVal -= scaleTime[len-1];
        }

        if(cTimePos < len-1)
            scale = interpolate(scaleTime[cTimePos], scaleTime[cTimePos+1], cTimeVal, scaleSignal[cTimePos], scaleSignal[cTimePos+1]);
        else
            scale = interpolate(scaleTime[cTimePos], scaleTime[0], cTimeVal, scaleSignal[cTimePos], scaleSignal[0]);

        //scale /= scaleDivider;    // WTF is this?
    }

    // pcout << "scale: " << scale << std::endl;
    
    if (flag >= OUTLET_RPA) {
        // setBoundaryDensity(*lattice, *boundingBox, 1.0);
        // setBoundaryDensity(*lattice, *boundingBox, PressureProfile3D<T,DESCRIPTOR>(&presArr, scale)); // For time dependent pressure boundary
                
    }
    else {
        setBoundaryVelocity(*lattice, *boundingBox, VelocityProfile3D<T,DESCRIPTOR>(&velArr, scale));
    }
}

OpeningHandler::~OpeningHandler()
{

}