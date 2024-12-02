#include "opening.h"


OpeningHandler::OpeningHandler(const unsigned short* flagAray, GeometryLabel flag_, OpeningType type_, T radius_lb, vec3d dirVec)
{
    flag = flag_;
    type = type_;
    hasScaleFunction = false;
    cTimePos = 0;
    cTimeVal = 0;
    cScale = 1.0;
    bcParameter = 1.0;

    vector<int> xCoord; vector<int> yCoord; vector<int> zCoord;
    
    pcout << "-> Configuring LBM nodes on the opening... " << std::endl;

    for(int i=0; i<Nx; i++)
        for(int j=0; j<Ny; j++)
            for(int k=0; k<Nz; k++) 
                if(flagAray[gT(i,j,k)]==flag) {
                    nodes.push_back({i,j,k});
                    xCoord.push_back(i); yCoord.push_back(j); zCoord.push_back(k);                
                }


    if(nodes.empty()) {
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

void OpeningHandler::printOpeningDetails(SimPar s)
{
    pcout << "---- Opening parameters -----" << std::endl;
    pcout << "-> Name: " << getName() << std::endl;
    pcout << "-> Geometry flag: " << flag << std::endl;
    pcout << "-> Opening type: " << type << std::endl;
    pcout << "-> Radius [m]: " << R * s.C_l << std::endl;
    pcout << "-> Center [lb]: " << center.x << " " << center.y << " " << center.z << std::endl;
    pcout << "-> Normal: " << direction.x << " " << direction.y << " " << direction.z << std::endl;
    pcout << "-> Area [lb]: " << nodes.size() << std::endl;
    pcout << "-> Scale function: " << hasScaleFunction << " length: " << scaleSignal.size() << std::endl;
    pcout << "-> BC parameter [lb]: " << bcParameter << std::endl;

    if(type == OPENING_VELOCITY) {
        T flowRate = 0.5 * bcParameter * nodes.size() * s.C_l * s.C_l * s.C_l / s.C_t; 
        pcout << "-> BC Q [m^3/s]: " << flowRate << std::endl;
        pcout << "-> BC Q [lbm]: " << 0.5 * bcParameter * nodes.size() << std::endl;

        T maxVel = bcParameter * s.C_l / s.C_t;
        pcout << "-> BC max vel. [m/s]: " << maxVel << std::endl;
        if(bcParameter > s.U_max_lb) {
            pcout << "-> *WARNING*: Maximum velocity exeeds expected maximum velocity: " << s.U_max_lb * s.C_l / s.C_t << std::endl;
        }
    }
    else if (type == OPENING_PRESSURE) {
        if(bcParameter < 0.5 || bcParameter > 1.5) {
            pcout << "-> *WARNING*: LB pressure is in an unstable regime: " << bcParameter << ". Consider decreasing the time-step size." << std::endl;
        }        
    }
}

// Set Q or p or other BC parameter, convert it to LBM units.
void OpeningHandler::setBCParameter(T bcParameter_, SimPar s)
{
    if(type == OPENING_VELOCITY){
        // [m/s] -> LB max velocity
        bcParameter = bcParameter_ * s.C_t / s.C_l;
    }
    else if (type == OPENING_PRESSURE) {
        // p -> LB density
        bcParameter = bcParameter_ / 3.0  / s.C_p + 1.0;   // 1.0 is defined as density for 0 pressure
    }
    else if (type == OPENING_MURRAY) {
        //VFR in LB units?
        bcParameter = bcParameter_; //* s.C_t/(pow(s.C_l,3));
    }
    else {
        // Murray of freeflow, nothing to be done
        return;
    }
}

// Return flow rate on the opening in LBM units of the predefined profile
T OpeningHandler::getProfileVelSum()
{
    T velSum = 0.0;

    for(auto const& v: nodes)
        velSum += velArr[v.x][v.y][v.z].norm();

    return velSum;
}

T OpeningHandler::getFlowRate(SimPar s)
{
    T velSum = 0.0;

    for (auto const &v : nodes)
        velSum += velArr[v.x][v.y][v.z].dot(direction);

    return velSum * bcParameter * cScale * s.C_l * s.C_l * s.C_l / s.C_t; //*getArea()/getArea()
}

// Get the current flowrate on a defined velocity boundary
T OpeningHandler::getScaledFlowRate()
{
    return bcParameter * 0.5 * getArea() * cScale;
}

// Scale the flow velocity array
void OpeningHandler::scaleFlowRate(T scale_)
{
    for(auto const& v: nodes)
        velArr[v.x][v.y][v.z] = velArr[v.x][v.y][v.z] * scale_;
}

// Scale the pressure array
void OpeningHandler::scalePressure(T scale_)
{
    for(auto const& v: nodes)
        presArr[v.x][v.y][v.z] = presArr[v.x][v.y][v.z] * scale_;
}

void OpeningHandler::normalizeFlowRate()
{
  T invVelSum = 1.0 / getProfileVelSum();

    for(auto const& v: nodes)
        velArr[v.x][v.y][v.z] = velArr[v.x][v.y][v.z] * invVelSum;
}

// Note the peak of the profile is v_norm == 1
void OpeningHandler::createPoiseauilleProfile()
{
    pcout << "-> Creating direction-corrected Pouseuille velocity profile on label: " << flag << std::endl;

    // Sanity check
    if ( !(boundingBox->x0 == boundingBox->x1 || boundingBox->y0 == boundingBox->y1 || boundingBox->z0 == boundingBox->z1) )
        pcout << "!!! ERROR: The opening is not parallel to any principal plane! This functionality is not implemented, the opening will not work properly!" << std::endl;
    if (type != OPENING_VELOCITY && type != OPENING_MURRAY && type != OUTLET_FREEFLOW)
        pcout << "WARNING! Setting velocity profile for a non-velocity opening! This will have no effect." << std::endl;

    // Paraboloid height
    T u_max = 1.0;

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

void OpeningHandler::createBluntVelocityProfile()
{
    pcout << "-> Creating blunt velocity profile on flag: " << flag << std::endl;

    // Sanity check
    if ( !(boundingBox->x0 == boundingBox->x1 || boundingBox->y0 == boundingBox->y1 || boundingBox->z0 == boundingBox->z1) )
        pcout << "!!! ERROR: The opening is not parallel to any principal plane! This functionality is not implemented, the opening will not work properly!" << std::endl;
    if (type != OPENING_VELOCITY && type != OPENING_MURRAY)
        pcout << "WARNING! Setting velocity profile for a non-velocity opening! This will have no effect." << std::endl;

    vec3d vel = direction.getNormal();

    for(auto const& v: nodes) {
        velArr[v.x][v.y][v.z].set(vel.x, vel.y, vel.z);
    }
}

void OpeningHandler::createConstantPressureProfile(T density)
{
    pcout << "-> Creating constant pressure on flag: " << flag << std::endl;

    if(type != OPENING_PRESSURE )
        pcout << "WARNING! Setting pressure profile for a non-pressure opening! This will have no effect." << std::endl;

    for(auto const& v: nodes) {
        presArr[v.x][v.y][v.z] = density;
    }
}

void OpeningHandler::loadScaleFunction(const string& fileName)
{
    if(type==OPENING_MURRAY || type==OUTLET_FREEFLOW){
        pcout << "-> *WARNING*: The opening type can't have a scale function. The scale function will be disregarded." << std::endl;
        hasScaleFunction = false;
        return;
    }

    pcout << "-> Loading and scale function: " << fileName << std::endl;

    plb_ifstream finSign(fileName.c_str());
    //istream &finSign = pfinSign.getOriginalStream();

    if(!finSign.good()) {
        pcout << "WARNING!!! Flow rate scale file " << fileName << " is not readable!" << std::endl;
        // hasScaleFunction = false;
    }

    T time, value;
    int Ns;

    finSign >> Ns;
    global::mpi().bCast(&Ns, 1);
    global::mpi().barrier();

    for(int i=0; i<Ns; i++) {
        finSign >> time; finSign >> value;
        global::mpi().bCast(&time, 1);
        global::mpi().bCast(&value, 1);
        global::mpi().barrier();
        scaleTime.push_back(time); scaleSignal.push_back(value);
    }

    finSign.close();

    hasScaleFunction = true;
    pcout << fileName << " loaded with " << scaleTime.size() << " data points." << std::endl;
}

void OpeningHandler::setBCType(MultiBlockLattice3D<T, DESCRIPTOR> *lattice)
{
    if (type == OPENING_VELOCITY || type == OPENING_MURRAY ) {
        OnLatticeBoundaryCondition3D<T, DESCRIPTOR> *bc = createLocalBoundaryCondition3D<T,DESCRIPTOR>();
        // OnLatticeBoundaryCondition3D<T, DESCRIPTOR> *bc = createZouHeBoundaryCondition3D<T,DESCRIPTOR>();
        bc->setVelocityConditionOnBlockBoundaries(*lattice, *boundingBox, boundary::dirichlet);
    }

    if(type == OPENING_PRESSURE) {
        // OnLatticeBoundaryCondition3D<T, DESCRIPTOR> *bc = createLocalBoundaryCondition3D<T,DESCRIPTOR>();
        // OnLatticeBoundaryCondition3D<T, DESCRIPTOR> *bc = createZouHeBoundaryCondition3D<T,DESCRIPTOR>();
        OnLatticeBoundaryCondition3D<T, DESCRIPTOR> *bc = createInterpBoundaryCondition3D<T,DESCRIPTOR>();
        bc->setPressureConditionOnBlockBoundaries(*lattice, *boundingBox, boundary::dirichlet);
    }

    if(type == OUTLET_FREEFLOW) {
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
    }
}

void OpeningHandler::setExternalVelocityProfile(MultiBlockLattice3D<T, DESCRIPTOR> *lattice, field3D &velocityArr)
{
    
    if(type==OPENING_VELOCITY || type==OPENING_MURRAY) {
        // Set the predetermined profile
        setBoundaryVelocity(*lattice, *boundingBox, VelocityProfile3D<T,DESCRIPTOR>(&velocityArr, 1.0));
    }
    else {
        pcout << "WARNING, opening " << flag << " is incorrectly addressed as velocity opening, while it is of type: " << type << std::endl;
    }
}

void OpeningHandler::setExternalPressureProfile(MultiBlockLattice3D<T, DESCRIPTOR> *lattice, scalar3D &pressureArr)
{
    if(type==OPENING_PRESSURE) {
        // Set predetermined pressure profile
        setBoundaryDensity(*lattice, *boundingBox, PressureProfile3D<T,DESCRIPTOR>(&pressureArr, 1.0));
    }
    else {
        pcout << "WARNING, opening " << flag << " incorrectly addressed as pressure opening, while it is of type:" << type << std::endl;
    }
}

void OpeningHandler::setScaledBoundaryProfile(MultiBlockLattice3D<T, DESCRIPTOR> *lattice, T scale = 1.0)
{   
    if(type == OPENING_VELOCITY || type == OPENING_MURRAY) 
        setBoundaryVelocity(*lattice, *boundingBox, VelocityProfile3D<T,DESCRIPTOR>(&velArr, scale));
    else if(type == OPENING_PRESSURE)
        setBoundaryDensity(*lattice, *boundingBox, PressureProfile3D<T,DESCRIPTOR>(&presArr, scale));
    else if(type == OUTLET_FREEFLOW)
        setBoundaryVelocity(*lattice, *boundingBox, VelocityProfile3D<T,DESCRIPTOR>(&velArr, scale));

    // TODO: errorhandling 'else'-case?

}

void OpeningHandler::progressTime(MultiBlockLattice3D<T, DESCRIPTOR> *lattice, T dt)
{
    // Do we need to change the scale value?
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
            cScale = interpolate(scaleTime[cTimePos], scaleTime[cTimePos+1], cTimeVal, scaleSignal[cTimePos], scaleSignal[cTimePos+1]);
        else
            cScale = interpolate(scaleTime[cTimePos], scaleTime[0], cTimeVal, scaleSignal[cTimePos], scaleSignal[0]);
    }

    // Apply the previously defined profile scaled with 'parameter' and the scale function if it exists.
    setScaledBoundaryProfile(lattice, cScale * bcParameter);

}

OpeningHandler::~OpeningHandler()
{

}