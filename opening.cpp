#include "opening.h"


OpeningHandler::OpeningHandler(unsigned short* flagAray, int flag_, T radius_lb, vec3d dirVec)
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

    pcout << "-> Getting bounding box, reduced radius, center... " << std::endl;
    boundingBox = new Box3D(*min_element(xCoord.begin(), xCoord.end()),
                            *max_element(xCoord.begin(), xCoord.end()),
                            *min_element(yCoord.begin(), yCoord.end()),
                            *max_element(yCoord.begin(), yCoord.end()),
                            *min_element(zCoord.begin(), zCoord.end()),
                            *max_element(zCoord.begin(), zCoord.end()) );
    

    R = radius_lb; //sqrt(nodes.size() / M_PI);

    center.x = accumulate( xCoord.begin(), xCoord.end(), 0.0) / xCoord.size();
    center.y = accumulate( yCoord.begin(), yCoord.end(), 0.0) / yCoord.size();
    center.z = accumulate( zCoord.begin(), zCoord.end(), 0.0) / zCoord.size();

    direction.x = dirVec.x; direction.y = dirVec.y; direction.z = dirVec.z;
}

void OpeningHandler::printOpeningDetails()
{
    pcout << "Opening parameters:" << std::endl;
    
    if(flag == INLET)
        pcout << "=> This is the inlet." << std::endl;
    if(flag == FIRST_OUTLET)
        pcout << "=> This is the smallest pressure outlet." << std::endl;

    pcout << "-> flag: " << flag << std::endl;
    pcout << "-> radius [lb]: " << R << std::endl;
    //pcout << "-> Q rate: " << Qr << std::endl;
    pcout << "-> Center [lb]: " << center.x << " " << center.y << " " << center.z << std::endl;
    pcout << "-> Normal: " << direction.x << " " << direction.y << " " << direction.z << std::endl;
    pcout << "-> Area [lb]: " << nodes.size() << std::endl;
    pcout << "-> Scale function: " << hasScaleFunction << " length: " << scaleSignal.size() << std::endl;

}

void OpeningHandler::createPoiseauilleProfile(T u_avg)
{
    // TODO: Implement.
}

void OpeningHandler::createBluntVelocityProfile(T u_avg)
{
    pcout << "-> Creating blunt velocity profile on flag: " << flag << std::endl;

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
    pcout << "-> Creating constant pressure on flag: " << flag << std::endl;

    if(flag < FIRST_OUTLET) {
        pcout << "WARNING! Calculating pressure profile for a non-pressure opening!" << std::endl;
    }

    for(auto const& v: nodes) {
        presArr[v.x][v.y][v.z] = density;
    }
}

void OpeningHandler::loadScaleFunction(string fileName)
{
    pcout << "-> Loading scale function: " << fileName << std::endl;

    plb_ifstream finSign(fileName.c_str());
    //istream &finSign = pfinSign.getOriginalStream();

    if(!finSign.good()) {
        pcout << "WARNING!!! Scale file " << fileName << " is not readable!" << std::endl;
        hasScaleFunction = false;
        return;
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


void OpeningHandler::setBC(MultiBlockLattice3D<T, DESCRIPTOR> *lattice, OnLatticeBoundaryCondition3D<T, DESCRIPTOR> *bc)
{
    if (flag == FIRST_OUTLET) {
        bc->setPressureConditionOnBlockBoundaries(*lattice, *boundingBox, boundary::dirichlet);
    }
    else {
        bc->setVelocityConditionOnBlockBoundaries(*lattice, *boundingBox, boundary::dirichlet);
    }
}

void OpeningHandler::setVelocityProfile(MultiBlockLattice3D<T, DESCRIPTOR> *lattice, OnLatticeBoundaryCondition3D<T, DESCRIPTOR> *bc, field3D &velocityArr)
{
    // For flowrate, additional measurements required after lattice initialization
    if(flag==INLET || flag > FIRST_OUTLET) {
        // Set given profile
        setBoundaryVelocity(*lattice, *boundingBox, VelocityProfile3D<T,DESCRIPTOR>(&velocityArr, 1.0));
    }
    else {
        pcout << "WARNING, opening " << flag << " incorrectly addressed as velocity opening (instead of pressure)!" << std::endl;
    }
}

void OpeningHandler::setPressureProfile(MultiBlockLattice3D<T, DESCRIPTOR> *lattice, OnLatticeBoundaryCondition3D<T, DESCRIPTOR> *bc, scalar3D &pressureArr)
{
    // For flowrate, additional measurements required after lattice initialization
    if(flag==FIRST_OUTLET) {
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

    if (flag != FIRST_OUTLET)
        setBoundaryVelocity(*lattice, *boundingBox, VelocityProfile3D<T,DESCRIPTOR>(&velArr, scale));
    else {
        setBoundaryDensity(*lattice, *boundingBox, 1.0);
        // setBoundaryDensity(*lattice, *boundingBox, PressureProfile3D<T,DESCRIPTOR>(&presArr, scale)); // For time dependent pressure boundary
    }
}

OpeningHandler::~OpeningHandler()
{

}