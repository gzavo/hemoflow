#### New fetures:

##### Pre-processing:

1. Added parallel processing to the voxelize() function.

2. Added the connected-component analysis method to detect the opening voxels.
   
   Aadvntage: very fast, execution time does not scale with problem size.
   
   Drawbacks: has problem detecting opening voxels on the boundary edge.

3. Added the feature to detect multiple inlets. However, hardcoded. One can change the following lines according to the specific geometry.
   
   Line 61 in `readCL.py` is used to define the inlets so that the direction of the inflow points inside.
   
   Line 216 ~ 249 in `detectOpenings.py` are used to paint the openings with different labels.

4. Seperated the opening information from the geometry flag. Now the opening information is stored in a npz file and the geometry flag is stored in a compressed HDF5 file.

##### HemoFlow:

1. Changed the XML structure, the details can be found in the `ExampleXML.xml`. 

2. Added the support of multiple inlets geometry. However, still hardcoded. Change the following lines according to the geometry.
   
   Line 247 ~ 268 in `hemoFlow.cpp` are used to define the boundary velocity/pressure profile and load the corresponding flowrate function.
   
   Line 731 ~ 740 in `hemoFlow.cpp` are used to read the flowrate functions.
   
   Line 52, 54, 198, 216, 228, 264 in `opening.cpp` are used to determine the pressure outlet.

3. Optimized the RAM utilization by MPI shared memory.

4. Added the Smagorinsky model to implement LES.

5. Improve the Output efficiency significantly by utilizing HDF5 and Xdmf3

#### TODO

##### Pre-processing:

1. Improve the efficiency of pre-processor, possibly needs a re-write in C or C++.

2. Improve the applicability of opening detection for multiple inlets and outlets, preferably with a graphicl interface.

3. Find a solution to deal with non-perpendicular openings.

##### HemoFlow:

1. Improve the applicability of opening detection, so that no hardcoded is needed.

2. Only read necessary geometry flag voxels to each MPI thread to prevent RAM waste.

3. Remove Xdmf3 since it is not actively maintained. Utilize the VTK HDF Reader ([VTK HDF Reader (kitware.com)](https://www.kitware.com/vtk-hdf-reader/)).
