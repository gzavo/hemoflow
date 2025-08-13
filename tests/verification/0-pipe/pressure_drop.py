import pyvista as pv
import os
import glob
import numpy as np
import argparse
import re

def mag(a: pv.pyvista_ndarray):
    "Returns the magnigude of an array of scalars/vectors."
    return np.sqrt((a * a).sum(axis=-1))
def calculate_mean(slice):
    slice['mean_density'] = slice['density [Pa]'] / slice['Area']
    slice['mean_velocity'] = slice['velocity_magnitude']/ slice['Area']
    slice['VFR']=slice['mean_velocity']*slice['Area']

def main(path):

    timesteps = glob.glob(os.path.join(path,"output/**/*.vti"), recursive=True)
    timesteps=sorted(timesteps)
    pointdata=pv.read(timesteps[-1])
    dx=pointdata.spacing[0]
    velocity_magnitude = mag(pointdata.point_data["velocity [m/s]"])
    pointdata.point_data["velocity_magnitude"] = velocity_magnitude
    voxeldata=pointdata.point_data_to_cell_data()
    threshold=voxeldata.threshold(value=0.0001,scalars='velocity_magnitude')
    inlet_slice=threshold.slice(origin=[voxeldata.bounds[0]+dx,0,0], normal=[1,0,0])
    outlet_slice=threshold.slice(origin=[voxeldata.bounds[1]-dx,0,0], normal=[1,0,0])
    inlet=inlet_slice.integrate_data()
    outlet=outlet_slice.integrate_data()
    calculate_mean(inlet)
    calculate_mean(outlet)
    pressure_drop=(inlet['mean_density']-outlet['mean_density'])[0]

    print(f"Pressure drop for pipe {pressure_drop}")
    print(f"Inlet VFR: {inlet['VFR'][0]:.4} m3/s")
    print(f"Outlet VFR: {outlet['VFR'][0]:.4} m3/s")
    


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("d", default="./", help="Path to directory", type=str)

    args = parser.parse_args()

    main(path=os.path.abspath(args.d))

