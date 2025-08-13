import pyvista as pv
import os
import glob
import numpy as np
import argparse

def mag(a: pv.pyvista_ndarray):
    "Returns the magnigude of an array of scalars/vectors."
    return np.sqrt((a * a).sum(axis=-1))
def calculate_mean(slice):
    slice['mean_density'] = slice['density [Pa]'] / slice['Area']
    slice['mean_velocity'] = slice['velocity_magnitude']/ slice['Area']
    slice['d']=np.sqrt(slice['Area']/np.pi*4)
    slice['VFR']=slice['mean_velocity']*slice['Area']

def main(path):
    timesteps = glob.glob(os.path.join(path,"output/**/*.vti"), recursive=True)
    timesteps=sorted(timesteps)
    time=0

    for timestep in timesteps:
        pointdata=pv.read(timestep)
        dx=pointdata.spacing[0]
        velocity_magnitude = mag(pointdata.point_data["velocity [m/s]"])
        pointdata.point_data["velocity_magnitude"] = velocity_magnitude
        voxeldata=pointdata.point_data_to_cell_data()

        threshold=voxeldata.threshold(value=0.001,scalars='velocity_magnitude')

        inlet_slice=threshold.slice(origin=[voxeldata.bounds[0]+dx,0,0], normal=[1,0,0])
        murray_outlet_slice=threshold.slice(origin=[voxeldata.bounds[1]-dx,0,0], normal=[1,0,0])
        small_outlet_slice=threshold.slice(origin=[0,0,voxeldata.bounds[5]-dx], normal=[0,0,1])
        small_outlet_slice['perp_velocity']=np.dot(small_outlet_slice['velocity [m/s]'],[0,0,1])
        inlet=inlet_slice.integrate_data()
        murray_outlet=murray_outlet_slice.integrate_data()
        small_outlet=small_outlet_slice.integrate_data()
        calculate_mean(inlet)
        calculate_mean(murray_outlet)
        calculate_mean(small_outlet)

        sum_vfr=inlet['mean_velocity']*inlet['Area']
        sum_d3=murray_outlet['d']**3+small_outlet['d']**3

        murray_outlet_theory=murray_outlet['d']**3*sum_vfr/sum_d3
        small_outlet_theory=small_outlet['d']**3*sum_vfr/sum_d3

        murray_diff=((murray_outlet['VFR']-murray_outlet_theory)/murray_outlet_theory)[0]
        small_diff=(((small_outlet['perp_velocity'])-small_outlet_theory)/small_outlet_theory)[0]
        print(f"Time:{time}")
        print(f"Inlet VFR:{inlet['VFR'][0]:.4}, Sum outlet VFR:{(murray_outlet['VFR'][0]+small_outlet['VFR'][0]):.4}")
        print(f"Murray outlet AVG velocity:{murray_outlet['mean_velocity'][0]:.4}, VFR: {murray_outlet['VFR'][0]:.4},Murray difference:{murray_diff:.2%}")
        print(f"Small outlet AVG velocity:{small_outlet['mean_velocity'][0]:.4}, VFR: {small_outlet['VFR'][0]:.4}, Murray difference:{small_diff:.2%}")
        
        time+=1

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("d", default="./", help="Path to directory", type=str)

    args = parser.parse_args()

    main(path=os.path.abspath(args.d))