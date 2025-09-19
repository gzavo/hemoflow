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

def make_box(origin, normal, size):
    # normal defines orientation, but Box requires bounds
    # so we make a local coordinate system aligned with the normal
    direction = np.array(normal) / np.linalg.norm(normal)
    # construct arbitrary perpendicular vectors for local axes
    if np.allclose(direction, [1,0,0]):
        perp1 = np.array([0,1,0])
    else:
        perp1 = np.cross(direction, [1,0,0])
        perp1 /= np.linalg.norm(perp1)
    perp2 = np.cross(direction, perp1)

    # compute corners of a cube around origin
    offsets = [
        direction*size, -direction*size,
        perp1*size, -perp1*size,
        perp2*size, -perp2*size,
    ]
    # get bounds
    points = [origin + o for o in offsets]
    mins = np.min(points, axis=0)
    maxs = np.max(points, axis=0)
    return [mins[0], maxs[0], mins[1], maxs[1], mins[2], maxs[2]]


def main(path):
    timesteps = glob.glob(os.path.join(path,"output/**/*.vti"), recursive=True)
    timesteps=sorted(timesteps)
    time=0
    
    compFile = np.load('./input/vox_murray_5M_c.npz')
    dx=compFile['dx']
    openingCenter=compFile['openingCenter']
    openingRadius=compFile['openingRadius']
    openingNormal=compFile['openingNormal']
    openingNormalizedQRatio=compFile['openingNormalizedQRatio']
    
    inlet_radius=openingRadius[0]
    murray_outlet1_radius=openingRadius[1]
    murray_outlet2_radius=openingRadius[2]
    pressure_outlet_radius=openingRadius[3]
    
    inlet_normal=openingNormal[0]
    murray_outlet1_normal=openingNormal[1]
    murray_outlet2_normal=openingNormal[2]
    pressure_outlet_normal=openingNormal[3]

    inlet_center=openingCenter[0]*dx[0]
    murray_outlet1_center=openingCenter[1]*dx[0]
    murray_outlet2_center=openingCenter[2]*dx[0]
    pressure_outlet_center=openingCenter[3]*dx[0]
    
    inlet_center[2] += 8*dx[0]
    murray_outlet1_center[0] -= 8*dx[0]
    murray_outlet2_center[2] -= 8*dx[0]
    pressure_outlet_center[1] += 8*dx[0]


    for timestep in timesteps:
        pointdata=pv.read(timestep)
        #dx=pointdata.spacing[0]
        velocity_magnitude = mag(pointdata.point_data["velocity [m/s]"])
        pointdata.point_data["velocity_magnitude"] = velocity_magnitude
        voxeldata=pointdata.point_data_to_cell_data()

        threshold=voxeldata.threshold(value=0.001,scalars='velocity_magnitude')
        
        inlet_box=make_box(inlet_center, inlet_normal, inlet_radius*3)
        inlet_local = threshold.clip_box(inlet_box, invert=False)
        murray1_box=make_box(murray_outlet1_center, murray_outlet1_normal, murray_outlet1_radius*3)
        murray1_local = threshold.clip_box(murray1_box, invert=False)
        murray2_box=make_box(murray_outlet2_center, murray_outlet2_normal, murray_outlet2_radius*3)
        murray2_local = threshold.clip_box(murray2_box, invert=False)
        pressure_box=make_box(pressure_outlet_center, pressure_outlet_normal, pressure_outlet_radius*3)
        pressure_local = threshold.clip_box(pressure_box, invert=False)
        
        inlet_slice=inlet_local.slice(origin=inlet_center, normal=inlet_normal)
        murray_outlet1_slice=murray1_local.slice(origin=murray_outlet1_center, normal=murray_outlet1_normal)
        murray_outlet2_slice=murray2_local.slice(origin=murray_outlet2_center, normal=murray_outlet2_normal)
        pressure_outlet_slice=pressure_local.slice(origin=pressure_outlet_center, normal=pressure_outlet_normal)
              
        inlet_slice['perp_velocity']=np.dot(inlet_slice['velocity [m/s]'], inlet_normal)
        murray_outlet1_slice['perp_velocity']=np.dot(murray_outlet1_slice['velocity [m/s]'], murray_outlet1_normal)
        murray_outlet2_slice['perp_velocity']=np.dot(murray_outlet2_slice['velocity [m/s]'], murray_outlet2_normal)
        pressure_outlet_slice['perp_velocity']=np.dot(pressure_outlet_slice['velocity [m/s]'], pressure_outlet_normal)
        
        inlet=inlet_slice.integrate_data()
        murray_outlet1=murray_outlet1_slice.integrate_data()
        murray_outlet2=murray_outlet2_slice.integrate_data()
        pressure_outlet=pressure_outlet_slice.integrate_data()
        #calculate_mean(inlet)
        #calculate_mean(murray_outlet1)
        #calculate_mean(murray_outlet2)
        #calculate_mean(pressure_outlet)
        
        sum_vfr=inlet['perp_velocity'][0]
        sum_outvfr=-1*(murray_outlet1['perp_velocity'][0]+murray_outlet2['perp_velocity'][0]+pressure_outlet['perp_velocity'][0])
        sum_d3=(2*murray_outlet1_radius)**3+(2*murray_outlet2_radius)**3+(2*pressure_outlet_radius)**3

        murray_outlet1_theory=((2*murray_outlet1_radius)**3)*sum_vfr/sum_d3
        murray_outlet2_theory=(2*murray_outlet2_radius)**3*sum_vfr/sum_d3
        pressure_outlet_theory=(2*pressure_outlet_radius)**3*sum_vfr/sum_d3

        murray1_diff=((-1*(murray_outlet1['perp_velocity'])-murray_outlet1_theory)/murray_outlet1_theory)
        murray2_diff=((-1*(murray_outlet2['perp_velocity'])-murray_outlet2_theory)/murray_outlet2_theory)
        pressure_diff=((-1*(pressure_outlet['perp_velocity'])-pressure_outlet_theory)/pressure_outlet_theory)
        
        print(f"Time:{time}")
        print(f"Inlet VFR:{inlet['perp_velocity'][0]:.4}, Sum outlet VFR:{sum_outvfr:.4}, Diff. VFR %:{(sum_vfr-sum_outvfr)/sum_vfr:.2%}")
        print(f"Murray outlet1 VFR: {-1*(murray_outlet1['perp_velocity'][0]):.4}, Murray1 difference:{murray1_diff[0]:.2%}")
        print(f"Murray outlet2 VFR: {-1*(murray_outlet2['perp_velocity'][0]):.4}, Murray2 difference:{murray2_diff[0]:.2%}")
        print(f"Pressure outlet1 VFR: {-1*(pressure_outlet['perp_velocity'][0]):.4}, Pressure difference:{pressure_diff[0]:.2%}")
        #print(f"Small outlet AVG velocity:{small_outlet['mean_velocity'][0]:.4}, VFR: {small_outlet['VFR'][0]:.4}, Murray difference:{small_diff:.2%}")
        
        time+=1

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("d", default="./", help="Path to directory", type=str)

    args = parser.parse_args()

    main(path=os.path.abspath(args.d))