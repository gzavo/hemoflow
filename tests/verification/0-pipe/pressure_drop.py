import pyvista as pv
import os
import glob
import numpy as np
import argparse
import matplotlib.pyplot as plt

rho_blood=1055
nu_blood=3.33e-6

def mag(a: pv.pyvista_ndarray):
    "Returns the magnigude of an array of scalars/vectors."
    return np.sqrt((a * a).sum(axis=-1))


def calculate_mean(slice):
    slice["mean_density"] = slice["density [Pa]"] / slice["Area"]
    slice["mean_velocity"] = slice["velocity_magnitude"] / slice["Area"]
    slice["VFR"] = slice["mean_velocity"] * slice["Area"]


def poiseuille_1d(y, D=1.0, mu=1.0, dp=1.0, L=1.0):
    """
    Analytical 1D Poiseuille velocity profile between parallel plates.

    Parameters
    ----------
    y : float or ndarray
        Vertical coordinate(s), 0 <= y <= H.
    H : float
        Channel height.
    mu : float
        Dynamic viscosity.
    G : float
        Pressure gradient magnitude (dp/dx = -G).

    Returns
    -------
    u : float or ndarray
        Velocity at position y.
    """
    G=dp/L
    u=(G / (4.0 * mu)) * y * (D - y)
    return u


def main(path):

    timesteps = glob.glob(os.path.join(path, "output/**/*.vti"), recursive=True)
    timesteps = sorted(timesteps)
    pointdata = pv.read(timesteps[-1])
    print(timesteps[-1])
    dx = pointdata.spacing[0]
    velocity_magnitude = mag(pointdata.point_data["velocity [m/s]"])
    pointdata.point_data["velocity_magnitude"] = velocity_magnitude
    voxeldata = pointdata.point_data_to_cell_data()
    threshold = voxeldata.threshold(value=0.0001, scalars="velocity_magnitude")
    x_bounds = voxeldata.bounds[1] - voxeldata.bounds[0]
    y_bounds = voxeldata.bounds[3] - voxeldata.bounds[2]
    inlet_slice = threshold.slice(origin=[x_bounds / 4, 0, 0], normal=[1, 0, 0])
    outlet_slice = threshold.slice(origin=[3 / 4 * x_bounds, 0, 0], normal=[1, 0, 0])
    inlet = inlet_slice.integrate_data()
    outlet = outlet_slice.integrate_data()
    calculate_mean(inlet)
    calculate_mean(outlet)
    pressure_drop = (inlet["mean_density"] - outlet["mean_density"])[0]

    
    rho_blood=1055
    nu_blood=3.22e-6
    L_dp=x_bounds/2
    D_pipe=y_bounds

    dp_Poiseuille=(8*nu_blood*rho_blood*L_dp*outlet['VFR'][0])/(3.14*(D_pipe/2)**4)
    dp_diff=abs(dp_Poiseuille-pressure_drop)/dp_Poiseuille

    velocity_profile=threshold.sample_over_line(pointa=[3 / 4 * x_bounds,y_bounds/2,0],pointb=[3 / 4 * x_bounds,y_bounds/2,y_bounds])
    analyitical_profile=poiseuille_1d(velocity_profile.points[:,2],D_pipe,mu=nu_blood*rho_blood,dp=dp_Poiseuille, L=L_dp)
    fig, ax = plt.subplots()
    ax.plot(velocity_profile.points[:,2], velocity_profile["velocity_magnitude"])
    ax.plot(velocity_profile.points[:,2], analyitical_profile)
    ax.set(xlabel='Z coordinate', ylabel='Velocity magnitude (m/s)')
    ax.grid()

    fig.savefig("profile_comp.png")


    print(f"Pipe diameter {D_pipe}")
    print(f"Distance between measurement {L_dp}")
    print(f"Pressure drop from Poiseuille equation (without stent) {dp_Poiseuille}")
    print(f"Pressure drop for pipe {pressure_drop}")
    print(f"Relative difference between analyitcal and numerical {dp_diff*100:.4}%")
    print(f"Inlet VFR: {inlet['VFR'][0]:.4} m3/s")
    print(f"Outlet VFR: {outlet['VFR'][0]:.4} m3/s")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("d", default="./", help="Path to directory", type=str)

    args = parser.parse_args()

    main(path=os.path.abspath(args.d))
