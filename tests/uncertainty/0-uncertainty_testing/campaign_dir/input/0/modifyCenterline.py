import os
import vmtk
from vmtk import pypes
import sys
from pathlib import Path
import glob

# YOU NEED VMTK FOR TIHS TO WORK

def cl_geom(path):
    surface = Path(path).resolve().as_posix()
    myArgs = (
        "vmtkcenterlines -ifile "
        + str(surface)
        + " -endpoints 1 --pipe vmtkcenterlineresampling -length 0.15 --pipe vmtkcenterlinegeometry -smoothing 1 -factor 1 -iterations 75 -outputsmoothed 0 -ofile "
        + str(surface).replace("_Srf_prep.stl", "_Ctl_prep.vtp")
    )
    pypes.PypeRun(myArgs)


def remesh_surf_cap_then_ctl(path):
    surface = Path(path).resolve().as_posix()
    myArgs = (
        "vmtksurfaceremeshing -ifile "
        + str(surface)
        + " -area 0.025 -iterations 5 --pipe vmtksurfacecapper -triangle 1 -method centerpoint -ofile "
        + str(surface).replace("_Srf_prep.stl", "_Srf_prep_remeshed.stl")
        + " --pipe vmtkcenterlines -endpoints 1 --pipe vmtkcenterlineresampling -length 0.15 --pipe vmtkcenterlinegeometry -smoothing 1 -factor 1 -iterations 75 -outputsmoothed 0 -ofile "
        + str(surface).replace("_Srf_prep.stl", "_Ctl_prep.vtp")
    )
    pypes.PypeRun(myArgs)


def calculate_curvature_torsion(path):
    surface = Path(path).resolve().as_posix()
    myArgs = (
        "vmtkcenterlinegeometry -ifile "
        + str(surface)
        + " -smoothing 1 -factor 1 -iterations 75 -outputsmoothed 0 -ofile "
        + str(surface)
    )
    pypes.PypeRun(myArgs)


def main(dir, ext):
    files = glob.glob(
        os.path.join(dir, "**", ext),
        recursive=True,
    )
    print(files)
    for file in files:
        calculate_curvature_torsion(file)


# cl_geom(sys.argv[1])
# remesh_surf_cap_then_ctl(sys.argv[1])
calculate_curvature_torsion(sys.argv[1])
# main(sys.argv[1], sys.argv[2])
