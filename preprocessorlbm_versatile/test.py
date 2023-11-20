import numpy as np
import h5py
import time

wrong = np.load("/home/kwang/VascuTreatCFD/data/fontan_mesh/1e7/Fontan_1e7_ops.npz")

wrong = dict(wrong)

print(wrong["openingTangent"])

# print(" ")

# wrong["openingTangent"][2][:] = -wrong["openingTangent"][2][:]
# wrong["openingTangent"][3][:] = -wrong["openingTangent"][3][:]
# wrong["openingTangent"][4][:] = -wrong["openingTangent"][4][:]
# wrong["openingTangent"][6][:] = -wrong["openingTangent"][6][:]

# np.savez("/home/kwang/VascuTreatCFD/data/fontan_mesh/1e7/Fontan_1e7_ops.npz", **wrong)

# print(wrong["openingTangent"])

print(" ")
correct = np.load("/home/kwang/VascuTreatCFD/data/fontan_mesh/1e10/Fontan_1e10_ops.npz")

print(correct["openingTangent"])