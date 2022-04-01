import numpy as np
from utils import formatToPrint, readFile, computeL1Error, computeL2Error, computeLInfError

# parameters
Lt = 300.
Lx = 50.
Ly = 50.

def computeErrors(uRef, uVector, norm = 1):
    global Lt, Lx, Ly

    ntRef = uRef.shape[0] - 1
    ntVector = [u.shape[0] - 1 for u in uVector]
    stepVector = [ntRef // nt for nt in ntVector]
    DeltaVector = [uRef[::step, :, :] - u for step, u in zip(stepVector, uVector)]

    if norm==1:
        return [computeL1Error(Delta) for Delta in DeltaVector]
    elif norm==2:
        L2Errors = []
        for Delta in DeltaVector:
            nt, nx, ny = Delta.shape
            dt = Lt / (nt - 1)
            dx = Lx / (nx - 1)
            dy = Ly / (ny - 1)
            L2Errors.append(computeL2Error(Delta, dt, dx, dy))
        return L2Errors
    return [computeLInfError(Delta) for Delta in DeltaVector]

# reference u
ntRef = 2400 * 2 ** 5
nxRef = 50; nyRef = 50
filenameRef = "../src/sol_{}_{}_{}.txt".format(nxRef, nyRef, ntRef)
uRef = readFile(filenameRef)

# computed u 
ntRange = [2400 * 2 ** i for i in range(4)]
nx = 50; ny = 50
filenames = ["../src/sol_{}_{}_{}.txt".format(nx, ny, nt) for nt in ntRange]
uVector = [readFile(filename) for filename in filenames]

# L1 errors
print("== L1 ERRORS ==")
L1Errors = computeErrors(uRef, uVector, norm = 1)
print(formatToPrint(L1Errors))
estimatedL1Rate = [np.log2(L1Errors[i] / L1Errors[i + 1]) for i in range(len(L1Errors) - 1)]
print(formatToPrint(estimatedL1Rate))

# L2 errors
print("\n== L2 ERRORS ==")
L2Errors = computeErrors(uRef, uVector, norm = 2)
print(formatToPrint(L2Errors))
estimatedL2Rate = [np.log2(L2Errors[i] / L2Errors[i + 1]) for i in range(len(L2Errors) - 1)]
print(formatToPrint(estimatedL2Rate))

# LInf errors
print("\n== LInf ERRORS ==")
LInfErrors = computeErrors(uRef, uVector, norm = "inf")
print(formatToPrint(LInfErrors))
estimatedLInfRate = [np.log2(LInfErrors[i] / LInfErrors[i + 1]) for i in range(len(LInfErrors) - 1)]
print(formatToPrint(estimatedLInfRate))