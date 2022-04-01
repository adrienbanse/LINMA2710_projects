import numpy as np
from utils import formatToPrint, readFile, computeL1Error, computeL2Error, computeLInfError

# parameters
Lt = 300.
Lx = 50.
Ly = 50.

def computeErrors(uRef, uVector, x = False, norm = 1):
    global Lt, Ly, Lx

    if x:
        nxRef = uRef.shape[1] - 1
        nxVector = [u.shape[1] - 1 for u in uVector]
        stepVector = [nxRef // nx for nx in nxVector]
        DeltaVector = [u - uRef[:, ::step, :] for u, step in zip(uVector, stepVector)]
    else:
        nyRef = uRef.shape[2] - 1
        nyVector = [u.shape[2] - 1 for u in uVector]
        stepVector = [nyRef // ny for ny in nyVector]
        DeltaVector = [uRef[:, :, ::step] - u for u, step in zip(uVector, stepVector)]

    if norm == 1:
        return [computeL1Error(Delta) for Delta in DeltaVector]
    elif norm == 2:
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
ntRef = 2400 * 2 ** 3
nxRef = 8 * 2 ** 4; nyRef = 8 * 2 ** 4
filenameRef = "../src/sol_{}_{}_{}.txt".format(nxRef, nyRef, ntRef)
uRef = readFile(filenameRef)

##################
### X ANALYSIS ###
##################

# computed u
nxRange = [8 * 2 ** i for i in range(4)]
nt = ntRef; ny = 128
filenames = ["../src/sol_{}_{}_{}.txt".format(nx, ny, nt) for nx in nxRange]
uVector = [readFile(filename) for filename in filenames]

print("== L1 ERRORS for x ==")
L1Errors = computeErrors(uRef, uVector, x = True, norm = 1)
print(formatToPrint(L1Errors))
estimatedL1Rate = [np.log2(L1Errors[i] / L1Errors[i + 1]) for i in range(len(L1Errors) - 1)]
print(formatToPrint(estimatedL1Rate))

print("\n== L2 ERRORS for x ==")
L2Errors = computeErrors(uRef, uVector, x = True, norm = 2)
print(formatToPrint(L2Errors))
estimatedL2Rate = [np.log2(L2Errors[i] / L2Errors[i + 1]) for i in range(len(L2Errors) - 1)]
print(formatToPrint(estimatedL2Rate))

print("\n== LInf ERRORS for x ==")
LInfErrors = computeErrors(uRef, uVector, x = True, norm = "Inf")
print(formatToPrint(LInfErrors))
estimatedLInfRate = [np.log2(LInfErrors[i] / LInfErrors[i + 1]) for i in range(len(LInfErrors) - 1)]
print(formatToPrint(estimatedLInfRate))

##################
### Y ANALYSIS ###
##################

nyRange = [8 * 2 ** i for i in range(4)]
nt = ntRef; nx = 128
filenames = ["../src/sol_{}_{}_{}.txt".format(nx, ny, nt) for ny in nyRange]
uVector = [readFile(filename) for filename in filenames]

print("\n== L1 ERRORS for y ==")
L1Errors = computeErrors(uRef, uVector, norm = 1)
print(formatToPrint(L1Errors))
estimatedL1Rate = [np.log2(L1Errors[i] / L1Errors[i + 1]) for i in range(len(L1Errors) - 1)]
print(formatToPrint(estimatedL1Rate))

print("\n== L2 ERRORS for y ==")
L2Errors = computeErrors(uRef, uVector, norm = 2)
print(formatToPrint(L2Errors))
estimatedL2Rate = [np.log2(L2Errors[i] / L2Errors[i + 1]) for i in range(len(L2Errors) - 1)]
print(formatToPrint(estimatedL2Rate))

print("\n== LInf ERRORS for y ==")
LInfErrors = computeErrors(uRef, uVector, norm = "Inf")
print(formatToPrint(LInfErrors))
estimatedLInfRate = [np.log2(LInfErrors[i] / LInfErrors[i + 1]) for i in range(len(LInfErrors) - 1)]
print(formatToPrint(estimatedLInfRate))