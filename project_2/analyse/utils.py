import pandas as pd
import numpy as np

def computeL1Error(Delta):
    return np.mean(np.abs(Delta))

def computeL2Error(Delta, dt, dx, dy):
    return np.sqrt(dt * dx * dy * np.sum(np.square(Delta)))

def computeLInfError(Delta):
    return np.max(np.abs(Delta))

def formatToPrint(values):
    return [f"{x:.3e}" for x in values]

def readFile(filename):
    f = open(filename, 'r')
    nt, nx, ny = list(map(int, f.readlines(1)[0][:-1].split(" ")))
    toSkip = 1
    u = pd.read_csv(filename, sep = " ", header = None, skiprows = toSkip).values
    return np.reshape(u, (nt, nx, ny))

# stability issues
def isStable(nt, nx, ny, alpha = 2., Lt = 300., Lx = 50., Ly =50.):
    dt = Lt / nt
    dx = Lx / nx
    dy = Ly / ny
    bound = 1 / (2 * alpha) * dx * dx * dy * dy / (dx * dx + dy * dy)
    if dt <= bound:
        print("Problem is stable : {} <= {}".format(dt, bound))
    else:
        print("Problem is instable : {} > {}".format(dt, bound))

toTest = [
    (2400, 50, 50),
    (4800, 50, 50),
    (9600, 50, 50),
    (19200, 50, 50),
    (76800, 50, 50),
    (19200, 8, 128),
    (19200, 16, 128),
    (19200, 32, 128),
    (19200, 64, 128),
    (19200, 128, 8),
    (19200, 128, 16),
    (19200, 128, 32),
    (19200, 128, 64),
    (19200, 128, 128)
]

for (nt, nx, ny) in toTest:
    isStable(nt, nx, ny)