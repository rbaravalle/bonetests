#TWO WAY PAINTING ALGORITHM
import numpy as np
from num_branches import compute_branches

Sx = 10
Sz = 10

data = np.random.rand(Sz,Sx)

data = 1.0*(data > 0.5)

connected = np.zeros((Sz,Sx))
connected2 = np.zeros((Sz,Sx))

print "DATA: "
print data
print

def Az(matrix):
    summ = 0.0

    for z in range(Sz):
        branches = compute_branches(matrix[z])
        Mz = len(branches)
        factor = 1.0 / Mz * Mz
        summ += factor * np.sum(1.0 / branches)

    return summ

def Ax(matrix):
    summ = 0.0

    for x in range(Sx):
        branches = compute_branches(matrix[:,z])
        Mx = len(branches)
        factor = 1.0 / Mx * Mx
        summ += factor * np.sum(1.0 / branches)

    return summ

# mark the first row as connected

# up -> down
connected[0] = data[0]


for z in range(1,Sz):
    for x in range(Sx):
        if (x > 0 and connected[z-1][x-1] == 1) \
        or connected[z-1][x] == 1  \
        or (x < Sx-1 and connected[z-1][x+1] == 1):
            # CONNECTED
            if data[z,x] ==1:
                connected[z,x] = 1

# down -> up
connected2[Sz-1] = data[Sz-1]


for z in range(Sz-2,-1,-1):
    for x in range(Sx):
        if (x > 0 and connected2[z+1][x-1] == 1) \
        or connected2[z+1][x] == 1 \
        or (x < Sz-1 and connected2[z+1][x+1] == 1):
            # CONNECTED
            if data[z,x] == 1:
                connected2[z,x] = 1

print connected
print ""
print connected2

loaded = np.logical_and((connected > 0), (connected2 > 0))
loaded = 1.0 * loaded

print "Loaded: "
print loaded

print Az(loaded), Ax(loaded)
