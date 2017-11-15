#TWO WAY PAINTING ALGORITHM
import numpy as np
from num_branches import compute_branches

Sx = 20
Sz = 20

data = np.random.rand(Sx,Sz)

data = 1.0*(data > 0.5)

connected = np.zeros((Sx,Sz))
connected2 = np.zeros((Sx,Sz))

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


for x in range(1,Sx):
    for z in range(Sz):
        if (z > 0 and connected[x-1][z-1] == 1) \
        or connected[x-1][z] == 1  \
        or (z < Sz-1 and connected[x-1][z+1] == 1):
            # CONNECTED
            if data[x,z] ==1:
                connected[x,z] = 1

# down -> up
connected2[Sx-1] = data[Sx-1]


for x in range(Sx-2,-1,-1):
    for z in range(Sz):
        if (z > 0 and connected2[x+1][z-1] == 1) \
        or connected2[x+1][z] == 1 \
        or (z < Sz-1 and connected2[x+1][z+1] == 1):
            # CONNECTED
            if data[x,z] == 1:
                connected2[x,z] = 1

print connected
print ""
print connected2

final = np.logical_and((connected > 0), (connected2 > 0))
final = 1.0 * final

print "FINAL: "
print final

print Az(final), Ax(final)
