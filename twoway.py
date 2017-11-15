#TWO WAY PAINTING ALGORITHM
import numpy as np
from num_branches import compute_branches

Sx = 10
Sz = 10

def twoway(data):

    connected = np.zeros((Sz,Sx))
    connected2 = np.zeros((Sz,Sx))


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

    #print connected
    #print ""
    #print connected2

    loaded = np.logical_and((connected > 0), (connected2 > 0))
    loaded = 1.0 * loaded

    return loaded


# Weinkammer 2004

# Matrix loaded architectural information
def Az(matrix):
    summ = 0.0

    for z in range(Sz):
        branches = compute_branches(matrix[z])
        Mz = len(branches)
        if Mz > 0:
            factor = 1.0 / Mz * Mz
            summ += factor * np.sum(1.0 / branches)

    return summ

# Matrix loaded architectural information
def Ax(matrix):
    summ = 0.0

    for x in range(Sx):
        branches = compute_branches(matrix[:,z])
        Mx = len(branches)
        if Mx > 0:
            factor = 1.0 / Mx * Mx
            summ += factor * np.sum(1.0 / branches)

    return summ

def compute_len_branch_x(matrix, z, x):
    res = 1 # this position
    i = x-1
    
    while i > 0 and matrix[z,i] > 0:
        res += 1
        i -= 1

    i = x+1

    while i < Sx and matrix[z,i] > 0:
        res += 1
        i += 1

    return res

def compute_len_branch_z(matrix, z, x):
    res = 1 # this position
    i = z-1
    
    while z > 0 and matrix[i, x] > 0:
        res += 1
        i -= 1

    i = z+1

    while i < Sz and matrix[i, x] > 0:
        res += 1
        i += 1

    return res

# maximal pressure p at a given position (in matrix matrix)
# (mechanical stimulus the bone cells respond to)
def p(matrix, z, x):

    # no bone, no pressure
    if matrix[z,x] == 0:
        return 0

    Fz = 9 # MISSING
    k = 0.9 # paper
    a = 1 # MISSING

    Ax_v = Ax(matrix)
    Az_v = Az(matrix)
    loaded = twoway(matrix)

    branches_x = compute_branches(matrix[:,x])
    branches_z = compute_branches(matrix[z])

    Mx = len(branches_x)
    Mz = len(branches_z)
    Njx = compute_len_branch_x(matrix, z, x)
    Njz = compute_len_branch_z(matrix, z, x)

    factor = Fz / (Ax_v + k*k*Az_v)

    v1 = k * Az_v / (Mx * a * Njx)
    v2 = Ax_v / (Mz * a * Njz)

    maxim = max(v1, v2)

    result = round(factor * maxim, 2)

    return result


data = np.random.rand(Sz,Sx)

data = 1.0*(data > 0.5)

print "DATA: "
print data
print

loaded = twoway(data)

print "Loaded: "
print loaded

P = np.zeros((Sz, Sx))

for z in range(Sz):
    for x in range(Sx):
        P[z,x] = p(loaded, z, x)

print "p"
print P
