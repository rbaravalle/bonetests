# Weinkammer 2004 Stochastic lattice model for bone remodeling and aging
import numpy as np
from num_branches import compute_branches
import Image

Sx = 40
Sz = 40
void_fraction_start = 0.95
steps = 30
pminus = 0.007
k = 0.9
Fz = 9 # MISSING
a = 1 # MISSING

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
        branches = compute_branches(matrix[:,x])
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
    
    while i > 0 and matrix[i, x] > 0:
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

# return possible neighbors of z,x in matrix
def get_neighbors(matrix, z, x):
    neighbors = []

    if z < 0 or x < 0 or z > matrix.shape[0] or x > matrix.shape[1]:
        print "Invalid call to get_neighbors!. Shape:", matrix.shape, z, x
        exit()


    if z-1 > 0:
        if x-1 > 0:
            neighbors.append([z-1, x-1])
            neighbors.append([z, x-1])

        neighbors.append([z-1, x])

        if x+1 < matrix.shape[1]:
            neighbors.append([z-1, x+1])
            neighbors.append([z, x+1])

    if z+1 < matrix.shape[0]:
        neighbors.append([z+1, x])
        if x-1 > 0:
            neighbors.append([z+1, x-1])
        if x+1 < matrix.shape[1]:
            neighbors.append([z+1, x+1])


    return neighbors


# probability of forming new material at pos matrix[z,x]
def Pplus(matrix, p_matrix, z, x):
    alpha = 0.003
    beta = 0.1

    neighbors_pos = get_neighbors(matrix, z, x)
    n = len(neighbors_pos)

    # pressure at current location
    pc = p_matrix[z,x]

    sum_gi = 0
    for i in neighbors_pos:
        # pressure at neighbor location
        pi = p_matrix[i[0],i[1]]

        sum_gi += 0 if pi < pc else float(pi) / pc

    return n * alpha + beta * sum_gi


# generate matrix with given void fraction
def generate_matrix(void_fraction, Sz, Sx):
    
    data = (np.random.rand(Sz,Sx) > 1 - void_fraction).astype(np.float32)

    return data

def simulate(loaded, steps):
    for t in range(steps):

        if t % 5 == 0:
            print t

        P = np.zeros((Sz, Sx))

        # compute pressures
        for z in range(Sz):
            for x in range(Sx):
                P[z,x] = p(loaded, z, x)


        # probabilities of new material
        pp = np.zeros((Sz, Sx))

        for z in range(Sz):
            for x in range(Sx):
                if P[z,x] > 0:
                    pp[z,x] = Pplus(loaded, P, z, x)

        # use probabilities to put or remove material
        for z in range(Sz):
            for x in range(Sx):
                if pp[z,x] > np.random.random():
                    loaded[z,x] = 1
                
                if pminus > np.random.random():
                    loaded[z,x] = 0


    return loaded
        

def save_img(arr, filename):
    im = Image.fromarray(255*arr)
    if im.mode != 'RGB':
        im = im.convert('RGB')
    im.save(filename)
    print "Saved", filename

def main():
    data = generate_matrix(void_fraction_start, Sz, Sx)

    print "DATA: "
    print data
    print

    # run two way algorithm
    loaded_orig = twoway(data)
    loaded = loaded_orig

    save_img(loaded, "original.png")

    print "Starting..."
    #print loaded

    final = simulate(loaded, steps)

    print "Finished"

    save_img(final, "result_bone.png")

if __name__ == '__main__':
    main()
