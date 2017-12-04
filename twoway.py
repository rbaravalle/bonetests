# Weinkammer 2004 Stochastic lattice model for bone remodeling and aging
import numpy as np
from scipy.ndimage.measurements import label
import Image
import argparse
import os

void_fraction_start = 0.95
pminus = 0.007
k = 0.9
alpha = 0.003
beta = 0.1

a = 1  # MISSING
Fz = 30 # MISSING
Lx = 0 # WILL BE DEFINED

ps_c = 2.5 # MEAN PAPER



def compute_amount_branches(array):
    return label(array)[1]

def compute_branches(array):
    branches, amount = label(array)
    
    lens = np.zeros(amount)

    for i in range(1,amount+1):
        lens[i-1] = len(np.where(branches == i)[0])

    return lens

def create_dir_if_not_exists(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

def twoway(data, Sz, Sx):

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
            or (x < Sx-1 and connected2[z+1][x+1] == 1):
                # CONNECTED
                if data[z,x] == 1:
                    connected2[z,x] = 1

    loaded = 1.0*np.logical_and((connected > 0), (connected2 > 0))

    return loaded


# Matrix loaded architectural information
def Az(matrix):
    summ = 0.0

    for z in range(matrix.shape[0]):
        branches = compute_branches(matrix[z])
        Mz = len(branches)
        if Mz > 0:
            factor = 1.0 / Mz * Mz
            summ += factor * np.sum(1.0 / branches)

    return summ

# Matrix loaded architectural information
def Ax(matrix):
    summ = 0.0

    for x in range(matrix.shape[1]):
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

    while i < matrix.shape[1] and matrix[z,i] > 0:
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

    while i < matrix.shape[0] and matrix[i, x] > 0:
        res += 1
        i += 1

    return res

# maximal pressure p at a given position (in matrix matrix)
# (mechanical stimulus the bone cells respond to)
def p(matrix, z, x, Mx, Mz, Ax_v, Az_v):

    # no bone, no pressure
    if matrix[z,x] == 0:
        return 0


    Njx = compute_len_branch_x(matrix, z, x)
    Njz = compute_len_branch_z(matrix, z, x)

    factor = Fz / (Ax_v + k*k*Az_v)

    v1 = k * Az_v / (Mx * a * Njx)
    v2 = Ax_v / (Mz * a * Njz)

    maxim = max(v1, v2)

    result = round(factor * maxim, 2)

    return result

# return possible neighbors of z,x in matrix
def get_4_neighbors(matrix, z, x):
    neighbors = []

    if z < 0 or x < 0 or z > matrix.shape[0] or x > matrix.shape[1]:
        print "Invalid call to get_neighbors!. Shape:", matrix.shape, z, x
        exit()

    if x-1 > 0:
        neighbors.append([z, x-1])

    if x+1 < matrix.shape[1]:
        neighbors.append([z, x+1])

    if z-1 > 0:
        neighbors.append([z-1, x])

    if z+1 < matrix.shape[0]:
        neighbors.append([z+1, x])


    return neighbors

# return possible neighbors of z,x in matrix
def get_8_neighbors(matrix, z, x):
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
def Pplus(matrix, p_matrix, z, x, pc):

    neighbors_pos = get_4_neighbors(matrix, z, x)

    n = len(neighbors_pos)

    sum_gi = 0.0
    for i in neighbors_pos:
        # pressure at neighbor location
        pi = p_matrix[i[0],i[1]]

        sum_gi += 0.0 if pi < pc else float(pi) / pc

    return n * alpha + beta * sum_gi


# generate matrix with given void fraction
def generate_matrix(void_fraction, Sz, Sx):
    
    data = (np.random.rand(Sz,Sx) > 1 - void_fraction).astype(np.float32)

    return data

def simulate(loaded, steps, Sz, Sx, str_id, out_dir, pc):
    for t in range(steps):

        if t % 5 == 0:
            save_img(loaded, out_dir+'/bone_out'+str_id+'_iteration_'+str(t)+'.png')

        P = np.zeros((Sz, Sx))

        Mz = np.zeros(Sz)
        Mx = np.zeros(Sx)
        for z in range(Sz):
            Mz[z] = compute_amount_branches(loaded[z])
        for x in range(Sx):
            Mx[x] = compute_amount_branches(loaded[:,x])


        Ax_v = Ax(loaded)
        Az_v = Az(loaded)

        # compute pressures
        for z in range(Sz):
            for x in range(Sx):
                P[z,x] = p(loaded, z, x, Mz[z], Mx[x], Ax_v, Az_v)


        # probabilities of new material
        pp = np.zeros((Sz, Sx))

        for z in range(Sz):
            for x in range(Sx):
                if P[z,x] > 0:
                    pp[z,x] = Pplus(loaded, P, z, x, pc)


        # use probabilities to put or remove material
        loaded = 1*np.logical_or(loaded, (pp > np.random.random((Sz, Sx))))
        loaded -= 1*(pminus > np.random.random((Sz, Sx)))
        loaded = np.maximum(loaded, np.zeros((Sz, Sx)))
        #loaded = loaded*1

        loaded = twoway(loaded, Sz, Sx)


    return loaded
        

def save_img(arr, filename):
    im = Image.fromarray(255*arr)
    if im.mode != 'RGB':
        im = im.convert('RGB')
    im.save(filename)
    print "Saved", filename

def main():
    parser = argparse.ArgumentParser(description='Run bone simulation')
    parser.add_argument("-s", dest="steps", type=int, required=True, nargs=1, help="Number of simulation steps")
    parser.add_argument("-sz", dest="Sz", type=int, required=True, nargs=1, help="Size in z direction")
    parser.add_argument("-sx", dest="Sx", type=int, required=True, nargs=1, help="Size in x direction")


    
    args = parser.parse_args()

    str_id = str(args.steps[0])+'_'+str(args.Sz[0])+'_'+str(args.Sx[0])+'_'

    data = generate_matrix(void_fraction_start, args.Sz[0], args.Sx[0])

    # run two way algorithm
    loaded_orig = twoway(data,args.Sz[0], args.Sx[0])
    loaded = loaded_orig

    Lx = float(a * args.Sx[0])
    pc = ps_c/(Fz/Lx)
    print "PC:", pc

    out_dir = 'output'
    create_dir_if_not_exists(out_dir)

    print "Starting..."

    final = simulate(loaded, args.steps[0], args.Sz[0], args.Sx[0], str_id, out_dir, pc)

    print "Finished"


    save_img(final, out_dir+'/bone_out'+str_id+'.png')

if __name__ == '__main__':
    main()
