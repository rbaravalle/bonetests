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
a = 1.0
pad = 5
min_void_fraction = 0.4

Fz = 30.0 # MISSING
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


def paint_x(data, op, from_v, to_v, step_v):
    array = np.zeros(data.shape)
    array[:,from_v] = data[:,from_v]
    for x in range(from_v, to_v, step_v):
        for z in range(data.shape[0]):   
            if array[z][x] and data[z][op(x,1)]:
                array[z][op(x,1)] = 1

                # additionally, neighbors
                i = 1
                while z+i < data.shape[0] and data[z+i][op(x,1)]:
                    array[z+i][op(x,1)] = 1
                    i+=1
                    
                i = 1
                while z-i > 0 and data[z-i][op(x,1)]:
                    array[z-i][op(x,1)] = 1
                    i+=1
    return array

def paint(data, op, from_v, to_v, step_v):
    array = np.zeros(data.shape)
    array[from_v] = data[from_v]
    for z in range(from_v, to_v, step_v):
        for x in range(data.shape[1]):
            if array[z][x] and data[op(z,1)][x]:
                array[op(z,1)][x] = 1

                # additionally, neighbors
                i = 1
                while x+i < data.shape[1] and data[op(z,1)][x+i]:
                    array[op(z,1)][x+i] = 1
                    i+=1
                    
                i = 1
                while x-i > 0 and data[op(z,1)][x-i]:
                    array[op(z,1)][x-i] = 1
                    i+=1
    return array

def twoway_new(data):

    # bottom -> top
    connected = paint_x(data, np.add, 0, data.shape[0]-1, 1)

    connected2 = paint_x(connected, np.subtract, data.shape[0]-1, -1, -1)

    return connected2

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

    if Mx == 0 or Mz == 0:
        return 0


    Njx = compute_len_branch_x(matrix, z, x)
    Njz = compute_len_branch_z(matrix, z, x)

    if Njx == 0 or Njz == 0 or Az_v == 0 or Ax_v == 0:
        return 0

    factor = Fz / (Ax_v + k*k*Az_v)

    v1 = k * Az_v / (Mx * a * Njx)
    v2 = Ax_v / (Mz * a * Njz)

    maxim = max(v1, v2)

    result = factor * maxim

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


def random_values(Sz, Sx):
    r = np.random.random((Sz,Sx))
    return r 

# generate matrix with given void fraction
def generate_matrix(void_fraction, Sz, Sx):
    data = np.ones((Sz, Sx)).astype(np.float32)
    data[pad:Sz-pad, pad:Sx-pad] = (random_values(Sz-2*pad,Sx-2*pad) > 1.0 - void_fraction).astype(np.float32)

    return data


def compute_void_fraction(loaded):
    amount = loaded.shape[0]*loaded.shape[1]
    summ = np.sum(loaded)
    return summ / amount

def simulate(loaded, steps, Sz, Sx, str_id, out_dir, pc):

    for t in range(steps):

        if t % 2 == 0:
            save_img(loaded, out_dir+'/bone_out'+str_id+'_iteration_'+str(t)+'.png')
            vf = compute_void_fraction(loaded)
            print vf
            if vf < min_void_fraction:
                print "Wrong simulation"
                exit()

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

        print np.max(P[21:, 21:])

        # probabilities of new material
        pp = np.zeros((Sz, Sx))

        for z in range(Sz):
            for x in range(Sx):
                if P[z,x] > 0:
                    pp[z,x] = Pplus(loaded, P, z, x, pc)




        pminus_m = np.ones((Sz, Sx))
        pminus_m[pad:Sz-pad, pad:Sx-pad] = random_values(Sz-2*pad, Sx-2*pad)
        # use probabilities to put or remove material
        loaded -= 1.0*(pminus > pminus_m)
        loaded = np.maximum(loaded, np.zeros((Sz, Sx)))
        loaded = 1.0*np.logical_or(loaded, (pp > random_values(Sz, Sx)))


        loaded = twoway_new(loaded)


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


    out_dir = 'output'
    create_dir_if_not_exists(out_dir)
    save_img(data, out_dir+"/original.png")

    # run two way algorithm
    loaded = twoway_new(data)

    save_img(loaded, out_dir+"/twoway.png")


    Lx = float(a * args.Sx[0])
    pc = ps_c*Fz/Lx
    print "PC:", pc


    print "Starting..."

    final = simulate(loaded, args.steps[0], args.Sz[0], args.Sx[0], str_id, out_dir, pc)

    print "Finished"


    save_img(final, out_dir+'/bone_out'+str_id+'.png')

if __name__ == '__main__':
    main()
