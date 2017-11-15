import numpy as np

#array = np.random.randint(2, size=(10,1))


# given an array, compute the number of branches
# i.e. amount of connected components
# also compute its lengths
def compute_branches(array):
    lens = []

    # prev holds the state for the next iteration
    prev = 0

    # start
    if array[0] == 1:
        prev = 1
        # add a component with length 1
        lens.append(1)

    for i in range(1, len(array)):
        # traversing a branch
        if array[i] == 1 and prev == 1:
            lens[len(lens)-1] += 1

        # starting a new branch
        if array[i] == 1 and prev == 0:
            lens.append(1)
 
        prev = array[i]


    return np.array(lens)

#print array
#branches = compute_branches(array)
#print branches, len(branches)
