#TWO WAY PAINTING ALGORITHM
import numpy as np

X = 10
Z = 10

data = np.random.rand(X,Z)

data = 1.0*(data > 0.5)

connected = np.zeros((X,Z))
connected2 = np.zeros((X,Z))

print "DATA: "
print data
print

# mark the first row as connected

# up -> down
connected[0] = data[0]


for x in range(1,X):
    for z in range(Z):
        if (z > 0 and connected[x-1][z-1] == 1) \
        or connected[x-1][z] == 1  \
        or (z < Z-1 and connected[x-1][z+1] == 1):
            # CONNECTED
            if data[x,z] ==1:
                connected[x,z] = 1

# down -> up
connected2[X-1] = data[X-1]


for x in range(X-2,-1,-1):
    for z in range(Z):
        if (z > 0 and connected2[x+1][z-1] == 1) \
        or connected2[x+1][z] == 1 \
        or (z < Z-1 and connected2[x+1][z+1] == 1):
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
