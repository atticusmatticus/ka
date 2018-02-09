import numpy as np
import sys

initial = float(sys.argv[1])
final = float(sys.argv[2])
dx = float(sys.argv[3])

n_pts = int( (final - initial)/dx )+1
data = np.zeros((n_pts,n_pts), dtype=float)
for i in range(n_pts):
    data[i,0] = initial + i*dx
    #data[i,1] = initial + i*dx
    data[i,1] = 1./(initial + i*dx + 1)

#print data

out = open('test.dat','w')
for i in range(n_pts):
    out.write('%f %f\n' %(data[i,0], data[i,1]))

out.close
