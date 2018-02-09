import numpy as np
from bisect import bisect_left

data = np.loadtxt('test.dat')
bin_centers = data[:,0]
energy = data[:,1]
delta_x = bin_centers[1] - bin_centers[0]

# calculate integral of a bounded region using Simpson's Rule.
def integrateSimpsons(r1, r2, bin_centers, function):
    Area = 0.

    r1_index = bisect_left(bin_centers,r1)## returns the index of the 'Coordinate' point in bin_centers that is closest to 'r_critical'
    r2_index = bisect_left(bin_centers,r2)## returns the index of the 'Coordinate' point in bin_centers that is closest to 'r_critical'
    index_range = r2_index - r1_index ## total number of elements in the monomer range. This will need to be an odd number for Simpson's rule.
    print 'r1 and r1_index: ',r1, r1_index
    print 'r2 and r2_index: ',r2, r2_index
    print 'bin_centers[r1_index]: ',bin_centers[r1_index]
    print 'bin_centers[r2_index]: ',bin_centers[r2_index]
    #print 'index_range = ',index_range

    #crit_bin_index_mon = index_range - 1 ## last index value of monomer_range ie. [0,monomer_range) == [0,crit_bin_index_mon]. This will need to be an even number for Simpson's rule.

    range_mod = index_range % 2 ## == 0 if divisible by 2. =/= 0 otherwise. For Simpson's Rule
    if range_mod:## This statement will execute if monomer_mod =/= 0
        index_range -= 1 ## minus instead of plus because our array doesnt have any more points in it to add.
    #    print 'Shifted the last index by 1 for Simpson\'s Rule: %s' %(crit_bin_index_mon)
    #else:## This statement will execute if monomer_mod == 0
    #    print 'Already even # of monomer subdivisions for Simpson\'s Rule: %s' %(crit_bin_index_mon)

    count = 0
    for i in range(r1_index, r2_index+1):
        #print count, i
        if count == 0 or count == (r2_index - r1_index):
            Area += function[i]
            #print 'first or last: ',i ,function[i]
        elif count%2 == 0:
            Area += 2*function[i]
            #print 'even: ',i, function[i]
        else:
            Area += 4*function[i]
            #print 'odd: ',i ,function[i]
        count += 1

    Area *= delta_x/3.
    
    return Area;


# Main Program

r1 = 2.
r2 = 3.

area = integrateSimpsons(r1, r2, bin_centers, energy)
print 'Area: ',area
