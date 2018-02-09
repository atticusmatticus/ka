## USAGE:
## python ka.py [pmf_file] [r_critical_value]

import numpy as np
import math
import sys
import scipy.constants as scc
from bisect import bisect_left
wham_file = sys.argv[1]
r_critical = float(sys.argv[2])

kT = 0.5921868
max_dist_5mM = 87.2516
bin_centers = []
free_energy = []
fe_err = []
C_values = []
with open(wham_file,'r') as f:
    for line in f:
        temp = line.split()## split each element of line into its own element of a list
        if temp[0] == '#Coor' or temp[0] == '#Window':
            continue
            ## just skip to next line if the line starts with #Coor or #Window, ie. the start of the pmf and pmf_vert listings
        elif temp[1] == 'inf' or temp[2] == 'nan':
            continue
            ## if the free energy of a line is infinite (ie. the second column = 'inf'), skip that line as well
        elif temp[0][0] != '#':## if the first character of the first element of a line is not '#' then do the following
            bin_centers.append(float(temp[0]))## add the first column (coordinate value) to this array
            free_energy.append(float(temp[1]))## add the second column (Free energy) to this array
            fe_err.append(float(temp[2]))## add the third column (+/- in FE) to this array
        # GRABBING DATA TO BE USED TO CREATE THE UNSTICKED FE PLOTS; THIS ELSE STATEMENT GRABS THE CONSTANT VALUES THAT ARE USED TO ALIGN THE WINDOWS WITH EACH OTHER;                                       
        else:## the only values left are the vertical offset values at the end of the pmf file
            C_values.append(float(temp[1]))## grap the second column (FE offset) and add it to this array

###### Turn these lists into arrays
bin_centers = np.array(bin_centers)
free_energy = np.array(free_energy)
delta_x = (bin_centers[-1] - bin_centers[0])/len(bin_centers) ## step size in Angstroms

fe_err = np.array(fe_err)
C_values = np.array(C_values)

min_index = np.argmin(free_energy)## return the index of the first occurrence of the global minimum in the 'free_energy' array and set it as min_index
max_index = len(bin_centers) - 1## set the total number of Coordinate points as the 'max_index'
PMF_max = free_energy[max_index]## set the FE value at the largest Coordinate distance and set it to 'PMF_max' (NOT THE MAXIMUM FE)!!!


print "Gw = %s kCal/mol, %s kJ/mol" %(-PMF_max-1.2*math.log(bin_centers[max_index]/bin_centers[min_index]),(-PMF_max-1.2*math.log(bin_centers[max_index]/bin_centers[min_index]))*4.184)## print the volume adjusted delta G between the PMF minimum value at min_index and the PMF value at max_index and set the FE value at PMF_max to be zero.

R_max = bin_centers[max_index]## coordinate distance value at max_index

delta_PMF = -kT*np.log(R_max*R_max) - PMF_max ## Why the 4*pi? 4*pi should be in every baseline correction but it only ends up shifting the entire PMF curve by log(4*pi) and since the FE scale is relative it doesn't matter here. Baseline entire pmf so that PMF_max is at: FE = 0.


## VOLUME CORRECT FREE ENERGY
for i in range(len(free_energy)):
    free_energy[i] += 2.0*kT*np.log(bin_centers[i])

free_energy += delta_PMF ## vertically shift every FE value of the PMF by the volume corrected value of the PMF at long R, (delta_PMF).

# add distance points till it gets to max_dist_5mM
while True:
    bin_centers = np.append(bin_centers, bin_centers[-1]+delta_x)
    free_energy = np.append(free_energy, 0.)
    if bin_centers[-1] >= max_dist_5mM:
        break
max_index = len(bin_centers) - 1 ## set the total number of Coordinate points as the 'max_index'

######## PLOT SHIFTED PMF
test_out = np.column_stack((bin_centers,free_energy))
np.savetxt('free_energy.out',test_out)## free_energy values are now volume corrected and shfited to zero at long R. They have been shifted by the delta_PMF which is itself volume corrected.
########

exp_free_energy = np.exp(-free_energy/kT)## make Boltzmann factor for each free_energy. g(r) = exp( - PMF(r) / kT ). So exp_free_energy is g(r). This makes sense because the PMF is the energy as a function of COM distance, thus if you put the PMF in a boltzmann factor you'll get a weighting function based on the energy as a function of COM displacement, ie g(r).

######## PLOT g(r)
test_out = np.column_stack((bin_centers,exp_free_energy))
np.savetxt('exp_free_energy.out',test_out)
########


######################################################### Integral the way it's written in De Jong
## Simpson's Rule requires an even number of subdivisions! (crit_bin_index must be even!)
## Since the index goes from 0 --> n, EVEN values of n indicate an EVEN number of subdivisions.

####################################################################################################
## NUMERATOR in De Jong eqn 44. ("dimer" radii) 
dimer_int = 0 ## initialize numerator integral area
#dimer_range = [] ## initialize dimer integration range array

crit_bin_index = bisect_left(bin_centers,r_critical)## returns the index of the 'Coordinate' point in bin_centers that is closest to 'r_critical'

#print 'crit_bin_index = %s' %crit_bin_index

crit_mod = crit_bin_index%2 ## == 0 if divisible by 2. =/= 0 otherwise. For Simpson's Rule
#print 'crit_mod = %s' %crit_mod
if crit_mod:## This statement will execute if crit_mod =/= 0
    crit_bin_index += 1
#    print 'Shifted your choice of r_critical by 1 index for Simpson\'s Rule: %s' %(crit_bin_index)
#else:## This statement will execute if crit_mod == 0
#    print 'Already even # of dimer subdivisions for Simpson\'s Rule: %s' %(crit_bin_index)

print 'R_critical = %s Angstroms' %(bin_centers[crit_bin_index])

#for i in range(crit_bin_index + 1):
#    dimer_range.append(float(bin_centers[i]))## Show distance values that correspond to "dimer" range for double check.

#print crit_bin_index, len(dimer_range)

for i in range(crit_bin_index + 1):
    if i == 0 or i == (crit_bin_index):
        dimer_int += exp_free_energy[i] * bin_centers[i]**2
        #print 'executed'
    elif i%2 == 0:
        dimer_int += 2*exp_free_energy[i] * bin_centers[i]**2
    else:
        dimer_int += 4*exp_free_energy[i] * bin_centers[i]**2

dimer_int *= delta_x/3.
print 'numerator_int: ',dimer_int

####################################################################################################
## DENOMINATOR in De Jong eqn 44. ("monomer" radii) 
monomer_int = 0

monomer_range = (max_index + 1) - crit_bin_index ## total number of elements in the monomer range. This will need to be an odd number for Simpson's rule.
#print max_index, crit_bin_index
#print 'total monomer_range = %s' %monomer_range

crit_bin_index_mon = monomer_range - 1 ## last index value of monomer_range ie. [0,monomer_range) == [0,crit_bin_index_mon]. This will need to be an even number for Simpson's rule.

monomer_mod = crit_bin_index_mon%2 ## == 0 if divisible by 2. =/= 0 otherwise. For Simpson's Rule
if monomer_mod:## This statement will execute if monomer_mod =/= 0
    crit_bin_index_mon -= 1 ## minus instead of plus because our array doesnt have any more points in it to add.
#    print 'Shifted the last index by 1 for Simpson\'s Rule: %s' %(crit_bin_index_mon)
#else:## This statement will execute if monomer_mod == 0
#    print 'Already even # of monomer subdivisions for Simpson\'s Rule: %s' %(crit_bin_index_mon)

print 'Maximum R integrated = %s Angstroms' %(bin_centers[crit_bin_index + crit_bin_index_mon])## Not double counting r_critical index. crit_bin_index is inclusive of the r_critical index and monomer_range is exclusive of it.

it = 0
for i in range(crit_bin_index, max_index + 1):
    #print it, i
    if it == 0 or it == ((max_index + 1) - (crit_bin_index) - 1):
        #print 'executed'
        monomer_int += exp_free_energy[i] * bin_centers[i]**2
    elif it%2 == 0:
        monomer_int += 2*exp_free_energy[i] * bin_centers[i]**2
    else:
        monomer_int += 4*exp_free_energy[i] * bin_centers[i]**2
    it += 1

monomer_int *= delta_x/3.
print 'denominator_int: ',monomer_int

#coeff = ( (bin_centers[max_index])**3 )/((bin_centers[crit_bin_index])**3)
coeff = (4*np.pi/3.) * ( (bin_centers[max_index])**3 )/( 1.66 )
print "Coefficient: ",coeff
dA = - kT * np.log( coeff * (dimer_int / monomer_int) )
Ka = np.exp(-dA / kT)

#########################################################


print "Ka = %e" %(Ka)
Ga = -kT*np.log(Ka)
print "Ga = %s" %(Ga)
