from inout.messag import *
from conf.conf import queue_name
qname=queue_name()

long_string = """
option                2

## Initial population
fragment_1            b9.xyz
fragment_2            B

calculator            mopac
qsys                  local

#GLOBAL SIMILARITY DISCRIMINATION OPTIONS
sym_algorithm_opt     0
sym_normalize_opt     0
disc_unconnected_opt  0

## Calculations and discrimination parameters (stages)
zero_energy_difference  1
number_of_stages        1
percent_of_convergence  100.0
num_of_attemps_to_conv  2
similarity_tolerance    0.95
maximum_energy_allowed  95.0
zpe_energy_correction   0
disc_neg_frequencies    0
normal_termination      0

---MOPAC1---
CHARGE=0 DOUBLET MNDO GNORM=0.0 xyz PRTXYZ

All coordinates are Cartesian
COORDINATES
---MOPAC1---""".format(queue=qname)

bilfile='INPUT.txt'
exfile = open(bilfile, "w")
exfile.write(growpal_short_string) 
exfile.write(welcome_growpal)
exfile.write(long_string)
exfile.close() 
