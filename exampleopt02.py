from inout.messag import *
from conf.conf import queue_name
qname=queue_name()

long_string = """
option                2

## Initial population
fragment_1            b7.xyz
fragment_2            B

## Discrimination
similarity_tol        0.95
disc_unconnected_opt  off
energy_range          100.0

## Queue system and Computer resources
nprocshared           4
memory_in_gb          2
queue                 qgpu
njobs                 1
timesleep             2.0

## Discrimination/Calculation parameters
number_of_stages        1
walltime                01:00:00
percent_of_convergence  100.0
no_attempts_opt         2
zpe_energy_correction   0
disc_neg_frequencies    0
normal_termination      0

---GAUSSIAN1---
# PBE1PBE 6-31G(d)  OPT=(CARTESIAN,MAXCYCLES=512) SCF=(XQC,MAXCYCLES=512)

Comment: GrowPAl

0   1
COORDINATES
---GAUSSIAN1---""".format(queue=qname)

bilfile='INPUT.txt'
exfile = open(bilfile, "w")
exfile.write(growpal_short_string)
exfile.write(welcome_growpal)
exfile.write(long_string)
exfile.close()
