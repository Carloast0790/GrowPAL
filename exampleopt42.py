from inout.messag import *
from conf.conf import queue_name
qname=queue_name()

long_string = """
#GLOMOS 1.0

option                2

## Initial population
fragment_1            Mo20.xyz
fragment_2            Mo

## Discrimination
similarity_tol        0.95
disc_unconnected_opt  off
energy_range          100.0

## Calculation parameters
zero_energy_difference  1
number_of_stages        1
latt_space              15.0

## Queue system and Computer resources
calculator              gulp
qsys                    local

---GULP1---
opti conj nosymmetry conv
switch_minimiser bfgs gnorm 0.001
vectors
LATTICEVECTORS
frac
COORDINATES
space
1
species
Mo  0.0
lennard epsilon
Mo Mo 1.000 3.000 0.0 12.0
maxcyc 850
---GULP1---

#### Please find the path to gulp executable and change it
---GULP.CONF---
exe_gulp=path/gulp
---GULP.CONF---""".format(queue=qname)

bilfile='INPUT.txt'
exfile = open(bilfile, "w")
exfile.write(growpal_short_string) 
exfile.write(welcome_growpal)
exfile.write(long_string)
exfile.close()
