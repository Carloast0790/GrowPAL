from discriminate.usrp import kick_similar_molecules, molin_sim_molref
from inout.getbilparam import get_a_float, get_a_str
from growpal.gen_addatom import make_many_rand
#------------------------------------------------------------------------------------------
tol_gen = get_a_float('similarity_tol', 0.95)
log_file = get_a_str('output_file','glomos_out.txt')
#------------------------------------------------------------------------------------------
def make_molecules_check_added(molist, symbol):
    molist_out = make_many_rand(molist, symbol)
    fopen = open(log_file,'a')
    print('-------------------------------------------------------------------',file=fopen)
    print('-------------------------------------------------------------------',file=fopen)
    print('After additions we have '+str(len(molist_out))+' structures',file=fopen)
    molist_out = kick_similar_molecules(molist_out, tol_gen,1)
    print('After discrimination we have '+str(len(molist_out))+' structures',file=fopen)
    print('-------------------------------------------------------------------',file=fopen)
    print('-------------------------------------------------------------------\n',file=fopen)
    fopen.close()
    return molist_out

#------------------------------------------------------------------------------------------
def run_sample():
    from utils.libmoleculas import readxyzs, writexyzs
    m = readxyzs('mo22.xyz')
    l = make_molecules_check_added(m,'Mo')
    writexyzs(l,'out.xyz')
# run_sample()
