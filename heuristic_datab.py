import sys
import os.path
from inout.readbil       import get_bilfile
from inout.getbilparam   import get_a_int, get_a_float, get_a_str, get_int_list, get_str_list
from inout.messag        import *
from utils.libmoleculas  import readxyzs, writexyzs, rename_molecule, sort_by_energy, sort_by_stoichiometry, align_all_inertia_axis_x
from discriminate.energy import cutter_energy
from discriminate.usrp   import kick_similar_molecules
from discriminate.connectivity import cutter_nonconnectedmoleculargraph
from growpal.popgen_initadded import make_molecules_check_added
from gaussian.multiplicity import multiplicity_compared, get_charge_and_multi
#====================================================================
bilfile = get_bilfile()
log_file = get_a_str('output_file','glomos_out.txt')
nofstages = get_a_int('number_of_stages', 2)
emax =  get_a_float('energy_range', 99.0)
sim = get_a_float('similarity_tol', 0.96)
ncnt = get_a_str('disc_unconnected_opt', 'off')
zed = get_a_int('zero_energy_difference', 0)
dbop = get_a_int('option', 2)
frag_1 = get_a_str('fragment_1', 'off')
frag_2 = get_a_str('fragment_2', 'off')
#====================================================================
fopen = open(log_file,'w')
pid = os.getpid()
print("Current pid          = %s:" %(pid),file=fopen)
fopen.write(growpal_short_string)
fopen.write(welcome_growpal)
print("-------------------------------------------------------------------", file=fopen)
print("Read input file      = %s:" %(bilfile),file=fopen)
print("Task                 = Growth Pattern algorithm", file=fopen)
print("Number of Stages     = %d" %(nofstages), file=fopen)
fopen.close()
#====================================================================
def build_population_0():
    initialfile = 'initial.xyz'
    fopen = open(log_file,'a')
    print("-------------------------------------------------------------------", file=fopen)
    print("-----------------------POPULATION  GENERATOR-----------------------", file=fopen)
    fopen.close()
    if not os.path.isfile(initialfile):
        if frag_1 != 'off' and frag_2 != 'off' and dbop == 2:
            data_set = readxyzs(frag_1)
            mol99 = make_molecules_check_added(data_set, frag_2)
            if flag == 'gaussian':
                chg_mlt = get_charge_and_multi('INPUT.txt')
                ml_check = multiplicity_compared(mol99[0],chg_mlt)
                if ml_check == 0:
                    fopen = open(log_file,'a')
                    print('\n',file=fopen)
                    print('------------------------------------------------------------',file=fopen)
                    print('The multiplicity seems to be incompatible, double check this',file=fopen)
                    print('------------------------------------------------------------',file=fopen)
                    fopen.close()
                    sys.exit() 
            moleculeout = align_all_inertia_axis_x(mol99)
            moleculeout = rename_molecule(mol99,'initial',4)
            writexyzs(moleculeout, initialfile)
        else:
            fopen = open(log_file,'a')
            print('----------------------------------------------------------------------------------------',file=fopen)
            print('There seems to be a problem with the fragments or the selected option, double check this',file=fopen)
            print('----------------------------------------------------------------------------------------',file=fopen)
            fopen.close()
            sys.exit()
    else:
        fopen = open(log_file,'a')
        print("%s exist .... we take it" %(initialfile), file=fopen)
        fopen.close()
        moleculeout=readxyzs(initialfile)
    return moleculeout
#====================================================================
def run_calculator(moleculelist, stage=0):
    folder='stage'+str(stage)+'/'
    fopen = open(log_file,'a')
    print("-------------------------------------------------------------------", file=fopen)
    print("------------------- LOCAL OPTIMIZATION: STAGE %d -------------------" %(stage), file=fopen)
    fopen.close()
    if flag=='gaussian':
        blockg='gaussian'+str(stage)
        from gaussian.calculator_all import calculator_gaussian_all_check
        moleculeout=calculator_gaussian_all_check(moleculelist, folder, blockg, stage)
    if flag=='vasppara':
        from vasp.calculator_all import calculator_vasp_all_check
        from utils.libmoleculas  import sort_by_stoichiometry
        mol00=sort_by_stoichiometry(moleculelist)
        moleculeout=calculator_vasp_all_check(mol00, folder, stage)
    if flag=='mopac':
        blockg='mopac'+str(stage)
        from mopac.calculator_all import calculator_mopac_all_check
        moleculeout=calculator_mopac_all_check(moleculelist, folder, blockg, stage)
    if flag=='dftb':
        from utils.libmoleculas import sort_by_stoichiometry
        from dftb.calculator_all import calculator_dftb_all_check
        moleculelist=sort_by_stoichiometry(moleculelist)
        blockfiles=get_str_list('specific_dftb_file', ['body1.txt','body2.txt'])
        blockname=blockfiles[stage]
        moleculeout=calculator_dftb_all_check(moleculelist, folder, blockname, stage)
    if flag=='gulp':
        from gulp.calculator_all import calculator_gulp_all_check
        from utils.libmoleculas  import sort_by_stoichiometry
        mol00=sort_by_stoichiometry(moleculelist)
        blockname='gulp'+str(stage)
        moleculeout=calculator_gulp_all_check(mol00, folder, blockname, stage)
    return moleculeout
#====================================================================
def run_discrimination(moleculelist, stage=0):
    fopen = open(log_file,'a')
    print("-------------------------------------------------------------------", file=fopen)
    print("----------------------DISCRIMINATION  PROCESS----------------------", file=fopen)
    fopen.close()
    if flag == 'gaussian':
        from gaussian.get_geometry import cutter_normterm_gaussian
        list_nt = get_int_list('normal_termination',[1,2])
        nt = list_nt[stage]
        moleculelist = cutter_normterm_gaussian(moleculelist,folder,nt)
    if ncnt == 'on':
        moleculelist = cutter_nonconnectedmoleculargraph(moleculelist,0)
    mol01 = cutter_energy(moleculelist, emax)
    if flag == 'gaussian':
        from gaussian.get_geometry import cutter_negfreq_gaussian
        listdnf = get_int_list('disc_neg_frequencies',[0,1])
        dscfr = listdnf[stage]
        if dscfr == int(1): mol01 = cutter_negfreq_gaussian(mol01,folder)
    mol02 = kick_similar_molecules(mol01, sim, 0)
    if len(mol01) == len(mol02):
        logfile = open(log_file,'a')
        print("ZERO elements discriminated by Similarity", file=logfile)
        logfile.close()
    moleculeout = sort_by_energy(mol02, zed)
    fopen = open(log_file,'a')
    print("\nAfter the discrimination process we have %d molecules" %(len(moleculeout)), file=fopen)
    fopen.close()
    return moleculeout
#====================================================================
def display_mol_info(moleculein, flag, stage=0):
    fopen = open(log_file,'a')
    print("-------------------------------------------------------------------", file=fopen)
    print("------------------------- SUMMARY STAGE %d -------------------------" %(stage), file=fopen)
    fopen.close()
    molzz = sort_by_energy(moleculein,1)
    fopen = open(log_file,'a')
    if flag == 'gulp':
        print("Number File-Name   Energy (eV)     Delta-E  NT", file=fopen)
    else:
        print("Number File-Name  Energy (kcal/mol) Delta-E  NT", file=fopen)
    fopen.close()
    for ii, imol in enumerate(molzz):
        deltae=imol.e - molzz[0].e
        nt = imol.c[0]
        fopen = open(log_file,'a')
        jj = str(ii+1).zfill(5)
        print("%5s %11s %14.9f %10.6f %d" %(jj, imol.i, imol.e, deltae, nt), file=fopen)
        fopen.close()
#====================================================================
flag = get_a_str('calculator','gaussian')
mol00 = build_population_0()
for stage in range(1,nofstages+1):
    basenm = 'stage'+str(stage)
    folder = basenm+'/'
    moluu = rename_molecule(mol00, basenm, 4)
    mol00 = run_calculator(moluu, stage)
    mol00 = run_discrimination(mol00, stage)
    if flag == 'gaussian':
        folder = 'stage'+str(stage)+'/'
        from gaussian.get_geometry import display_info_gaussian
        display_info_gaussian(mol00, folder, fopen)
    else:
        display_mol_info(mol00, flag, stage)
    writexyzs(mol00, basenm+'.xyz', 1)
fopen = open(log_file,'a')
print ("-------------------------------------------------------------", file=fopen)
print ("GLOMOS HAS FINISHED SUCCESSFULLY", file=fopen)
print ("-------------------------------------------------------------", file=fopen)
fopen.close()
# exit()


