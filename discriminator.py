import spglib
import numpy as np
from vasp.libperiodicos import readposcars, writeposcars
from vasp.libperiodicos import direct2cartesian, cartesian2direct
from utils.atomic import get_atomic_number
from inout.getbilparam import get_a_int, get_a_str, get_a_float

log_file = get_a_str('output_file','glomos_out.txt')
min_dE = get_a_float('min_energy_difference',0.00004)
min_dV = get_a_float('min_volume_difference',0.0001) 
#----------------------------------------------------------------------------------------------------------
def find_spacegroup(xtal_in):
	'''
	This function uses spglib to find the space group of the crystal structures
	
	in: xtal_in (Molecule), the structure of interest
	out: spc_grp (str), the symbol for the space group of the structure
	'''
	lattice = xtal_in.m
	atms_position = []
	atms_specie = []
	for a in xtal_in.atoms:
		x,y,z = cartesian2direct(a.xc,a.yc,a.zc,xtal_in.m)
		x,y,z = round(x,2),round(y,2),round(z,2)
		atms_position.append([x,y,z])
		n = get_atomic_number(a.s)
		atms_specie.append(n)
	atms_specie.sort()
	x = (lattice,atms_position,atms_specie)
	spc_grp = spglib.get_spacegroup(x)
	return spc_grp

#----------------------------------------------------------------------------------------------------------
def get_volume(xtal_in):
	i,j,k = xtal_in.m[0],xtal_in.m[1],xtal_in.m[2]
	jxk = np.cross(j,k)
	vol = np.dot(i,jxk)
	return vol

#----------------------------------------------------------------------------------------------------------
def get_fingerprint(xtal_in):
	'''
	This function builts the fingerprint of each structure. According to Zurek et.al. this comprises a list
	like this: [energy,spacial_group,volume]. 
	'''
	sg = find_spacegroup(xtal_in)
	v = get_volume(xtal_in)
	e = xtal_in.e
	fp = [e,sg,v]
	return fp

#----------------------------------------------------------------------------------------------------------
def compare_fingerprints(xtal_a,xtal_b,vol_restr):
	t_atma, t_atmb = len(xtal_a.atoms), len(xtal_b.atoms)
	fp_a = get_fingerprint(xtal_a)
	fp_b = get_fingerprint(xtal_b)
	# get the energy difference per atom of the two structures
	edif = abs(fp_a[0] - fp_b[0])/t_atma
	# get the volume difference too
	vdif = abs(fp_a[2] - fp_b[2])/t_atma
	sym_a, sym_b = fp_a[1],fp_b[1]
	equal = True
	if edif <= min_dE:
		# dE is too small, considere them equal
		return equal, edif, vdif, sym_a, sym_b
	elif vdif <= min_dV and vol_restr == False:
		# dE is large enough but dV is too small and there is no contrain on volume imposed by the user, we considere them equal
		return equal, edif, vdif, sym_a, sym_b
	else:
		# dH and dV are large enough, we considere them differet 
		equal = False
		return equal, edif, vdif, sym_a, sym_b

#---------------------------------------------------------------------------------------------------------
def discriminate_calculated(xtal_list, vol_restr):
	xtalist_out = xtal_list.copy()
	l_list = len(xtal_list)-1
	fopen = open(log_file,'a')
	print('-------------------------------------------------------------------',file=fopen)
	print('------------------- DISCRIMINATION Generation ---------------------',file=fopen)
	fopen.close()
	for i in range(l_list):
		str_a = xtal_list[i]
		for j in range(i+1,l_list):
			str_b = xtal_list[j]
			e, dE, dV, sym_a, sym_b = compare_fingerprints(str_a,str_b,vol_restr)
			if e == True:
				if str_b in xtalist_out:
					fopen = open(log_file,'a')
					print('%s sym %s discriminated, too similar to %s sym %s, dE = %.5f dV = %.5f' %(str_b.i,sym_b,str_a.i,sym_a,dE,dV),file=fopen)
					fopen.close()
					xtalist_out.remove(str_b)
	if len(xtalist_out) == len(xtal_list):
		fopen = open(log_file,'a')
		print('No similar structures found!',file=fopen)
		fopen.close()
	return xtalist_out

#---------------------------------------------------------------------------------------------------------
def discriminate_calculated_vs_pool(calulation_list, pool_list,vol_restr):
	fopen = open(log_file,'a')
	print('-------------------------------------------------------------------',file=fopen)
	print('------------------ DISCRIMINATION Gen vs Pool ---------------------',file=fopen)
	fopen.close()
	lista = calulation_list.copy()
	listb = pool_list.copy()
	xtalist_out = calulation_list.copy()
	for xtala in lista:
		div = len(listb)//2
		fh = listb[:div]
		sh = listb[div:]
		init = True
		while init == True:
			for xtalb in fh:
				equal, dE, dV, sym_a, sym_b = compare_fingerprints(xtala,xtalb,vol_restr)
				if equal == True:
					fopen = open(log_file,'a')
					print('%s sym %s discriminated, too similar to %s sym %s, dE = %.5f dV = %.5f' %(xtala.i,sym_a,xtalb.i,sym_b,dE,dV),file=fopen)
					fopen.close()
					init = False
					xtalist_out.remove(xtala)
					break
			if init == True:
				div = len(sh)//2
				fh = sh.copy()
				sh = fh[div:]
				fh = fh[:div]
				if len(fh) == 0:
					init = False
	if len(xtalist_out) == len(calulation_list):
		fopen = open(log_file,'a')
		print('No similar structures found!',file=fopen)
		fopen.close()
	return xtalist_out

# from vasp.libperiodicos import readposcars
# a = readposcars('initial000.vasp')
# b = readposcars('initial001.vasp')
# discriminate_calculated_vs_pool(a,b)

# def discriminate_list_list(xtal_lista, xtal_listb):
# 	lista = xtal_lista.copy()
# 	listb = xtal_listb.copy()
# 	for ela in lista:
# 		div = len(listb)//2
# 		fh = listb[:div]
# 		sh = listb[div:]
# 		init = False
# 		while init = False:
# 			if ela in fh:
# 				break
# 			else:
# 				div = len(sh)//2
# 				fh = sh.copy()
# 				sh = fh[div:]
# 				fh = fh[:div]
# 				if len(fh) == 0:
# 					break
