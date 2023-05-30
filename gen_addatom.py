import numpy as np
from utils.libmoleculas import Molecule, Atom, get_min_binding_distance, copymol, atom_molecule_dist, translate_to_cm, get_covalent_radius
from inout.getbilparam   import get_a_str, get_a_int, get_a_float, get_str_list
from discriminate.connectivity import conectmx
from pymatgen.core import Molecule as MoleculePMG
from pymatgen.symmetry.analyzer import PointGroupAnalyzer

log_file = get_a_str('output_file','glomos_out.txt')
add_type = get_a_int('addition_type',0)
l_add_type = get_str_list('addition_types',['a'])
crit_size = get_a_float('critical_size',20)
sym_rest = get_a_str('sym_restrict','on')
num_iso  = get_a_int('take_isomers',5)
#------------------------------------------------------------------------------------------
def closest_atm_to_cm(moleculein):
    '''
    This function finds the atom that is closest to the center of mass.
    in : moleculein (Molecule); The molecule from wich the atom will be obtained
    out: label (int); The position of the atom in the molecule.atoms list
    '''
    from utils.libmoleculas import center_off_mass
    cm = center_off_mass(moleculein)
    xcm, ycm, zcm = cm[0],cm[1],cm[2]
    d = 1000
    c = 0
    for iatm in moleculein.atoms:
        x,y,z = iatm.xc, iatm.yc, iatm.zc
        xv,yv,zv = (xcm-x),(ycm-y),(zcm-z)
        aux_d = np.sqrt(xv*xv + yv*yv + zv*zv)
        if aux_d < d:
            d = aux_d
            label = c
        c = c + 1
    return label

#------------------------------------------------------------------------------------------
def neighbor_finder(adj_mtx):
    '''
    This function builds a dictionary that contains the neighbors of all atoms
    
    in: adj_mtx (list); a list of list elements that represent the adjacency matrix 
        of the molecule as a graph.
    out: dict_neig (dict); a dictionary containing as keys the position of the atom 
        in the molecule.atoms list, and as values a list of all its neighbors.
        ({0:[1,3,7,11], 1:[0,7],...})
    '''
    dict_neig = {}
    conta = 0
    for i in adj_mtx:
        neig = []
        contb = 0
        for j in i:
            if j == 1:
                neig.append(contb)
            contb = contb + 1
        dict_neig[conta] = neig
        conta = conta + 1
    return dict_neig

#------------------------------------------------------------------------------------------
def moleculeglomos2pymatgen(moleculeglomos):
    coords, atoms=[],[]
    for iatom in moleculeglomos.atoms:
        atoms.append(iatom.s)
        coords.append([iatom.xc, iatom.yc, iatom.zc])
    moleculepymatgen=MoleculePMG(atoms,coords)
    return moleculepymatgen

#------------------------------------------------------------------------------------------
def inequivalent_finder(moleculeglomos, tolerance=0.3, eigen_tolerance=0.008, matrix_tolerance=0.1):
    molpmg=moleculeglomos2pymatgen(moleculeglomos)
    molx=PointGroupAnalyzer(molpmg, tolerance, eigen_tolerance, matrix_tolerance)
    dictionary=molx.get_equivalent_atoms()
    dict,inequiv, equiv=dictionary['eq_sets'], [], []
    for key in dict:
        listu=list(dict[key])
        listu.sort()
        inequiv.append(listu[0])
        b=listu[1:] if len(listu) > 1 else []
        equiv=equiv+b
    inequiv.sort()
    equiv.sort()
    return inequiv, equiv

#------------------------------------------------------------------------------------------
def triangle_finder(neighbors_dict,atom):
    '''
    This function finds the triangles formed by the binded atoms in a molecule
    in: neighbors_dict (Dict), the keys are position of the atom in the Molecule.atoms list, 
        and the values are a list of all their neighbors
        atom (int), the atom whose triangles are to be found
    out: triangles (list), a list of all 3 atoms that form each triangle 
    '''
    triangles = []
    # For this atom find its neighbors
    neiga = neighbors_dict[atom]
    for n in neiga:
        # For every neighbor find its own neighbors
        neigb = neighbors_dict[n]
        for n2 in neigb:
            # Check if this neighbor is in common with the original and if so save it
            if n2 in neiga:
                auxlist = [atom,n,n2]
                auxlist.sort()
                if auxlist in triangles:
                    continue
                triangles.append(auxlist)
    return triangles

#------------------------------------------------------------------------------------------
def radial_addition(original_molecule, atom_list, species):
    '''
    This function gets a list of atoms, extracts their position vector and radially, to a 
    safe distance, adds one extra atom of the required species to the molecule.
    in: original_molecule (Molecule), the molecule on which the new atoms will be added
        atoms_list (list[int]), a list of the atoms upon which the new atom will be added
        species (str), the chemical symbol of the atom that is going to be added
    out: molist_out (list[Molecule]), a list of all the new molecules
    '''
    molist_out = []
    dist_min = get_min_binding_distance(original_molecule)
    for atm in atom_list:
        cp_mol = copymol(original_molecule)
        # find the atom in the molecule
        atm_aux = cp_mol.atoms[atm]
        # and find its position vector
        p_vect = np.array([atm_aux.xc,atm_aux.yc,atm_aux.zc])
        c = 0.1
        add_vect = np.array([p_vect[0]+c,p_vect[1]+c,p_vect[2]+c])
        n_atm = Atom(species,add_vect[0],add_vect[1],add_vect[2])
        # check that the distance is at least the smallest of the molecule
        d = atom_molecule_dist(n_atm,cp_mol)
        factor = 1.0
        # and correct it if it isn't
        while d < dist_min:
            add_vect = add_vect*factor
            n_atm = Atom(species,add_vect[0],add_vect[1],add_vect[2])
            d = atom_molecule_dist(n_atm,cp_mol)
            factor = factor + 0.0001
        # add this new atom to the molecule and the molecule to the out_list
        cp_mol.add_atom(n_atm)
        cp_mol.c = 'radial'
        molist_out.append(cp_mol)
    return molist_out

#------------------------------------------------------------------------------------------
def bridge_addition(original_molecule,atom_list,neighbor_dict,species):
    '''
    This function adds one atom of the required species in between two of the neighboring atoms 
    of the ones contained in the atom_list.
    in: original_molecule (Molecule), the molecule on which the new atoms will be added
        atoms_list (list[int]), the atoms that will be used to add the new one
        neighbor_dict (Dict,{atm:[neigbors]}), a dictionary with all the neighbors of the atom_list
                      elements
        species (str), the chemical symbol of the atom that is going to be added
    out: molist_out (list[Molecule]); A list of all the new molecules
    '''
    molist_out = []
    visited_neig = []
    dist_min = get_min_binding_distance(original_molecule)
    for a in atom_list:
        aux_atm = original_molecule.atoms[a]
        # build a vector that points into that atom
        a_vect = np.array([aux_atm.xc,aux_atm.yc,aux_atm.zc])
        # and find all of its neighbors
        neig = neighbor_dict[a]
        for n in neig:
            al = [a,n]
            al.sort() 
             # for each of its neighbors, check if it has been visited already before
            if al in visited_neig:
            # if al in visited_neig or n not in atom_list:
                continue
            # if it hasn't, add it to the list and build a copy of the original molecule    
            visited_neig.append(al)
            cp_mol = copymol(original_molecule)
            # find the corresponding neighboring Atom, its radius and xyz coordinates
            aux_atmb = cp_mol.atoms[n]
            # build its position vector
            b_vect = np.array([aux_atmb.xc,aux_atmb.yc,aux_atmb.zc])
            # find the vector that points right in between the original atom and its neighbor and normalize it
            p_vect = (a_vect + b_vect) / 2
            add_vect = np.array([p_vect[0],p_vect[1],p_vect[2]])
            # and add in there the new atom
            n_atm = Atom(species,add_vect[0],add_vect[1],add_vect[2])
            d = atom_molecule_dist(n_atm,cp_mol)
            factor = 1 
            # but check if this new position is not overlapped with its surroundings
            while d < dist_min:
                add_vect = add_vect*factor
                n_atm = Atom(species,add_vect[0],add_vect[1],add_vect[2])
                d = atom_molecule_dist(n_atm,cp_mol)
                factor = factor + 0.0001
            # finally, add this atom to the molecule and add the molecule to the molist_out
            cp_mol.add_atom(n_atm)
            cp_mol.c = 'bridge'
            molist_out.append(cp_mol)
    return molist_out
#------------------------------------------------------------------------------------------
def hollow_addition(original_molecule,atom_list,neighbors_dict,species):
    '''
    This function adds one atom of the required species in the triangle formed by its neighbors.
    in: original_molecule (Molecule), the molecule on which the new atoms will be added
        atoms_list (list[int]), the atoms that will be used to add the new one
        neighbor_dict (Dict,{atm:[neigbors]}), a dictionary with all the neighbors of the atom_list
                      elements
        species (str), the chemical symbol of the atom that is going to be added
    out: molist_out (list[Molecule]), a list of all the new molecules
    '''
    molist_out = []
    visited_t = []
    for a in atom_list:
        # find the triangles composed by each atom in the list
        triangles = triangle_finder(neighbors_dict,a)
        c = 1
        for t in triangles:
            if t in visited_t:
                continue
            add_mol = copymol(original_molecule)
            # get all the three vectors and its middle point
            atm1 = add_mol.atoms[t[0]]
            v1 = np.array([atm1.xc,atm1.yc,atm1.zc])
            atm2 = add_mol.atoms[t[1]]
            v2 = np.array([atm2.xc,atm2.yc,atm2.zc])        
            atm3 = add_mol.atoms[t[2]]
            v3 = np.array([atm3.xc,atm3.yc,atm3.zc])
            mid_vect = (v1 + v2 + v3) / 3
            # to ensure that the new atom is above the triangle, find its plane
            c_vect1 = v2-v1
            c_vect2 = v3-v1
            # find a perpendicular unit vector to this plane
            c_vect = np.cross(c_vect1,c_vect2)
            n = np.linalg.norm(c_vect)
            c_vect = c_vect / n
            dot_vect = np.dot(c_vect,mid_vect)
            if dot_vect < 0:
                c_vect = (-1) * c_vect
            add_vect = (mid_vect + c_vect) * 1.1
            # add the new atom in there and add this molecule to the molist_out
            add_atom = Atom(species,add_vect[0],add_vect[1],add_vect[2])
            add_mol.add_atom(add_atom)
            add_mol.c = 'hollow'
            molist_out.append(add_mol)
            visited_t.append(t)
    return molist_out

#------------------------------------------------------------------------------------------
def focused_expansion(molin,vector):
    molout = copymol(molin)
    mag_p = np.linalg.norm(vector)
    for a in molout.atoms:
        if a.xf == 0:
            continue
        a_vect = np.array([a.xc,a.yc,a.zc])
        r = a_vect - vector
        delta = np.linalg.norm(r)
        rnorm = r/delta
        # print(delta)
        edel = np.exp(-delta)*1.5
        # print(edel,'\n')
        a_vect = a_vect+(rnorm*edel)
        a.xc, a.yc, a.zc = a_vect[0],a_vect[1],a_vect[2]
    return molout

def hollow_expanded_addition(original_molecule,atom_list,neighbors_dict,species):
    '''
    This function adds one atom of the required species in the triangle formed by its neighbors.
    in: original_molecule (Molecule), the molecule on which the new atoms will be added
        atoms_list (list[int]), the atoms that will be used to add the new one
        neighbor_dict (Dict,{atm:[neigbors]}), a dictionary with all the neighbors of the atom_list
                      elements
        species (str), the chemical symbol of the atom that is going to be added
    out: molist_out (list[Molecule]), a list of all the new molecules
    '''
    molist_out = []
    visited_t = []
    for a in atom_list:
        # find the triangles composed by each atom in the list
        triangles = triangle_finder(neighbors_dict,a)
        for t in triangles:
            if t in visited_t:
                continue
            add_mol = copymol(original_molecule)
            # get all the three vectors and its middle point
            atm1 = add_mol.atoms[t[0]]
            v1 = np.array([atm1.xc,atm1.yc,atm1.zc])
            # ta1 = Atom('H',v1[0],v1[1],v1[2])
            # add_mol.add_atom(ta1)
            atm2 = add_mol.atoms[t[1]]
            v2 = np.array([atm2.xc,atm2.yc,atm2.zc])
            # ta2 = Atom('H',v2[0],v2[1],v2[2])
            # add_mol.add_atom(ta2)        
            atm3 = add_mol.atoms[t[2]]
            v3 = np.array([atm3.xc,atm3.yc,atm3.zc])
            # ta3 = Atom('H',v3[0],v3[1],v3[2])
            # add_mol.add_atom(ta3)
            mid_vect = (v1 + v2 + v3) / 3
            # ta = Atom('H',mid_vect[0],mid_vect[1],mid_vect[2])
            # add_mol.add_atom(ta)
            # add_atom = Atom(species,mid_vect[0],mid_vect[1],mid_vect[2])
            # add_atom.xf = 0
            # add_mol.add_atom(add_atom)
            exp_mol = focused_expansion(add_mol,mid_vect)
            add_atom = Atom(species,mid_vect[0],mid_vect[1],mid_vect[2])
            add_atom.xf = 0
            exp_mol.add_atom(add_atom)
            # add_atom = Atom('H',mid_vect[0],mid_vect[1],mid_vect[2])
            # exp_mol.add_atom(add_atom)
            exp_mol.c = 'interstitial'
            molist_out.append(exp_mol)
            visited_t.append(t)
    return molist_out

#------------------------------------------------------------------------------------------
def make_many_rand(molist_in, species):
    '''
    This function masters all the above to add a new atom to each of the inequivalent atoms that belong
    to the outter layer of the molecule. It uses two types of addition: a.- Directly above of the atoms and 
    b.- in between two adjacent neighbors.
    in: molist_in (list[Molecule]); all the molecules that will receive an extra atom
        species (str); the species of the atom to be additioned
    out: molist_out (list[Molecule]);  a list that contains all the molecules with the extra atom 
    '''
    molist_out = []
    req_molist = molist_in[:num_iso]
    for imol in req_molist:
        # for each Molecule in the list, make a copy and translate it to the center of mass
        org_mol = copymol(imol)
        org_mol = translate_to_cm(org_mol)
        # get its adjacency matrix
        mtx = conectmx(org_mol)
        # find a dict with all the neighbors
        all_neighbors = neighbor_finder(mtx)
        if sym_rest == 'on':
            # find the inequivalent atoms of the molecule
            inequivalent, equivalent = inequivalent_finder(org_mol)
            elegible_atoms = inequivalent.copy()
        else:
            elegible_atoms = []
            [elegible_atoms.append(i) for i in range(len(org_mol.atoms))]
        # make the additions of each type
        for at in l_add_type:
            if at == 'a':
                list_a = radial_addition(org_mol,elegible_atoms,species)
                molist_out.extend(list_a)
                list_b = bridge_addition(org_mol,elegible_atoms,all_neighbors,species)
                molist_out.extend(list_b)
                list_c = hollow_addition(org_mol,elegible_atoms,all_neighbors,species)
                molist_out.extend(list_c)
                list_d = hollow_expanded_addition(org_mol,elegible_atoms,all_neighbors,species)            
                molist_out.extend(list_d)
                break
            elif at == 'r':
                list_a = radial_addition(org_mol,elegible_atoms,species)
                molist_out.extend(list_a)
            elif at == 'b':
                list_b = bridge_addition(org_mol,elegible_atoms,all_neighbors,species)
                molist_out.extend(list_b)
            elif at == 'h':
                list_c = hollow_addition(org_mol,elegible_atoms,all_neighbors,species)
                molist_out.extend(list_c)
            elif at == 'i':
                list_d = hollow_expanded_addition(org_mol,elegible_atoms,all_neighbors,species)            
                molist_out.extend(list_d)
    fopen = open(log_file,'a')
    for i,m in enumerate(molist_out):
        si=str(i+1).zfill(5)
        print('stage1%s from %s growtype %s' %(si,m.i,m.c), file=fopen)
    fopen.close()
    return molist_out

#------------------------------------------------------------------------------------------
def run_sample():
    from utils.libmoleculas import readxyzs, writexyzs
    m = readxyzs('mo20.xyz')
    l = make_many_rand(m,'C')
    writexyzs(l,'tests_initial.xyz')
# run_sample()

#------------------------------------------------------------------------------------------
# def bfs_algorithm(moleculein, neighbors_dict):
#     '''
#     This algorithm labels iteratively the neighbors of the atoms which is closer to the
#     molecule's center of mass and returns a list with all the 'n' neighbors labeled.
#     in: moleculein (Molecule); The molecule on which the algorithm will be performed
#     out: out_list (list); list of list like [[closest atom],[first neighbors],...,[last neighbors]] 
#     '''
#     moleculeout = copymol(moleculein)
#     l = closest_atm_to_cm(moleculeout)
#     tup_list = []
#     # this list will have tuples like (node,label)
#     for i in neighbors_dict:
#         if i == l:
#             aux_tup = (i,0)
#         else:
#             aux_tup = (i,-1)
#         tup_list.append(aux_tup)
#     stop_crit, step_id = 0, 0
#     out_list, labeled = [[l]], [l]
#     while stop_crit == 0:
#         # finding the atoms whose neighbors will be found and labeled
#         atms_of_interest = []
#         for tup in tup_list:
#             if tup[1] == step_id:
#                 atms_of_interest.append(tup[0])
#         # finding all the non-repeated neighbors of all the previous atoms
#         non_repeated_neigh = []
#         for atm in atms_of_interest:
#             neighbors = neighbors_dict[atm]
#             for n in neighbors:
#                 if n in labeled:
#                     continue
#                 non_repeated_neigh.append(n)
#         non_repeated_neigh = list(set(non_repeated_neigh))
#         # adding these non repeated neighbors to the counted list
#         labeled.extend(non_repeated_neigh)
#         if non_repeated_neigh:
#             stop_crit = 0
#         else:
#             break
#         # adding the non-repeated neighbors to the out_list
#         out_list.append(non_repeated_neigh)
#         # updating the tuple list with the new label
#         step_id = step_id + 1
#         for lab in labeled:
#             if tup_list[lab][1] < 0:
#                 tup_list[lab] = (lab,step_id)
#     return out_list


#------------------------------------------------------------------------------------------
# def random_addition(original_molecule, atom_list, species):
#     '''
#     This function gets one molecule and randomly adds one atom into one of the existing ones 
#     to a safe distance, taking into consideration that this new atom is not one of the inequivalent ones.
#     in: original_molecule (Molecule), the molecule that will receive the addition
#         atom_list (list[int]); The equivalent atoms that will receive a random addition
#         species (str); The species of the atom to be added
#     out: molist_out (list[Molecule]); a list that contains all the molecules with the extra atom    
#     '''
#     from utils.general import randunitvector
#     molist_out = []
#     n_list = []
#     d_ref = get_min_binding_distance(original_molecule)
#     ra = get_covalent_radius(species)
#     for a in atom_list:
#         # create a copy of the original molecule, find the corresponding Atom and its position vector
#         cp_mol = copymol(original_molecule)
#         atm = cp_mol.atoms[a]
#         rb = get_covalent_radius(atm.s)
#         p_vect = np.array([atm.xc, atm.yc, atm.zc])
#         # build a random unit vector
#         n_vect = randunitvector() / (ra + rb)
#         d = 0
#         fact = 1.0
#         while d < d_ref:
#             fact = fact + 0.0001
#             n_vect = fact * n_vect
#             add_vect = p_vect + n_vect
#             try_atm = Atom(species,add_vect[0],add_vect[1],add_vect[2])
#             d = atom_molecule_dist(try_atm,cp_mol)
#         # finally add this atom to the molecule and the molecule to the molist_out
#         cp_mol.add_atom(try_atm)
#         cp_mol.c = 'random'
#         molist_out.append(cp_mol)
#     return molist_out