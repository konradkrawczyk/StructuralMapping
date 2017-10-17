#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
import math
import numpy

from prosci.common import *
from prosci.util.pdb import Atom, Pdb, intersectAtomTypes
from prosci.util.seq import find_aligned_residues
from prosci.util.gaps import length_ungapped
from prosci.util.residue import ResidueList



class NearestNeighbourSearcher:
    def __init__(self, atoms):
        from prosci.util.geohash import GeometricHash
        self.atoms = atoms
        self.gh = GeometricHash(get_coords(atoms))
    
    def get_neighbours(self, query, num=1, one_per_res=True, res_dist=1, max_dist=None):
        """
        NearestNeighbourSearcher.get_neighbours(query, num=1, one_per_res=True, res_dist=1) : returns a list of atom indeces of the nearest neighbours of the query
        
        query is an Atom object
        num is the number of neighbours to return
        one_per_res=True means we should only return one atom for each neibouring residue
        one_per_res=False means multiple hits to the same residue are allowed
        res_dist is the minimum distance, in the sequence, that a neighbour has to be away from query.
        If query is from a different structure, res_dist should be set to 0.
        """
        assert isinstance(query, Atom)
        
        query_coords = (query.x, query.y, query.z)
        
        def residue_filter(results, candidate):
          cand_atom = self.atoms[candidate[1]]
          if abs(query.ires - cand_atom.ires) < res_dist:
            return False
          if one_per_res:
            for sqd, i in results:
              r_atom = self.atoms[i]
              if r_atom.ires == cand_atom.ires and r_atom.chain == cand_atom.chain and r_atom.inscode == cand_atom.inscode:
                  return False
          return True
        
        return self.gh.get_neighbours(query_coords, num, max_dist, residue_filter)


#class PdbAligner:
    #def __init__(self, translation_vector, rotation_matrix, rmsd):
        #self.translation_vector = translation_vector
        #self.rotation_matrix = rotation_matrix
        #self.rmsd = rmsd
    
    
    #def transform(self, struc, invert=False):
        #### Apply changes to all atoms ###
        
        #coords = get_coords(struc)
        
        #translate(coords, translation_vector)
        #rotate(coords, rotation_matrix)
      
        #set_coords(struc, coords)


#def build_aligner(struc1, struc2, subset1, subset2):
    
    #### Get the translation vector, rotation matrix and RMSD value. For this, we only use the subset of aligned, corresponding atoms. ###
    
    #coords_sub1 = get_coords(subset1)
    #coords_sub2 = get_coords(subset2)
    
    #translation_vector1 = centroid(coords_sub1)
    #translation_vector2 = centroid(coords_sub2)
    
    #translate(coords_sub1, - translation_vector1)
    #translate(coords_sub2, - translation_vector2)
    
    #rotation_matrix, rmsd = optimal_rotation(coords_sub2, coords_sub1)
    
    #return PdbAligner(translation_vector2 - translation_vector1, rotation_matrix, rmsd)

def rmsd_static(coordA, coordB, atom_types=("N", "CA", "C", "O")):
    "Calculate RMSD between two corresponding set of residues, without moving anything. Expects Pdb objects."
    
    if atom_types is not None:
      coordA = coordA.get_atoms_by_type(atom_types)
      coordB = coordB.get_atoms_by_type(atom_types)
    
    residuesA = coordA.residues()
    residuesB = coordB.residues()
    assert len(residuesA)==len(residuesB)
    
    sumsqrdist=0.0
    count=0
    
    for i in xrange(len(residuesA)):
      resA, resB = intersectAtomTypes(residuesA[i], residuesB[i])
      for a, b in zip(resA, resB):
        assert a.atom == b.atom
        sumsqrdist+=(a.x-b.x)**2 + (a.y-b.y)**2 + (a.z-b.z)**2
        count+=1
    
    if count==0:
      return float("inf")
    
    return math.sqrt(float(sumsqrdist)/float(count))


def rmsd_backbone(residuesA, residuesB):
    "Calculate backbone RMSD between two corresponding set of residues, without moving anything. Expects ResidueList objects."
    
    assert len(residuesA)==len(residuesB)
    
    sumsqrdist=0.0
    count=0
    
    for resA, resB in zip(residuesA, residuesB):
      for a, b in ((resA.N, resB.N), (resA.CA, resB.CA), (resA.C, resB.C), (resA.O, resB.O)):
        if a is None or b is None:
          continue
        sumsqrdist+=(a.x-b.x)**2 + (a.y-b.y)**2 + (a.z-b.z)**2
        count+=1
    
    if count==0:
      return float("inf")
    
    return math.sqrt(float(sumsqrdist)/float(count))



def sqrdist(a, b):
    return (a.x-b.x)**2 + (a.y-b.y)**2 + (a.z-b.z)**2

def dist(a, b):
    return math.sqrt((a.x-b.x)**2 + (a.y-b.y)**2 + (a.z-b.z)**2)


def gdt_ts(coordA, coordB, length=None):
    "Calculate GDT between two corresponding set of residues, without moving anything."
    return average(gdt(coordA, coordB, (1.0, 2.0, 4.0, 8.0), length))


def gdt_ha(coordA, coordB, length=None):
    "Calculate GDT between two corresponding set of residues, without moving anything."
    return average(gdt(coordA, coordB, (0.5, 1.0, 2.0, 4.0), length))


def gdt(coordA, coordB, distance_cutoffs, length=None):
    "Calculate GDT between two corresponding sets of residues, without moving anything."
    
    assert tuple(distance_cutoffs) == tuple(sorted(set(distance_cutoffs))), str(distance_cutoffs)
    
    sqcutoffs = tuple([x**2 for x in distance_cutoffs])
    
    #atom_types=("CA",)
    #coordA = coordA.get_atoms_by_type(atom_types)
    #coordB = coordB.get_atoms_by_type(atom_types)
    
    residuesA = ResidueList(coordA)
    residuesB = ResidueList(coordB)
    assert len(residuesA)==len(residuesB)
    
    if length is None:
      length = len(residuesA)
    
    counts=[0] * len(distance_cutoffs)
    drange = range(len(distance_cutoffs))
    drange.reverse() # we want to start by checking the largest cut-offs, then move on to the smaller ones
    
    for a, b in zip(residuesA, residuesB):
      if a.CA is None or b.CA is None:
        continue
      d = sqrdist(a.CA, b.CA)
      for j in drange:
        if d >= sqcutoffs[j]:
          break
        counts[j] += 1
    
    for j in drange:
      counts[j] = float(counts[j]) / length
    
    return counts



def find_closest(pdb1_all, pdb2_all):
    if not len(pdb1_all) or not len(pdb2_all):
      raise ArgumentError("Pdbs must be non-empty. Lengths given: %d, %d" % (len(pdb1_all), len(pdb2_all)))
    pdb1 = pdb1_all.get_atoms_by_type("CA")
    pdb2 = pdb2_all.get_atoms_by_type("CA")
    assert len(pdb1) == len(pdb2)
    m=dist(pdb1[0], pdb2[0])
    ix = 0
    for i in xrange(1, len(pdb1)):
      d = dist(pdb1[i], pdb2[i])
      if d < m:
        m = d
        ix = i
    return pdb1[ix], pdb2[ix], m

    



def get_coords(atoms):
    coords = []
    for a in atoms:
        c = numpy.array([a.x, a.y, a.z])
        coords.append(c)
    return numpy.array(coords)
    

def set_coords(pdb, coords):
    assert coords.shape == (len(pdb), 3)
    for i,a in enumerate(pdb):
        a.x, a.y, a.z = coords[i]


def centroid(coords):
    return numpy.mean(coords, axis=0)


def translate(coords, vector):
    assert vector.shape == (3,)
    for c in coords:
        c += vector


def optimal_rotation(crds1, crds2):
    assert crds1.shape[1] == 3
    assert crds2.shape[1] == 3
    assert crds1.shape == crds2.shape, "optimal_rotation(): Shapes differ: %s  %s" %(crds1.shape, crds2.shape)

    n_vec = len(crds1)
    correlation_matrix = numpy.dot(numpy.transpose(crds1), crds2)
    v, s, w = numpy.linalg.svd(correlation_matrix)
    is_reflection = (numpy.linalg.det(v) * numpy.linalg.det(w)) < 0.0
    if is_reflection:
        s[-1] = - s[-1]
        v[:,-1] = -v[:,-1]
    
    E0 = sum(sum(crds1 * crds1)) + \
         sum(sum(crds2 * crds2))
    rmsd_sq = (E0 - 2.0*sum(s)) / float(n_vec)
    rmsd_sq = max([rmsd_sq, 0.0])
    rmsd = numpy.sqrt(rmsd_sq)
    return numpy.transpose(numpy.dot(v,w)), rmsd


#~ def vectors2rotation_matrix(old_direction, new_direction):
    #~ from math import cos, sin
    #~ unitvector_old = old_direction / numpy.linalg.norm(old_direction)
    #~ unitvector_new = new_direction / numpy.linalg.norm(new_direction)
    #~ 
    #~ xangle = angle([0.0, old_direction[1], old_direction[2]], [0.0, new_direction[1], new_direction[2]])
    #~ yangle = angle([old_direction[0], 0.0, old_direction[2]], [new_direction[0], 0.0, new_direction[2]])
    #~ zangle = angle([old_direction[0], old_direction[1], 0.0], [new_direction[0], new_direction[1], 0.0])
    #~ 
    #~ xrot = numpy.matrix([[1.0, 0.0, 0.0], [0.0, cos(xangle), -sin(xangle)], [0.0, sin(xangle), cos(xangle)]])
    #~ yrot = numpy.matrix([[cos(yangle), 0.0, sin(yangle)], [0.0, 1.0, 0.0], [-sin(yangle), 0.0, cos(yangle)]])
    #~ zrot = numpy.matrix([[cos(zangle), -sin(zangle), 0.0], [sin(zangle), cos(zangle), 0.0], [0.0, 0.0, 1.0]])
    #~ 
    #~ return xrot * yrot #* zrot

def skew_symmetric(vector):
  return numpy.array([[0.0, -vector[2], vector[1]],[vector[2], 0.0, -vector[0]],[-vector[1], vector[0], 0.0]])


def vectors2rotation_matrix(old_direction, new_direction):
    from math import cos, sin
    
    unitvector_old = old_direction / numpy.linalg.norm(old_direction)
    unitvector_new = new_direction / numpy.linalg.norm(new_direction)
    
    rotaxis = numpy.cross(unitvector_old, unitvector_new)
    rotaxis /= numpy.linalg.norm(rotaxis)
    rotangle = -angle(unitvector_old, unitvector_new)
    
    UxU = numpy.outer(rotaxis, rotaxis)
    rotmat = UxU + cos(rotangle) * ( numpy.array([[1,0,0],[0,1,0],[0,0,1]]) - UxU ) + sin(rotangle) * skew_symmetric(rotaxis)
    
    return rotmat


def axis_angle2rotation_matrix(axis, rotangle):
    from math import cos, sin
    
    rotaxis = axis / numpy.linalg.norm(axis)
    
    UxU = numpy.outer(rotaxis, rotaxis)
    rotmat = UxU + cos(rotangle) * ( numpy.array([[1,0,0],[0,1,0],[0,0,1]]) - UxU ) + sin(rotangle) * skew_symmetric(rotaxis)
    
    return rotmat


def rotate(coords, rotation_matrix):
    #assert rotation_matrix.shape == ...
    for i,c in enumerate(coords):
        coords[i] = numpy.dot(c, rotation_matrix)


def fit_line(atoms, atom_types=("N", "CA", "C", "O")):
    if isinstance(atoms, Pdb):
      a = atoms.get_atoms_by_type(atom_types)
      crds = get_coords(a)
    else:
      crds = atoms
    
    mean = centroid(crds)
    v, s, w = numpy.linalg.svd(crds - mean)
    is_reflection = (numpy.linalg.det(v) * numpy.linalg.det(w)) < 0.0
    if is_reflection:
      #s[-1] = - s[-1]
      v[:,-1] = -v[:,-1]
    
    best_fit_unit_vector = w[0] ## direction of first principle component
    length = numpy.linalg.norm(crds[-1] - crds[0])
    
    # Check if fitted line is going in the right direction. If not, turn it around.
    if numpy.linalg.norm(crds[-1] - crds[0] - best_fit_unit_vector) > length:
      best_fit_unit_vector = -best_fit_unit_vector
    
    linepoints = best_fit_unit_vector * [[-length/2], [length/2]]
    linepoints += mean
    
    return best_fit_unit_vector, linepoints
    
    

def rmsd_numpy(coords1, coords2):
    translation_vector1 = centroid(coords1)
    translation_vector2 = centroid(coords2)
    
    translate(coords1, - translation_vector1)
    translate(coords2, - translation_vector2)
    
    rotation_matrix, rmsd = optimal_rotation(coords2, coords1)
    return rmsd


def rmsd(struc1, struc2, subset1=None, subset2=None, modify=False):
    ### Get the translation vector, rotation matrix and RMSD value. For this, we only use the subset of aligned, corresponding atoms. ###
    
    if subset1 is None:
      subset1 = struc1
    if subset2 is None:
      subset2 = struc2

    coords_sub1 = get_coords(subset1)
    coords_sub2 = get_coords(subset2)
    
    translation_vector1 = centroid(coords_sub1)
    translation_vector2 = centroid(coords_sub2)
    
    translate(coords_sub1, - translation_vector1)
    translate(coords_sub2, - translation_vector2)
    
    rotation_matrix, rmsd = optimal_rotation(coords_sub2, coords_sub1)
    
    if not modify:
      return rmsd
    
    
    ### Apply changes to all atoms ###
    
    coords1 = get_coords(struc1)
    
    #print "Ex-Centroid1:     ", translation_vector1
    #print "Ex-Centroid2:     ", translation_vector2
    
    translate(coords1, - translation_vector1)
    rotate(coords1, rotation_matrix)
    translate(coords1, translation_vector2)

    set_coords(struc1, coords1)

    #print "Centroid1:        ", centroid(get_coords(subset1))
    #print "Centroid2:        ", centroid(get_coords(subset2))
    #print
    
    return rmsd


def rotate_translate(struc1, struc2, subset1, subset2):
    return rmsd(struc1, struc2, subset1, subset2, True)


def reduceToAlignable(struc1_allchains, struc2_allchains, seq1, seq2, subset1=None, subset2=None, atom_types=("N", "CA", "C", "O"), modify_structures=True):
    assert type(struc1_allchains) == type(struc2_allchains)
    assert type(subset1) == type(subset2)
    assert isinstance(struc1_allchains, Pdb)
    assert None == subset1 or isinstance(subset1, Pdb)
    
    if subset1 is None:
        subset1 = struc1_allchains
    if subset2 is None:
        subset2 = struc2_allchains
    if seq1 is None:
        seq1 = subset1.get_seq()
    if seq2 is None:
        seq2 = subset2.get_seq()

    if not (seq1 and seq2):
      raise ValueError("Need to have non-empty sequence to align proteins:\nseq1:%s\nseq2:%s\n"%(seq1, seq2))
    ## if structure has more than 1 chain, only use the first one
    #if subset1.chaincount() > 1:
    #  subset1 = subset1.get_first_chain()
    #if subset2.chaincount() > 1:
    #  subset2 = subset2.get_first_chain()
    
    
    subset1_resbounds = subset1.residue_boundaries()
    subset2_resbounds = subset2.residue_boundaries()
    
    # residue count, according to the structure data
    pdb1_rescount = len(subset1_resbounds)
    pdb2_rescount = len(subset2_resbounds)
    
    #print deGappify(seq1)
    #print deGappify(subset1.get_seq())
    #print deGappify(seq2)
    #print deGappify(subset2.get_seq())
    
    # Make sure the residue counts coincide in sequence and structure data
    #
    assert length_ungapped(seq1) == pdb1_rescount, "length_ungapped(seq1) = %d, pdb1_rescount = %d" % (length_ungapped(seq1), pdb1_rescount)
    assert length_ungapped(seq2) == pdb2_rescount, "length_ungapped(seq2) = %d, pdb2_rescount = %d" % (length_ungapped(seq2), pdb2_rescount)
    
    # Get the residue indeces of aligned residues
    #
    aligned_indeces1, aligned_indeces2 = find_aligned_residues(seq1, seq2)
    assert len(aligned_indeces1) == len(aligned_indeces2)
    
    if not aligned_indeces1:
        raise ParsingError("No aligned residues?")
    
    
    #
    # Get the subset of backbone atoms corresponding to the aligned residues
    #
    #
    
    #subset1_CA = subset1.get_CA()
    #subset2_CA = subset2.get_CA()

    aligned_pdb1 = Pdb(subset1, [])
    aligned_pdb2 = Pdb(subset2, [])
    
    for ix1, ix2 in zip(aligned_indeces1, aligned_indeces2):
      #residue1 = subset1.get_residue(subset1_CA[ix1])
      #residue2 = subset2.get_residue(subset2_CA[ix2])
      #residue1 = residue1.get_atoms_by_type(atom_types)
      #residue2 = residue2.get_atoms_by_type(atom_types)
      residue1 = subset1.get_atoms(slice=subset1_resbounds[ix1], atom_types=atom_types)
      residue2 = subset2.get_atoms(slice=subset2_resbounds[ix2], atom_types=atom_types)
      
      if len(atom_types) != len(residue1) or len(atom_types) != len(residue2):
        residue1, residue2 = intersectAtomTypes(residue1, residue2)
      
      assert len(residue1) == len(residue2)
      aligned_pdb1.append_atoms(residue1)
      aligned_pdb2.append_atoms(residue2)
    
    
    
    
    assert len(aligned_pdb1) == len(aligned_pdb2)
    
    return aligned_pdb1, aligned_pdb2


def superimpose(struc1_allchains, struc2_allchains, seq1=None, seq2=None, subset1=None, subset2=None, align_atoms=("N", "CA", "C", "O"), modify=True, modify_structures=True):
    "Rotates and translates struc1 onto struc2"
    if not modify or not modify_structures:
      modify = modify_structures = False
    try:
      aligned_pdb1, aligned_pdb2 = reduceToAlignable(struc1_allchains, struc2_allchains, seq1, seq2, subset1, subset2, atom_types=align_atoms)
      rmsd_value = rmsd(struc1_allchains, struc2_allchains, aligned_pdb1, aligned_pdb2, modify=modify)
    except:
      sys.stderr.write("Error superimposing structures: %s and %s"%(struc1_allchains.code, struc2_allchains.code))
      raise
    
    return rmsd_value


def compare_structures(struc1_allchains, struc2_allchains, seq1, seq2, subset1=None, subset2=None, realign=True, modify=False, normalise_by=None):
    coordA, coordB = reduceToAlignable(struc1_allchains, struc2_allchains, seq1, seq2, subset1, subset2, atom_types=("N", "CA", "C", "O"))
    
    scores={}
    
    if realign:
      #~ if not modify:
      coordA = coordA.deep_copy()
      coordB = coordB.deep_copy()
      rotate_translate(coordA, coordB, coordA, coordB) # this aligns the structure
    
    scores["backbone_rmsd"] = rmsd_static(coordA, coordB, atom_types=("N", "CA", "C", "O"))
    scores["ca_rmsd"] = rmsd_static(coordA, coordB, atom_types=("CA",))
    
    try:
      i = normalise_by.__iter__
    except AttributeError:
      gdtscores         = gdt(coordA, coordB, (0.5, 1.0, 2.0, 4.0, 8.0), length=normalise_by)
      for cutoff, score in zip((0.5, 1.0, 2.0, 4.0, 8.0), gdtscores):
        scores["gdt_%.1f"%cutoff] = score
      scores["gdt_ha"]  = average(gdtscores[:-1])
      scores["gdt_ts"]  = average(gdtscores[1:])
    else:
      for i,norm in enumerate(normalise_by):
        gdtscores         = gdt(coordA, coordB, (0.5, 1.0, 2.0, 4.0, 8.0), length=norm)
        scores["gdt_ha_%d"%i]  = average(gdtscores[:-1])
        scores["gdt_ts_%d"%i]  = average(gdtscores[1:])
        for cutoff, score in zip((0.5, 1.0, 2.0, 4.0, 8.0), gdtscores):
          scores["gdt_%.1f_%d"%(cutoff, i)] = score
    
    return scores




def averageCoordinates(fragScaffold, fragPrediction, invertWeights=False):
    """Averages the coordinates of the two Pdb arguments.
    
    Modifies the first object to reflect the result. By default, the first object is assumed to be part of the N-terminal chain, the second is part of the C-terminal chain. This can be reversed by setting invertWeights=True."""
    
    from prosci.util.pdb import intersectAtomTypes

    def transformCoord(P, S, D, l):
      # P = coordinate of prediction
      # S = coordinate of scaffold
      # D = distance (in residues) from (loop+anchor)-fragment end
      # l = anchor length
      return 1.0/(l+1) * (D*P + (l+1-D)*S)
    
    resS = fragScaffold.residues()
    resP = fragPrediction.residues()
    
    assert len(resS) == len(resP), "Anchor lengths differ in template and loop?? : %d, %d" % (len(resS), len(resP))
    
    for i,rS in enumerate(resS):
      rP = resP[i]
      
      if not rS.has_atom_type('CA') or not rP.has_atom_type('CA'):
        raise IllegalStateError("Residue must contain a CA atom:\nrS:\n"+str(rS)+"rP:\n"+str(rP))
        
      rS, rP = intersectAtomTypes(rS, rP)
        #raise IllegalStateError("Residue must contain only/all main-chain atoms (and optionally a CB):\nrS:\n"+str(rS)+"rP:\n"+str(rP))
      
      #print "##################################################################"
      #print "Residue", i
      #print
      #print "Scaffold:"
      #print rS
      #print "Prediction:"
      #print rP
      
      CA_vector=None
      CA_b=999
      for j in xrange(min(len(rS), len(rP))):
        aS = rS[j]
        aP = rP[j]
        assert aS.atom == aP.atom, "Atom order differs between residues. %s != %s" % (aS.atom, aP.atom)
        
        if invertWeights:
          D = len(resS)-i
        else:
          D = i+1
        
        if aS.atom == "CA":
          assert None == CA_vector, "More than one CA atom in residue??"
          CA_vector=[aS.x, aS.y, aS.z]
        
        #print "aS.x = transformCoord(%.3f, %.3f, %d, %d)" % (aP.x, aS.x, dist, len(resS))
        aS.x = transformCoord(aP.x, aS.x, D, len(resS))
        aS.y = transformCoord(aP.y, aS.y, D, len(resS))
        aS.z = transformCoord(aP.z, aS.z, D, len(resS))
        aS.b = transformCoord(aP.b, aS.b, D, len(resS)) # average B factors too
        
        if aS.atom == "CA":
          CA_vector[0] = aS.x - CA_vector[0]
          CA_vector[1] = aS.y - CA_vector[1]
          CA_vector[2] = aS.z - CA_vector[2]
          CA_b = aS.b
      
      if len(rS) > len(rS.parent):
        #assert rS[-1].atom == "CB", "Last atom in residue expected to be a CB atom. Instead '%s' found." % (rS[-1].atom)
        assert CA_vector is not None, "No CA atom in residue. Cannot use missing CA movement vector for unpaired atom transformation."
        for atm in rS.parent:
          if atm not in rS:
            atm.x += CA_vector[0]
            atm.y += CA_vector[1]
            atm.z += CA_vector[2]
            atm.b  = CA_b
      
      #print "Averaged:"
      #print rS
      #print "##################################################################"


def angle(crds1, crds2):
  return math.acos(numpy.dot(crds1, crds2)/(numpy.linalg.norm(crds1)*numpy.linalg.norm(crds2)))
  

def dihedral_angle(v1, v2, v3):
  n1 = numpy.cross(v1,v2)
  n2 = numpy.cross(v2,v3)
  
  phi = math.degrees(angle(n1,n2))
  
  if numpy.dot(v1,n2) < 0:
    return -phi
  return phi
  
