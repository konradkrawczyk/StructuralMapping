#!/usr/bin/python
# -*- coding: utf-8 -*-
#

import os
import re
import gzip
import sqlite3 as sql
import subprocess
from StringIO import StringIO
from math import degrees, atan2, pi, ceil # , sin, cos
from array import array
from glob import glob

import numpy

from prosci.util.pdb3d import dihedral_angle, dist, sqrdist, vectors2rotation_matrix, axis_angle2rotation_matrix
from prosci.util.residue import Residue, ResidueList
from prosci.util.residue import is_residue_consecutive as is_consecutive
from prosci.util.pdb import Atom, residueCode
from prosci.util.ali import Ali
from prosci.util.geohash import GeometricHash
from prosci.common import join

# Used implicitly in this file:
#
#from prosci.util.residue import Residue



# Database anchor hash granularity
#ANCHOR_GRANULARITY = 1.0

ANCHOR_LENGTH = 2

#STRUC_HASH_GRANULARITY = 0.3




class Data:
  pass


def get_max_loop_length(loop_length):
  return 3.8 * (1+loop_length)

# def is_consecutive(a, b, minlen=1.0, maxlen=2.0):
# #def is_consecutive(residues, a, b):
#     if a.C is None or b.N is None:
#       #raise RuntimeError("Cannot reliably test if two residues are consecutive, purely based on CA distance! Need N and C atoms!")
#       return False
#     
#     #d=dist(r.CA, q.CA)
#     #if (d < 3.6 or d > 4.0):
#     d = dist(a.C, b.N)
#     if (d < minlen or d > maxlen):
#       return False
#     return True


def get_loop_structure(dbdir, strucname, start, length):
  
  def parse(f):
    residues=[]
    r=[]
    i=-1
    for line in f:
      if line.startswith('>'):
        i+=1
        if i >= start+length:
          residues.append(Residue(r))
          break
        if r:
          residues.append(Residue(r))
          r=[]
      elif i>= start:
        r.append(Atom(line))
    return residues
  
  source = os.path.join(dbdir, "structures", strucname[1:3], strucname)
  f = None
  try:
    #p = subprocess.Popen(["zcat", source + ".atm.gz"], stdout = subprocess.PIPE)
    #f = StringIO(p.communicate()[0])
    #if p.returncode != 0:
    #  raise IOError()
    f = gzip.open(source + ".atm.gz")
    residues = parse(f)
  except IOError:
    f = open(source + ".atm")
    residues = parse(f)
  finally:
    if f is not None:
      f.close()
  
  return ResidueList(residues)


# # TODO: Doesn't work yet. Needs to get anchor residues too!
# def get_loop_structure_from_decoyname(dbdir, name):
#   "Parse the name of a loop decoy and get its structure from the database"
#   
#   strucname, startres, length = name.split("_")
#   startinscode = ""
#   if startres[-1:].isalpha():
#     startinscode = startres[-1:]
#     startres = startres[:-1]
#   startres = int(startres)
#   length = int(length) + 2*ANCHOR_LENGTH
#   
#   source = os.path.join(dbdir, "structures", strucname[1:3], strucname)
#   residues=[]
#   r=[]
#   i=-1
#   f = gzip.open(source + ".atm.gz")
#   skip = False
#   try:
#     for line in f:
#       if line.startswith('>'):
#         if i < 0:
#           skip = False
#           continue
#         i+=1
#         if i >= length:
#           residues.append(Residue(r))
#           break
#         if r:
#           residues.append(Residue(r))
#           r=[]
#       elif i < 0:
#         if not skip:
#           a = Atom(line)
#           if a.ires < startres or (a.ires == startres and  a.inscode < startinscode):
#             skip = True
#           else:
#             r.append(Atom(line))
#             i = 0
#       else:
#         r.append(Atom(line))
#     # No need to check if r has not yet been added, as loops never
#     # reach the end of the protein.
#   finally:
#     f.close()
#   return ResidueList(residues)


def get_native_contacts(dbdir, strucname, start, length):
  source = os.path.join(dbdir, "structures", strucname[1:3], strucname)
  a = Ali(source + ".tem")
  dihedrals = a[0]["FREAD dihedral class"].seq[start:start+length]
  contacts = a[0]["FREAD contact class"].seq[start:start+length]
  assert "?" not in dihedrals
  return contacts


def add_oxygens(residues, start=0, end=None, dO=1.23, force=False):
  if end is None:
    end = len(residues) - 1
  for i in xrange(start, end):
    # add missing mainchain oxygens, making them planar
    # with CA, C and the next residue's N
    #
    r = residues[i]
    if r is None:
      continue
    
    if r.O is None:
      q = residues[i+1]
      if q is None or q.N is None or r.C is None:
        continue
      
      if not force and not is_consecutive(r, q):
        continue
      
      O = r.C.copy()
      O.atom = "O"
      if r.C.element:
        O.element = "O"
      O.iatom = r.C.iatom+1

      vO_1 = r.C.xyz - r.CA.xyz
      vO_2 = r.C.xyz - q.N.xyz
      vO_3 = (vO_1 / numpy.linalg.norm(vO_1)) + (vO_2 / numpy.linalg.norm(vO_2))
      O.xyz = r.C.xyz + (vO_3 / numpy.linalg.norm(vO_3)) * dO

      r.O = O


def is_structure_ok(residues, start, end, max_b=float("inf")):
  for i in xrange(start, end):
    r1=residues[i]
    if r1.N is None or r1.CA is None or r1.C is None or r1.O is None or \
       r1.N.b > max_b or r1.CA.b > max_b or r1.C.b > max_b or r1.O.b > max_b:
      return False
  
  r1=residues[start]  
  for i in xrange(start+1, end):
    r2=residues[i]
    if (r2.CA.ires != r1.CA.ires + 1):
      return False
    if not is_consecutive(r1, r2):
      return False
    # d=dist(r1.CA, r2.CA)
    # if (d < 2.8 or d > 4.0):
    #   return False
    r1=r2
  
  return True



# dihed_phi   = zip([-180, -110,-110,-180,-180,-180,  20,   0,  0],
#                   [   0,    0,   0,-110,-110,   0, 140, 180,180])
# dihed_psi     = zip([ -90,  100,-180, 100,-180,  45, -40,-180, 80],
#                   [  45,  180, -90, 180, -90, 100,  80, -40,180])
# dihed_class =     [   0,    1,   1,   2,   2,   3,   4,   5,  5]

# dihed_phi   = zip([-180, -110,-110,-180,-180,-180,  20,   0],
#                   [   0,    0,   0,-110,-110,   0, 140, 180])
# dihed_psi   = zip([ -90,  100,-180, 100,-180,  80, -40,-180],
#                   [  80,  180, -90, 180, -90, 100,  80, 180])
# dihed_class =     [   0,    1,   1,   2,   2,   3,   4,   5]
# #                     A     B    B    C    C    D    E    F

dihed_phi   = zip([-180, -110,-110,-180,-180,-180,  20,   0],
                  [   0,    0,   0,-110,-110,   0, 140, 180])
dihed_psi   = zip([ -90,  100,-180, 100,-180,  45, -40,-180],
                  [  45,  180, -90, 180, -90, 100,  80, 180])
dihed_class =     [   0,    1,   1,   2,   2,   3,   4,   5]
#                     A     B    B    C    C    D    E    F


def get_dihedral_angle(residues, i):
  r1 = residues[i-1]
  r2 = residues[i]
  r3 = residues[i+1]
  phi = dihedral_angle(r2.N.xyz  - r1.C.xyz, r2.CA.xyz - r2.N.xyz, r2.C.xyz - r2.CA.xyz)
  psi = dihedral_angle(r2.CA.xyz - r2.N.xyz, r2.C.xyz - r2.CA.xyz, r3.N.xyz - r2.C.xyz)
  return (phi, psi)


def angle2class(phi, psi):
  for i in xrange(9):
    minphi,maxphi = dihed_phi[i]
    minpsi,maxpsi = dihed_psi[i]
    if phi <= maxphi and phi >= minphi and psi <= maxpsi and psi >= minpsi:
      return dihed_class[i]
  raise IllegalStateError("illegal dihedral angles : %f, %f" % (phi, psi))


def get_dihedral_class(residues, i):
  phi, psi = get_dihedral_angle(residues, i)
  return angle2class(phi, psi)


def describe_anchors(anchor_n, anchor_c, loop_length):
  max_actual_length = get_max_loop_length(loop_length)
  
  n1, n2 = anchor_n[-2:]
  c1, c2 = anchor_c[:2]
  
  # get CA coordinates of anchors
  coords = numpy.array([n1.CA.xyz, n2.CA.xyz, c1.CA.xyz, c2.CA.xyz, n2.O.xyz, c1.O.xyz])
  
  origin = numpy.copy(coords[1])
  
  # translate second CA of N-anchor onto origin
  for i,c in enumerate(coords):
    coords[i] -= origin
  
  #print coords
  
  # rotate first CA of N-anchor onto Y axis in the negative direction, such that
  # the actual loop goes off in the positive direction along Y axis
  #
#   if numpy.round(coords[0][0], 4) != 0 or numpy.round(coords[0][2], 4) != 0:
#     rotmat1 = vectors2rotation_matrix(coords[0], numpy.array([0.0,-1.0,0.0]))
#     for i,c in enumerate(coords):
#       coords[i] = numpy.dot(c, rotmat1)
#   else:
#     rotmat1 = None
  
#   print "Experimental rotation matrix:"
#   print vectors2rotation_matrix(numpy.array([-1.0,0.0,0.0]), numpy.array([1.0,0.0,0.0]))
#   print "Experimental rotation matrix2:"
#   exprotmat = axis_angle2rotation_matrix(numpy.array([0.0,1.0,0.0]), pi)
#   print numpy.dot(numpy.array([-1.0,0.0,0.0]), exprotmat)
#   raise RuntimeError()
  
  # rotate first CA of C-anchor onto X axis in the positive direction
  #
  if numpy.round(coords[2][1], 4) != 0 or numpy.round(coords[2][2], 4) != 0:
    rotmat1 = vectors2rotation_matrix(coords[2], numpy.array([1.0,0.0,0.0]))
    for i,c in enumerate(coords):
      coords[i] = numpy.dot(c, rotmat1)
  elif coords[2][0] < 0:
    rotmat1 = axis_angle2rotation_matrix((0.0,1.0,0.0), pi)
    for i,c in enumerate(coords):
      coords[i] = numpy.dot(c, rotmat1)
  else:
    rotmat1 = None
  
  
  #print coords
  
  # rotate first CA of N-anchor onto X-Y plane, in the negative Y direction, so that the loop always goes off into the positive Y direction.
  if numpy.round(coords[0][2], 4) != 0:
    rotmat2 = vectors2rotation_matrix(numpy.array([0.0, coords[0][1], coords[0][2]]), numpy.array([0.0,-1.0,0.0]))
    for i,c in enumerate(coords):
      coords[i] = numpy.dot(c, rotmat2)
  elif coords[0][1] > 0:
    rotmat2 = axis_angle2rotation_matrix((1.0,0.0,0.0), pi)
    for i,c in enumerate(coords):
      coords[i] = numpy.dot(c, rotmat2)
  else:
    rotmat2 = None
  
  #print coords
  
  coords = coords.round(4)
  
  try:
    assert (coords[1] == (0,0,0)).all()  # second N atom should be the origin
    assert coords[2][1] == 0   # first C atom should be on the X axis
    assert coords[2][2] == 0   # first C atom should be on the X axis
    assert coords[2][0] >  0   # first C atom should be be in the positive X range
    assert coords[0][2] == 0   # first N atom should be on the X-Y plane
    assert coords[0][1] <= 0   # first N atom should be on in the negative Y range
    assert not numpy.isnan(coords).any()   # must not have NaN in coords
  except:
    print "Anchor coordinates:"
    for c in coords:
      print c
      print
    raise
  
  d = Data()
  
  # save coordinates
  d.N1 = coords[0]
  # N2 = [0, 0, 0]   # N2 is always on the origin
  d.C1 = coords[2]
  d.C2 = coords[3]
  
  # compute loop stretch (can just use X coords of first CA of C-anchor, as it's on the X axis)
  stretch = min(1.0, coords[2][0] / max_actual_length)
  
  
  # get angle between N anchor and Y axis
  Nanchor_angle = atan2(-coords[0][0], -coords[0][1])
  
  # get rotation angle of N anchor's O atom around the N anchor's own axis
  #Nanchor_rotmat = numpy.array([[cos(Nanchor_angle), -sin(Nanchor_angle), 0], [sin(Nanchor_angle), cos(Nanchor_angle), 0], [0, 0, 1]]); # matrix rotates N anchor onto Y axis
  Nanchor_rotmat = axis_angle2rotation_matrix((0.0,0.0,1.0), Nanchor_angle)
  Nanchor_O = numpy.dot(coords[4], Nanchor_rotmat)
  Nanchor_twist = atan2(Nanchor_O[0], Nanchor_O[2])
  
  
  # get X-Y-Z unit vector from second C atom to first C atom
  Canchor_tilt = coords[2] - coords[3]
  Canchor_tilt /= numpy.linalg.norm(Canchor_tilt) # make unit vector
  Canchor_tilt = numpy.dot(Canchor_tilt, Nanchor_rotmat)
  
  
  # get rotation angle of C anchor's O atom around the C anchor's own axis
  Canchor_rotmat = vectors2rotation_matrix(coords[3] - coords[2], numpy.array([0.0,-1.0,0.0]))
  Canchor_O = numpy.dot(coords[5] - coords[2], Canchor_rotmat)
  Canchor_twist = atan2(Canchor_O[0], Canchor_O[2])
  
  d.stretch = stretch
  d.Nangle = degrees(Nanchor_angle)%180
  d.Ctilt = Canchor_tilt
  d.Ntwist = degrees(Nanchor_twist)%360
  d.Ctwist = degrees(Canchor_twist)%360
  
  return d, (origin, rotmat1, rotmat2)



# def hash_loop_anchor_numeric(anchor_desc):
#   return (
#     int(anchor_desc.N1[0] / ANCHOR_GRANULARITY),
#     int(anchor_desc.N1[1] / ANCHOR_GRANULARITY),
#     int(anchor_desc.C1[0] / ANCHOR_GRANULARITY),
#     int(anchor_desc.C2[0] / ANCHOR_GRANULARITY),
#     int(anchor_desc.C2[1] / ANCHOR_GRANULARITY),
#     int(anchor_desc.C2[2] / ANCHOR_GRANULARITY)
#     )
  

# def hash_loop_anchor_string(anchor_desc):
#   return hash_numeric_to_string(hash_loop_anchor_numeric(anchor_desc))



# def hash_numeric_to_string(numeric):
#   return join(".", numeric)

# def hash_string_to_numeric(string):
#   numeric=[]
#   x = string.split('.')
#   for y in x:
#     numeric.append(int(y))
#   return tuple(numeric)



# def hash_structure(coords, transform, bin_size=STRUC_HASH_GRANULARITY):
#   hash=array('c')
#   for i in xrange(len(coords)):
#     c = transform_xyz(coords[i], transform)
#     hash.extend(str(int((c[0] + bin_size/2) / bin_size)))
#     hash.append('.')
#     hash.extend(str(int((c[1] + bin_size/2) / bin_size)))
#     hash.append('.')
#     hash.extend(str(int((c[2] + bin_size/2) / bin_size)))
#     hash.append('.')
#   hash.pop()
#   return hash.tostring()


def transform_xyz(xyz, transform):
  c = xyz - transform[0]
  if transform[1] is not None:
    c = numpy.dot(c, transform[1])
  if transform[2] is not None:
    c = numpy.dot(c, transform[2])
  return c




###########################################################


# def score_decoy(seq, esss_vectors, ascii2index):
#   # Calculate ESSS
#   score=0
#   for i,x in enumerate(seq):
#     score += esss_vectors[i][ascii2index[ord(x)]]
#   return score


def get_db_path(dbdir, loop_length):
  return os.path.join(dbdir, "length%d.sqlite"%loop_length)


def get_db_file_dict(dbdir):
  prefix = os.path.join(dbdir, "length")
  suffix = ".sqlite"
  d={}
  for fname in glob(prefix+"*"+suffix):
    L = int(fname[len(prefix):-len(suffix)])
    d[L] = fname
  return d


def is_loop_length_in_database(dbdir, loop_length):
  return os.path.exists(get_db_path(dbdir, loop_length))


def get_min_loop_length_in_database(dbdir):
  dbfiles = glob(os.path.join(dbdir, "length*.sqlite"))
  if not dbfiles:
    raise RuntimeError("No database found in directory: "+dbdir)
  minlen=1000000
  for fname in dbfiles:
    length = os.path.splitext(os.path.basename(fname))[0].replace("length", "")
    minlen = min(minlen, int(length))
  return minlen



def iterate_database(dbdir, loop_length, esst, anchor, sequence, strict_rmsd_cutoff, score_cutoff):
  q = anchor
  dbfile = get_db_path(dbdir, loop_length)
  
  # If we cannot possibly meet the score_cutoff with a perfect match to this sequence, reduce cut-off to allow perfect matches to be found
  score_cutoff = min(score_cutoff, esst.get_perfect_score(sequence))
  
  #print "Score cut-off:", score_cutoff
  
  tables = esst.tables
  ascii2index = esst.ascii2index
  seqmap = tuple([ascii2index[ord(s)] for s in sequence])
  
  assert os.path.exists(dbfile), "Database file missing: "+dbfile
  conn = sql.connect(dbfile)
  
  prevdihed=["", [None]*len(seqmap)]
  def score_sequence(loopseq, dihed):
      if prevdihed[0] != dihed:
        prevdihed[0] = dihed
        for i,x in enumerate(dihed):
          prevdihed[1][i] = tables[int(x)][seqmap[i]]
      score=0
      for i,x in enumerate(loopseq):
        score += prevdihed[1][i][ascii2index[ord(x)]]  # Speed-optimised version
        #score += tables[int(dihed[i])][seqmap[i]][ascii2index[ord(x)]]  # Naive version
      return score
  
  
  # DEBUGGING: See if this function throws an exception
  score_sequence(sequence, "0"*len(sequence))
  
  
  conn.create_function("score_seq", 2, score_sequence)
  
  casep1 = q.C1[0]
  casep2 = numpy.linalg.norm(q.C2)
  casep3 = numpy.linalg.norm(q.N1 - q.C1)
  casep4 = numpy.linalg.norm(q.N1 - q.C2)
  min_casep1 = max(0, casep1 - strict_rmsd_cutoff*3)
  max_casep1 =        casep1 + strict_rmsd_cutoff*3
  min_casep2 = max(0, casep2 - strict_rmsd_cutoff*3)
  max_casep2 =        casep2 + strict_rmsd_cutoff*3
  min_casep3 = max(0, casep3 - strict_rmsd_cutoff*3)
  max_casep3 =        casep3 + strict_rmsd_cutoff*3
  min_casep4 = max(0, casep4 - strict_rmsd_cutoff*3)
  max_casep4 =        casep4 + strict_rmsd_cutoff*3
  
  try:
#      contacts,
    for row in conn.execute("""
    SELECT
      score_seq(sequence, dihedral) AS score,
      dihedral,
      sequence,
      pdbcode,
      start,
      ((casep1-?)*(casep1-?) + (casep2-?)*(casep2-?) + (casep3-?)*(casep3-?) + (casep4-?)*(casep4-?)) / 16 AS internal_rmsd_sq
    FROM
      loops
    WHERE
      casep1 BETWEEN ? AND ?
      AND casep2 BETWEEN ? AND ?
      AND casep3 BETWEEN ? AND ?
      AND casep4 BETWEEN ? AND ?
      AND internal_rmsd_sq <= ?
      AND score >= ?
    """, (casep1, casep1, casep2, casep2, casep3, casep3, casep4, casep4, min_casep1, max_casep1, min_casep2, max_casep2, min_casep3, max_casep3, min_casep4, max_casep4, strict_rmsd_cutoff**2, score_cutoff)):
      d = Data()
      d.score, d.dihedrals, d.seq, d.struc, d.start = row[:5]
      #d.score, d.dihedrals, d.seq, d.struc, d.start, d.native_contacts = row[:6]
      d.internal_rmsd = numpy.sqrt(row[-1])
#       d.casep1_diff = abs(row[5] - casep1)
#       d.casep2_diff = abs(row[6] - casep2)
#       d.casep3_diff = abs(row[7] - casep3)
#       d.casep4_diff = abs(row[8] - casep4)
      yield d
  finally:
    conn.close()



###########################################################



def relabel_loop(decoy_residues, sequence, startnum=0, endnum=0, endinscode="", chain="A", prevatom=None, nextatom=None):
  assert len(decoy_residues) == len(sequence)
  
  if prevatom:
    chain = prevatom.chain
    startnum = prevatom.ires + 1
  
  if nextatom:
    if prevatom:
      assert prevatom.chain == nextatom.chain
    else:
      chain = nextatom.chain
    
    if nextatom.inscode > "A":
      endinscode = chr(ord(nextatom.inscode) - 1)
    
    if nextatom.inscode:
      endnum = nextatom.ires
    else:
      endnum = nextatom.ires - 1
  
  inscode = "A"
  for i,r in enumerate(decoy_residues):
    r.set_type(residueCode(sequence[i]))
    if startnum:
      if startnum + i <= endnum:
        for a in r:
          a.ires = startnum + i
          a.inscode = ""
          a.chain = chain
      else:
        for a in r:
          a.ires = endnum
          a.inscode = inscode
          a.chain = chain
        inscode = chr(ord(inscode)+1)


vdw_radii = {'C':1.70, 'N':1.55,'O':1.52,'S':1.80,'P':1.80,'H':1.20}

def is_clash(gh, atoms, decoy, vdw_factor=0.7):
  if vdw_factor <= 0:
    return False
  
  for r in decoy:
    for a in r:
      if a.atom not in ("N", "CA", "C", "O", "CB"):
        continue
      nbr = gh.get_neighbours(a.xyz, 1, max_dist=3.6*vdw_factor)
      if nbr:
        ix, dst = nbr[0]
        try:
          r1 = vdw_radii[a.atom[0]]
        except KeyError:
          r1 = 1.7
        try:
          r2 = vdw_radii[atoms[ix].atom[0]]
        except KeyError:
          r2 = 1.7
        if dst < (r1 + r2) * vdw_factor:
          return True
  return False


def calculate_rmsd(a, b):
    assert len(a) == len(b)
    total = 0.0
    count = 0
    for i, ra in enumerate(a):
      rb = b[i]
      if ra.N and rb.N:
        total += numpy.linalg.norm(ra.N.xyz - rb.N.xyz)**2
        count += 1
      if ra.CA and rb.CA:
        total += numpy.linalg.norm(ra.CA.xyz - rb.CA.xyz)**2
        count += 1
      if ra.C and rb.C:
        total += numpy.linalg.norm(ra.C.xyz - rb.C.xyz)**2
        count += 1
      if ra.O and rb.O:
        total += numpy.linalg.norm(ra.O.xyz - rb.O.xyz)**2
        count += 1
    return numpy.sqrt(total / count)


def calculate_rmsdCA(a, b):
    assert len(a) == len(b)
    total = 0.0
    for i, ra in enumerate(a):
      rb = b[i]
      total += numpy.linalg.norm(ra.CA.xyz - rb.CA.xyz)**2
    return numpy.sqrt(total / len(a))


def norm_angle(a):
  if a > 180:
    a = 360 - a
  elif a <= -180:
    a += 360
  return a



def compare_hit_anchor(query, hit):
  d = Data()
  # N anchor twist difference
  d.ntwist_diff = abs(norm_angle(query.Ntwist - hit.Ntwist))
  # C anchor twist difference
  d.ctwist_diff = abs(norm_angle(query.Ctwist - hit.Ctwist))
  return d



def find_contacts(gh, residue, maxdist=6.0):
  atoms = gh.atoms
  
  # gh is a GeometricHash object
  # atoms is a list of Atom objects that corresponds to the co-ordinates used to initialise gh
  # residue is a Residue object to find contacts for
  if maxdist <= 0:
    return []
  
  a = residue.CB
  if a is None:
    a = residue.CA
  if a is None:
    return []
  
  results = []
  nbrs = gh.get_neighbours(a.xyz, 999, max_dist=maxdist)
  for n in nbrs:
    ix, dst = n
    a = atoms[ix]
    
    element = a.atom[0:1]
    if element.isdigit():
      element = a.atom[1:2]
    if element == "H" or element.isdigit() or not element:
      continue
    
    if a.chain != residue.chain or abs(a.ires - residue.ires) > 10:
      results.append(a)
  
  return results


def find_contacts_simple(residues, i, maxdist=6.0):
  maxsqrdist = maxdist ** 2
  a1 = residues[i].CB
  if a1 is None:
    a1 = residues[i].CA
  if a1 is None:
    return []
  
  results = []
  #for residue in residues[:i-10] + residues[i+11:]:
  for residue in residues:
    a2 = residue.CB
    if a2 is None:
      a2 = residue.CA
    if a2 is None:
      continue
    
    if a1.chain != a2.chain or abs(a1.ires - a2.ires) > 10:
      sqd = sqrdist(a1, a2)
      if sqd < maxsqrdist:
        results.append(a2)
  
  return results


def get_contact_class(thischain, contact_atoms):
  x = 0 # no contacts
  for a in contact_atoms:
    if a.chain != thischain:
      x |= 1 # inter-chain contact
    else:
      x |= 2 # intra-chain contact
  
  # x = 0 ... no contacts
  # x = 1 ... only inter-chain contacts
  # x = 2 ... only intra-chain contacts
  # x = 3 ... both intra and inter-chain contacts
  return str(x)



def make_contact_gh(protein):
  coords = []
  gh_atoms = []
  for chain in protein:
    for r in chain:
      a = r.CB
      if a is None:
        a = r.CA
      if a is None:
        continue
      coords.append(a.xyz)
      gh_atoms.append(a)
  
  for chain in protein.ligands:
    for r in chain:
      for a in r:
        element = a.atom[0:1]
        if element.isdigit():
          element = a.atom[1:2]
        if element == "H" or element.isdigit() or not element:
          continue
        coords.append(a.xyz)
        gh_atoms.append(a)
  
  gh = GeometricHash(numpy.array(coords))
  gh.atoms = gh_atoms
  
  return gh
