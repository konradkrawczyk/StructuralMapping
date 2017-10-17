#!/usr/bin/python
# -*- coding: utf-8 -*-
#

import sys
import os.path
import sqlite3 as sql
import gzip
import errno
import subprocess
from array import array
from glob import glob
from StringIO import StringIO

import numpy


from prosci.util.protein import Pdb, ResidueList, Protein, Residue, Atom
from prosci.util.ali import Ali
from prosci.shell import Params

from prosci.util.geohash import GeometricHash
from prosci.loops.loopmodel import ANCHOR_LENGTH, add_oxygens, get_dihedral_class, is_structure_ok, describe_anchors, make_contact_gh, find_contacts, get_contact_class


MIN_CHAIN_LENGTH = 40


# Taken from:
# http://stackoverflow.com/questions/600268/mkdir-p-functionality-in-python
def mkdir_p(path):
  try:
    os.makedirs(path)
  except OSError as exc:
    if exc.errno == errno.EEXIST:
      pass
    else:
      raise


def gzip_open(source, mode="rb"):
  if mode.startswith("r"):
    p = subprocess.Popen(["zcat", source], stdout = subprocess.PIPE)
    f = StringIO(p.communicate()[0])
    if p.returncode != 0:
      raise IOError()
  else:
    f = gzip.open(source, mode)
  return f


def write_structure(f, residues):
  for r in residues:
    f.write(">\n")
    f.write(str(r))


def read_structure(f):
  if isinstance(f, str):
    fnamestub = os.path.splitext(os.path.basename(f))[0]
    try:
      fil = gzip_open(f)
      fnamestub = os.path.splitext(fnamestub)[0]
    except IOError:
      fil = open(f)
    try:
      x = read_structure(fil)
      x.code = fnamestub
    finally:
      fil.close()
    return x
  
  residues=[]
  r=[]
  i=-1
  for line in f:
    if line.startswith('>'):
      if r:
        residues.append(Residue(r))
        r=[]
    else:
      r.append(Atom(line))
  if r:
    residues.append(Residue(r))
  return ResidueList(residues)


class FreadDB(object):
  
  
  def __init__(self, dbdir):
    self.dbdir = os.path.abspath(dbdir)
    self.strucdir = os.path.join(self.dbdir, "structures")
    #self.origstrucdir = os.path.join(self.dbdir, "original_structures")
    
    
    self.loop_lengths = []
    self.terminal_cutoff = 5
    #self.max_b = 100
    
    #self.structures_added = []
    
    self.verbose = False
    
    try:
      self.load_db_options()
    except IOError:
      pass
  
  
  
  
  def load_db_options(self, fname=None):
    if fname is None:
      fname = os.path.join(self.dbdir, "options.txt")
    
    f = open(fname)
    for line in f:
      key, value = line.split(None, 2)
      if key == "loop_lengths":
        loop_lengths = value.split(",")
        for i,x in enumerate(loop_lengths):
          loop_lengths[i] = int(x)
        self.loop_lengths = loop_lengths
      elif key in self.__dict__:
        t = type(self.__dict__[key])
        self.__dict__[key] = t(value)
      else:
        raise ParsingError("Unknown key in FREAD database options file: %s = %s" % (key, value))
    
    self.loop_lengths = self.get_indexed_loop_lengths()
  
  
  def save_db_options(self, fname=None):
    if fname is None:
      fname = os.path.join(self.dbdir, "options.txt")
    
    mkdir_p(self.dbdir)
    f = open(fname, "wb")
    f.write("terminal_cutoff %s\n" % (self.terminal_cutoff))
    #f.write("max_b %s\n" % ())
    f.close()
  
  
  
  
  def get_index_filename(self, loop_length):
    return os.path.join(self.dbdir, "length%d.sqlite"%loop_length)
  
  
  def get_index_filenames(self):
    return glob(os.path.join(self.dbdir, "length*.sqlite"))
  
  
  def get_indexed_loop_lengths(self):
    lengths = []
    for fname in self.get_index_filenames():
      L = os.path.splitext(os.path.basename(fname))[0].replace("length", "")
      lengths.append(int(L))
    lengths.sort()
    return lengths
  
  
  def get_structure_path_stub(self, code):
    return os.path.join(self.strucdir, code[1:3], code)
  
  
  def get_protein_complex(self, code):
    code = code[:4]
    p = Protein([])
    files = glob(self.get_structure_path_stub(code)+"*.gz")
    files.sort()
    for fname in files:
      f = gzip_open(fname)
      try:
        rl = read_structure(f)
        if p and p[-1].chain == rl.chain:
          p[-1].extend(rl)
        else:
          p.append(rl)
      finally:
        f.close()
    p.code = code
    return p
  
  
  def add_structure(self, pdb, code=None, chain_codes=None, skip_existing=False):
    p = Protein(pdb)
    if code is None:
      code = p.code
    code = code[:4].lower() + code[4:]
    
    if self.verbose:
      print "Adding structure: %s" % (code)
    
    for chain in p:
      add_oxygens(chain)
    
    
    #dest_orig = os.path.join(self.origstrucdir, code[1:3], code+".atm.gz")
    #mkdir_p(os.path.dirname(dest_orig))
    #
    #f = gzip.open(dest_orig, "wb")
    #try:
    #  f.write(str(p))
    #finally:
    #  f.close()
    
    
    mkdir_p(os.path.join(self.strucdir, code[1:3]))
    
    for chain in p:
      # Save chains we're not interested in as ligands
      if len(chain) >= MIN_CHAIN_LENGTH and (not chain_codes or chain.chain in chain_codes):
        suffix = "atm"
      else:
        #suffix = "lig"
        continue
      
      dest_model = self.get_structure_path_stub(code+chain.chain)+"."+suffix+".gz"
      if skip_existing and os.path.exists(dest_model):
        return
      
      f = gzip.open(dest_model, "wb", 6)
      try:
        write_structure(f, chain)
      finally:
        f.close()
    
    """
    for chain in p.ligands:
      dest_model = self.get_structure_path_stub(code+chain.chain)+".lig.gz"
      if skip_existing and os.path.exists(dest_model):
        return
      
      f = gzip.open(dest_model, "wb", 6)
      try:
        write_structure(f, chain)
      finally:
        f.close()
    """
    
    # Inter-residue contacts
    gh = make_contact_gh(p)
    
    for chain in p:
      # Skip chains we're not interested in
      if len(chain) < MIN_CHAIN_LENGTH or (chain_codes and (chain.chain not in chain_codes)):
        continue
      
      annotation_file_destination = self.get_structure_path_stub(code+chain.chain)+".tem"
      
      sequence = chain.get_seq()
      assert len(sequence) == len(chain)
      
      # Dihedral angles. Also serves as a mask - residues marked '?' are excluded from the database.
      dihedrals = array('c', "?")
      for i in xrange(1, len(chain)-1):
        if not is_structure_ok(chain, i-1, i+2) or sequence[i] == 'X':
          x = '?'
        else:
          try:
            x = str(get_dihedral_class(chain, i))
          except ValueError:
            # Math domain error, if computing angles between overlayed atoms
            x = '?'
        dihedrals.append(x)
      dihedrals.append("?")
      dihedrals = dihedrals.tostring()
      assert len(dihedrals) == len(chain)
      
      contacts = array('c')
      for i, r in enumerate(chain):
        contact_atoms = find_contacts(gh, r, maxdist=6.0)
        x = get_contact_class(r.chain, contact_atoms)
        contacts.append(x)
      contacts = contacts.tostring()
      
      f = open(annotation_file_destination, "wb")
      try:
        f.write(">%s\nsequence\n%s\n>%s\nFREAD dihedral class\n%s\n>%s\nFREAD contact class\n%s\n"%(code+chain.chain, sequence, code+chain.chain, dihedrals, code+chain.chain, contacts))
      finally:
        f.close()
  
  
  def build_index(self, loop_length):
    terminal_cutoff = self.terminal_cutoff
    #max_b = self.max_b
    
    
    mkdir_p(self.dbdir)
    
    structure_files = glob(os.path.join(self.strucdir, "*", "*.atm.gz"))
    structure_files.sort()
    
    if not structure_files:
      return
    
    dbfile = self.get_index_filename(loop_length) + ".new"
    if os.path.exists(dbfile):
      os.remove(dbfile)
    
    conn = sql.connect(dbfile)
    cursor = conn.cursor()
    
    #cursor.execute('''CREATE TABLE loops (dihedral text, sequence text, pdbcode text, start int, casep1 real, casep2 real, casep3 real, casep4 real, contact text)''')
    cursor.execute('''CREATE TABLE loops (dihedral text, sequence text, pdbcode text, start int, casep1 real, casep2 real, casep3 real, casep4 real)''')
    
    done = False
    proteins_not_annotated = []
    
    try:
      for ifile, fname in enumerate(structure_files):
        annotation_file = fname[:-len(".atm.gz")] + ".tem"
        
        print "%d/%d : %s" % (ifile+1, len(structure_files)-1, fname)
        
        residues = read_structure(fname)
        
        if len(residues) < 5:
          print "Length of residues < 5"
          continue
        
        strucname = residues.code
        sequence = residues.get_seq()
        assert len(sequence) == len(residues)
        
        try:
          annotation = Ali(annotation_file)
        except IOError:
          print "Annotation not found for entry %s" % (residues.code)
          proteins_not_annotated.append(residues.code)
          continue
          #p = self.get_protein_complex(strucname)
          #self.add_structure(p)
          #annotation = Ali(annotation_file)
          
        
        annotseq = annotation[0].master.seq
        dihedrals = annotation[0]["FREAD dihedral class"].seq
        contacts = annotation[0]["FREAD contact class"].seq
        assert sequence == annotseq
        
        
        totallen = loop_length + 2*ANCHOR_LENGTH
        
        # Loop over all possible fragments within the protein
        for start in xrange(terminal_cutoff, len(residues)-terminal_cutoff-totallen):
          end = start + totallen
          loopstart = start+ANCHOR_LENGTH
          loopend = end-ANCHOR_LENGTH
          
          loopdihedrals = dihedrals[loopstart:loopend]
          
          # Skip bad loops
          if '?' in dihedrals[start:end]:
            continue
          
          loopseq = sequence[loopstart:loopend]
          
          loopcontacts = contacts[loopstart:loopend]
          
  #         struc_coords = struc_coords_list[start:end]
          
  #         b_factor = 0.0
  #         for i in xrange(loopstart, loopend):
  #           b_factor += residues[i].CA.b
  #         b_factor /= loop_length
          
          # Descsribe anchor
          anchor_description, transform = describe_anchors(residues[start:loopstart], residues[loopend:end], loop_length)
          
          casep1 = anchor_description.C1[0]
          casep2 = numpy.linalg.norm(anchor_description.C2)
          casep3 = numpy.linalg.norm(anchor_description.N1 - anchor_description.C1)
          casep4 = numpy.linalg.norm(anchor_description.N1 - anchor_description.C2)
          
          #cursor.execute("INSERT INTO loops VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)", (loopdihedrals, loopseq, strucname, start, casep1, casep2, casep3, casep4, loopcontacts))
          cursor.execute("INSERT INTO loops VALUES (?, ?, ?, ?, ?, ?, ?, ?)", (loopdihedrals, loopseq, strucname, start, casep1, casep2, casep3, casep4))
        
      done = True
    finally:
      if done:
        if self.verbose:
          print "Done creating data for loop length %d. Now sorting database for optimal run speed." % (loop_length)
        conn.executescript("""
        CREATE TABLE loops2 AS SELECT * FROM loops ORDER BY dihedral;
        DROP TABLE loops;
        ALTER TABLE loops2 RENAME TO loops;
        VACUUM;
        """)
        conn.commit()
        conn.close()
        if os.path.exists(dbfile[:-4]):
          os.remove(dbfile[:-4])
        os.rename(dbfile, dbfile[:-4])
        if self.verbose:
          print "Done with loop length %d." % (loop_length)
      else:
        conn.close()
      
      if proteins_not_annotated:
        print ""
        print "These proteins lack annotation files:"
        for code in proteins_not_annotated:
          print code
