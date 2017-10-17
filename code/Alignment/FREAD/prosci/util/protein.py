#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Protein class represents a protein, as defined by a PDB file.
#
# Author: Sebastian Kelm
# Created: 19/04/2012
#
from prosci.util.pdb import Pdb, Atom
from prosci.util.residue import ResidueList, Residue


class Protein(list):
  "A list of ResidueList objects, each of which represents one PDB chain."
  
  def __init__(self, pdb, code=None):
    "Takes a Pdb object or a list of ResidueLists as the first argument. Otherwise passes the argument to Pdb() first, which can deal with filenames. Second argument is an optional short decription (usually a PDB code)."
    self.code = ""
    self.ligands = []
    
    if not pdb:
      pass
    elif not isinstance(pdb, Pdb):
      if isinstance(pdb[0], ResidueList):
        for reslist in pdb:
          self.append(reslist)
        self.code = pdb[0].code
      else:
        pdb = Pdb(pdb)
    
    if isinstance(pdb, Pdb):
      self.extend(ResidueList(pdb).split_chains())
      self.code = pdb.code
      self.ligands = Protein(ResidueList(pdb.ligands).split_chains())
    
    if code is not None:
      self.code = code
  
  def __repr__(self):
    return "Protein(%s)"%(list.__repr__(self))
  
  def __str__(self):
    out=""
    for r in self:
      if r:
        out += str(r)
    return out
  
  def __getslice__(self, start=None, end=None):
    return Protein(list.__getslice__(self, start, end))
  
  
  def get_chain(self, chaincode):
    "Get a particular chain by its chain code"
    for chain in self:
      if chain.chain == chaincode:
        return chain
    for chain in self.ligands:
      if chain.chain == chaincode:
        return chain
    raise KeyError("Chain %s not found in Protein object with code '%s'" % (chaincode, self.code))
  
  
  def get_seq(self):
    "Returns a FASTA string representation of the sequences."
    seq = ""
    for chain in self:
      if chain is not None:
        seq += chain.get_seq()+"/"
    return seq[:-1]
  
  
  def get_coords(self, atomfilter=lambda atom: True, ligands=False):
    "Get a numpy array of coordinates"
    coords=[]
    for chain in self:
      for r in chain:
        if r is not None:
          for a in r:
            if atomfilter(a):
              coords.append(a.xyz)
    if ligands:
      for chain in self.ligands:
        for r in chain:
          if r is not None:
            for a in r:
              if atomfilter(a):
                coords.append(a.xyz)
    return numpy.array(coords)
  
  
  def to_pdb(self, atomfilter=lambda atom: True):
    "Returns a Pdb object containing the Atom objects within this object"
    p = Pdb(self.code, [])
    for chain in self:
      for res in chain:
        for atm in res:
          if atomfilter(atm):
            p.data.append(atm)
    for chain in self.ligands:
      for res in chain:
        for atm in res:
          if atomfilter(atm):
            p.ligands.append(atm)
    return p
  
  
  def iter_atoms(self, ligands=False):
    for chain in self:
      for r in chain:
        for a in r:
          yield a
    if ligands:
      for chain in self.ligands:
        for r in chain:
          for a in r:
            yield a
  
  
  def to_residuelist(self, ligands=False):
    "Returns a single ResidueList object"
    rl = ResidueList([], self.code)
    for chain in self:
      rl.extend(chain)
    if ligands:
      for chain in self.ligands:
        rl.extend(chain)
    return rl
  
  
  def renumber_atoms(self):
    "Renumbers all atoms, starting from 1 upwards"
    i = 1
    for chain in self:
      for r in chain:
        for a in r:
          a.iatom = i
          i += 1
    for chain in self.ligands:
      for r in chain:
        for a in r:
          a.iatom = i
          i += 1
  
  
  def remove_residues(self, func, ligands=False):
    i=0
    while i < len(self):
      self[i].remove_residues(func)
      if not self[i]:
        del self[i]
      else:
        i += 1
    
    i=0
    while i < len(self.ligands):
      self.ligands[i].remove_residues(func)
      if not self.ligands[i]:
        del self.ligands[i]
      else:
        i += 1
  
  
  def iter_backbone(self):
    "Iterate over all backbone atoms that aren't None"
    for chain in self:
      for r in chain:
        for a in r.iter_backbone():
          yield a
  
  
  def find_residue(self, ires, inscode="", chain=""):
    if isinstance(ires, Residue):
      inscode = ires.inscode
      chain = ires.chain
      ires = ires.ires
    
    for ichain, chainresidues in enumerate(self):
      if chainresidues.chain == chain:
        return (ichain, chainresidues.find_residue(ires, inscode, chain))
    
    return (-1, -1)

