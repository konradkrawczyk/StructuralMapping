#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

import sys
import os
import time
import itertools

from heapq import heappush, heappop

# Compatibility with Python < 2.6
#
try:
  from heapq import heappushpop
except ImportError:
  def heappushpop(heap, item):
    heappush(heap, item)
    return heappop(heap)

import numpy

from prosci.util.protein import Pdb, ResidueList, Protein
from prosci.util.pdb3d import dist
from prosci.util.geohash import GeometricHash

from prosci.loops import ccd
from prosci.loops.loopmodel import ANCHOR_LENGTH, get_max_loop_length, describe_anchors, iterate_database, get_loop_structure, add_oxygens, is_loop_length_in_database, get_min_loop_length_in_database, relabel_loop, is_clash, calculate_rmsd, is_consecutive, make_contact_gh, find_contacts, get_contact_class, find_contacts_simple, get_native_contacts
from prosci.loops.esst import load_esst


# Define consecutive residues very loosely, for purposes of checking the framework
def is_framework_consecutive(a, b):
  # Set max distance for consecutive C and N atoms to 2.8 Angstroms
  # ... which is very loose
  #
  return is_consecutive(a, b, maxlen_sq=7.84)



def _heapadd(heap, item, maxitems):
  if len(heap) < maxitems:
    heappush(heap, item)
  else:
    heappushpop(heap, item)



class FreadError(RuntimeError):
  pass

class BadInputError(FreadError):
  pass

class UnsupportedLoopLengthError(BadInputError):
  pass

class LoopStretchError(BadInputError):
  pass

class NonConsecutiveAnchorError(BadInputError):
  pass


class Fread(object):
  def __init__(self, db=None, subst_table_path=None, score_cutoff=25, open_rmsd_cutoff=1.0, closed_rmsd_cutoff=0.3, vdw_factor=0.7, close=True, verbose=False, errstream=sys.stderr, meld=True, max_melding_rmsd=1.0, nostruc=False, mutate=False, ccd_target_rmsd=0.15, ccd_iterations=5000, max_decoys=100, first_decoys=0, extension_size=sys.maxint, extension_minimum=0, calculate_contacts=False, contact_distance=6.0, contact_identity=0.8):
    self.score_cutoff = score_cutoff
    self.open_rmsd_cutoff = open_rmsd_cutoff
    self.closed_rmsd_cutoff = closed_rmsd_cutoff
    self.vdw_factor = vdw_factor
    self.close = close
    self.verbose = verbose
    self.errstream = errstream
    self.meld = meld
    self.max_melding_rmsd = max_melding_rmsd
    self.min_database_loop_length = 1
    self.no_structure_output = nostruc
    self.mutate = mutate
    self.ccd_target_rmsd = ccd_target_rmsd
    self.ccd_iterations = ccd_iterations
    self.max_decoys = max_decoys
    self.first_decoys = first_decoys
    self.extension_size = extension_size
    self.extension_minimum = extension_minimum
    self.min_b_factor = 0.0

    self.calculate_contacts = calculate_contacts
    self.contact_distance = contact_distance
    self.contact_identity = contact_identity
    
    self.warnings = []
    
    if db:
      self.set_db(db)
    
    self.set_subst_tables(subst_table_path)
    
    self.results = []
    self.counter = itertools.count(1)
    self.idecoy = 0
  
  def set_options(self, opts):
    for o in opts:
      if o in self.__dict__ and not o.startswith('_'):
        self.__dict__[o] = opts[o]
  
  def set_silent(self):
    #self.verbose = False
    self.errstream = open(os.devnull, 'w')
  
  def note(self, message, force=False):
    if force or self.verbose:
      self.errstream.write(str(message)+"\n")
  
  def warn(self, exception, force=False):
    if isinstance(exception, str):
      exception = FreadError(exception)
    self.note("WARNING: "+str(exception), force)
    self.warnings.append(exception)
    return exception
    
  def die(self, exception, force=True):
    raise self.warn(exception, force)
  
  
  def set_subst_tables(self, subst_table_path):
    self.subst_tables = load_esst(subst_table_path)
  
  def set_db(self, db):
    self.min_database_loop_length = get_min_loop_length_in_database(db)
    self.db = db
  
  def set_structure(self, structure):
    self.clear_results()
    p = Pdb(structure)
    self.ligand_atoms = p.ligands
    self.residues = ResidueList(p)
    add_oxygens(self.residues)
  
  
  def get_loop_index(self, residue_number, inscode="", chain=None):
    if inscode:
      inscode = inscode.strip()
    if chain:
      chain = chain.strip()
    
    if isinstance(residue_number, str) and residue_number[-1].isalpha():
        ires = int(residue_number[:-1])
        inscode = residue_number[-1]
    else:
        ires = int(residue_number)
    
    start_of_loop=None
    for i,r in enumerate(self.residues):
      if r.ires == ires and r.inscode == inscode and (not chain or chain == r.chain):
          start_of_loop = i+1
          break
    
    return start_of_loop
  
  
  def get_structure_sequence(self, startindex=0, endindex=None):
    return self.residues[startindex:endindex].get_seq()
  
  
  def set_loop(self, start_of_loop, loop_sequence, resnum=False, chain=None, structure_known=None):
    self.clear_results()
    if chain is not None:
      chain = chain.strip()
    
    # We've been given a residue number instead of a start index
    if resnum:
      if isinstance(start_of_loop, str):
        if start_of_loop[-1].isalpha():
          ires = int(start_of_loop[:-1])
          inscode = start_of_loop[-1]
        else:
          ires = int(start_of_loop)
          inscode = ""
      else:
        ires = int(start_of_loop)
        inscode = ""
      
      found=False
      for i,r in enumerate(self.residues):
        if r.ires == ires and r.inscode == inscode and (not chain or chain == r.chain):
            start_of_loop = i+1
            found = True
            self.note("Located residue left of loop to be modelled: %d%s. Loop starting index: %d. Sequence surrounding start of loop: %s|%s" % (r.ires, r.inscode, start_of_loop, self.residues[max(0,start_of_loop-8):start_of_loop].get_seq(), self.residues[start_of_loop:start_of_loop+8].get_seq()))
            break
      
      if not found:
        msg = "Residue before start of loop (residue number %d%s" % (ires, inscode)
        if chain:
          msg += ", chain '%s'" % chain
        msg += ") not found in query structure"
        self.die(BadInputError(msg))
    
    if start_of_loop < ANCHOR_LENGTH or start_of_loop > len(self.residues) - ANCHOR_LENGTH:
      self.die(BadInputError("Cannot model loop closer than %d residues to the terminus (need a complete anchor on both sides)." % (ANCHOR_LENGTH)))
    
    # Assuming no coordinates are present in file, start and end of the loop are the same
    end_of_loop = start_of_loop
    
    if structure_known is not None:
      # We were told if loop coordinates are present in input file.
      # Skip auto-detection.
      #
      self.loop_structure_is_known = structure_known
      if structure_known:
        self.note("User tells me loop structure is present.")
      else:
        self.note("User tells me loop structure is not present.")
    # Auto-detect if loop coordinates are present in input file.
    else:
      # Are the coordinates are in the file?
      self.loop_structure_is_known = is_framework_consecutive(self.residues[start_of_loop-1],
                                                    self.residues[start_of_loop])
      self.note("Is loop structure present: "+str(self.loop_structure_is_known))
    
    # If we have the native loop structure, adjust end_of_loop.
    # Also do some additional sanity checks
    #
    if self.loop_structure_is_known:
      end_of_loop += len(loop_sequence)
      strucseq = ""
      for i,r in enumerate(self.residues[start_of_loop:end_of_loop]):
          strucseq += r.get_seq()
      if (not self.mutate) and (loop_sequence != strucseq):
        self.die(BadInputError("Residues differ between sequence and structure input: %s, %s" % (loop_sequence, strucseq)))
      self.note("Loop sequence of given ATOM co-ordinates: %s" % (strucseq))
      del strucseq
    
    if end_of_loop != start_of_loop:
      if len(loop_sequence) != end_of_loop - start_of_loop:
        self.die(BadInputError("Loop co-ordinates present in input, but number of residues (%d) does not match length of input loop sequence (%d)." % (end_of_loop - start_of_loop, len(loop_sequence))))
      if end_of_loop > len(self.residues) - ANCHOR_LENGTH:
        self.die(BadInputError("Cannot model loop closer than %d residues to the terminus (need a complete anchor on both sides)." % (ANCHOR_LENGTH)))
    
    
    # Set ourselves up for loop extension
    #
    self.seq = self.loop_sequence = loop_sequence
    self.s = self.start_of_loop = start_of_loop
    self.e = self.end_of_loop = end_of_loop
    
    self.verify_anchors()
    
    while len(self.seq) < self.min_database_loop_length:
      try:
        self.verify_stretch()
      except LoopStretchError:
        pass
      self.extend_loop()
    
    while True:
      try:
        self.verify_stretch()
        break
      except LoopStretchError as e:
        try:
          self.extend_loop()
        except BadInputError:
          self.die(LoopStretchError(str(e)+" Cannot extend gap further."))
    

    
  
  def extend_loop(self):
      """Extend loop by one residue on each side. Raises NonConsecutiveAnchorError() or UnsupportedLoopLengthError() upon failure.
      """
      start_of_loop = self.start_of_loop
      end_of_loop = self.end_of_loop
      s = self.s - 1
      e = self.e + 1
      
      length = (start_of_loop - s) + (e - end_of_loop) + len(self.loop_sequence)

      if s<ANCHOR_LENGTH or not is_framework_consecutive(self.residues[s-ANCHOR_LENGTH], self.residues[s-ANCHOR_LENGTH+1]) or e>len(self.residues)-ANCHOR_LENGTH or not is_framework_consecutive(self.residues[e+ANCHOR_LENGTH-2], self.residues[e+ANCHOR_LENGTH-1]):
        self.die(NonConsecutiveAnchorError("Cannot extend loop to length %d, due to gaps in the query structure or proximity to the termini"%(length)))
      
      
      if not is_loop_length_in_database(self.db, length):
        self.die(UnsupportedLoopLengthError("Cannot extend loop to length %d, due to database limitations"%(length)))
      
      seq = ""
      for i in xrange(s, start_of_loop):
        seq += self.residues[i].get_seq()
      seq += self.loop_sequence
      for i in xrange(end_of_loop, e):
        seq += self.residues[i].get_seq()
      
      assert len(seq) == length, str([s, start_of_loop, end_of_loop, e, length, len(seq), seq, dbdir, pdb_file, start_of_loop, loop_sequence])
      
      self.s = s
      self.e = e
      self.seq = seq
      
      self.note("Extending loop to length %d" % (length))
  
  
  def verify_anchors(self):
    # Ensure anchors are consecutive stretches of amino acids
    for x in xrange(self.s-ANCHOR_LENGTH+1, self.s):
      if not is_framework_consecutive(self.residues[x-1], self.residues[x]):
        self.die(NonConsecutiveAnchorError("Anchor residues not consecutive in framework structure: residue index %s, %s"%(self.s, self.residues.code)))
    for x in xrange(self.e, self.e+ANCHOR_LENGTH-1):
      if not is_framework_consecutive(self.residues[x], self.residues[x+1]):
        self.die(NonConsecutiveAnchorError("Anchor residues not consecutive in framework structure: residue index %s, %s"%(self.s, self.residues.code)))
  
  def verify_stretch(self):
    # Ensure anchors are not too far apart for loop to stretch gap
    anchor_distance = dist(self.residues[self.s-1].C, self.residues[self.e].N)
    if anchor_distance > get_max_loop_length(len(self.seq)) * 1.05:
      self.die(LoopStretchError("Loop span (%.2f Angstrom) too large to be closed using %d residues. Trying to extend gap." % (anchor_distance, len(self.seq))))
  
  def verify_length(self):
    # Ensure loop is not too long for the database    
    if not is_loop_length_in_database(self.db, len(self.seq)):
      self.die(UnsupportedLoopLengthError("Cannot model loop of length %d, due to database limitations"%(len(self.seq))))
  
  def verify(self):
    self.verify_anchors()
    self.verify_stretch()
    self.verify_length()
  
  
  def model(self, top=None, stop_after=None, f_rank_decoy=None, f_stop_search=None, f_filter=None):
    """Model loop using the FREAD algorithm. This is the method handling most of the work. Raises UnsupportedLoopLengthError if loop length is not supported.
    """
    if top is None:
      top = self.max_decoys
    
    if stop_after is None:
      stop_after = self.first_decoys
    
    if top <= 0:
      top = sys.maxint
    
    if stop_after > 0 and self.idecoy >= stop_after:
      return self.results
    
    if not f_rank_decoy:
      f_rank_decoy = FREAD_RANKING
    if not f_stop_search:
      f_stop_search = lambda x: False
    if not f_filter:
      f_filter = lambda x: True
    
    while len(self.seq) < self.min_database_loop_length:
      self.extend_loop()
    
    self.verify()
    
    meld_anchors = self.meld  # and (self.open_rmsd_cutoff <= self.max_melding_rmsd)
    close_loop = self.close
    
    verbose = self.verbose
    residues = self.residues
    start_of_loop = self.s
    end_of_loop = self.e
    loop_sequence = self.seq
    
    start = start_of_loop - ANCHOR_LENGTH
    end = end_of_loop + ANCHOR_LENGTH
    loop_length = len(loop_sequence)
    total_length = loop_length + 2*ANCHOR_LENGTH
    
    self.note("Loop region: N(%4d, %4d)  C(%4d,%4d)  Loop(%3d,'%s')" % (start, start_of_loop, end_of_loop, end, loop_length, loop_sequence))
    
    
    ############################################################################
    
                # Get query anchors and prepare for clash checking #
    
    print "Getting anchor coordinates."
    # Get anchor coordinates
    #
    anchor_N = residues[start:start_of_loop]
    anchor_C = residues[end_of_loop:end]
    
    
    # Describe anchors in query structure
    #
    anchor_description, query_transform = describe_anchors(anchor_N, anchor_C, loop_length)
    
    print "Building a geometric hash"
    # Build a GeometricHash of the query structure (without the loop region), for
    # clash checking
    #
    coords = []
    gh_atoms = []
    for r in residues[:start]+residues[end:]:
      for a in r:
        #if a.atom in ("N", "CA", "C", "O", "CB"):
        coords.append(a.xyz)
        gh_atoms.append(a)
    gh = GeometricHash(numpy.array(coords))
    del coords
    
    print "Calculating inter-residue contacts"
    # Inter-residue contacts
    if self.calculate_contacts:
      p = Protein((residues[:start]+residues[end:]).split_chains())
      p.ligands = Protein(ResidueList(self.ligand_atoms).split_chains())
      contact_gh, contact_gh_atoms = make_contact_gh(p)
    
    
    ############################################################################
    
                                # Search the database #
    
    results = [] # Heap Queue of the top-ranking decoys
    for decoy in self.results:
      
      _heapadd(results, (f_rank_decoy(decoy), decoy.idecoy, decoy), top)
    
    
    for decoy in iterate_database(self.db, loop_length, self.subst_tables, anchor_description, loop_sequence, self.open_rmsd_cutoff, self.score_cutoff):
        # Structure-based steps: loop closure, clash check, more filtering #
        print decoy.score,decoy.struc,decoy.seq
	
        if len(results) >= top:
          if decoy.internal_rmsd > results[0][-1].anchor_rmsd_open:
            continue
        
        # Retrieve loop structure from database
        decoy_residues = get_loop_structure(self.db, decoy.struc, decoy.start, total_length)
        
        assert len(decoy_residues) == total_length
        
        
        # Superimpose anchors and check anchor RMSD before starting
        anchor_rmsd_open = ccd.superimpose(decoy_residues[:ANCHOR_LENGTH]+decoy_residues[-ANCHOR_LENGTH:], anchor_N+anchor_C, decoy_residues)
        
        if anchor_rmsd_open > self.open_rmsd_cutoff:
          self.note("%s_%d_%d : Anchor RMSD too large: %.3f"%(decoy.struc, decoy.start, loop_length, anchor_rmsd_open))
          continue
        
        if len(results) >= top:
          if anchor_rmsd_open > results[0][-1].anchor_rmsd_open:
            continue
        
        
        # Save start residue number of loop, in database structure
        decoy.startres = decoy_residues[ANCHOR_LENGTH].ires
        decoy.startinscode = decoy_residues[ANCHOR_LENGTH].inscode
        
        # Relabel residues and discard non-matching atoms
        relabel_loop(decoy_residues[ANCHOR_LENGTH:-ANCHOR_LENGTH], loop_sequence, prevatom=anchor_N[-1].CA, nextatom=anchor_C[0].CA)
        
        
        if self.loop_structure_is_known:
          loop_rmsd_open = calculate_rmsd(decoy_residues[ANCHOR_LENGTH:-ANCHOR_LENGTH], residues[start_of_loop:end_of_loop])
        
        
        
        # Are we allowed to meld this decoy?
        #
        meld_this_decoy = meld_anchors and (anchor_rmsd_open <= self.max_melding_rmsd) #and not self.no_structure_output
        close_this_decoy = close_loop and anchor_rmsd_open > self.ccd_target_rmsd and (not meld_anchors or anchor_rmsd_open > self.max_melding_rmsd)
        
        #if not self.no_structure_output or close_this_decoy:
        decoy.nanchor = anchor_N.deep_copy()
        decoy.canchor = anchor_C.deep_copy()
        
        if meld_this_decoy:
          meld(decoy.nanchor, decoy_residues[:ANCHOR_LENGTH])
          meld(decoy.canchor, decoy_residues[-ANCHOR_LENGTH:], invertWeights=True)
        
        if not close_this_decoy:
          anchor_rmsd_closed = anchor_rmsd_open
          iterations = 0
        else:
          # Loop closure
          anchor_rmsd_closed, iterations = ccd.close_loop(decoy.nanchor, decoy_residues, decoy.canchor, target_rmsd=self.ccd_target_rmsd, iterations=self.ccd_iterations)
          if anchor_rmsd_closed > self.closed_rmsd_cutoff:
            self.note("Failed to close loop: %s_%d_%d"%(decoy.struc, decoy.start, loop_length))
            continue
        
        # Cut off the decoy loop's anchors
        decoy_residues = decoy_residues[ANCHOR_LENGTH:-ANCHOR_LENGTH]
        
        
        #if not self.no_structure_output:
        # Restore main chain oxygens, which got lost during the melding/closing procedure
        decoy_residues[-1].O = None # This oxygen is wrong anyway, so delete it
        add_oxygens(decoy.nanchor+decoy_residues[:1], force=True)
        add_oxygens(decoy_residues[-1:]+decoy.canchor+residues[end:end+1], force=True)
        
        
        # Clash check
        is_clashing = is_clash(gh, gh_atoms, decoy_residues, self.vdw_factor)
        if is_clashing:
          self.note("Clash detected in decoy: %s_%d_%d"%(decoy.struc, decoy.start, loop_length))
          continue
        
        
        if self.loop_structure_is_known:
          decoy.loop_rmsd_open = loop_rmsd_open
          if not iterations:
            decoy.loop_rmsd_closed = loop_rmsd_open
          else:
            decoy.loop_rmsd_closed = calculate_rmsd(decoy_residues, residues[start_of_loop:end_of_loop])
        
        decoy.length = loop_length
        decoy.anchor_rmsd_open = anchor_rmsd_open
        decoy.iterations = iterations
       
        #if not self.no_structure_output:
        decoy.loop = decoy_residues
        
        
        if self.calculate_contacts:
          contacts = ""
          for i, r in enumerate(decoy_residues):
            contact_atoms = find_contacts(contact_gh, r, maxdist=self.contact_distance)
            contact_atoms.extend(find_contacts_simple(decoy.nanchor+decoy_residues+decoy.canchor, i, maxdist=self.contact_distance))
            contacts += get_contact_class(r.chain, contact_atoms)
          
          decoy.native_contacts = get_native_contacts(self.db, decoy.struc, decoy.start, total_length)
          
          id = 0
          for x,y in zip(decoy.native_contacts, contacts):
            if x == y:
              id += 1
          decoy.contact_identity = float(id) / len(contacts)
          decoy.contacts = contacts
          
          if (self.contact_identity > 0) and (decoy.contact_identity < self.contact_identity):
            continue
        
        
        # Get per-residue substitution scores
        #
        tables = self.subst_tables.tables
        ascii2index = self.subst_tables.ascii2index
        seqmap = tuple([ascii2index[ord(s)] for s in loop_sequence])
        dihed = decoy.dihedrals
        score_list = []
        for i,x in enumerate(decoy.seq):
          score_list.append(tables[int(dihed[i])][seqmap[i]][ascii2index[ord(x)]])
        assert sum(score_list) == decoy.score
        decoy.residue_scores = score_list
        
        # Save per-residue scores in the occupancy column and
        # the total score in the B factor column
        #
        #
        for i,r in enumerate(decoy.loop):
          for a in r:
            a.b = self.min_b_factor + max(0, 30 - decoy.score)
            # Average score over a 3 residue window...
            # because that's how MEDELLER does it
            sc = score_list[max(0,i-1):i+2]
            a.occup = float(sum(sc))/len(sc)
        
        
        
        
        if not f_filter(decoy):
          self.note("User-defined filter is excluding decoy: %s_%d_%d"%(decoy.struc, decoy.start, loop_length))
          continue
        
        self.idecoy = decoy.idecoy = next(self.counter)
        # Save result
        _heapadd(results, (f_rank_decoy(decoy), decoy.idecoy, decoy), top)
        
        if stop_after > 0 and self.idecoy >= stop_after:
          break
        
        if f_stop_search(decoy):
          break
    
    self.note("%d decoys found"%(len(results)))
    
    self.results = []
    while results:
      r = heappop(results)
      self.results.append(r[-1])
    self.results.reverse()
    return self.results
  
  
  def clear_results(self):
    self.results = []
    self.warnings = []
    self.counter = itertools.count(1)
    self.idecoy = 0
    
  
  def write_summary(self, outstream=sys.stdout, write_decoy_sequence=False):
    """Write summary information to specified output stream.
    """
    for decoy in self.results:
      outstream.write("%s_%d%s_%d\t%d\t%.3f\t%.3f\t%d\t%s" % (decoy.struc, decoy.startres, decoy.startinscode, decoy.length, decoy.score, decoy.internal_rmsd, decoy.anchor_rmsd_open, decoy.iterations,decoy.seq))
      if write_decoy_sequence:
        outstream.write("\t%s" % (decoy.seq))
      if self.calculate_contacts:
        outstream.write("\t%s" % (decoy.contacts))
        outstream.write("\t%s" % (decoy.native_contacts))
        outstream.write("\t%.3f" % (decoy.contact_identity))
      if self.loop_structure_is_known:
        outstream.write("\t%.3f\t%.3f" % (decoy.loop_rmsd_open, decoy.loop_rmsd_closed))
      outstream.write("\n")
  
  
  def write_decoy_structures(self, out_dir, top=0, suffix=".loop.atm", idfilter=[]):
    """Write decoy structure files (PDB format) to specified directory.
    """
    if self.no_structure_output:
      self.die("Cannot write decoy structures as structure output is disabled")

    if out_dir:
      if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    
    for i, decoy in enumerate(self.results):
      if top > 0 and i >= top:
        break
      
      decoyname = "%s_%d%s_%d"%(decoy.struc, decoy.startres, decoy.startinscode, decoy.length)
      
      if idfilter and (decoyname not in idfilter):
        continue
      
      decoyfile = os.path.join(out_dir, decoyname+suffix)
      f = open(decoyfile, "w")
      try:
        for r in decoy.nanchor:
          f.write(str(r))
        for r in decoy.loop:
          f.write(str(r))
        for r in decoy.canchor:
          f.write(str(r))
      finally:
        f.close()
  
  
  def assemble_model(self, decoy):
    # Start and end of loop region in model
    start = self.s - ANCHOR_LENGTH
    end = self.e + ANCHOR_LENGTH
    
    if isinstance(decoy, int):
      decoy = self.results[decoy]
      if self.no_structure_output:
        self.die("Cannot write decoy structures as structure output is disabled")
    
    model = ResidueList([])
    for i in xrange(start):
      model.append(self.residues[i].copy())
    model.extend(decoy.nanchor)
    model.extend(decoy.loop)
    model.extend(decoy.canchor)
    for i in xrange(end, len(self.residues)):
      model.append(self.residues[i].copy())
    
    # Replace main chain oxygen before the loop, to ensure the peptide bond there is planar
    model[start-1].O = None
    add_oxygens(model, start=start-1, end=start, force=True)
    
    return model
  
  
  def write_model_structures(self, out_dir, top=0, suffix=".model.atm", idfilter=[]):
    """Write model structure files (PDB format), including decoy residues, to specified directory.
    """
    
    if self.no_structure_output:
      self.die("Cannot write model structures as structure output is disabled")
    
    if out_dir:
      if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    else:
      out_dir = ""
    
    for i, decoy in enumerate(self.results):
      if top > 0 and i >= top:
        break
      
      decoyname = "%s_%d%s_%d"%(decoy.struc, decoy.startres, decoy.startinscode, decoy.length)
      
      if idfilter and (decoyname not in idfilter):
        continue
      
      decoyfile = os.path.join(out_dir, decoyname+suffix)
      
      model = self.assemble_model(decoy)
      
      f = open(decoyfile, "w")
      try:
        f.write(str(model))
      finally:
        f.close()


  def open_errstream(self, filename):
      if filename == "-":
        self.errstream = sys.stderr
      else:
        self.errstream = open(filename, "w")


  def close_errstream(self):
      if self.errstream != sys.stderr:
        self.errstream.close()
        self.set_silent()
  
  
  def automodel_loop(self, start_of_loop, loop_sequence,
                      loopsuffix = "",
                      dbdir = '.',
                      strucdir = 'decoys',
                      summary_file_name = 'summary',
                      write_decoys = True,
                      write_models = True,
                      write_decoy_sequence = False,
                      resnum = False,
                      chain = None,
                      structure_known = None,
                      f_rank_decoy = None,
                      f_stop_search = None,
                      f_filter = None,
                      idfilter = [],
                      **kwargs):
      """Search the given database and return a list of decoys for the given loop.
      """
      
      # Set options defined in the Fread constructor
      self.set_options(kwargs)
      
      if self.extension_size < 0:
        max_extension = sys.maxint
      else:
        max_extension = self.extension_size
      
      if strucdir is not None:
        strucdir = strucdir.strip()
      if (not strucdir) or (not write_decoys and not write_models):
          strucdir = None
          self.no_structure_output = True
      
      for db_group in dbdir.split(":"): # try these databases until we find something
        databases = db_group.split("|") # combine results from all these databases
        self.set_db(databases[0])
        
        try:
	 
          self.set_loop(start_of_loop, loop_sequence.upper(), resnum=resnum, chain=chain, structure_known=structure_known)
        except FreadError:
          break
       
        # Extend loop if extension_minimum is set
        #
        for i in xrange(self.extension_minimum):
          try:
            self.extend_loop()
          except FreadError:
            break
        
        # Try to model, and extend loop if nothing was found
        anyresults = False
        for db in databases:
	  
          self.set_db(db)
          try:
	    print "First attempt"
            anyresults |= bool(self.model(f_rank_decoy=f_rank_decoy, f_stop_search=f_stop_search, f_filter=f_filter))
          except FreadError:
            break
        
        for i in xrange(max_extension):
          if anyresults:
            break
          try:
            self.extend_loop()
          except FreadError:
            break
          for db in databases:
            self.set_db(db)
	    print "Second attempt"
            anyresults |= bool(self.model(f_rank_decoy=f_rank_decoy, f_stop_search=f_stop_search, f_filter=f_filter))
        
        if anyresults:
          break
      
      # Write results to STDOUT or a file, if one was specified
      if summary_file_name:
        if summary_file_name == "-":
          self.write_summary(sys.stdout, write_decoy_sequence=write_decoy_sequence)
        else:
          root, ext = os.path.splitext(summary_file_name)
          outstream = open(root + loopsuffix + ext, "w")
          self.write_summary(outstream, write_decoy_sequence=write_decoy_sequence)
          outstream.close()
      
      if not self.no_structure_output:
        if write_decoys:
          self.write_decoy_structures(strucdir + loopsuffix, suffix=".loop.pdb", idfilter = idfilter)
        if write_models:
          self.write_model_structures(strucdir + loopsuffix, suffix=".model.pdb", idfilter = idfilter)
      
      return self.results


def meld(fragScaffold, fragPrediction, invertWeights=False):
    """Averages the coordinates of the two Pdb arguments. Move the first object's co-ordinates onto the averaged position.
    
    By default, the first object is assumed to be part of the N-terminal fragment, the second is part of the C-terminal fragment. This can be reversed by setting invertWeights=True."""
    
    resS = fragScaffold
    resP = fragPrediction
    L = len(resS)
    
    assert len(resS) == len(resP)
    
    def averageCoord(P, S):
      # P = coordinate of prediction
      # S = coordinate of scaffold
      # D = distance (in residues) from (loop+anchor)-fragment end
      # L = anchor length
      return 1.0/(L+1) * (D*P + (L+1-D)*S)
    
    for i, (rS, rP) in enumerate(zip(resS, resP)):
      if invertWeights:
        D = len(resS)-i
      else:
        D = i+1
      
      newN = averageCoord(rP.N.xyz, rS.N.xyz)
      newCA = averageCoord(rP.CA.xyz, rS.CA.xyz)
      newC = averageCoord(rP.C.xyz, rS.C.xyz)
      
      T_from, T_to, rotmat = ccd.get_rotmat([rS.N.xyz, rS.CA.xyz, rS.C.xyz], [newN, newCA, newC])
      
      rS.O = None # Remove main chain oxygen - can regenerate this later
      
      for a in rS:
        a.xyz = numpy.dot(a.xyz - T_from, rotmat) + T_to
      

def FREAD_RANKING(decoy):
  """Returns a comparable object (a tuple of scores) used for ranking a given decoy.
  
  Bigger values will be ranked higher."""
  return (-decoy.anchor_rmsd_open, decoy.score, -decoy.internal_rmsd, -decoy.iterations)


def FREAD_RANKING_BY_ESSS(decoy):
  """Returns a comparable object (a tuple of scores) used for ranking a given decoy.
  
  Bigger values will be ranked higher."""
  return (decoy.score, -decoy.anchor_rmsd_open, -decoy.internal_rmsd, -decoy.iterations)
