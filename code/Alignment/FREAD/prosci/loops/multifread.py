#!/usr/bin/env python
import sys
import os
import difflib
import copy
from itertools import product, count
from multiprocessing import Pool, cpu_count
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

from . import fread
from .loopmodel import vdw_radii  # is_clash
from ..util.geohash import GeometricHash
from ..util.residue import ResidueList


class MultiFread(fread.Fread):
    """
    """
    def __init__(self, struc=None, seq=None, **kwargs):
        fread.Fread.__init__(self, **kwargs)
        
        self.decoys = None
        self.scored_decoy_combos = None
        self.combocounter = count(1)
        
        if seq:
            self.set_target_sequence(seq)

        if struc and seq:
            self.set_structure(struc)
            self.set_loops(self.sequence, self.residues)
    
    
    def clear_results_multi(self):
      self.clear_results()
      self.decoys = None
      self.scored_decoy_combos = None
      self.combocounter = count(1)
    
    
    def set_target_sequence(self, seq):
        try:
            f = open(seq)
            self.sequence = f.read().strip('\n')
            f.close()
        except:
            self.sequence = seq


    def set_loops(self, ref_seq, struc):
        d = difflib.Differ()
        diff = list(d.compare(ref_seq, struc.get_seq()))

        self.template_model = ResidueList([])
        self.loops = []
        loop_start = -1
        n = 0

        for i, c, in enumerate(diff): 
            if c[0] != '-':
                try:
                    self.template_model.append(self.residues[n])
                except IndexError:
                    self.die(fread.BadInputError("Sequence cannot be mapped to "
                                                +"structure. Sequence has more upper "
                                                +"case letters than residues present "
                                                +"in structure.\n"))
                if loop_start >= 0:
                    self.loops.append((loop_start, i, n))
                    loop_start = -1
                n += 1
            else:
                self.template_model.append(None)
                if loop_start < 0:
                    loop_start = i


    def build_model(self, decoys):
        """Return structure with decoys inserted into missing loop regions."""
        out_model = self.template_model.deep_copy()

        for (start, end, gappystart), decoy in zip(self.loops, decoys):
            if not decoy:
                continue
            
            # Insert co-ordinates of the top-ranking loop decoy into model
            offset = (decoy.length + start - end) / 2
            for j,i in enumerate(xrange(start-offset, 
                                        start-offset+decoy.length)):
                res = out_model[i] = decoy.loop[j]
                minb = (decoy.nanchor[-1].CA.b + decoy.canchor[0].CA.b) / 2.0
                for atom in res:
                    atom.b = max(minb, atom.b)
            
            # Replace anchor co-ordinates
            for j,i in enumerate(xrange(start-offset-len(decoy.nanchor), 
                                        start-offset)):
                out_model[i] = decoy.nanchor[j]
            
            for j,i in enumerate(xrange(start-offset+decoy.length, 
                              start-offset+decoy.length+len(decoy.canchor))):
                out_model[i] = decoy.canchor[j]

        return out_model


    def get_missing_loops(self, **kwargs):
        """Model each missing loop in the target structure using FREAD.

    Keywords arguments:
    db : path 
        Path to database. 
    esst : path
        Specify file containing environment-specific substitution tables 
        (defualt is ./esst.txt).
    strucdir : str
        Directory in which to write decoy structures, suffixed with XX, the
        number of the loop.
    summary : str
        Write search to specified file, suffixed with 'XX.table', where XX 
        is the number of the loop.
    nostrucdir : bool
        Disables structure output if True (other than final model).
    verbose : bool
        If True, print status messages to the file specified by messages.
    messages : str
        Redirect status messages to the specified file.
    max_decoys : int
        Keep only the specified number of top-ranked decoys for each loop.
        Default=100
    first_decoys : int
        Stop the search after finding the specified number of decoys for 
        each loop. Note that this is different from max_decoys.
    open_rmsd : float
        Sets the anchor RMSD cut-off.
    score : int
        Sets the substitution score cut-off. If the target sequence cannot 
        possibly attain the cut-off score, the cut-off is lowered to allow
        perfect matches to be found.
    nomeld : bool
        Do not meld loop anchors; this means all loops are closed instead, 
        unless noclose=True.
    noclose : bool
        Do not close loops. This means that if the open_rmsd is set too high, 
        decoys that do not fit the model's anchors will be returned.
    closed_rmsd : float
        Set the anchor RMSD cut-off after loop closure (defualt is 0.3). Note
        that the target anchor RMSD during loop closure is 0.15, so setting 
        this cut-off lower than 0.15 will result in almost zero coverage.
    max_melding_rmsd : float
        Set the RMSD cut-off that decides when to use melding and when to use
        loop closure.
    vdw_factor : float
        Set the factor for Van der Waal's radii. If a clash is detected, 
        the loop is discarded. A value </= disables clash checking.
    extension_max :
        If no decoys are found, extend the loop by one residue on each side 
        and repeat the search. This option specifies the maximum number of 
        extensions to try.
    extension_min :
        Extend loop INT times before the first modelling attempt. This 
        does not count towards the limit imposed by extension_max.
    debug : bool 
        Disable the multiprocessing model (for use with the debugger).

        """
        
        debug = kwargs.get('debug', False) # Disable multiprocessing for 
                                           # debugging.
        
        self.close = not kwargs.get('noclose', not self.close)
        self.meld = not kwargs.get('nomeld', not self.meld)
        self.verbose = kwargs.get('verbose', self.verbose)
        
        self.score_cutoff = kwargs.get('score', self.score_cutoff)
        self.open_rmsd_cutoff = kwargs.get('open_rmsd', self.open_rmsd_cutoff)
        self.closed_rmsd_cutoff = kwargs.get('closed_rmsd', self.closed_rmsd_cutoff)
        self.vdw_factor = kwargs.get('vdw_factor', self.vdw_factor)
        self.max_melding_rmsd = kwargs.get('max_melding_rmsd', self.max_melding_rmsd)
        
        self.no_structure_output = kwargs.get('nostrucdir', self.no_structure_output)
        
        self.min_b_factor = kwargs.get('min_b_factor', self.min_b_factor)
        
        try:
            self.subst_tables
        except:
            self.set_subst_tables(kwargs.get('esst', 'esst.txt'))
        
        max_decoys = kwargs.get('max_decoys', 100)
        first_decoys = kwargs.get('first_decoys', 1000)
        
        options = {
                  "dbdir" : kwargs.get('db', '.'),
                  "strucdir" : kwargs.get('strucdir', 'decoys'),
                  "summary_file_name" : kwargs.get('summary', 'summary.table'),
                  "max_decoys" : max_decoys,
                  "first_decoys" : first_decoys,
                  "extension_size" : kwargs.get('extension_max', 1000000),
                  "extension_minimum" : kwargs.get('extension_min', 0),
                  "f_rank_decoy" : kwargs.get('f_rank_decoy', None),
                  "f_filter" : kwargs.get('f_filter', None)
                  }
        
        args = [] # Pack args into a list of tuples for use with pool.map
        
        for i, (start, end, input_struc_start) in enumerate(self.loops):
            loop_sequence = self.sequence[start:end]
            loop_suffix = '%02d' %i
            args.append((copy.copy(self), (input_struc_start, loop_sequence, loop_suffix), options))
        
        if debug:
            fread_copies = map(_get_loop_decoys_debug, args)
        else:
            # Pass each loop to a worker process
            j = kwargs.get('pool_size', 1)
            if j <= 0:
                if cpu_count() > len(self.loops):
                    j = len(self.loops)
                else:
                    j = cpu_count()
            pool = Pool(j)
            fread_copies = pool.map_async(_get_loop_decoys, args).get()
            for frd in fread_copies:
              assert frd is not self
        
        self.decoys = []
        for frd in fread_copies:
          for warning in frd.warnings:
            self.warn(warning)
          self.decoys.append(frd.results)
    
    
    def model_structure(self, *args, **kwargs):
        """Build and output the top_models highest ranked models.
        """
        
        top_models = kwargs.get('top_models', 1)
        no_file_out = kwargs.get('nofileout', False)
        out_file_name = kwargs.get('outfile', 'multifread')

        vdw_factor = self.vdw_factor
        
        self.open_errstream(kwargs.get('messages', '-'))
        
        
        if not os.path.splitext(out_file_name)[1]:
            out_file_name += '.pdb'
        
        if not self.sequence:
            try:
                self.set_target_sequence(args[1])
            except:
                try:
                    self.set_target_sequence(kwargs.get('seq'))
                except:
                    self.die(fread.BadInputError('Error: No target sequence specified.'))
                    
        if not self.residues:
            try:
                self.set_structure(args[0])
            except:
                try:
                    self.set_structure(kwargs.get('struc'))
                except:
                    self.die(fread.BadInputError('Error: No target structure specified.'))
        
        
        if self.decoys is None:
            self.scored_decoy_combos = None
            self.get_missing_loops(**kwargs)
        
        if self.scored_decoy_combos is None:
            # Add the option of deleting a loop (and a placeholder for non-predicted loops)
            for d in self.decoys:
              if None not in d:
                d.append(None)
            
            #print [len(x) for x in self.decoys]
            
            
            def is_clash(gh, atoms, decoyatoms):
              if vdw_factor <= 0:
                return False
              
              for a in decoyatoms:
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
            
            
            pair_clashing = {}
            
            def check_clash(I, J, K, L):
              key = (I, J, K, L)
              if key in pair_clashing:
                return pair_clashing[key]
              
              d1 = self.decoys[I][J]
              d2 = self.decoys[K][L]
              
              if (d1 is None) or (d2 is None):
                pair_clashing[key] = False
                return False
              
              try:
                d1.ghatoms
              except:
                atoms  = []
                for r in d1.loop:
                  for a in (r.N, r.CA, r.C, r.O, r.CB):
                    if a is not None:
                      atoms.append(a)
                d1.ghatoms = atoms
              
              try:
                d2.ghatoms
              except:
                atoms  = []
                for r in d2.loop:
                  for a in (r.N, r.CA, r.C, r.O, r.CB):
                    if a is not None:
                      atoms.append(a)
                d2.ghatoms = atoms
              
              try:
                d1.gh
              except:
                try:
                  pair_clashing[key] = is_clash(d2.gh, d2.ghatoms, d1.ghatoms)
                  return pair_clashing[key]
                except:
                  d1.gh = GeometricHash(numpy.array([a.xyz for a in d1.ghatoms]))
              pair_clashing[key] = is_clash(d1.gh, d1.ghatoms, d2.ghatoms)
              return pair_clashing[key]
            
            
            def score_combo(combo):
              combo_rmsd = 0.0
              for i, j in enumerate(combo):
                decoy = self.decoys[i][j]
                if decoy is not None:
                  combo_rmsd += decoy.anchor_rmsd_open**2
                else:
                  combo_rmsd += self.open_rmsd_cutoff**2
              return combo_rmsd # / len(combo)
            
            
            combo_heap = []
            combo_history = set([])
            
            def clash_search(combo):
              #print
              #print combo
              combo = tuple(combo)
              # Don't try combinations more than once
              if combo in combo_history:
                #print combo, "is in the history", combo_history
                return
              combo_history.add(combo)
              
              score = score_combo(combo)
              
              # Stop search if we have enough decoys and only have worse results to go
              if len(combo_heap) >= top_models and score >= -combo_heap[0][0]:
                #print "Heap:", [x[0] for x in combo_heap]
                #print score, ">=", -combo_heap[0][0]
                return
              
              clashes = [False]*len(self.decoys)
              # Re-calculate clashes
              for i in xrange(len(combo)-1):
                for k in xrange(i+1, len(combo)):
                  if check_clash(i, combo[i], k, combo[k]):
                    clashes[i] = clashes[k] = True
              
              #print "clashes:", clashes
              
              # Accept non-clashing combinations
              if True not in clashes:
                _heapadd(combo_heap, (-score, next(self.combocounter), [self.decoys[i][j] for i,j in enumerate(combo)]), top_models)
                if len(combo_heap) < top_models:
                  # Modify combination and keep searching
                  for i in xrange(len(combo)):
                    cmb = list(combo)
                    cmb[i] += 1
                    if len(self.decoys[i]) > cmb[i]:
                      clash_search(cmb)
              else:
                # Modify clashing combinations
                for i in xrange(len(combo)):
                  if clashes[i]:
                    cmb = list(combo)
                    cmb[i] += 1
                    if len(self.decoys[i]) > cmb[i]:
                      clash_search(cmb)
            
            
            clash_search([0] * len(self.decoys))
            
            # Get a sorted list of (score, decoy_combination) tuples
            #
            self.scored_decoy_combos = []
            while combo_heap:
              c = heappop(combo_heap)
              self.scored_decoy_combos.append((-c[0], c[-1]))
            self.scored_decoy_combos.reverse()
        
        results = []
        for score, combo in self.scored_decoy_combos:
          model = self.build_model(combo)
          results.append(model)
        
        if not no_file_out:
          for i, model in enumerate(results):
            root, ext = os.path.splitext(out_file_name)
            f = open('%s%04d%s' %(root, i+1, ext), 'w')
            f.write(str(model))
            f.close()
        
        self.note('Finished')
        
        self.close_errstream()
        self.errstream = sys.stderr
        
        return results
    
    
#     def is_decoy_combo_clashing(self, decoys):
#         """Check whether each loop decoy clashes with other loops in the given set of loops"""
#         if self.vdw_factor <= 0:
#           return False
#         
#         def backbone(a):
#           return a.atom in ("N", "CA", "C", "O", "CB")
#         
#         for decoy in decoys:
#             if not decoy:
#                 continue
#             
#             # Skip all None decoys
#             # Make a list of all residues/atoms in all decoys apart from the current one
#             #
#             alldecoys = ResidueList([])
#             for d in decoys:
#               if d is None or d == decoy:
#                 continue
#               alldecoys.extend(d.loop)
#             
#             gh = GeometricHash(alldecoys.get_coords(atomfilter = backbone))
#             allatoms = alldecoys.to_pdb(atomfilter = backbone)
#             #gh = GeometricHash(alldecoys.get_coords())
#             #allatoms = alldecoys.to_pdb()
#             
#             if is_clash(gh, allatoms, decoy.loop, self.vdw_factor):
#               return True
#         
#         return False
# 
# 
# def score_decoy_combo(combo):
#         """Compute the combined anchor RMSD for the given combination of decoys.
#         
#         If r_i is the anchor RMSD for loop i then the combined anchor RMSD is
#         calculated by sqrt(sum(r_i**2)).
#         """
#         
#         if not combo:
#           return 0.0
#         
#         combo_rmsd = 0.0
#         for decoy in combo:
#           if decoy is not None:
#             combo_rmsd += decoy.anchor_rmsd_open**2
#         
#         return combo_rmsd / len(combo)


def _get_loop_decoys(args):
    """Wrapper for Fread.automodel_loop(), to be used with map().
    
    args[0] is a the MultiFread object, args[1] are the sequential
    arguments for Fread.automodel_loop(), args[2] are the keyword
    arguments for Fread.automodel_loop().
    
    Returns a MultiFread object that is a copy of the input object,
    which also contains the modelling results, warnings, etc.
    You can access the results and warnings like so:
    
    fread_objects = map(_get_loop_decoys, argslist)
    for obj in fread_objects:
      for decoy in obj.results:
        # do something
      for warning in obj.warnings:
        # do something
    
    """
    frd, seqargs, kwargs = args
    
    frd.set_silent()
    frd.automodel_loop(*seqargs, write_models=False, **kwargs)
    
    return frd

def _get_loop_decoys_debug(args):
    """Same as _get_loop_decoys, but does not silence Fread."""
    frd, seqargs, kwargs = args
    #frd = copy.deepcopy(frd)
    #frd.errstream =args[0].errstream
    #frd.set_silent()
    frd.automodel_loop(*seqargs, write_models=False, **kwargs)
    return frd


def _heapadd(heap, item, maxitems):
  if len(heap) < maxitems:
    heappush(heap, item)
    return True
  else:
    return heappushpop(heap, item) != item
