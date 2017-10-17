#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# Author: Sebastian Kelm
# Created: 25/06/2009
#





import os
import re
import copy
from array import array

from prosci.common import *
from prosci.util.gaps import gappify, deGappify, length_ungapped, isGap, Gappifier


def _entryIsMaster(entry):
  t = entry.type
  return t == "structure" or t == "sequence"



class Ali(object):

    class InvalidEntryGroupError(IndexError):
      pass
      
    class InvalidEntryError(IndexError):
      pass
    
    class EntryGroup(object):
    
        class MasterEntry(object):
          "Descriptor to access an EntryGroup's master entry"
          def __get__(self, obj, cls=None):
              return obj.getMasterEntry()
          def __set__(self, obj, val): 
              raise AttributeError("Illegal operation: cannot set attribute")
          def __delete__(self, obj): 
              raise AttributeError("Illegal operation: cannot delete attribute")
        
        
        master = MasterEntry()

        class Code(object):
          "Descriptor to access an EntryGroup's code"
          def __get__(self, obj, cls=None):
              return obj[0].code
          def __set__(self, obj, val):
#               raise raise AttributeError("Illegal operation: cannot set attribute")
              for entry in obj:
                entry.code = val
          def __delete__(self, obj): 
              raise AttributeError("Illegal operation: cannot delete attribute")
        
        code = Code()
        
        
        def __init__(self, *args):
            self.entries = []
            if len(args)==1 and not isinstance(args[0], Ali.Entry):
              args = args[0]
            for e in args:
              self.add(e)
        
        def __iter__(self):
            return self.entries.__iter__()

        def __getitem__(self, i):
            if isinstance(i, int):
              return self.entries[i]
            
            if i.startswith("structure"):
              i = "structure"
            elif i.startswith("sequence"):
              i = "sequence"
            
            if self.entries[0].type == i:
              return self.entries[0]
            
            ix = binarySearch(self.entries, i, 1)
            if ix >= 0:
              return self.entries[ix]
            
            raise Ali.InvalidEntryError("Ali.EntryGroup.__getitem__() : Entry not found: %s" % (str(i)))
        
        def __delitem__(self, i):
            if isinstance(i, int) or isinstance(i, slice):
              del self.entries[i]
              return
            
            if i.startswith("structure"):
              i = "structure"
            elif i.startswith("sequence"):
              i = "sequence"
            
            if self.entries[0].type == i:
              del self.entries[0]
              return
            
            ix = binarySearch(self.entries, i, 1)
            if ix >= 0:
              del self.entries[ix]
              return
            
            raise Ali.InvalidEntryError("Ali.EntryGroup.__getitem__() : Entry not found: %s" % (str(i)))
        
        def __getslice__(self, start=None, end=None):
            return Ali.EntryGroup(self.entries[start:end])

        def __len__(self):
            return len(self.entries)
        
        def __str__(self):
            return join("", self.entries)
        
        def toFastaString(self):
            "String representation of the object, in FASTA format (master sequences only)."
            return self.master.toFastaString()
        
        def __repr__(self):
            return "Ali.EntryGroup(%s)" % (repr(self.entries))
       
        def __cmp__(self, other):
            return cmp(self.code, other.code)

        def add(self, entry, replace=False):
            assert isinstance(entry, Ali.Entry)

            if len(self) == 0:
              self.entries.append(entry)
              return

            assert entry.code == self.code
            if entry.isMaster():
              if self.hasMasterEntry():
                if replace:
                  self.entries[0] = entry
                else:
                  raise ValueError("Duplicate entry in Ali.EntryGroup : %s" % (entry.desc))
              else:
                self.entries.insert(0, entry)
            else:
              i = binarySearch(self.entries, entry, 1)
              if i >= 0:
                if replace:
                  self.entries[i] = entry
                else:
                  raise ValueError("Duplicate entry in Ali.EntryGroup : %s" % (entry.desc))
              else:
                self.entries.insert( -(1+i), entry )

        def getCode(self):
            if not len(self):
              return None
            return self.entries[0].code
        
        def hasMasterEntry(self):
            return len(self.entries) > 0 and _entryIsMaster(self.entries[0])
        
        def getMasterEntry(self):
            assert self.hasMasterEntry()
            return self.entries[0]
        
        def getMasterType(self):
            return self.master.type
        
        def make_gaps_consistent(self):
            mainseq = self.master.seq
            for i in xrange(1,len(self)):
              entry = self[i]
              seq = array('c', entry.seq)
              for i,c in enumerate(mainseq):
                if isGap(c):
                  seq[i] = '-'
                elif isGap(seq[i]):
                  seq[i] = '?'
              entry.seq = seq.tostring()
        
        def gappify(self, force=True, gappifier=None):
            """gappify() : Gappify (align) entries within the same EntryGroup"""
            
            if gappifier is not None:
              for entry in self:
                entry.seq = gappifier.gappify(entry.seq)
              return
            
            template = self.master.seq
            
            if force:
              diff = True
            else:
              diff = False
              for i in xrange(1,len(self)):
                lcur = len(self[i].seq)
                if lcur != len(template):
                  if length_ungapped(template) != length_ungapped(self[i].seq):
                    raise IllegalStateError("Ungapped lengths differ between master and slave entries in same EntryGroup!")
                  diff = True
                  break
            
            if diff:
              for i in xrange(1,len(self)):
                entry = self[i]
                newseq = gappify(template, deGappify(entry.seq))
                if not newseq:
                  raise ValueError("Gappification of entry group '%s' failed!" % (entry.code))
                entry.seq = newseq
        
        def align(self, other, degappify=False):
            """align(other, degappify=False) : Align this EntryGroup to another using a sequence alignment of the master entries.
            
            If degappify==True, will first remove all gaps from sequences before realigning.
            Calls gappify() on both EntryGroups after performing the alignment of the master sequences.
            """
            self.master.align(other.master, degappify)
            self.gappify()
            other.gappify()
        
        def delete_types(self, typenames):
          entries = []
          for e in self:
            if e.isMaster() or (e not in typenames):
              entries.append(e)
          self.entries = entries
        
        def keep_types(self, typenames):
          entries = []
          for e in self:
            if e.isMaster() or (e in typenames):
              entries.append(e)
          self.entries = entries
    
    
    class Entry(object):
    
        class Type(object):
          "Descriptor to access an Entry's type"
          def __get__(self, obj, cls=None):
              return obj.getType()
          def __set__(self, obj, val):
              raise RuntimeError("Illegal operation: cannot set an Entry's type")
          def __delete__(self, obj): 
              raise AttributeError("Illegal operation: cannot delete attribute")
        
        type = Type()
        
        
        def __init__(self, *args):
            if len(args) == 1:
              args = args[0]
            assert len(args) in (2,3)
            
            if len(args) == 2 or not args[2]:
              # FASTA format
              self.code = args[0]
              self.desc = "sequence"
              self.seq  = args[1]
            else:
              # PIR format ( .ali .tem )
              self.code, self.desc, self.seq = args
            
            if self.code and self.code[0] == '>':
              self.code = self.code[1:]
            if self.seq and self.seq[-1] == '*':
              self.seq = self.seq[:-1]
        
        # Make Entry object seem like a list
        #
        def __iter__(self):
            yield self.code
            yield self.desc
            yield self.seq
        
        def __getitem__(self, i):
            if   0==i:
                return self.code
            elif 1==i:
                return self.desc
            elif 2==i:
                return self.seq
            raise IndexError("list index out of range")
            
        def __len__(self):
            return 3
        
        def __str__(self):
            return ">%s\n%s\n%s*\n" % (self.code, self.desc, self.seq)
        
        def toFastaString(self):
            "String representation of the object, in FASTA format."
            return ">%s\n%s\n"%(self.code, self.seq)

        
        def __repr__(self):
            return "Ali.Entry(%s,%s,%s)" % (repr(self.code), repr(self.desc), repr(self.seq))
        
        def __cmp__(self, other):
            if isinstance(other, Ali.Entry):
              return cmp(self.desc, other.desc)
            else:
              return cmp(self.desc, str(other))

        def getCode(self):
            return self.code

        def getDesc(self):
            return self.desc
        
        def getStructureFilename(self, searchdirs=None):
            def resolve(fname):
              if searchdirs:
                for dir in searchdirs:
                  if os.path.exists("%s/%s"%(dir, fname)):
                    return os.path.abspath("%s/%s"%(dir, fname))
              if os.path.exists(fname):
                return os.path.abspath(fname)
              else:
                return None
            
            if not self.isMaster():
              raise ValueError("Entry is not a master entry (type 'structure' or 'sequence'). Actual type: %s"%(self.type))
            code = self.code
            if code.startswith("P1;"):
              code=code[3:]
            codelist = [code]
            fields = self.getDesc().split(":")
            if len(fields)>1:
              codelist.append(fields[1])
            
            for code in codelist:
              for ext in ["pdb", "atm", "ent"]:
                fname = resolve("%s.%s"%(code, ext))
                if fname:
                  return fname
              
            fname = resolve("pdb%s.ent"%(code))
            if fname:
              return fname
            
            raise NotFoundError("Could not find structure file for entry with code '%s' in directories %s. Example filename: %s"%(self.code, str(searchdirs), codelist[0]+".pdb"))

        def getType(self):
            if self.desc.startswith("structure"):
              return "structure"
            elif self.desc.startswith("sequence"):
              return "sequence"
            else:
              return self.desc
        
        def isMaster(self):
            t = self.type
            return t == "structure" or t == "sequence"
        
        def get_seq(self):
            return self.seq
        
        def align(self, other, degappify=False):
            """align(other, degappify=False) : Align this Entry to another using a sequence alignment.
            
            If degappify==True, will first remove all gaps from sequences before realigning.
            """
            from prosci.util.seq import align as seqalign
            assert self.code != other.code
            
            if degappify:
              self.seq = deGappify(self.seq)
              other.seq = deGappify(other.seq)
            
            aligned = Ali(seqalign(self.toFastaString()+other.toFastaString()), fasta_mode=True)
            assert len(aligned) == 2
            assert deGappify(aligned[self.code].master.seq) == deGappify(self.seq)
            assert deGappify(aligned[other.code].master.seq) == deGappify(other.seq)
            self.seq  = aligned[self.code].master.seq
            other.seq = aligned[other.code].master.seq



  # BEGINNING OF ALI CLASS FUNCTIONS

    def __init__(self, inputdata, parent=None, ignore_duplicates=False, fasta_mode=False):
        self.parent      = None # points to a parent Ali instance, if available
        self.entrygroups = []
        
        if parent is not None and not isinstance(parent, Ali):
            raise ArgumentError("Parent must be an instance of the Ali class!")
        
        if inputdata:
          if isinstance(inputdata, Ali) or \
            ((isinstance(inputdata, list) or isinstance(inputdata, tuple)) and isinstance(inputdata[0], Ali.EntryGroup)):
              self.parent = parent
              if None != inputdata:
                  #self.entrygroups = [Ali.EntryGroup([Ali.Entry(e) for e in eg]) for eg in inputdata]
                  self.add(inputdata)
          elif parent == None:
              self.parse_file(inputdata, ignore_duplicates, fasta_mode)
        
    
    # Make Ali object seem like a list
    #
    def __iter__(self):
        return self.entrygroups.__iter__()
    
    def __getitem__(self, i):
        if isinstance(i, int):
          return self.entrygroups[i]
        for e in self.entrygroups:
          if e.code == i:
            return e
        raise Ali.InvalidEntryGroupError("Ali.__getitem__() : Ali.EntryGroup not found: %s" % (str(i)))
        
    def __len__(self):
        return len(self.entrygroups)
    
    def __str__(self):
        return join("", self.entrygroups)
    
    def __repr__(self):
        return "Ali(%s, %s)" % (repr(self.entrygroups), repr(self.parent))
    
    def __cmp__(self, other):
        if type(other) == type(None):
          return 1;
        return cmp(self.entrygroups, other.entrygroups)
    
    
    def add(self, entrygroups, merge_duplicates=False, replace_duplicates=False):
        entrygroups = copy.deepcopy(entrygroups)
        if isinstance(entrygroups, Ali.EntryGroup):
          entrygroups = [entrygroups]
        
        for eg in entrygroups:
          assert isinstance(eg, Ali.EntryGroup)
          assert eg.code != None
          
          if not self.has_entry(eg.code):
            self.entrygroups.append(eg)
          else:
            if not merge_duplicates:
              raise ValueError("Duplicate Ali.EntryGroup in Ali: %s" % (eg.code))
            
            oldeg = self[eg.code]
            oldMaster = oldeg.master.seq
            newMaster = eg.master.seq
            #if deGappify(oldMaster) == deGappify(newMaster):
            #  for e in eg:
            #    if not e.isMaster():
            #      oldeg.add(e, replace=replace_duplicates)
            #  oldeg.gappify()
            #else:
            if True:
              if oldMaster != newMaster:
                # Sequences differ... need to align old and new masters
                merger = AlignmentMerger(oldMaster, newMaster)
                
                # Gappify all old EntryGroups
                for egx in self:
                  if merger.is_ungapped_length1(egx.master.seq):
                    for e in egx:
                      if merger.is_ungapped_length1(e.seq):
                        e.seq = merger.gappify1(e.seq)
                        
                # Gappify all new EntryGroups
                for egx in entrygroups:
                  if merger.is_ungapped_length2(egx.master.seq):
                    for e in egx:
                      if merger.is_ungapped_length2(e.seq):
                        e.seq = merger.gappify2(e.seq)
              
              # Add entries of the current new entrygroup to oldeg
              for e in eg:
                if e.isMaster():
                  if e.seq == oldeg.master.seq:
                    continue
                  n=1
                  skip=False
                  while True:
                    desc = "alternative sequence %d" % (n)
                    if desc in oldeg or desc in eg:
                      if desc in oldeg and oldeg[desc].seq == e.seq:
                        skip=True
                        break
                      if desc in eg and eg[desc].seq == e.seq:
                        skip=True
                        break
                      n += 1
                      continue
                    e = copy.deepcopy(e)
                    e.desc = desc
                    break
                  if skip:
                    continue
                oldeg.add(e)
              
    
    def has_entry(self, code):
      try:
        self[code]
        return True
      except Ali.InvalidEntryGroupError:
        return False
    
    
    def has_equivalent_entry(self, entry):
        assert isinstance(entry, Ali.Entry)
        for eg in self:
          if eg.code == entry.code:
            for e in eg:
              if e.type == entry.type or (e.isMaster() and entry.isMaster()):
                return True
            break
        return False
    
    
    def add_entry(self, entry, ignore_duplicates=False, replace_duplicates=False):
        """add_entry(entry, ignore_duplicates=False) : Returns True if the entry was added, False otherwise"""
        assert isinstance(entry, Ali.Entry)
        for eg in self:
          if eg.code == entry.code:
            try:
              eg.add(entry, replace=replace_duplicates)
            except ValueError:
              if ignore_duplicates:
                return False
              raise
            return True
        self.add(Ali.EntryGroup(entry))
        return True
    
    
    def parse_file(self, fileobject, ignore_duplicates=False, fasta_mode=False):
        if not fileobject:
            return
        
        doClose = False
        
        if isinstance(fileobject, str):
          fileobject = fileobject.split("\n")
          if len(fileobject) == 1:
            fileobject = file(fileobject[0])
            doClose = True
        
        title = ''
        desc  = ''
        seq   = ''
        
        isAli   = False
        isFasta = False
        
        try:
          
          for line in fileobject:
              line = line.strip()
              if not line or line.startswith('*'):
                  continue
              elif line.startswith('>'):
                  if len(line) == 1:
                    raise ParsingError("ERROR: Stray '>' in sequence file!")
                  if title:
                      if seq:
                        isAli = True
                        if isFasta:
                          raise ParsingError("Sequence file must not mix ALI and FASTA formats:\n"+str(fileobject))
                      else:
                        isFasta = True
                        if isAli:
                          raise ParsingError("Sequence file must not mix ALI and FASTA formats:\n"+str(fileobject))
                      self.add_entry(Ali.Entry(title, desc, seq), ignore_duplicates)
                  title = line
                  desc = ''
                  seq = ''
              elif not desc:
                  desc = line
              elif fasta_mode:
                if desc.startswith("sequence") or desc.startswith("structure"):
                  desc = ""
                desc += line
              else:
                  seq += line
          
          if title and desc:
              if seq:
                isAli = True
                if isFasta:
                  raise ParsingError("Sequence file must not mix ALI and FASTA formats:\n"+str(fileobject))
              else:
                isFasta = True
                if isAli:
                  raise ParsingError("Sequence file must not mix ALI and FASTA formats:\n"+str(fileobject))
              self.add_entry(Ali.Entry(title, desc, seq), ignore_duplicates)
        
        finally:
          if doClose:
            fileobject.close()
    
    
    def gappify(self, gappifier=None):
      """gappify([gappifier]) : Gappify (align) entries within the same EntryGroup"""
      for eg in self:
        eg.gappify(gappifier=gappifier)
    
    
    def keep_columns(self, selection):
        """keep_columns(selection) : remove all columns not in the selection.
        
        accepts a list of booleans or a list of integer indeces.
        Assumes that all EntryGroups are aligned to each other."""
        
        if not isinstance(selection[0], bool):
          numsel = selection
          selection = [False] * len(self[0][0].seq)
          for i in numsel:
            selection[i] = True
        
        for eg in self:
          for e in eg:
            assert len(selection) == len(e.seq)
            newseq = array('c')
            for i,c in enumerate(e.seq):
              if selection[i]:
                newseq.append(c)
            e.seq = newseq.tostring()
        
              
        
    
    def remove_gaps(self, template_seq=None):
        """remove_gaps([template_seq]) : Removes all alignment columns containing gaps.
        
        Only applicable if sequences are aligned (at least within each EntryGroup).
        If template_seq is given, all entries must be aligned to template_seq.
        """
        # Note that using this function will unalign entries with different codes!
        
        if template_seq:
          assert isinstance(template_seq[0], str)
        
        for eg in self:
            code = eg.code
            
            if template_seq:
                template = template_seq
            else:
                if not eg.hasMasterEntry():
                  continue
                template = eg[0].seq
            
            seqs = []
            for entry in eg:
                assert len(entry.seq) == len(template)
                seqs.append(array('c', entry.seq))
            
            iend = len(template)
            ingap  = False
            for i in xrange(len(template), 0, -1):
                c = template[i-1]
                #print i,c,ingap
                if ingap:
                    if not isGap(c):
                        ingap = False
                        for seq in seqs:
                            del seq[i:iend]
                else:
                    if isGap(c):
                        ingap = True
                        iend  = i
            if ingap:
                #print "removing final gap"
                for seq in seqs:
                    del seq[:iend]
            
            for i,entry in enumerate(eg):
                entry.seq = seqs[i].tostring()
    
    
    def rename_entries(self, fromname, toname):
        
        if fromname[0] == ">":
          fromname = fromname[1:]
        if toname[0] == ">":
          toname = toname[1:]
        
        renamed = 0
        
        eg = self[fromname]

        for entry in eg:
           entry.code = toname
        
        return len(eg)
    
    
    def get_entries_by_type(self, etype):
        output = []
        for eg in self:
          try:
            output.append(eg[etype])
          except Ali.InvalidEntryError:
            pass
        return output

    
    def find_entry(self, code, etype):
        if code[0] == '>':
          code = code[1:]
        return self[code][etype]
    

    def find_entry_struc(self, code):
        return self.find_entry(code, "structure")
    
    
    def write_to_file(self, outfile, format="%s\n%s\n%s*"):
        """
        Write the entire Ali object to a specified file.
        
        write_to_file(outfile, format="%s\n%s\n%s*")
        
        Outfile can be a filename or a file object.
        The specified format is applied to each entry.
        """
        
        f = None
        doclose = True
        try:
          outfile.mode
          f = outfile
          doclose = False
        except AttributeError:
          f = open(outfile, 'w')
          doclose = True

        
        for entry in self:
          f.write(format % tuple(entry))
          f.write("\n")
        
        if doclose:
          f.close()
    
    
    def toFastaString(self):
        "String representation of the Ali object, in FASTA format."
        a=[]
        for eg in self:
          a.append(eg.toFastaString())
        return "".join(a)


    def align(self, degappify=False):
        """align(degappify=False) : Align all EntryGroups in a multiple sequence alignment, using the master entries.
        
        If degappify==True, will first remove all gaps from sequences before realigning.
        Calls gappify() on all EntryGroups after performing the alignment.
        """
        from prosci.util.seq import align as seqalign
        
        if degappify:
          self.remove_gaps()
        
        aligned = Ali(seqalign(self.toFastaString()), fasta_mode=True)
        
        for eg_self, eg_aligned in zip(self, aligned):
          assert deGappify(eg_aligned.master.seq) == deGappify(eg_self.master.seq)
          eg_self.master.seq  = eg_aligned.master.seq
          eg_self.gappify()
    
    
    def unalign(self):
        """unalign() : Remove all gaps from all sequences.
        """
        for eg in self:
          e = eg.master
          e.seq = deGappify(e.seq)
          eg.gappify()


    def delete_types(self, typenames):
      for eg in self:
        eg.delete_types(typenames)
    
    def keep_types(self, typenames):
      for eg in self:
        eg.keep_types(typenames)




class AlignmentMerger(object):
  def __init__(self, seq1, seq2, gaplist='-'):
    self.gaplist = gaplist
    self.seqlen1 = len(seq1)
    self.seqlen2 = len(seq2)
    
    seq1_aligned = deGappify(seq1)
    seq2_aligned = deGappify(seq2)
    
    if seq1_aligned != seq2_aligned:
      alignment = Ali(">seq1\nsequence\n%s\n>seq2\nsequence\n%s\n"%(seq1, seq2))
      alignment.align() # This removes all previously existing gaps
      
      seq1_aligned = alignment["seq1"].master.seq
      seq2_aligned = alignment["seq2"].master.seq
      
      assert len(seq1_aligned) == len(seq2_aligned)
    
    
    # We now have 4 sequences:
    # seq1 with initial gaps
    # seq1 with new gaps (seq1_aligned)
    # seq2 with new gaps (seq2_aligned)
    # seq2 with initial gaps
    #
    # We want to map seq1_aligned back to seq1 but keep any gaps we've introduced.
    # We want to map seq2_aligned back to seq2 but keep any gaps we've introduced.
    #
    # We need to do all this in such a way that every change we make to one sequence, we
    # also make to the other, to keep them the same length.
    
    # Calculate the mapping of seq1_aligned to seq1
    self.gfr1_aligned2output = Gappifier(seq1_aligned, seq1, gaplist)
    
    # Calculate the mapping of seq2_aligned to seq2
    self.gfr2_aligned2output = Gappifier(seq2_aligned, seq2, gaplist)
    
    # Calculate merge instructions
    #
    mergerules1 = []
    mergerules2 = []
    n1=0
    n2=0
    while n1 < len(self.gfr1_aligned2output.insertions) or n2 < len(self.gfr2_aligned2output.insertions):
      while n1 < len(self.gfr1_aligned2output.insertions):
        a = self.gfr1_aligned2output.insertions[n1]
        if a >= 0:
          break
        mergerules1.append(-1)
        mergerules2.append(-1)
        n1 += 1
      
      while n2 < len(self.gfr2_aligned2output.insertions):
        b = self.gfr2_aligned2output.insertions[n2]
        if b >= 0:
          break
        mergerules1.append(-1)
        mergerules2.append(-1)
        n2 += 1
      
      if n1 >= len(self.gfr1_aligned2output.insertions) or n2 >= len(self.gfr2_aligned2output.insertions):
        assert n1 >= len(self.gfr1_aligned2output.insertions) and n2 >= len(self.gfr2_aligned2output.insertions)
        break
      
      mergerules1.append(a)
      mergerules2.append(b)
      n1 += 1
      n2 += 1
    
    
    self.gfr1_aligned2output.insertions = mergerules1
    self.gfr2_aligned2output.insertions = mergerules2
    
    seq1_final = self.gfr1_aligned2output.gappify(seq1_aligned)
    seq2_final = self.gfr2_aligned2output.gappify(seq2_aligned)
    
    """
    print seq1
    print seq2
    print seq1_aligned
    print seq2_aligned
    print seq1_final
    print seq2_final
    """
    
    self.gfr1_input2output = Gappifier(seq1, seq1_final)
    self.gfr2_input2output = Gappifier(seq2, seq2_final)
  
  
  def is_ungapped_length1(self, seq):
    return len(seq) == self.seqlen1

  def is_ungapped_length2(self, seq):
    return len(seq) == self.seqlen2
  
  def gappify1(self, seq, gap=None):
    "AlignmentMerger.gappify1(seq[, gap]) : gappify seq"
    return self.gfr1_input2output.gappify(seq, gap)
  
  def gappify2(self, seq, gap=None):
    "AlignmentMerger.gappify2(seq[, gap]) : gappify seq"
    return self.gfr2_input2output.gappify(seq, gap)
  
  def deGappify1(self, seq):
    "AlignmentMerger.deGappify1(seq) : degappify seq"
    return self.gfr1_input2output.deGappify(seq)
  
  def deGappify2(self, seq):
    "AlignmentMerger.deGappify2(seq) : degappify seq"
    return self.gfr2_input2output.deGappify(seq)

