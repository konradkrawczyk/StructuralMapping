#!/usr/bin/python
# -*- coding: utf-8 -*-

from prosci.util.gaps import isGap

try:
  from prosci.util.clustalo import align as align_clustalo
  align = align_clustalo
except ImportError:
  pass

try:
  from prosci.util.muscle import align as align_muscle
  align = align_muscle
except ImportError:
  pass


def __stringRemoveIndeces(source, indeces):
    #indeces.sort()
    a = ""
    n=0
    for i in indeces:
      a += source[n:i]
      n = i+1
    a += source[n:len(source)]
    return a


def find_aligned_residues(seq1, seq2, equiv=None, eqchars=":."):
    assert len(seq1) == len(seq2)
    assert None == equiv or len(seq1) == len(equiv)
    
    numbering_list1 = []
    numbering_list2 = []
    
    n1=0
    n2=0
    if None == equiv:
      for i in xrange(len(seq1)):
          if not isGap(seq1[i]) and not isGap(seq2[i]):
              numbering_list1.append(n1)
              numbering_list2.append(n2)
              n1 += 1
              n2 += 1
          elif not isGap(seq1[i]):
              n1 += 1
          elif not isGap(seq2[i]):
              n2 += 1    
    else:
      for i, e in enumerate(equiv):
          if equiv[i] in eqchars:
              numbering_list1.append(n1)
              numbering_list2.append(n2)
          if not isGap(seq1[i]):
              n1 += 1
          if not isGap(seq2[i]):
              n2 += 1

    return numbering_list1, numbering_list2


def __get_pid_counts(seq1, seq2):
  assert len(seq1) == len(seq2), "\n>seq1\n%s\n>seq2\n%s\n"%(seq1, seq2)
  
  l1 = 0   # ungapped length of seq1
  l2 = 0   # ungapped length of seq2
  cov = 0  # number of aligned pairs (absolute coverage)
  idn = 0  # number of identical pairs (absolute identity)
  
  for i in xrange(len(seq1)):
    c1 = seq1[i]
    c2 = seq2[i]
    if (not isGap(c1)):
      l1 += 1
      if (not isGap(c2)):
        l2 += 1
        cov += 1
        if c1 == c2:
          idn += 1
    elif (not isGap(c2)):
      l2 += 1
  
  return (l1, l2, cov, idn)



def pid(seq1, seq2, mode="shorter"):
  l1, l2, cov, idn = __get_pid_counts(seq1, seq2)
  sh = min(l1, l2) # length of the shorter sequence
  ln = max(l1, l2) # length of the longer sequence
  
  # Note: len(seq1) == len(seq2) == alignment length
  # This is distinct from l1 and l2, which are the ungapped lengths of the two sequences.
  
  if mode == "shorter":
    return (float(idn) / sh, float(cov) / sh, sh, l1, l2)
  elif mode == "longer":
    return (float(idn) / ln, float(cov) / ln, ln, l1, l2)
  elif mode == "first":
    return (float(idn) / l1, float(cov) / l1, l1, l1, l2)
  elif mode == "second":
    return (float(idn) / l2, float(cov) / l2, l2, l1, l2)
  elif mode == "aligned":
    return (float(idn) / cov, 1.0, cov, l1, l2)
  elif mode == "all":
    return (float(idn) / len(seq1), float(cov) / len(seq1), len(seq1), l1, l2)
  raise ArgumentError("Illegal mode parameter: '%s'" %mode)


def consensus(seqs):
  cons=[]
  totals=[]
  length = len(seqs[0])
  for s in seqs:
    assert len(s) == length, "Input sequences must be aligned (of the same length)"
  for i in xrange(length):
    d = {}
    n=0
    for j in xrange(len(seqs)):
      c = seqs[j][i]
      if isGap(c):
        continue
      n+=1
      if c not in d:
        d[c]=1
      else:
        d[c]+=1
    residues=[]
    for res in sorted(d):
      residues.append((res, d[res]))
    cons.append(residues)
    totals.append(n)
  return cons, totals


def remove_gapped_columns(seqs, master=0):
  # Keep only those columns where master sequence has no gap
  #
  residuesToDelete = []
  n = 0
  for i,c in enumerate(seqs[master]):
    if isGap(c):
      residuesToDelete.append(i)

  for i in xrange(len(seqs)):
    assert len(seqs[i]) == len(seqs[-1]), "Sequences not aligned"
    seqs[i] = __stringRemoveIndeces(seqs[i], residuesToDelete)


