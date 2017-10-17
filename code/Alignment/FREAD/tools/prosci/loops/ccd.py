#!/usr/bin/python
# -*- coding: utf-8 -*-
#

import sys
import os
import subprocess

import numpy


class CCDException(RuntimeError):
  pass



def build_ccd_input(n, loop, c):
  txt = ""
  txt += "%.4f %.4f %.4f\n"%(n.N.x, n.N.y, n.N.z)
  txt += "%.4f %.4f %.4f\n"%(n.CA.x, n.CA.y, n.CA.z)
  txt += "%.4f %.4f %.4f\n"%(n.C.x, n.C.y, n.C.z)
  txt += "#\n"
  
  for a in loop:
    txt += "%.4f %.4f %.4f\n"%(a.N.x, a.N.y, a.N.z)
    txt += "%.4f %.4f %.4f\n"%(a.CA.x, a.CA.y, a.CA.z)
    txt += "%.4f %.4f %.4f\n"%(a.C.x, a.C.y, a.C.z)
  
  txt += "#\n"
  txt += "%.4f %.4f %.4f\n"%(c.N.x, c.N.y, c.N.z)
  txt += "%.4f %.4f %.4f\n"%(c.CA.x, c.CA.y, c.CA.z)
  txt += "%.4f %.4f %.4f\n"%(c.C.x, c.C.y, c.C.z)
  
  return txt


def parse_ccd_output(txt):
  coords = []
  lines = txt.splitlines()
  for line in lines:
    fields = line.split()
    coords.append(numpy.array([float(fields[0]), float(fields[1]), float(fields[2])]))
  return coords


def get_rotmat(cfrom, cto):
    cfrom = numpy.array(cfrom)
    cto = numpy.array(cto)
    
    assert len(cfrom) == len(cto)
    
    
    T_to = numpy.mean(cto, axis=0)
    T_from = numpy.mean(cfrom, axis=0)
    
    for i,c in enumerate(cto):
       cto[i] = c - T_to
    for i,c in enumerate(cfrom):
       cfrom[i] = c - T_from
    
    
    correlation_matrix = numpy.dot(cto.transpose(), cfrom)
    v, s, w = numpy.linalg.svd(correlation_matrix)
    
    is_reflection = (numpy.linalg.det(v) * numpy.linalg.det(w)) < 0.0
    if is_reflection:
        #s[-1] = - s[-1]
        v[:,-1] = -v[:,-1]
    
    rotation_matrix = numpy.dot(v,w).transpose()
    
    
    return T_from, T_to, rotation_matrix


def get_rotmat_and_rmsd(cfrom, cto):
    cfrom = numpy.array(cfrom)
    cto = numpy.array(cto)
    
    assert len(cfrom) == len(cto)
    
    
    T_to = numpy.mean(cto, axis=0)
    T_from = numpy.mean(cfrom, axis=0)
    
    for i,c in enumerate(cto):
       cto[i] = c - T_to
    for i,c in enumerate(cfrom):
       cfrom[i] = c - T_from
    
    
    correlation_matrix = numpy.dot(cto.transpose(), cfrom)
    v, s, w = numpy.linalg.svd(correlation_matrix)
    
    is_reflection = (numpy.linalg.det(v) * numpy.linalg.det(w)) < 0.0
    if is_reflection:
        s[-1] = - s[-1]
        v[:,-1] = -v[:,-1]
    
    rotation_matrix = numpy.dot(v,w).transpose()
    
    E0 = sum(sum(cfrom * cfrom)) + sum(sum(cto * cto))
    rmsd_sq = (E0 - 2.0*sum(s)) / float(len(cfrom))
    rmsd = numpy.sqrt(max(rmsd_sq, 0.0))
    
    return T_from, T_to, rotation_matrix, rmsd



def coords2residues(coords, loop):
    assert len(coords) == len(loop)*3
    
    for i in xrange(1, len(loop)-1):
      #   calculate translation/rotation for original main chain coords onto new coords
      #   rotate all old coords into new position
      #   handle O atoms separately, making them planar with CA, C and the next residue's N
      
      r = loop[i]
      j = i*3

      dO = numpy.linalg.norm(r.O.xyz - r.C.xyz)
      
      cfrom = [r.N.xyz.copy(), r.CA.xyz.copy(), r.C.xyz.copy()]
      cto = [coords[j], coords[j+1], coords[j+2]]
      T_from, T_to, rotmat = get_rotmat(cfrom, cto)
      
      for a in r:
        if a.atom != "O":
          a.xyz = numpy.dot(a.xyz - T_from, rotmat) + T_to
      
      vO_1 = r.C.xyz - r.CA.xyz
      vO_2 = r.C.xyz - numpy.array(coords[j+3]).squeeze()
      vO_3 = (vO_1 / numpy.linalg.norm(vO_1)) + (vO_2 / numpy.linalg.norm(vO_2))
      
      r.O.xyz = r.C.xyz + (vO_3 / numpy.linalg.norm(vO_3)) * dO


def superimpose(rfrom, rto, loop):
  cfrom = []
  cto   = []
  for r,q in zip(rfrom, rto):
    cfrom.append(r.N.xyz)
    cfrom.append(r.CA.xyz)
    cfrom.append(r.C.xyz)
    cto.append(q.N.xyz)
    cto.append(q.CA.xyz)
    cto.append(q.C.xyz)
    if r.O is not None and q.O is not None:
      cfrom.append(r.O.xyz)
      cto.append(q.O.xyz)
  
  T_from, T_to, rotmat, rmsd = get_rotmat_and_rmsd(cfrom, cto)
  
  for r in loop:
    for a in r:
      a.xyz = numpy.dot(a.xyz - T_from, rotmat) + T_to
  
  return rmsd


def close_loop(n, loop, c, target_rmsd=0.15, iterations=5000):
  "needs one residue per terminal, and loop must contain the corresponding residues as well"
  
  assert len(n) == len(c)
  anchor_length = len(n)
  
  #rmsd_before = min_rmsd((loop[0], loop[-1]), (n[-1], c[0]))
  #rmsd_before = superimpose(loop[:anchor_length]+loop[-anchor_length:], n+c, loop)
  
  loop = loop[anchor_length-1:-anchor_length+1]
  n = n[-1]
  c = c[0]
  
  input = build_ccd_input(n, loop, c)
  
  #p = subprocess.Popen("%s/cpp/ccd"%(os.path.dirname(sys.argv[0])), stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
  p = subprocess.Popen("%s/pyfread_cpp/ccd %f %d"%(os.path.dirname(sys.argv[0]), target_rmsd, iterations), stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
  
  output, messages = p.communicate(input)
  error_code = p.wait()
  
  messages = messages.strip()
  
  if error_code > 1:
    e = CCDException("Error code %d: %s"%(error_code, messages))
    e.returncode = error_code
    raise e
  
  try:
    iter, rmsd_after = messages.split()
  except ValueError:
    e = CCDException("Unexpected CCD output: %s"%(messages))
    e.returncode = error_code
    raise e
  
  iter = int(iter)
  rmsd_after = float(rmsd_after)

  coords = parse_ccd_output(output)
  coords2residues(coords, loop)
  
  return (rmsd_after, iter)
