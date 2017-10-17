#
#  geohash.py
#  
#
#  Created by Sebastian Kelm on 18/09/2010.
#  Copyright (c) 2010 Sebastian Kelm. All rights reserved.
#

import sys
from bisect import bisect

import numpy as np


def sqrdist(a, b):
    return np.sum((a-b)**2)


class GeometricHash:
    def __init__(self, coords, cell_size=None):
        """GeometricHash(coords [, cell_size]) : Put the given coordinates on a latice, to allow fast nearest-neighbour searches.
        
        coords is a numpy.array of 3d coordinates
        cell_size determines the granularity of the lattice and hence the speed of
        the nearest-neighbour search. By default, this is set automatically, based
        on the distance between the first few successive atoms in the structure.
        """
        # Guess the optimal cell_size. This is fairly arbitrary, but it works well
        # enough.
        # This procedure allows good access speeds, independently of whether the
        # input structure is all-atom, backbone-only or CA-only.
        #
        if cell_size is None:
            cell_size=10.0
            for i in xrange(min(10, len(coords))):
                cell_size=min(cell_size, np.sqrt(sum((coords[i]-coords[i+1])**2)))
            cell_size = max(1.0, cell_size*1.27)

        self.cell_size = cell_size
        self.sqr_cell_size = cell_size**2
        self.coords = coords
        self.max = self.coords.max(0)
        self.min = self.coords.min(0)
        self.offset = np.floor(self.min) - self.cell_size/2
        
        self.origin = self.index((0.0, 0.0, 0.0))
        
        self.dimensions=self.index(self.max)+1
        
        # Create empty 3d array
        self.data = []
        for i in xrange(self.dimensions[0]):
            col = []
            self.data.append(col)
            for j in xrange(self.dimensions[1]):
                col.append([None]*self.dimensions[2])
        
        # Populate the array
        # (values are held in lists [the 4th dimension], so
        # multiple atoms at one location are allowed)
        #
        for i, c in enumerate(coords):
            self.add(self.index(c), i)
        
    def index(self, c):
        "GeometricHash.index(c) : return the index at which the given coordinates would be located in this GeometricHash"
        return np.array(np.floor((c-self.offset)/self.cell_size), dtype=int)
    
    def get(self, ix):
        "GeometricHash.get(index) : get value at the specified index"
        try:
            a = self.data[ix[0]][ix[1]][ix[2]]
        except IndexError:
            #print >>sys.stderr, str(ix), str(self.dimensions)
            #raise
            return None
        
        if ix[0]<0 or ix[1]<0 or ix[2]<0:
            return None
            #raise IndexError("Index out of bounds: index = %s, dimensions = %s"%(str(ix), str(self.dimensions)))
        
        return a
    
    def occupied(self, ix):
        "GeometricHash.occupied(index) : check if a given bin exists and is occupied by at least one atom"
        try:
          if self.data[ix[0]][ix[1]][ix[2]]:
            if ix[0]<0 or ix[1]<0 or ix[2]<0:
              return False
            return True
        except IndexError:
          pass
        return False
    
    def add(self, ix, value):
        L = self.data[ix[0]][ix[1]][ix[2]]
        if L is None:
            L = []
            self.data[ix[0]][ix[1]][ix[2]] = L
        L.append(value)
    
    def get_neighbours(self, c, num, max_dist=None, filter=None):
        """GeometricHash.get_neighbours(c, num, [max_dist], [filter]): find the num closest neighbours (in space) to c
        
           If c is an integer, it is treated as an Atom index in the coordinate list. In this case, c
           itself is filtered out of the neighbour list.
           If c is a coordinate (x,y,z), then it is assumed that c is not present in the original
           coordinate list and the neighbours are returned unfiltered.
           max_dist is the maximum distance in space any neighbour can be away from the query.
           filter is a function of the form:
              def filter(results, candidate)
              where results is the list of already found neighbours, candidate is the new candidate neighbour.
              The function should return True if candidate should be saved in the results, False otherwise.
              candidate and each entry in results is a tuple of the form (square_distance, atom_index).
        """
        if isinstance(c, int):
            query_index = c
            c = self.coords[c]
        else:
            query_index = None
        
        if max_dist is not None:
          max_dist = max_dist**2
        
        result=[]
        
        # Procedure to save and rank neighbours
        # Only ever keeps the closest "num" neighbours found so far
        #
        def save(idx):
            #print "saving", idx
            L = self.get(idx) # list of atoms in grid cell 'idx'
            if L is not None:
              for i in L:
                val = (sqrdist(self.coords[i], c), i)
                if not max_dist or val[0] < max_dist:
                  index = bisect(result, val)
                  if index < num:
                    if filter is None or filter(result, val):
                      result.insert(index, val)
                      if len(result) > num:
                          result.pop()
        
        # Save neighbours in the same lattice cell as the query
        # Filter out the query itself, if we know its index
        #
        ix = self.index(c)
        if query_index is None:
            save(ix)
        else:
            num+=1
            save(ix)
            num-=1
            for j, (i, sqd) in enumerate(result):
                if  i == query_index:
                    assert sqd == 0.0
                    result.pop(j)
                    break
            if len(result) > num:
                result.pop()
        
        
        # Walk through neighbouring cells in cubic lattice, to find spatial neighbours
        #
        d=0
        while len(result) < num or result[-1][0] > self.sqr_cell_size * d * d:
            d+=1
            
            if max_dist and self.sqr_cell_size * d * d > max_dist:
                break
            
            #print ix, d
            if (ix+d >= self.dimensions).all() and (ix-d < 0).all():
                break
            
            for x in xrange(max(0, ix[0]-d), min(self.dimensions[0], ix[0]+d+1)):
                
                # Two opposing sides of a cube
                for y in xrange(max(0, ix[1]-d), min(self.dimensions[1], ix[1]+d+1)):
                    if ix[2]-d > 0:
                        save((x, y, ix[2]-d))
                    if ix[2]+d < self.dimensions[2]:
                        save((x, y, ix[2]+d))
                
                # Two more sides
                # skipping the top & bottom rows on each side
                for z in xrange(max(0, ix[2]-d+1), min(self.dimensions[2], ix[2]+d)):
                    if ix[1]-d > 0:
                        save((x, ix[1]-d, z))
                    if ix[1]+d < self.dimensions[1]:
                        save((x, ix[1]+d, z))
                    
            for y in xrange(max(0, ix[1]-d+1), min(self.dimensions[1], ix[1]+d)):
                
                # The last two sides of the cube
                # skipping top & bottom rows, left- & right-most columns
                for z in xrange(max(0, ix[2]-d+1), min(self.dimensions[2], ix[2]+d)):
                    if ix[0]-d > 0:
                        save((ix[0]-d, y, z))
                    if ix[0]+d < self.dimensions[0]:
                        save((ix[0]+d, y, z))
        
        
        return [(x[1], np.sqrt(x[0])) for x in result]
