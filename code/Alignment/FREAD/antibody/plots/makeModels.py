#!/usr/bin/env python
import pdb
import os
import numpy as np
from prosci.loops.multifread import MultiFread
from prosci.util.residue import ResidueList
from tools import getFitCommands, getRMSD

dbdir = '../../../db_sol5'
subst_tables = '../../esst.txt'

strucdir = 'fw'
seqdir = 'seq'

strucs = sorted(os.listdir(strucdir), key = lambda x: x[:-4])
seqs = sorted(os.listdir(seqdir), key = lambda x: x[:-4])

strucseqs = zip(strucs, seqs)

class DataAggregator(MultiFread):
    def __init__(self, struc, seq, ref):
        MultiFread.__init__(self, struc, seq)
        self.refStruc = ref 
    def get_decoy_combos(self):
        try:
            self.scored_decoy_combos
        except:
            self.score_decoy_combos()

        combos = list(self.scored_decoy_combos)
        self.top_decoy_combos = []
        while True:
            if len(self.top_decoy_combos) >= 500:
                break

            combo = combos.pop(combos.index(min(combos, 
                                                key = lambda x: x[0])))[1]
            clashing = data.clash_check(combo)

            if not clashing:
                self.top_decoy_combos.append(combo)

    def compile_data(self, code):
        try:
            self.top_decoy_combos
        except:
            self.get_decoy_combos()

        self.data = []

        for i, combo in enumerate(self.top_decoy_combos):
            model = self.build_model(combo)
            outfile = os.path.join('models', code, 'model%i.pdb'%i)
            f = open(outfile, 'w')
            f.write(str(model))
            f.close()

            self.data.append(Datum(self.refStruc, combo, outfile, self.loops))



class Datum():
    def __init__(self, refStruc, decoys, model, loops):
        self.decoys = decoys
        self.model = model
        
        self.cdr_zones = [(i+1, i+1, j-i) for (i, j, k) in self.loops]

        self.cdr_rmsd = []
        for i in xrange(6):
            self.cdr_rmsd.append(getRMSD(refStruc, self.model, 
                                         zones=self.cdr_zones[i]))

        self.rmsd(getRMSD(refStruc, self.model))

if __name__ == '__main__':
    for s in strucseqs:

        code = s[0][:-4]

        whole = 'whole/%s.pdb'%code
        data = DataAggregator(os.path.join(strucdir, s[0]), 
                              os.path.join(seqdir, s[1]), whole)
        data.get_missing_loops(esst=subst_tables, 
                             db=dbdir,
                             outfile=os.path.join('models', 
                                                  '%s_model.pdb'%code),
                             top_decoys=5,
                             summary=os.path.join('summary', 
                                                  '%s_summary'%code),
                             nostrucdir=True)
                             
        
        data.compile_data(code)
        pdb.set_trace()
        print 'Done'
