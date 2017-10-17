#!/usr/bin/env python

import os
import numpy as np
import matplotlib.pyplot as plt

from subprocess import call

datafiles = [f for f in os.listdir('data') if f[-4:] == '.txt']

for file in datafiles:
    code = file[:4]
    data = np.loadtxt(os.path.join('data', file))
    anchor_rmsd = data[:, 0]
    fw_rmsd = data[:, 1]
    rmsd = data[:, 2]

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(data[:, 0], data[:, 1], 'bo')
    ax.set_title('RMSD vs framework for %s'%code)
    ax.set_xlabel('Combined anchor RMSD')
    ax.set_ylabel('RMSD')
    ax.set_ylim(1, 4) 

    out_file = os.path.join('figures', 'rmsd_fw_plot_%s.eps'%code)
    fig.savefig(out_file)
     
    os.system('epstopdf %s'%out_file)

    fig = plt.figure()
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(data[:, 0], data[:, 2], 'bo')
    ax.set_title('RMSD vs whole structure for %s'%code)
    ax.set_xlabel('Combined anchor RMSD')
    ax.set_ylabel('RMSD')
    ax.set_ylim(0, 2) 

    out_file = os.path.join('figures', 'rmsd_plot_%s.eps'%code)
    fig.savefig(out_file)
     
    os.system('epstopdf %s'%out_file)

os.system('rm figures/*.eps')
