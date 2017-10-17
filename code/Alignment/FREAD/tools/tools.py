import pdb
import difflib
import subprocess
from prosci.util.residue import ResidueList

def fasta2seq(fasta, ref_struc):
    """Return a sequence of the chains in fasta, referencing ref_struc."""
    s = open(ref_struc).readlines()
    f = open(fasta).readlines()
    
    for line in s:
        if line.split()[0] == 'ATOM':
            firstChain = line.split()[4]
            break
        else:
            continue
    
    
    if firstChain == 'H':
        seq = f[1].strip('\n') + f[3].strip('\n')
    if firstChain == 'L':
        seq = f[3].strip('\n') + f[1].strip('\n')
   
    return seq 

def getFitCommands(ref=None, mobile=None, zones=None):
    """Get matching zones and output commands for ProFit RMSD calculation.
    
    Determine fit commands for calculating RMSD with ProFit. Structures
    ref and model can either be .pdb files or ResidueList objects.

    To specify zones manually, each zone should be in the form (i, j, n), 
    where i and j are the residue indices at which the zone starts in ref
    and mobile, respectively, and n is the number of residues in the zone.
    """
    if not zones:
        ref, model = ResidueList(ref).get_seq(), ResidueList(mobile).get_seq()
        s = difflib.SequenceMatcher(None, ref, model)
        zones = s.get_matching_blocks()[:-1]

    if not isinstance(zones, list):
        zones = [zones]

    cmd = 'ATOMS CA,O,N,C \nNUMBER SEQUENTIAL \n'
    for z in zones:
        (i, j, n) = z
        cmd += 'ZONE %i-%i:%i-%i \n' %(i+1, i+n+1, j+1, j+n+1)
    cmd += 'FIT\n'

    return cmd

def getRMSD(ref, mobile, zones=None, fitCommands=None):
    """Calculate the RMSD between ref and mobile structures using ProFit."""

    p = subprocess.Popen(['profit', ref, mobile], stdin=subprocess.PIPE,
                         stdout=subprocess.PIPE)

    if not fitCommands:
        fitCommands = getFitCommands(ref, mobile, zones)

    fit = p.communicate(fitCommands)[0]

    return float(fit.split('\n')[-2].split()[-1])
