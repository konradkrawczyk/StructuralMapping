#! /usr/bin/env python
# -*- coding: utf-8 -*-
#

import sys
import os
import datetime
from optparse import OptionParser
where_am_i = os.path.dirname(os.path.realpath(__file__))

import numpy # Just to test if we have numpy installed.

# Make sure we can import stuff from this file's directory
sys.path.append(os.path.abspath(os.path.dirname(sys.argv[0])))
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(os.path.dirname(sys.argv[0])), "lib-python")))
#sys.path.extend(os.environ['PATH'].split(':'))

from prosci.loops.fread import Fread

scriptpath = os.path.dirname(sys.argv[0])

VERSION = "PyFREAD 1.02 (build "+str(datetime.date.fromtimestamp(os.path.getmtime(sys.argv[0])))+")"


longVersion = """
                              +---------+
+-----------------------------+ PyFREAD +------------------------------+
                              +---------+

          PyFREAD : HIGH-SPEED FRAGMENT-BASED LOOP MODELLING
        copyright 2011-2012 Sebastian Kelm <kelm@stats.ox.ac.uk>

              Version: {version}

      The "prosci" Python library, used extensively by PyFREAD, is
        copyright 2007-2012 Sebastian Kelm

      PyFREAD's implementation of Cyclic Co-ordinate Descent is
        copyright 2011 Yoonjoo Choi        (prototype)
        copyright 2011 Sebastian Kelm      (C++ implementation)

      CCD for protein loop closure was first described in
        Canutescu AA, Dunbrack RL (2003) Protein Sci., 12(5):963–972

      PyFREAD is based on the FREAD algorithm described in
        Choi Y, Deane CM (2010) Proteins, 78(6):1431-40
        Deane CM, Blundell TL (2001) Protein Sci., 10(3):599-612

                          All rights reserved.
                This program is free for academic use.
           Commercial use requires the author's permission.

          For feedback, bug reports and questions contact:
                Sebastian Kelm <kelm@stats.ox.ac.uk>

""".format(version=VERSION)


USAGE="""
pyfread [OPTIONS] <db_dir> <pdb_file> <start_residue> <loop_sequence>

<db_dir>
    Root directory of PyFREAD's loop database

<pdb_file>
    Query PDB file with or without missing bits. PyFREAD detects
    automatically if your specified loop is present in the file. These
    coordinates are then ignored, except for reporting each decoy's RMSD
    to this structure.

<start_residue>
    The number (and insertion code) of the residue BEFORE the beginning of
    the loop. This residue is always present in the input model. If the
    residue has an insertion code, <start_residue> can have the format 51A,
    where 51 is the residue number and A is the insertion code.
    The meaning of this argument can be changed to a Python-style array index
    using the --arrayindex option.

<loop_sequence>
    Literally, the amino acid sequence of the loop, e.g. "LTPVA"

"""

outputHelp = """
                             +----------+
+----------------------------+  OUTPUT  +------------------------------+
                             +----------+

PyFREAD outputs a summary of the database hits it found. This is written
to STDOUT by default, but can be redirected to a file using the option
--summary FILENAME.

In addition, the 3D structure of each database hit is returned, but only
if the option --strucdir DIRNAME is specified. All structure files will
be placed in the directory DIRNAME.

The result summary provides a ranking of database loops, with the first
loop decoy being the top prediction. Decoy names refer to the filename
of the loop's 3D structure.

This summary might look a little like this:

1jgyM_254_8	61	0.037	0.095	0
2gmrM_251_8	50	0.063	0.122	0
2j8dM_254_8	57	0.033	0.147	0
2uxmM_254_8	39	0.069	0.151	0
1dv6S_252_8	33	0.036	0.156	0

Here's what each column means:

1 : decoy_name : filename of the loop's 3D structure
2 : score : environment-specific substitution score
3 : iRMSD : "internal" RMSD of C-alpha separations of loop anchors
4 : RMSD : all-backbone-atom RMSD of loop anchors
5 : CCD : number of CCD iterations needed to close the loop
          (always 0 unless the --close option was specified)

NOTE: If the input structure already contains (native) co-ordinates for
      the query loop, two more columns appear in the results summary:

6 : all-backbone-atom RMSD of the loop, before closure
7 : all-backbone-atom RMSD of the loop, after closure
    (These two numbers will be identical unless --close was specified)

The summary is sorted by: RMSD, -score, iRMSD, CCD
You could recreate the same sort order with this linux command:
  sort -k4n -k2nr -k3n -k5n FILENAME

+----------------------------------------------------------------------+
"""

def run_fread(dbdir,pdb_file,start_of_loop,loop_sequence,chain,output_file):
	
    
 
    
    loop_sequence = loop_sequence.strip().upper()
    
    frd = Fread()
    
    frd.open_errstream('-')
    # {    'strucdir': 'decoys', 'vdw_factor': 0.7, 'messages': '-', 'summary': 'summary.table', 'write_decoys': False, 'idfilter': '', 'first_decoys': 3}>

    frd.score_cutoff = 25#Default value.
    frd.open_rmsd_cutoff = 1.0#Defulat value
    frd.closed_rmsd_cutoff = 0.3#Default value
    frd.vdw_factor = 0.7#Default valueclo
    frd.close = True#Default
    frd.meld = True#Default
    frd.verbose = False
    frd.max_melding_rmsd = 1.0#Default
    frd.mutate = False #Default value
    frd.calculate_contacts = False#Default
    frd.contact_identity = 0.8#Default
    
    frd.set_subst_tables(os.path.join(where_am_i,'esst.txt'))#Default
    
    frd.set_structure(pdb_file)
   
    return frd.automodel_loop(start_of_loop, loop_sequence,
        loopsuffix = "",#Default
        resnum = True,#TODO this is NOT as in defualts, but I do not see why we would always got for 'array numbering' rather than native array numbering.
        chain = chain,
        structure_known = None,#TODO I am not sure if this is default.
        dbdir = dbdir,
        strucdir = 'decoys',#Default
        summary_file_name = output_file,#Not default
        max_decoys = 100,#Defualt
        first_decoys = 100,#Default
        extension_size = 1000000,#Default
        extension_minimum = 0,#Default
        write_decoys =False,#Default
        write_models = False,#Default
        write_decoy_sequence = False,#Default
        idfilter = "")#default 
  

def main():
    
    parser = OptionParser(version=VERSION, usage=USAGE)
    
    #allowed=("arrayindex", "verbose", "noclose", "nomeld", "version", "nostrucdir", "nodecoys", "nomodels", "loop_sequence", "chain", "mutate")
    
    #withargument=("strucdir", "esst", "summary", "messages", "score", "open_rmsd", "closed_rmsd", "vdw_factor", "extension", "extension_min", "max_decoys", "first_decoys", "max_melding_rmsd", "chain", "structure_known")
    
    parser.add_option("--copyright", dest="copyright", action="store_true",
                      help="Display copyright message and exit.")

    parser.add_option("--output_help", dest="output_help", action="store_true",
                      help="Display a description of the program's output, then exit.")
    

    parser.add_option("--arrayindex", dest="resnum", action="store_false",
                      default=True,
                      help="Change the meaning of the <start_residue> argument to refer to an array index instead of a residue number. Array indeces start from 0 and number all residues within the model continuously, ignoring gaps entirely. The number should refer to the index of the first loop residue (which is usually absent from the input model). In Python terms, this is a slice index pointing to the boundary between the N-terminal loop anchor and the loop itself.")

    parser.add_option("--chain", dest="chain",
                      help="Look for <start_residue> in the chain specified by CHAIN (a single-character chain identifier). Only has an effect if --arrayindex is not specified. If --chain is not specified, the first residue with the matching residue number and insertion code will be used.")

    parser.add_option("--mutate", dest="mutate", action="store_true", default=False,
                      help="Allow the input loop sequence to differ from the residues in the input model. Only has an effect if the input model already contains 3D co-ordinates for the loop to be modelled.")

    parser.add_option("--structure_known", dest="structure_known", type="int",
                      help="Specifying this option disables auto-detection of whether the loop's structure is present in the input model. It can have a value of 0 or 1.")

    parser.add_option("--esst", dest="esst",
                      default=os.path.join(scriptpath, 'esst.txt'),
                      help="Load environment-specific substitution tables from specified file, rather than from the default location ('esst.txt' in pyfread directory).")

    parser.add_option("--strucdir", dest="strucdir", default="decoys",
                      help="Create DIRECTORY and write all decoy structures to files therein. Setting DIRECTORY to the empty string has the same effect as --nostrucdir.")

    parser.add_option("--nostrucdir", dest="strucdir", action="store_const", const=None,
                      help="Disables structure output. Overrides --strucdir. Same as --nomodels --nodecoys.")

    parser.add_option("--nomodels", dest="write_models", action="store_false",
                      default=True,
                      help="Disables output of 'completed' models, i.e. given framework+decoy.")

    parser.add_option("--nodecoys", dest="write_decoys", action="store_false",
                      default=True,
                      help="Disables output of loop decoy structures.")

    parser.add_option("--summary", dest="summary", default="summary.table", 
                      help="Write search summary to specified file. A value of '-' will direct output to STDOUT.")

    parser.add_option("--loopsuffix", dest="loopsuffix", default="", 
                      help="Add SUFFIX to the end of the paths specified by --strucdir DIRECTORY and --summary FILENAME (before the file extension).")
                      
    parser.add_option("--loop_sequence", dest="write_decoy_sequence", default=False,
                      action="store_true", help="Include each database fragment's sequence in the result summary.")

    parser.add_option("-v", "--verbose", dest="verbose", default=False,
                      action="store_true", help="Print status messages to the file specified by --messages.")

    parser.add_option("--messages", dest="messages", default="-",
                      help="Redirect status messages to the specified file.")

    parser.add_option("--max_decoys", dest="max_decoys",
                      type="int", default=100,
                      help="Specify the maximum number of decoys to retain for each gap in the input structure. Higher numbers mean longer running times.")

    parser.add_option("--first_decoys", dest="first_decoys", type="int", 
                      default=0, 
                      help=" Stops the search after finding INT decoys matching the search criteria. This means that the results will not necessarily be the best possible ones present in the database, but will result in significant speed-ups when cut-offs are loose.")

    parser.add_option("--open_rmsd", dest="open_rmsd", type="float", 
                      default=1.0,
                      help="Set the anchor RMSD cut-off to FLOAT.")

    parser.add_option("--score", dest="score", type="int",
                      default=25,
                      help="Set substitution score cut-off to specified integer. If the given sequence cannot possibly attain this cut-off score, the cut-off is lowered to allow perfect matches to be found.")

    parser.add_option("--nomeld", dest="meld", default=True, action="store_false",
                      help="Do not meld loop anchors. This means all loops are closed instead, unless the --noclose option is given.")

    parser.add_option("--noclose", dest="close", default=True, action="store_false", 
                      help="Do not close loops. This means if you set --open_rmsd too high, you will get loops that do not fit the input model's anchors.")

    parser.add_option("--closed_rmsd", dest="closed_rmsd", type="float", default=0.3,
                      help="Set the anchor RMSD cut-off after loop closure to specified float (default is 0.3). Note that the target anchor RMSD during loop closure is 0.15. Setting this cut-off lower than 0.15 will thus result in almost zero coverage.")

    parser.add_option("--max_melding_rmsd", dest="max_melding_rmsd", type="float",
                      default=1.0, 
                      help="Set the RMSD cut-off that decides when to use melding and when to use loop closure.")

    parser.add_option("--vdw_factor", dest="vdw_factor", type="float", 
                      default=0.7,
                      help="Set factor for Van der Waal's radii to FLOAT. If a clash is detected, the loop is discarded. A value <= 0 disables clash checking.")

    parser.add_option("--extension_max", "--extension", dest="extension_max", type="int", 
                      default=1000000,
                      help="If no hits are found, the loop is extended by one residue on each end and the search repeated. INT is the maximum number of times this extension should happen before giving up. A value of INT <= 0 disables loop extension.")

    parser.add_option("--extension_min", dest="extension_min", type="int",
                      default=0,
                      help="Extend loop INT times before the first modelling attempt. Do not even try searching the database for matches shorter than that. This does not count towards the limit imposed by --extension.")
    
    parser.add_option("--contacts", dest="contacts", default=False, action="store_true",
                      help="Calculate residue contacts for each loop within the input structure and compare them to the fragment's native contacts.")
    
    parser.add_option("--contact_identity", dest="contact_identity", type="float",
                      default=0.8,
                      help="Minimum fraction of matching contacts needed to accept a loop.")
    
    parser.add_option("--idfilter", dest="idfilter", type="string",
                      default="",
                      help="FREAD summary file (or list of FREAD loop IDs separated by newlines) listing the only loops to be accepted in this FREAD run. This is useful for first running FREAD without structure output, deciding on some loops of interest, and then re-running FREAD to obtain the 3D models of only a few interesting loops.")
    
    (opts, args) = parser.parse_args()
    
    if opts.copyright:
        print longVersion
        sys.exit(1)
    
    if len(args) != 4:
        print longVersion
        parser.print_help()
        print outputHelp
        sys.exit(1)
    
    dbdir, pdb_file, start_of_loop, loop_sequence = args
    
    loop_sequence = loop_sequence.strip().upper()
    
    if opts.structure_known is not None:
      opts.structure_known = bool(opts.structure_known)
    
    idfilter = []
    if opts.idfilter:
      f = open(opts.idfilter)
      for line in f:
        fields = line.split(None, 1)
        if fields:
          idfilter.append(fields[0])
      f.close()
    
    
    frd = Fread()
    
    frd.open_errstream(opts.messages)
    
    frd.score_cutoff = opts.score
    frd.open_rmsd_cutoff = opts.open_rmsd
    frd.closed_rmsd_cutoff = opts.closed_rmsd
    frd.vdw_factor = opts.vdw_factor
    frd.close = opts.close
    frd.meld = opts.meld
    frd.verbose = opts.verbose
    frd.max_melding_rmsd = opts.max_melding_rmsd
    frd.mutate = opts.mutate
    frd.calculate_contacts = opts.contacts
    frd.contact_identity = opts.contact_identity
    
    frd.set_subst_tables(opts.esst)
    
    frd.set_structure(pdb_file)
    
   
    return frd.automodel_loop(start_of_loop, loop_sequence,
        loopsuffix = opts.loopsuffix,
        resnum = opts.resnum,
        chain = opts.chain,
        structure_known = opts.structure_known,
        dbdir = dbdir,
        strucdir = opts.strucdir,
        summary_file_name = opts.summary,
        max_decoys = opts.max_decoys,
        first_decoys = opts.first_decoys,
        extension_size = opts.extension_max,
        extension_minimum = opts.extension_min,
        write_decoys = opts.write_decoys,
        write_models = opts.write_models,
        write_decoy_sequence = opts.write_decoy_sequence,
        idfilter = idfilter)
    
    
    


if __name__ == "__main__":
  main()
