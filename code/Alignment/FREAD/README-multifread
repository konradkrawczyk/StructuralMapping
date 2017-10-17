
Summary
-------
MultiFREAD is a multi-loop modeling algorithm. The algorithm uses FREAD to 
fetch sets of decoys for individual loops of an incomplete input structure, 
then ranks possible combinations of the decoys by their combined anchor RMSD.
The algorithm returns N structures, where N is a number specified by the user,
and each structure "i" is the input structure with the i'th best combo of
decoys with no detected atomic clashes. 


Compatibility
-------------

multifread is not compatible with python of version <2.6 or </=3.0.


Notes on usage
--------------
    
    -- A full listing of command-line options can be viewed by executing 
                
                        ./multifread --help

    This shows all possible command line options, e.g. RMSD cutoff, sequence
    score cutoff, etc.

    -- Setting the --max_decoys option higher than 10 will result in an 
    extremely long run time, as the number of score calculations performed
    is < O(n^M), where M is the number of loops being modeled and n is 
    max_decoys. However, even using 5 for a structure with six loops will
    result in 15625 combinations of decoys, so there will likely be plenty
    of combinations without atomic clashes.


USAGE
-----

    ./multifread [OPTIONS] /database/path /incomplete/struc/path complete_seq

    Example:

    ./multifread --esst ../esst.txt \
                 --nostrucdir \
                 --max_decoys=5 \
                 --messages messages.log \
                 ../../db_sol5 ./1A0Q.pdb ./1A0Q.seq 

    -- OR --
 
    ./multifread --esst ../esst.txt \
                 --nostrucdir \
                 --max_decoys=5 \
                 --messages messages.log \
                 ../../db_sol5 ./1A0Q.pdb IELTQSPSSLSASLGGKVTITCkasqdikkyigWYQHKPGKQPRLLIHytstllpGIPSRFRGSGSGRDYSFSISNLEPEDIATYYClqyynlrtFGGGTKLEIKRAVQLQESDAELVKPGASVKISCKASgytftdhVIHWVKQKPEQGLEWIGYIspgngdIKYNEKFKGKATLTADKSSSTAYMQLNSLTSEDSAVYLCKRgyyvdyWGQGTTLTVSS

                

    The above example takes the structure 1A0Q.pdb and determines the missing
    loop regions by comparing the sequence of the structure to the complete 
    sequence in 1A0Q.seq. The maximum number of decoys to get for each loop is 
    set to 5. 

    The environment-specific substitution tables are in esst.txt, and the 
    database in which to search for decoys is set to ../../db_sol5. 


