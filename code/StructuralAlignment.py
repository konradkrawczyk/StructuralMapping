import os,pickle
from os import listdir
from os.path import join
#Aligning full variable region and frameworks.
from Alignment.Align import align_sequences
#Aligning CDRs.
from Alignment.LoopAlignment import perform_loop_alignment,extract_cdrs
from Common.Common import get_sequence,numbered_datasets_location,structural_map_location,number_sequence
from DataManagement.SAbDab import structural_reference
import json

#Find the best sequence identity and the corresponding template -- this is chiefly used to get region annotatinos like full sequence, framework etc. 
#CDRs can be done here BUT it is better to use a more nuanced approach using FREAD which is written in get_best_cdr_match
def get_best_match(query,structures,region=None):
	curr_best_sid = 0
	curr_best_pdb = 0		
	j = 0
	for struc in structures:
		
		#Compare heavy to heavy and light to light...
		if query[1] == structures[struc][1]:
			j+=1
			#print struc
			sid = align_sequences(query[0],structures[struc][0],region=region)
			
			if sid> curr_best_sid:
				curr_best_sid = sid
				curr_best_pdb = struc
	#Results object	
	results = {'best_sid':curr_best_sid,'best_pdb':curr_best_pdb}
	return results

#Given anchors, employ FREAD to get the best CDR-templates.
def get_best_cdr_match(query,fread_template):
	#pickle.dump(query,open('sample_sequence.pckl','w'))
	cdrs = extract_cdrs(query[0])
	pdb_template = fread_template[0:4]
	pdb_chain = fread_template[4]

	cdr_results = {}
	for cdr in cdrs:
		cdr_results[cdr] = perform_loop_alignment(cdr,pdb_template,pdb_chain,cdrs[cdr])
	return cdr_results
	
#given a numbered anarci sequence and reference structures, get the best structural hits.
def align_single_sequence(query, structures):
	#Get best full-sequence match.
	full_results = get_best_match(query,structures)
	#Get best framework match. False region stands for 'not CDR' ergo framework
	frame_results = get_best_match(query,structures,region=[False])
	#Get best CDR match -- full region.
	cdr_results = get_best_match(query,structures,region=['H1','H2','H3','L1','L2','L3'])
	#Get the nuanced-fread based hits for each CDR in turn.
	fread_results = None
			
	if 'best_pdb' in full_results:
		try:
			fread_results = get_best_cdr_match(query,full_results['best_pdb'])
		except:
			pass
	return full_results,frame_results,cdr_results,fread_results



if __name__ == '__main__':
	
	import sys
	cmd = sys.argv[1]

	
	if cmd == 'count_structures':#Count how many structures in the currrent release.
		strucs = structural_reference()
		print "We have ",len(strucs)

	
	#Processing a single sequence
	if cmd == 'process_single':
		
		sequence = sys.argv[2]
		#Number the sequence
		#Format of the query: (numbered_sequence,chain)
		query = number_sequence(sequence)
		
		#Load the structural reference.

		strucs = structural_reference()

		print align_single_sequence(query, strucs)

		

	
