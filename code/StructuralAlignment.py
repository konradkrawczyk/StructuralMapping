import os,pickle
from os import listdir
from os.path import join
#Aligning full variable region and frameworks.
from Alignment.Align import align_sequences
#Aligning CDRs.
from Alignment.LoopAlignment import perform_loop_alignment,extract_cdrs
from Common.Common import get_sequence,numbered_datasets_location,structural_map_location
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
	pickle.dump(query,open('sample_sequence.pckl','w'))
	cdrs = extract_cdrs(query[0])
	pdb_template = fread_template[0:4]
	pdb_chain = fread_template[4]
	cdr_results = {}
	for cdr in cdrs:
		cdr_results[cdr] = perform_loop_alignment(cdr,pdb_template,pdb_chain,cdrs[cdr])
	return cdr_results
	
		
#Go through all the chunks in the 
def perform_comparison(sample_name,structures,start,end):
	source = join(numbered_datasets_location,sample_name)
	i = 0
	for chunk in sorted(listdir(source)):
		
		#for parallelization
		i+=1
		if i<start or i>end:
			continue  
		d = pickle.load(open(join(source,chunk)))
		print "Number of seqs in chunk",len(d)
		chunk_results = {}
		for query in d:
			query_name = query
			query = d[query]
			
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
			print fread_results
			
			#print get_sequence(query[0])
			#print 'FULL',full_results
			#print 'Frame',frame_results
			#print 'CDR',cdr_results
			alignment_results = {'full':full_results,'frame':frame_results,'cdr':cdr_results,'fread':fread_results}
			chunk_results[query_name] = alignment_results
		print "Dumping structural map results for chunk",chunk
		#Json-format the results so they can be easily interpreted by different pieces of software.
		js = json.dumps(chunk_results)
		f = open(join(structural_map_location,sample_name,chunk+'.json'),'w')
		f.write(js)
		f.close()

if __name__ == '__main__':
	
	import sys
	cmd = sys.argv[1]

	#Usage: python StructuralAlignment.py parallel [experiment name] number-of-cpus
	if cmd == 'parallel':
		exp_name = sys.argv[2]
		ncpus = float(sys.argv[3])
		#Figure out how many chunks there are to process.
		nchunks = len(listdir(join(numbered_datasets_location,exp_name)))
		dchunk = int(float(nchunks)/ncpus)
		done_chunks = 0
		while done_chunks<nchunks:
			print 'python StructuralAlignment.py structurallymap ',exp_name,' ',done_chunks,' ',done_chunks+dchunk,' &'
			done_chunks+=dchunk

	if cmd == 'structurallymap':
		
		strucs = structural_reference()
		
		exp_name = sys.argv[2]
		start = int(sys.argv[3])
		end =int(sys.argv[4])
		results_directory = join(structural_map_location,exp_name)
		#Check if directory to save results exists.
		if not os.path.exists(results_directory):
			os.mkdir(results_directory)
		print "Got",len(strucs),'structures to compare against.'
	
		perform_comparison(exp_name,strucs,start,end)
	
	if cmd == 'debug':
		strucs = structural_reference()
		
		exp_name = 'sample'
		start = 1
		end =2
		results_directory = join(structural_map_location,exp_name)
		#Check if directory to save results exists.
		if not os.path.exists(results_directory):
			os.mkdir(results_directory)
		print "Got",len(strucs),'structures to compare against.'
	
		perform_comparison(exp_name,strucs,start,end)

	#Load structural reference.

	#For each sequence in the set to be numbered, compare.
		
	#dump results.
	
