from Common.Common import bulk_number,number_sequence,structural_map_location
from DataManagement.DataHandling import save_data
from StructuralAlignment import  align_single_sequence
from DataManagement.SAbDab import structural_reference
import json,pprint
from os.path import join
from os import listdir
import os
import pickle
#How many numbered sequences to hold in a single file.
chunk_size = 100

#Given raw sequences, number them and save.
def number_and_save(data,experiment,chunk="random_chunk_name"):

	chunk_results = {}

	for sequence_name in data:
		
		sequence = data[sequence_name]
		query = number_sequence(sequence)
		full_results,frame_results,cdr_results,_fread_results = align_single_sequence(query, strucs)
		fread_results = {}
		#For fread, we only want the top prediction to assess if there is anything there.
		try:
			for cdr in _fread_results:
				fread_results[cdr] = _fread_results[cdr][0][1]
			
		except:
			fread_results = None

		alignment_results = {'full':full_results,'frame':frame_results,'cdr':cdr_results,'fread':fread_results}
		chunk_results[sequence_name] = alignment_results
	
	
	print "Dumping structural map results for chunk",chunk
	#Json-format the results so they can be easily interpreted by different pieces of software.
	js = json.dumps(chunk_results)
	f = open(join(structural_map_location,experiment,chunk+'.json'),'w')
	f.write(js)
	f.close()

	#numbered_results = bulk_number(data)
	#Save the chunk
	#save_data(experiment,numbered_results,filename=filename)

#Update the antibody portion of the DB.
def update_sabdab():

	from DataManagement.SAbDab import fetch_sabdab_seqs
	
	antibodies = fetch_sabdab_seqs()
	#Translate info format suitbale for numbering -- format the name as PDBchainHchainL
	heavies = {}
	lights = {}
	prog = 0
	for ab in antibodies:
		
		prog+=1
		print "[DataProcessing.py] Renumbering SABDAB ",prog,'out of ',len(antibodies),'antibodies'
		if ab['H']!=None:
			heavy_nm = ab['pdb']+ab['pdb-h']
			heavies[heavy_nm] = ab['H']
		if ab['L']!=None:
			light_nm = ab['pdb']+ab['pdb-l']
			lights[light_nm] = ab['L']
		#Number lights.
		if len(lights) > chunk_size:
			number_and_save(lights,'sabdab')
			lights = {}
		#Number lights.
		if len(heavies) > chunk_size:
			number_and_save(heavies,'sabdab')
			heavies = {}

	#See if we have some leftover sequences.
	#Number lights.
	if len(lights) > 0:
		number_and_save(lights,'sabdab')
	#Number lights.
	if len(heavies) > 0:
		number_and_save(heavies,'sabdab')

#Parse out the sequences from a fasta file, number and save them
def parse_fasta_and_number(experiment_name,fasta_location,start,finish):
	sequences = {}
	curr_name = ""
	
	#How far did we get through the file.
	prog = 0
	
	#used for checkpointing,seeing when the calcaultion stopped.
	curr_prog = 0
	#Find the latest if any.
	if os.path.exists(join(structural_map_location,experiment_name)):
		
		for fname in listdir(join(structural_map_location,experiment_name)):
			fname = int(fname.replace('.json',''))
			print fname,start,finish
			if fname<start or finish<fname:
				continue
			else:
				if fname>start and finish>=fname and fname>curr_prog:
					curr_prog = fname
	
	for line in open(fasta_location):
		if 'QVQLVQSGPGLVKPSQTLSLTCAISGDIVSSNNAAWNWIRQSPSRGLEWLGRTYYRSKWYNDYAVSVKSRITINPDTSKNQFSLQLNSVTPEDTAVYYCVRDNLSTGHREFDYWGQGPLVTVSS' in line:
			print prog
		if '>' in line:#Count progress only by the sequence line, to make sure we start with the next one.
			prog+=1
		if prog<start or prog<curr_prog:
			continue
		if prog>=finish:
			break
		#print "[DataProcessing.py] Custom dataset parsing ",prog,'sequences read... '
		line = line.strip()
		if '>' in line:
			curr_name = line.replace('>','')
		else:
			
			sequences[curr_name] = line
		if len(sequences)> chunk_size:
			number_and_save(sequences,experiment_name,str(min(prog,finish)))

			sequences = {}
	
	#See if we got any leftovers
	if len(sequences)> 0:
		number_and_save(sequences,experiment_name,str(min(prog,finish)))
	
if __name__ == '__main__':

	import sys
	command = sys.argv[1]
	
	#Update the sabdab datbaase
	#Usage: python DataProcessing.py update_sabdab
	#Saves pickled numbered files into [datadirectory]/numbered/sabdab/
	if command == 'update_sabdab':
		update_sabdab()
	#Create a numbered dataset from fasta file.
	#Usage: python DataProcessing.py number_dataset [experiment_name] [fasta location]
		#Saves pickled numbered files into [datadirectory]/numbered/[experiment_name]/
	if command == 'process_dataset':

		exp_name = sys.argv[2]
		fasta_location = sys.argv[3]

		#For parallelization.
		start = 0
		finish = 10000000000
		if len(sys.argv) >5:
			start = sys.argv[4]
			finish = sys.argv[5]
		
		#Load the structural reference.
		strucs = structural_reference()


		results_directory = join(structural_map_location,exp_name)
		#Check if directory to save results exists.
		if not os.path.exists(results_directory):
			os.mkdir(results_directory)

		parse_fasta_and_number(exp_name,fasta_location,int(start),int(finish))

	#Create scripts for parallel processing
	#python DataProcessing.py parallel [exp_name] [fasta_location] number-of-cpus number-of-seqs
	if command == 'parallel':
		experiment_name = sys.argv[2]
		fasta_location = sys.argv[3]
		
		ncpus = float(sys.argv[4])
		nseqs = float(sys.argv[5])
		dseqs = int(float(nseqs)/float(ncpus))
		
		done = 0
		while done< nseqs:
			
			print 'python DataProcessing.py process_dataset ',experiment_name,' ',fasta_location,' ',done,done+dseqs,' &'
			done+=dseqs
	

	#
		 
	#DEUBGGING
