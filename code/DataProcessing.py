from Common.Common import bulk_number
from DataManagement.DataHandling import save_data
import json,pprint
#How many numbered sequences to hold in a single file.
chunk_size = 100

#Given raw sequences, number them and save.
def number_and_save(data,experiment):
	
	numbered_results = bulk_number(data)
	#Save the chunk
	save_data(experiment,numbered_results)

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
	prog = 0
	for line in open(fasta_location):
		prog+=1
		if prog<start or prog>finish:
			continue
		print "[DataProcessing.py] Custom dataset parsing ",prog,'sequences read... '
		line = line.strip()
		if '>' in line:
			curr_name = line.replace('>','')
		else:
			sequences[curr_name] = line
		if len(sequences)> chunk_size:
			number_and_save(sequences,experiment_name)
			sequences = {}
	#See if we got any leftovers
	if len(sequences)> 0:
		number_and_save(sequences,experiment_name)
	
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
	if command == 'number_dataset':
		experiment_name = sys.argv[2]
		fasta_location = sys.argv[3]
		start = sys.argv[4]
		finish = sys.argv[5]
		parse_fasta_and_number(experiment_name,fasta_location,int(start),int(finish))
	if command == 'parallel':
		experiment_name = sys.argv[2]
		fasta_location = sys.argv[3]
		
		ncpus = float(sys.argv[4])
		nseqs = float(sys.argv[5])
		dseqs = int(float(nseqs)/float(ncpus))
		
		done = 0
		while done< nseqs:
			
			print 'python DataProcessing.py number_dataset ',experiment_name,' ',fasta_location,' ',done,done+dseqs,' &'
			done+=dseqs
		 
		 
	#DEUBGGING
