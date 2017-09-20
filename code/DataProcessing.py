from Common.Common import bulk_number
from DataManagement.DataHandling import save_data
import json,pprint
#How many numbered sequences to hold in a single file.
chunk_size = 100

#Update the antibody portion of the DB.
def update_sabdab():
	from DataManagement.SAbDab import fetch_sabdab_seqs
	
	antibodies = fetch_sabdab_seqs()
	#Translate info format suitbale for numbering -- format the name as PDBchainHchainL
	heavies = {}
	lights = {}
	numbered = {}
	for ab in antibodies:
		if ab['H']!=None:
			heavy_nm = ab['pdb']+ab['pdb-h']
			heavies[heavy_nm] = ab['H']
		if ab['L']!=None:
			light_nm = ab['pdb']+ab['pdb-l']
			lights[light_nm] = ab['L']
		#Number lights.
		if len(lights) > chunk_size:
			print "[DataProcessing.py] Currently numbered sequences",len(numbered),'out of',len(antibodies)
			numbered_results = bulk_number(lights)
			#Save the chunk
			save_data('sabdab',numbered_results)
			for sequence in numbered_results:
				numbered[sequence] = numbered_results[sequence]
	#TODO See if we have some leftover sequences.
	
if __name__ == '__main__':

	import sys
	command = sys.argv[1]

	
	#Update the sabdab datbaase
	if command == 'update_sabdab':
		 update_sabdab()
		 
		 
	#DEUBGGING
