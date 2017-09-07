from Common.Common import bulk_number

#How many numbered sequences to hold in a single file.
chunk_size = 1

#####Dealing with processed data#############

#Receives a dictionary of sequences to be saved mapped to numbering.
def save_data():
	pass

#Update the antibody portion of the DB.
def update_sabdab():
	from DataManagement.SAbDab import fetch_sabdab_seqs
	
	antibodies = fetch_sabdab_seqs()
	#Translate info format suitbale for numbering -- format the name as PDBchainHchainL
	heavies = {}
	lights = {}
	for ab in antibodies:
		heavy_nm = ab['pdb']+ab['pdb-h']
		light_nm = ab['pdb']+ab['pdb-h']
		heavies[heavy_nm] = ab['H']
		lights[light_nm] = ab['H']
		
		if len(lights) > chunk_size:
			print "[SABDAB Update] Bulk Numbering..."
			numbered = bulk_number(lights)
	#See if we have some leftover sequences.
	
if __name__ == '__main__':

	import sys
	command = sys.argv[1]

	#Update the sabdab datbaase
	if command == 'update_sabdab':
		 update_sabdab()
		 
		 
	#DEUBGGING
