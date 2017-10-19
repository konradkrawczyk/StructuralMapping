from FREAD.pyfread_api import run_fread

#For fread, one residue BEFORE the start of the loop.
loop_starts = {'L1':23,'L2':49,'L3':88,'H1':25,'H2':51,'H3':94}

#Perform fread search for the given sequence.
#loop -- H1, H2 etc.
#template_pdb -- e.g. 12e8
#templaet_chain -- e.g. P
#sequence -- e.g. 'DPEIGD'
def perform_loop_alignment(loop,template_pdb,template_chain,sequence):
	
	#database location
	db = '../../data/fread_db/db_CDR'+loop

	#Template location.
	template = '../../data/structures/'+template_pdb+'/structure/chothia/'+template_pdb+template_chain+'_no_cdrs.pdb'

	
	
	res = run_fread(db,template,loop_starts[loop],sequence,template_chain,'summary.txt')

	for decoy in res:
		print decoy.struc, decoy.startres, decoy.startinscode, decoy.length, decoy.score,decoy.seq

if __name__ == '__main__':

	loop = 'H2'
	sequence = 'DPEIGD'
	perform_loop_alignment(loop,'12e8','P',sequence)
	
