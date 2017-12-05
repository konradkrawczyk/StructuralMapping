from FREAD.pyfread_api import run_fread
from Common.Common import is_CDR

#For fread, one residue BEFORE the start of the loop.
loop_starts = {'L1':23,'L2':49,'L3':88,'H1':25,'H2':51,'H3':94}

#Get the CDR sequences from the numbered object.
def extract_cdrs(numbered):
	cdrs = {}
	for chid in sorted(numbered):
		
		chothia = chid[2]+str(chid[0])
		cdr_code = is_CDR(chothia)#H1, L3 if CDR, False if freamework
		
		if cdr_code != False:
			#means we are a CDR
			aa = numbered[chid][0]
			if cdr_code not in cdrs:
				cdrs[cdr_code] = ""
			cdrs[cdr_code]+=aa
	return cdrs			

#Perform fread search for the given sequence.
#loop -- H1, H2 etc.
#template_pdb -- e.g. 12e8
#templaet_chain -- e.g. P
#sequence -- e.g. 'DPEIGD'
#Result: json-formatted best match {'seq': u'DPEIGD', 'str': u'12e8P', 'scr': 45} or None
def perform_loop_alignment(loop,template_pdb,template_chain,sequence):
	
	#database location
	db = '../data/fread_db/db_CDR'+loop

	#Template location.
	template = '../data/structures/'+template_pdb+template_chain+'_no_cdrs.pdb'
	print template
	results = run_fread(db,template,loop_starts[loop],sequence,template_chain,'')

	#Find the highest score one.
	best_match = {'str':None,'seq':None,'scr':None,'qu':sequence}
	curr_max = -1
	for decoy in results:
		#Debugging...
		#print decoy.struc, decoy.startres, decoy.startinscode, decoy.length, decoy.score,decoy.seq
		if decoy.score> curr_max:
			best_match = {'str':decoy.struc,'seq':decoy.seq,'scr':decoy.score,'qu':sequence}
			curr_max = decoy.score
	return best_match

#Strictly for debugging constituent funcionality.
if __name__ == '__main__':

	import sys
	cmd = sys.argv[1]
	if cmd == 'fetch_loop':
		loop = 'H2'
		sequence = 'DPEIGD'
		print perform_loop_alignment(loop,'12e8','P',sequence)
	if cmd == 'fetch_cdrs':
		import pickle
		query = pickle.load(open('sample_sequence.pckl'))
		print extract_cdrs(query[0])
	
