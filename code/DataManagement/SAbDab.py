import pprint
from Common.Common import numbered_datasets_location
import os
from os.path import join
from os import listdir
import pickle
####################
##Sequence Retrieval
####################

#Everything is one level in relation to root.
#numbered_datasets_location = numbered_datasets_location

#For now this is a dummy before a suitable API is available.
def fetch_sabdab_seqs():
	print "[SABDAB.py] Fetching sequences"
	from ABDB import database
	structures = []
	
	for pdb in database:
		pdb_details = database.fetch(pdb) # get the details of a pdb in the database
		raw_seqs =  pdb_details.get_raw_sequences()
		single_fab = {'pdb':pdb,'pdb-h':None,'pdb-l':None,'H':None,'L':None}
		for fab_details in pdb_details.get_fabs():
			
			if fab_details.VH!='NA':
				single_fab['pdb-h'] = fab_details.VH
				single_fab['H'] = raw_seqs[fab_details.VH]
			if fab_details.VL!='NA':
				single_fab['pdb-l'] = fab_details.VL
				single_fab['L'] = raw_seqs[fab_details.VL]
			
		structures.append(single_fab)

	return structures

#This can be loaded in whole.
def structural_reference():
	source = join(numbered_datasets_location,'sabdab')
	all_strucs = {}
	for chunk in listdir(source):
		print "Loading Structures chunk",chunk
		d = pickle.load(open(join(source,chunk)))
		for pdb in d:
			all_strucs[pdb] = d[pdb]
		
	return all_strucs

if __name__ == '__main__':
	strucs = fetch_sabdab_seqs()
	pprint.pprint(strucs)
