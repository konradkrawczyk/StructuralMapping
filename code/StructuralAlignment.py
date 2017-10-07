import os,pickle
from os import listdir
from os.path import join
from Alignment.Align import align_sequences
from Common.Common import get_sequence

#This can be loaded in whole.
def structural_reference():
	source = '../data/numbered/sabdab/'
	all_strucs = {}
	for chunk in listdir(source):
		print "Loading",chunk
		d = pickle.load(open(join(source,chunk)))
		for pdb in d:
			all_strucs[pdb] = d[pdb]
		
	return all_strucs

#Go through all the chunks in the 
def perform_comparison(sample_name,structures):
	source = join('../data/numbered/',sample_name)
	i = 0
	for chunk in listdir(source):
		d = pickle.load(open(join(source,chunk)))
		
		for seq in d:
			curr_best = 0
			
			j = 0
			for struc in structures:
				#COMpare heavy to heavy and light to light...
				if d[seq][1] == structures[struc][1]:
					j+=1
					#print struc
					sid = align_sequences(d[seq][0],structures[struc][0])
					if sid> curr_best:
						curr_best = sid
				print get_sequence(d[seq][0])
			print seq,curr_best
			print "Checked against",j,'strucs'
			quit()

#Results folder -- holds data on how many sequences aligned how well.
results_folder = '../data/results/'
if __name__ == '__main__':
	
	strucs = structural_reference()
	print "Got",len(strucs),'structures to compare against.'
	perform_comparison('KELLY_HEPB_BOOSTER',strucs)
		

	#Load structural reference.

	#For each sequence in the set to be numbered, compare.
		
	#dump results.
	
