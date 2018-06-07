from os.path import join
import pickle
import os
from Common.Common import gene_db,list_file_paths
from Plotting import plot_histogram

#Load the pickled file as a dictionary.
def load_dictionary(fname):
	dic = {}
	d = pickle.load(open(fname))
	for elem in d:
		dic[str(elem[0][0])+elem[0][1]] = elem[1]
	return dic


#Given vregion and jregion, get me the germline annotation - CHothia numbered.
def fetch_germline_alignment(vregion,jregion):
	fname = join(gene_db,vregion+'_'+jregion)
	if os.path.exists(fname):
		return pickle.load(open(fname))
	else:
		print "We do not have the gene",vregion,jregion,'in our germline database'
		return {}
			
if __name__ == '__main__':

	
	import sys
	cmd = sys.argv[1]

	
	if cmd == 'sequence_identities':
		which_chain = sys.argv[2]
		done = []
		hist_data = []
		
		for gene_1 in list_file_paths(gene_db):
			if 'D' in gene_1:
				continue
			if 'IG'+which_chain not in gene_1:
				continue
			g1_data = load_dictionary(gene_1)
			done.append(gene_1)
			print len(done)
			#if len(done) > 5:
			#	break
			for gene_2 in list_file_paths(gene_db):	
				if 'D' in gene_2:
					continue
				if gene_2 in done:
					continue
				if 'IG'+which_chain not in gene_2:
					continue
				if gene_1 == gene_2:
					continue
				
				g2_data = load_dictionary(gene_2)
				matches = 0
				for g1_ch in g1_data:
					try:
						if g1_data[g1_ch] == g2_data[g1_ch]:
							matches+=1
					except KeyError:
						pass
				matches = int(100*(float(matches)/float(len(g1_data))))
				hist_data.append(matches)
				if matches>90:
					print gene_1,gene_2
		plot_histogram(hist_data,which_chain+"_vj.png",xlabel='x',ylabel='y',title='')
				
			
		

	
	
