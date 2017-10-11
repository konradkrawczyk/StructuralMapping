#Module for aligning the sequences.
import pickle

#Attempt to align query to template.
#query therefore can be a subset of a longer template and will have 100% sid.
#Region stands for where to perform the alignment.
#False - framework
#H1 - CDR H3
#L2 - CDR L2
#etc...
def align_sequences(query, template,region=None):
	
	#Number of residues that align.
	aligned = 0
	#NUmber of residues that we actually consider (ie regions like framework)
	considered = 0
	for elem in sorted(query):
		#Check if this is the portion of the antibody that we are interested in. (region=None means we look at the entire sequence)
		if region== None or query[elem][1] in region:
			considered+=1
		else:#Not interested in this region.
			continue
		if elem in template:

			if template[elem][0] == query[elem][0]:
				#print "Aligning",elem
				#print "T:",template[elem]
				#print "Q:",query[elem]
				aligned+=1
		
				
	#If there were no residues that we could look at, ignore.
	if considered == 0:
		return 0
	sid = int((float(aligned)/float(considered))*100)
	#print aligned,'/',considered
	#print "sid",sid
	#print '======'
	
	return sid

#Strictly for debugging.
if __name__ == '__main__':

	#Test aligning
	sqs = pickle.load(open('../../data/numbered/sabdab/704057578313450684'))
	
	align_sequences(sqs['2dqiH'][0],sqs['4z5rS'][0],region="H3")

	pass
