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
	
	aligned = 0
	for elem in sorted(query):
		
		if elem in template:
			
			if template[elem][0] == query[elem][0] and (region== None or query[elem][1] == region):
				#print "Aligning",elem
				#print "T:",template[elem]
				#print "Q:",query[elem]
				aligned+=1
	
	sid = int((float(aligned)/float(len(query)))*100)
	print aligned,'/',len(query)
	print "sid",sid
	print '======'
	
	return sid

#Strictly for debugging.
if __name__ == '__main__':

	#Test aligning
	sqs = pickle.load(open('../../data/numbered/sabdab/704057578313450684'))
	
	align_sequences(sqs['2dqiH'][0],sqs['4z5rS'][0],region="H3")

	pass
