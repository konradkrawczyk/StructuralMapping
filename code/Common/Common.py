####Common functions & locaions##########
import random,os
#################
#Common functions
#################

#Generates a random filename
def random_filename():
	filename = str(int(1000000*random.random()))+str(int(1000000*random.random()))+str(int(1000000*random.random()))
	return filename
	


#################
#Common locations
#################

#Where the numbered_datasets are stored
numbered_datasets_location = '../data/numbered'
structural_map_location = '../data/structuralmap'


############
#CDR matters
############

definitions = {		
		
		"chothia" : {"L1" : ["L24", "L25", "L26", "L27", "L28", "L29", "L30", "L31", "L32", "L33", "L34"], 
				  "L2" : ["L50", "L51", "L52", "L53", "L54", "L55", "L56"],
				  "L3" : ["L89", "L90", "L91", "L92", "L93", "L94", "L95", "L96", "L97"], 
				  "H1" : ["H26", "H27", "H28", "H29", "H30", "H31", "H32"], 
				  "H2" : ["H52", "H53", "H54", "H55", "H56"] ,
				  "H3" : ["H95", "H96", "H97", "H98", "H99", "H100", "H101", "H102"]} 
		}


#Is this a CDR? If so tell which.
def is_CDR(res,deff='chothia'):
	for CDR in definitions[deff]:
		if res in definitions[deff][CDR]:
			return CDR

	return False


####################
#Numbering by ANARCI
####################

#ANARCI-- number.
from anarci import anarci

#BUlk number sequences. 
#Input: dictionary from sequences ids to raw sequences
#Output: dictionary from sequence ids to numbered dictionaries of chothia ids to aa
def bulk_number(sequences):
	print "[Common.py] Bulk Numbering ",len(sequences),'sequences...'
	numbered_sequences = {}
	for sequence_id in sequences:
		numbered_sequences[sequence_id] = number_sequence(sequences[sequence_id])
		
	return numbered_sequences

#From a mapping from chothia ids to amino acids, get the full sequence.
def get_sequence(s):
	sequence = ""
	for chid in sorted(s):
		
		sequence += s[chid][0]
	return sequence	

#Given a raw sequence, perform anarci annotation.
#Input: raw sequence
#output dictionary from chothia ids to amino acids and chothia cdr annotations
def number_sequence(query_seq):
	#Number the query sequence
    try:
	res = anarci([('q',query_seq)],scheme='chothia',assign_germline=True)
    except:
        return None,None
    #Numbering failed.
    if res[0][0] == None:
	return None,None
    query= res[0][0][0][0]
	
    chain =  res[1][0][0]['chain_type']
	
    #Translate the query into the format we want
    #[chid,insert,chain] = (aa,is_CDR)
    _query = {}
    for sid in sorted(query):
		
		chid = (sid[0][0],sid[0][1].replace(' ',''),chain)
		if sid[1] == '-':
			continue
		_query[chid] = (sid[1],is_CDR(chain+str(sid[0][0])))
	
    query = _query

    return query,chain
	
###Strictly for testing purposes.	
if __name__ == '__main__':

	#Test bulk numbering
	seqs = {
				'1' :'DIKMTQSPSSMYASLGERVTITCKASQDIRKYLNWYQQKPWKSPKTLIYYATSLADGVPSRFSGSGSGQDYSLTISSLES',
				'2' :'DIKMTQSPSSMYASLGERVTITCKASQDIRKYLNWYQQKPWKSPKTLIYYATSLADGVPSRFSGSGSGQDYSLTISSLES',
				'3' :'DIKMTQSPSSMYASLGERVTITCKASQDIRKYLNWYQQKPWKSPKTLIYYATSLADGVPSRFSGSGSGQDYSLTISSLES',
				'4' :'DIKMTQSPSSMYASLGERVTITCKASQDIRKYLNWYQQKPWKSPKTLIYYATSLADGVPSRFSGSGSGQDYSLTISSLES'
		   }
	print bulk_number(seqs)
	
	
	#Test numbering
	print number_sequence('DIKMTQSPSSMYASLGERVTITCKASQDIRKYLNWYQQKPWKSPKTLIYYATSLADGVPSRFSGSGSGQDYSLTISSLES')

