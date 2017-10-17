#Handling some of the tasks related to pre-processing sabdab data before structural mapping.
from Common.Common import structures_location,is_CDR

import os
from os.path import join
from StructuralAlignment import structural_reference
from DataManagement.SAbDab import structural_reference

#All paths are relative wrt root directory.
structures_location = '../'+structures_location

#Given a PDb file -- remove the CDRs from it (according to Chothia definitions)
#HL_chain -- indicate if the chain is heavy (H) or light (L)

def constrain_structure(HL_chain,pdb_chain,pdb):

	#Chothia structures are stored in their special folder.
	pdb_location = join(structures_location,pdb,'structure','chothia',pdb+'.pdb')
	#Our output file is in the same chothia location but called differently.
	output_location = join(structures_location,pdb,'structure','chothia',pdb+pdb_chain+'_no_cdrs.pdb')
	#Do not redo computational effort.
	if os.path.exists(output_location):
		return
	output_file = open(output_location,'w')
	for line in open(pdb_location):
		if line[0:4]!='ATOM':
			continue
		line = line.strip()
		#print "Chain:",line[21]
		chain = line[21]
		if chain!=pdb_chain:
			continue
		#print line
		#print "Sid:",line[22:28].replace(" ","")
		sid = line[22:26].replace(" ","")
				
		sid = str(int(sid))#Sanity check that we are not picking up insertions -- is_CDR function does not want these.
		
		chid = HL_chain+sid

		if is_CDR(chid) == False:
			output_file.write(line+'\n')

	output_file.close()

		#Check if this residue is a Chothia CDR.

		

if __name__ == '__main__':
	#constrain_structure('H','B','1ahw')
	strucs = structural_reference()
	
	for s in strucs:
		pdb = s
		pdb_chain = s[4]
		target_chain = strucs[s][1]
		if target_chain == 'K':
			target_chain = 'L'
		pdb = pdb[0:4]
		print pdb,pdb_chain,target_chain
		try:
			constrain_structure(target_chain,pdb_chain,pdb)
		except:
			pass
