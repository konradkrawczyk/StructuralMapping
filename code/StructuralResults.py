##This file deals with compiling the results of structurally mapped sequencing experiments.
import json,os,pprint
from os import listdir
from os.path import join
from Common.Common import structural_map_location,aggregates_location

#From a single map file, create an aggreagate of sequence identitities.
def process_single_map(loc):
	
	#cdr, frame and full results count the number of percentages (e.g. how many times we got sequence identity of 83%)
	#fread holds the cdrs:
	#h1->
	#--->sequence
	#---->#total,pdb_found,#modelable (as in certain cases because of anchors we might not have been able to model?)	
	results = {'cdr':{},'frame':{},'full':{},'fread':{}}
	

	data = json.load(open(loc))

	for sequence in data:
		
		for region in data[sequence]:
			if region=='fread':
						
				if data[sequence][region]!=None:
					for cdr in data[sequence][region]:
						fread_data = data[sequence][region][cdr]
						
						cdr_sequence = fread_data['qu']
						if cdr not in results['fread']:
							results['fread'][cdr] = {}
						if cdr_sequence not in results['fread'][cdr]:
							results['fread'][cdr][cdr_sequence] = {'total':0,'pdb':False,'model':0}
						results['fread'][cdr][cdr_sequence]['total']+=1
						if fread_data['str'] !=None:
							results['fread'][cdr][cdr_sequence]['model']+=1
							#If the structure we found is identical to the query, say we have a pdb hit.
							if fread_data['qu'] == fread_data['seq']:
								results['fread'][cdr][cdr_sequence]['pdb'] = True
						
				
			else:
				sid = data[sequence][region]['best_sid']
				if sid not in results[region]:
					results[region][sid] = 0
				results[region][sid]+=1
	return results

#Given a dict from regions to sid list, fold the results of the new such list into the current one
def merge_aggregates(current,new):

	for region in new:
		if region == 'fread':
			for cdr in new[region]:
				if cdr not in current[region]:
					current[region][cdr] = {}
				for sequence in new[region][cdr]:
					if sequence not in current[region][cdr]:
						current[region][cdr][sequence] = {'total':0,'pdb':False,'model':0}
					current[region][cdr][sequence]['total']+=new[region][cdr][sequence]['total']
					current[region][cdr][sequence]['model']+=new[region][cdr][sequence]['model']
					if new[region][cdr][sequence]['pdb'] == True:
						current[region][cdr][sequence]['pdb'] = True
					
		else:

			for sid in new[region]:
				if sid not in current[region]:
					current[region][sid] = 0
				current[region][sid]+=new[region][sid]
	return current

#Process the results of a single experiment
def process_experiment(exp_name):
	results_location = join(structural_map_location,exp_name)
	i = 0
	aggregate = {'cdr':{},'frame':{},'full':{}}

	for f in sorted(listdir(results_location)):
		i+=1
		print i
		
		results = process_single_map(join(results_location,f))
		
		#Fold results into the global aggregate.
		aggregate = merge_aggregates(aggregate,results)
	
	#Save aggregate.
	aggregate_target = join(aggregates_location,exp_name)
	if not os.path.exists(aggregate_target):
		os.mkdir(aggregate_target)
	f = open(join(aggregate_target,'aggregate.json'),'w')
	f.write(json.dumps(aggregate))
	f.close()
if __name__ == '__main__':
	
	import sys
	
	cmd = sys.argv[1]

	if cmd == 'create_aggregate':#Usage: python StructuralResults.py create_aggregate [experiment name]
		#Create aggregate given experiment name
		process_experiment(sys.argv[2])
	
	if cmd == 'debug_single':#Debugging processing single maps
		exp_name = 'sample'
		results_location = join(structural_map_location,exp_name)
		for f in sorted(listdir(results_location)):
			print process_single_map(join(results_location,f))
	if cmd == 'debug_multiple':#Debugging processing multiple maps
		exp_name = 'sample'
		results_location = join(structural_map_location,exp_name)
		i = 0
		aggregate = {'cdr':{},'frame':{},'full':{},'fread':{}}

		for f in sorted(listdir(results_location)):
			i+=1
			print i
		
			results = process_single_map(join(results_location,f))
		
			#Fold results into the global aggregate.
			aggregate = merge_aggregates(aggregate,results)
	
		#Save aggregate.
		aggregate_target = join(aggregates_location,exp_name)
		if not os.path.exists(aggregate_target):
			os.mkdir(aggregate_target)
		f = open(join(aggregate_target,'aggregate.json'),'w')
		f.write(json.dumps(aggregate))
		f.close()
	

	pass
