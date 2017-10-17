##This file deals with compiling the results of structurally mapped studies.
import json,os,pprint
from os import listdir
from os.path import join
from Common.Common import structural_map_location,aggregates_location

#From a single map file, create an aggreagate of sequence identitities.
def process_single_map(loc):
	
	results = {'cdr':{},'frame':{},'full':{}}

	data = json.load(open(loc))

	for sequence in data:
		for region in data[sequence]:
			sid = data[sequence][region]['best_sid']
			if sid not in results[region]:
				results[region][sid] = 0
			results[region][sid]+=1
	return results

#Given a dict from regions to sid list, fold the results of the new such list into the current one
def merge_aggregates(current,new):
	for region in new:
		for sid in new[region]:
			if sid not in current[region]:
				current[region][sid] = 0
			current[region][sid]+=new[region][sid]
	return current
	

if __name__ == '__main__':


	exp_name = 'KELLY_HEPB_BOOSTER'
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
	

	pass
