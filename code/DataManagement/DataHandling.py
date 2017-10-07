from Common.Common import numbered_datasets_location,random_filename
import os,pickle
from os.path import join

#####Dealing with processed data#############

#Receives a dictionary of sequences to be saved mapped to numbering.
def save_data(experiment_name,data):
	#See if experiment_name folder exists.
	target_path = join(numbered_datasets_location,experiment_name)
	if not os.path.exists(target_path):
		os.mkdir(target_path)
	
	#Random name for the chunk
	filename = random_filename()
	while os.path.exists(join(target_path,filename)):
		filename = random_filename()
	
	#Save the data as pickle.
	print "[DataHandling.py] Dumping chunk",filename,'for experiment',experiment_name
	pickle.dump(data,open(join(target_path,filename),'w'))

if __name__ == '__main__':
	import sys
	command = sys.argv[1]

	if command == 'saving':
		experiment_name  =  'sample data'
		save_data(experiment_name,{'random':'data'})
		 
