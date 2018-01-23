###This file servers to plot the results of the aggregates.
from Common.Common import aggregates_location
import json,pprint,os
from os.path import join
import pickle
import numpy as np

import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib

#Where to look for aggregate files.
from Common.Common import aggregates_location

###############
#FullSequence##
###############

colors = ['green','blue','magenta']

def plot_rmsds(rmsds,sids,rmsderr,labels,title,ub,lb):

	N = len(labels)
	t = np.arange(N)
	bar_width = 0.6

	fig, ax1 = plt.subplots()
	
	ax2 = ax1.twinx()
	j = 0
	for exp_name in sids:
		
		for i in range(0,len(t)):
			print exp_name,colors[j]
			ax2.bar(t[i], sids[exp_name][1][i], bar_width,color='magenta',alpha=0.3)#,alpha=color_range)
			#print i,rmsds[i],color_range
		j+=1
	ax2.set_ylabel('Number of Representatives', color='black')
	for tl in ax2.get_yticklabels():
	    tl.set_color('k')
	plt.title(title)
	plt.xticks(t+bar_width/2 , labels)


	#ax1.plot(t+bar_width/2, rmsds,color='r')
	ax1.errorbar(t+bar_width/2, rmsds, yerr=rmsderr,marker='s',mec='g',elinewidth=0,mfc='k',ms=0)
	ax1.set_xlabel('(%) Sequence identity')
	# Make the y-axis label and tick labels match the line color.
	ax1.set_ylabel('Mean RMSD (with standard deviation)', color='b')
	for tl in ax1.get_yticklabels():
	    tl.set_color('b')
	ax1.yaxis.grid()

	#ax1.grid()
	#plt.show()
	plt.savefig("plots/"+title.replace(' ','')+".png")

def get_coverage(exp_name,region,start_range):
	
	#Load data.
	data = json.load(open(join(aggregates_location,exp_name,'aggregate.json')))
	
	#Prepare coverage.
	coverage = dict()
	for i in range(0,101):
		coverage[i] = 0
	above_80 = 0
	total = 0
	for i in data[region]:
		if int(i)>= 80:
			above_80+=data[region][i]
		total+=data[region][i]
		coverage[int(i)]+=data[region][i]

	coverage_final = []
	for i in range(start_range,101):
		if i in coverage:
			coverage_final.append((i,coverage[i]))
		else:
			coverage_final.append((i,0))
	plotdata =  zip(*sorted(coverage_final))

	print "Total",total
	print "Above 80",above_80,float(above_80)/float(total)

	return plotdata




#######
#CDRs##
#######

#For input order is #CDR1, #CDR2 and #CDR3
def bar_plot(can_model=(200, 350, 300),cant_model=(200, 300, 400),pdb=(20, 30, 40),exp_name="noname",redundant="non-redundant"):
	import numpy as np
	import matplotlib.pyplot as plt

	#For aesthetic purposes
	the_alpha = 0.4

	#We always have three CDRs
	N = 3
	
	
	ind = np.arange(N)  # the x locations for the groups
	width = 0.3       # the width of the bars
	
	fig, ax = plt.subplots(figsize=(9, 6))
	#can model
	rects1 = ax.bar(ind, can_model, width, color='g',alpha=the_alpha)
	#cant model
	rects2 = ax.bar(ind + width, cant_model, width, color='r',alpha=the_alpha)
	#PDB
	rects3 = ax.bar(ind + width+width, pdb, width, color='b',alpha=the_alpha)

	# add some text for labels, title and axes ticks
	ax.set_ylabel('#Loops')
	ax.set_xlabel('CDR')
	

		
	ax.set_title('How many '+redundant+' CDRs from '+exp_name+' can we model?')
	ax.set_xticks(ind + width / 2)
	if exp_name == 'UCB_L':#I will burn in hell.
		ax.set_xticklabels(('L1', 'L2', 'L3'))
	
	else:
		ax.set_xticklabels(('H1', 'H2', 'H3'))
	
	ax.legend((rects1[0], rects2[0],rects3[0]), ('Can model', 'Cant model','In PDB'))

	def autolabel(rects):
	    """
	    Attach a text label above each bar displaying its height
	    """
	    for rect in rects:
		height = rect.get_height()
		ax.text(rect.get_x() + rect.get_width()/2., 1.01*height,
		        '%d' % int(height),
		        ha='center', va='bottom')

	autolabel(rects1)
	autolabel(rects2)
	autolabel(rects3)
	#plt.show()
	plt.savefig('plots/cdrs_'+redundant+'_'+exp_name+'.png')

#Fetch the cdr data for a single experiment.
#Counting in a non-redundant fashion.
def fetch_cdr_data(exp_name,redundant=False):
	
	data =json.load(open(join(aggregates_location,exp_name,'aggregate.json')))
	#Statistics if we can find templates or not.
	plot_can_model = []
	plot_cant_model = []
	plot_in_pdb = []
	
	#Plotting ESS scores.
	#Format cdr(h1)->length(6)->ess(26)->number
	scores = {}
	for cdr in sorted(data['fread']):
		print 'Fetching data for ',cdr
		if cdr not in scores:
			scores[cdr] = {}
		
		#Found in PDB
		in_pdb = 0
		#Can we model?
		can_model = 0
		#Can't we model?
		cant_model = 0
		#In the PDB
		in_pdb = 0
		#How many sequences in total?
		total = 0
		for sequence in data['fread'][cdr]:


			if redundant == True:
				total+=data['fread'][cdr][sequence]['total']
			else:
				total+=1
			if data['fread'][cdr][sequence]['pdb'] == True:
				
				if redundant == True:
					in_pdb+=data['fread'][cdr][sequence]['total']
				else:
					in_pdb+=1
				continue
			#Check if we have a template and ess score for it is in line with length predictions.
			if data['fread'][cdr][sequence]['model'] >0:
				if len(sequence)<13 and data['fread'][cdr][sequence]['score']>=25:
					if redundant == True:#If we can model even one properly -- consider it as model as we know the shape of the loop.
						can_model+=data['fread'][cdr][sequence]['total']
					else:
						can_model+=1
				elif len(sequence)<17 and data['fread'][cdr][sequence]['score']>=40:
					if redundant == True:#If we can model even one properly -- consider it as model as we know the shape of the loop.
						can_model+=data['fread'][cdr][sequence]['total']
					else:
						can_model+=1
				elif len(sequence)>16 and data['fread'][cdr][sequence]['score']>=55:
					if redundant == True:#If we can model even one properly -- consider it as model as we know the shape of the loop.
						can_model+=data['fread'][cdr][sequence]['total']
					else:
						can_model+=1
				else:
					if redundant == True:#If we can model even one properly -- consider it as model as we know the shape of the loop.
						cant_model+=data['fread'][cdr][sequence]['total']
					else:
						cant_model+=1
			else:
				if redundant == True:#If we can model even one properly -- consider it as model as we know the shape of the loop.
					cant_model+=data['fread'][cdr][sequence]['total']
				else:
					cant_model+=1

			#Deal with sequence lengths.
			sequence_length = len(sequence)
			if sequence_length not in scores[cdr]:
				scores[cdr][sequence_length]= {}

			#print cdr,sequence,data['fread'][cdr][sequence]
		plot_can_model.append((cdr,can_model))
		plot_cant_model.append((cdr,cant_model))
		plot_in_pdb.append((cdr,in_pdb))
	
	can_model = zip(*plot_can_model)[1]
	cant_model = zip(*plot_cant_model)[1]
	in_pdb = zip(*plot_in_pdb)[1]
	return can_model,cant_model,in_pdb

#Plotting entire regions, full sequence, frame, three cdrs.
def plot_region(chain,region,experiment):
	
	#Formatting the title correctly.
	proper_region = 'full variable'
	if region == '_framework':
		proper_region = 'Framework only'
	
	
	start_range = 50

	#exp_name = 'KELLY_HEPB_BOOSTER'
	experiments = [experiment]
	plotdatas = {}
	for exp_name in experiments:
		if region =='_framework':
			plotdatas[exp_name] = get_coverage(exp_name,'frame',start_range)
		elif region=='':
			plotdatas[exp_name] = get_coverage(exp_name,'full',start_range)

	print plotdatas

	labels = []
	
	rmsd_means = []
	rmsd_err = []

	rmsds = pickle.load(open('rmsd_serial/rmsd'+region+'_serial'+chain+'.txt'))
	for perc in range(start_range,101):
		
		if perc in rmsds:
			summ = 0
			for rmsd in rmsds[perc]:
				summ+=rmsd
			rmsd_means.append((perc,summ/float(len(rmsds[perc]))))
			rmsd_err.append((perc,np.std(rmsds[perc])))
		else:
			rmsd_err.append((perc,0))
			rmsd_means.append((perc,0))
		
		if perc % 10 == 0:
			labels.append((perc,str(perc)))
		else:
			labels.append((perc,''))

	

	plot_rmsd = zip(*rmsd_means)# pickle.load(open(join('rmsd/plot'+region+'_serial'+chain+'.txt')))
	
	

	plot_rmsd_err = zip(*rmsd_err)# pickle.load(open(join('rmsd/plot'+region+'_serial'+chain+'.txt')))
	labels = zip(*labels)
	
	plot_rmsds(plot_rmsd[1],plotdatas,plot_rmsd_err[1],labels[1],'Coverage of experiment '+experiment+' on region '+proper_region,1.0,0.6)

	#pickle.dump(plot_rmsd,open(join('rmsd/plot'+region+'_serial'+chain+'.txt'),'w'))
	#pickle.dump(plot_rmsd_err,open(join('rmsd/error'+region+'_serial'+chain+'.txt'),'w'))
	#Display as histogram.

def plot_es_scores(exp_name,cdr):
	
	data = json.load(open(join(aggregates_location,exp_name,'aggregate.json')))

	lengths = {}

	for sequence in data['fread'][cdr]:

		l = len(sequence)
		if l not in lengths:
			lengths[l] = []
		
		
		lengths[l].append(int(data['fread'][cdr][sequence]['score']))
		
	sorted_lengths = []
	sorted_len_nums = []
	for l in lengths:
		if len(lengths[l])>1:
			sorted_lengths.append((l,lengths[l]))
			sorted_len_nums.append((l,len(lengths[l])))
	
	sorted_lengths = sorted(sorted_lengths)

	labels,data = zip(*sorted_lengths)
	labels,len_coverage = zip(*sorted_len_nums)
	data_to_plot = []
	for elem in data:
		data_to_plot.append(elem)

	# Create a figure instance
	fig = plt.figure(1, figsize=(9, 6))

	# Create an axes instance
	ax = fig.add_subplot(111)
	ax.set_title('ESS score profile for '+cdr+', experiment '+exp_name)
	ax2 = ax.twinx()
	ax.plot([0,len(labels)+1],[25,25],'--g',label='FREAD score cutoff')
	ax.set_ylabel('FREAD score')
	ax.set_xlabel(cdr+' length (Chothia)')	
	ax2.plot(labels,len_coverage,'bo-',alpha=0.7,label="Distribution of CDR lengths")
	ax2.set_ylabel('Number of Loops')
	# Create the boxplot
	bp = ax.boxplot(data_to_plot)
	ax.set_xticklabels(labels)
	plt.legend()
	#plt.show()
	#quit()
	# Save the figure
	fig.savefig(join('plots','ESS_'+exp_name+cdr+'.png'), bbox_inches='tight')
	plt.close('all')

if __name__ == '__main__':

	import sys
	cmd = sys.argv[1]
	exp_name = sys.argv[2]
	if cmd == 'redundancies':
		can_model_r,cant_model_r,in_pdb_r = fetch_cdr_data(exp_name,redundant=True)
		can_model,cant_model,in_pdb = fetch_cdr_data(exp_name,redundant=False)
		print "Redundant"
		print can_model_r,cant_model_r,in_pdb_r
		print "nonredudnant"
		print can_model,cant_model,in_pdb
		for i in range(0,3):
			print 'CDR',i
			redu=can_model_r[i]+cant_model_r[i]+in_pdb_r[i]
			uniq=can_model[i]+cant_model[i]+in_pdb[i]
			print redu,uniq
			print 'Redundancy', float(uniq)/float(redu)
			print "In PDB R",float(in_pdb_r[i])/float(redu)
			print "In PDB U",float(in_pdb[i])/float(uniq)
		
	if cmd == 'frame':
		plot_region(sys.argv[3],'_framework',exp_name)
	if cmd == 'full':
		print "Plotting full",exp_name
		plot_region(sys.argv[3],'',exp_name)
	if cmd == 'cdr':#USAGE: python PlotResults.py cdr [exp name]  redundant
		can_model,cant_model,in_pdb = fetch_cdr_data(exp_name,redundant=(sys.argv[3]=='redundant'))
		bar_plot(can_model,cant_model,in_pdb,exp_name=exp_name,redundant=sys.argv[3])
	if cmd == 'ess':#Plot ess scores against lengths.
		plot_es_scores(exp_name,'H1')
		plot_es_scores(exp_name,'H2')
		plot_es_scores(exp_name,'H3')
		
	
	#plot_region('H','','KELLY_MENINGO')
	#Plot ESS scores.
	

	#Plotting cdrs
	
	
