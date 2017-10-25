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
		print exp_name
		for i in range(0,len(t)):
		
			print exp_name,colors[j]
			ax2.bar(t[i], sids[exp_name][1][i], bar_width,color=colors[j],alpha=0.3)#,alpha=color_range)
			#print i,rmsds[i],color_range
		j+=1
	ax2.set_ylabel('Number of Representatives', color='k')
	for tl in ax2.get_yticklabels():
	    tl.set_color('k')
	plt.title(title)
	plt.xticks(t+bar_width/2 , labels)


	ax1.plot(t+bar_width/2, rmsds,color='r')
	ax1.set_xlabel('(%) Sequence identity')
	# Make the y-axis label and tick labels match the line color.
	ax1.set_ylabel('Mean RMSD (with standard deviation)', color='b')
	for tl in ax1.get_yticklabels():
	    tl.set_color('b')

	print "Showing"
	ax1.grid()
	plt.show()
	#plt.savefig("/homes/krawczyk/Downloads/Figures/"+title+".png")

def get_coverage(exp_name,region):
	
	#Load data.
	data = json.load(open(join(aggregates_location,exp_name,'aggregate.json')))
	
	#Prepare coverage.
	coverage = dict()
	for i in range(0,101):
		coverage[i] = 0
	
	for i in data[region]:
		coverage[int(i)]+=data[region][i]
	coverage_final = []
	for i in range(start_range,101):
		if i in coverage:
			coverage_final.append((i,coverage[i]))
		else:
			coverage_final.append((i,0))
	plotdata =  zip(*sorted(coverage_final))
	return plotdata

#Plotting entire regions, full sequence, frame, three cdrs.
def plot_region():
	chain = 'H'
	region = '_framework'

	start_range = 50

	#exp_name = 'KELLY_HEPB_BOOSTER'
	experiments = ['KELLY_HEPB_BOOSTER']
	plotdatas = {}
	for exp_name in experiments:
		plotdatas[exp_name] = get_coverage(exp_name,'frame')

	labels = []
	
	#rmsd_means = []
	#rmsd_err = []
	
	for perc in range(start_range,101):
		'''
		if perc in rmsds:
			summ = 0
			for rmsd in rmsds[perc]:
				summ+=rmsd
			rmsd_means.append((perc,summ/float(len(rmsds[perc]))))
			rmsd_err.append((perc,np.std(rmsds[perc])))
		else:
			rmsd_err.append((perc,0))
			rmsd_means.append((perc,0))
		'''
		if perc % 10 == 0:
			labels.append((perc,str(perc)))
		else:
			labels.append((perc,''))
	plot_rmsd = pickle.load(open(join('rmsd/plot'+region+'_serial'+chain+'.txt')))
	
	plot_rmsd_err = pickle.load(open(join('rmsd/plot'+region+'_serial'+chain+'.txt')))
	labels = zip(*labels)
	
	plot_rmsds(plot_rmsd[1],plotdatas,plot_rmsd_err[1],labels[1],'bla',1.0,0.6)

	

	#pickle.dump(plot_rmsd,open(join('rmsd/plot'+region+'_serial'+chain+'.txt'),'w'))
	#pickle.dump(plot_rmsd_err,open(join('rmsd/error'+region+'_serial'+chain+'.txt'),'w'))
	

	#Display as histogram.


#######
#CDRs##
#######

def bar_plot():
	import numpy as np
	import matplotlib.pyplot as plt

	#For aesthetic purposes
	the_alpha = 0.4

	#We always have three CDRs
	N = 3
	#CDR1, CDR2 and CDR3
	can_model = (200, 350, 300)
	cant_model = (200, 300, 400)
	pdb = (20, 30, 40)
	
	ind = np.arange(N)  # the x locations for the groups
	width = 0.3       # the width of the bars

	fig, ax = plt.subplots()
	#can model
	rects1 = ax.bar(ind, can_model, width, color='g',alpha=the_alpha)
	#cant model
	rects2 = ax.bar(ind + width, cant_model, width, color='r',alpha=the_alpha)
	#PDB
	rects3 = ax.bar(ind + width+width, pdb, width, color='b',alpha=the_alpha)


	# add some text for labels, title and axes ticks
	ax.set_ylabel('#Loops')
	ax.set_xlabel('CDR')
	ax.set_title('How many loops can we model?')
	ax.set_xticks(ind + width / 2)
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

	plt.show()

if __name__ == '__main__':

	bar_plot()
	
	
