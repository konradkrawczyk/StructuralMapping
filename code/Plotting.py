
import matplotlib.pyplot as plt


#Generic code to plot histograms.
def plot_histogram(data,output,xlabel='x',ylabel='y',title=''):
	fig, ax = plt.subplots(figsize=(9, 6))
	ax.hist(data,100,color='g',alpha=0.6)
	ax.set_xlabel(xlabel)
	ax.set_ylabel(ylabel)
	ax.set_title(title)
	if output!=None:
		print "Saving to"
		plt.savefig(output)
	else:
		plt.show()
	plt.close()

if __name__ == '__main__':
	pass
