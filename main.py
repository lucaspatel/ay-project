#!
'''
@author: lucas
'''

### import
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import fisher_exact as fish
from scipy.stats import chi2_contingency as chi

import os, gzip, re
import pickle, random, time
import notifier

def natural_key(string_):
	return [int(s) if s.isdigit() else s for s in re.split(r'(\d+)', string_)]

BASE_SUBTRACT = 20000
BASE_DIVIDE = 40000
SAMPLE_SIZE = 1000

output = open("../results/output.txt",'w')
pvals = []

def load(f1, f2, chrNum):

	hash12 = {}
	sums1 = {}
	sums2 = {}

	startt = time.time()
	'''
	for line in f1:
		if line.startswith("chr%s" % chrNum):
			nums = re.findall('\d+', line)
			x = (int(nums[1])-BASE_SUBTRACT)/BASE_DIVIDE
			y = (int(nums[3])-BASE_SUBTRACT)/BASE_DIVIDE
			z = int(nums[4])

			if x not in hash12:
				hash12[x] = {}
			if y not in hash12[x]:
				hash12[x][y] = [0,0]
			hash12[x][y][0] = z

			if x not in sums1:
				sums1[x] = 0
			sums1[x] += int(nums[4])

			if y not in sums1:
				sums1[y] = 0
			sums1[y] += int(nums[4])	

	for line in f2:
		if line.startswith("chr%s" % chrNum):
			nums = re.findall('\d+', line)
			x = (int(nums[1])-BASE_SUBTRACT)/BASE_DIVIDE
			y = (int(nums[3])-BASE_SUBTRACT)/BASE_DIVIDE
			z = int(nums[4])

			if x not in hash12:
				hash12[x] = {}
			if y not in hash12[x]:
				hash12[x][y] = [0,0]
			hash12[x][y][1] = z

			if x not in sums2:
				sums2[x] = 0
			sums2[x] += int(nums[4])

			if y not in sums2:
				sums2[y] = 0
			sums2[y] += int(nums[4])
	'''
	endt = time.time()
	print ("Loading chromosome %s took %f" % (chrNum,endt - startt))

	pickle.dump(hash12, open("hash12",'wb'))
	pickle.dump(sums1, open("sums1",'wb'))
	pickle.dump(sums2, open("sums2",'wb'))
	
def compute(matrixname, chrNum):

	hash12 = pickle.load(open("hash12",'rb'))
	sums1 = pickle.load(open("sums1",'rb'))
	sums2 = pickle.load(open("sums2",'rb'))

	p = []
	sumtime = 0

	#takes longer each iteration, ask Dr. Ay
	for x in hash12:
		if x not in sums1 or x not in sums2:
			continue
		startt = time.time()
		for y in hash12[x]:
			chi2, pval, dof, exp = chi([[sums1[x],sums2[x]],[hash12[x][y][0],hash12[x][y][1]]])
			p.append(pval)
		endt = time.time()
		#print ("Index\t%s\ttime\t%f\tnumberOfTests\t%d" % (x,endt - startt,len(hash12[x])))
		sumtime += (endt-startt)

	print ("Computing chromsome %s took %s" % (chrNum, sumtime))

	return plothist(p, len(p), matrixname, chrNum)

def plothist(pvals, count, matrixname, chrNum):
	n, bins = np.histogram(pvals, bins=100)
	pvals = [1]
	plt.hist(pvals, bins=100, weights=np.zeros_like(pvals) + 100. / len(pvals))
	plt.hist(pvals, bins=100, weights=np.zeros_like(pvals) + 100. / len(pvals), histtype='step', cumulative=True)

	plt.title("P-Values for Locus Interactions in %s" % matrixname)
	plt.xlabel("Bins")
	plt.ylabel("Frequency")
	plt.savefig("../results/plots/%schr%s_plot.png" % (matrixname, chrNum))

	#Compute score
	total = 0
	score = 0
	for x in n:
		total += x
		score += (.01*(1.0*total/len(pvals)))
	print "Reproducibiity score for chromosome %s is %f" % (chrNum, 1.0-score)
	output.write("Score for chromosome %s is %f \n" % (chrNum, 1.0-score))
	return (1.0-score)

def main():
	files = sorted(os.listdir("../data/"), key=natural_key)
	for f in files:
		if f.startswith('.'):
			files.remove(f) #remove pesky .DS_Store, implement more elegant solution later
	fileiter = iter(files)
	
	for x,y in zip(fileiter, fileiter):
		scores = []
		f1 = gzip.open("../data/"+x, 'r')
		f2 = gzip.open("../data/"+y, 'r')
		output.write("%s: \n" % x[:-12])

		for i in range(1,23):
			try:
				load(f1,f2,i)
				scores.append(compute(x[:-12], i))
			except Exception as e: #for my OSX notifications
				notifier.promptError("Script failed",str(e))
				quit()

		try:
			load(f1,f2,"X")
			scores.append(compute(x[:-12], "X"))
			load(f1,f2,"Y")
			scores.append(compute(x[:-12], "Y"))
		except Exception as e: #for my OSX notifications
			notifier.promptError("Script failed",str(e))
			quit()

		output.write("Average reproducibility score is %f \n" % (sum(scores)/len(scores)))

		f1.close()
		f2.close()
	
	notifier.promptSuccess("Script success!", "No errors")

if __name__ == "__main__":
    main()
