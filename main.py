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


BASE_SUBTRACT = 20000
BASE_DIVIDE = 40000
SAMPLE_SIZE = 1000

#paths now relative
chrNum = 16
chrNumString = "chr%d" % chrNum

path1 = "../data/Matrix%dA.int.bed.gz" % chrNum
f1 = gzip.open(path1, 'r')
path2 = "../data/Matrix%dB.int.bed.gz" % chrNum
f2 = gzip.open(path2, 'r')

def load():
	print "Loading..."

	hash12 = {}
	sums1 = {}
	sums2 = {}

	startt = time.time()

	for line in f1:
		if line.startswith(chrNumString):
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
		if line.startswith(chrNumString):
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

	endt = time.time()
	print ("Loading took %f" % (endt - startt))

	pickle.dump(hash12, open("hash12",'wb'))
	pickle.dump(sums1, open("sums1",'wb'))
	pickle.dump(sums2, open("sums2",'wb'))

	f1.close()
	f2.close()

	print "Pickled!"
	
def compute():

	hash12 = pickle.load(open("hash12",'rb'))
	sums1 = pickle.load(open("sums1",'rb'))
	sums2 = pickle.load(open("sums2",'rb'))

	p = []
	sumtime = 0

	print "Computing..."
	for x in hash12:
		# need to check whether x is in sums1 and also sums2 --ferhat
		if x not in sums1 or x not in sums2:
			continue
		startt = time.time()
		for y in hash12[x]:
			chi2, pval, dof, exp = chi([[sums1[x],sums2[x]],[hash12[x][y][0],hash12[x][y][1]]])
			p.append(pval)
		endt = time.time()
		#print ("Index\t%s\ttime\t%f\tnumberOfTests\t%d" % (x,endt - startt,len(hash12[x])))
		sumtime += (endt-startt)
		# uncomment the below line if you want to test things on a few rows only -- ferhat
		#if x>100: break
	print ("Total time\t%s" % (sumtime))

	plothist(p, len(p), 'chiHistogram.png')

def plothist(pvals, count, name):
	n, bins = np.histogram(pvals, bins=100)
	plt.hist(pvals, bins=100, weights=np.zeros_like(pvals) + 100. / len(pvals))
	plt.hist(pvals, bins=100, weights=np.zeros_like(pvals) + 100. / len(pvals), histtype='step', cumulative=True)

	plt.title("P-Values for Locus Interactions in " + chrNumString)
	plt.xlabel("Bins")
	plt.ylabel("Frequency")
	plt.show()
	plt.savefig(name)

	#Compute score
	total = 0
	score = 0
	for x in n:
		total += x
		score += (.01*(1.0*total/len(pvals)))
	print "Reproducibiity score is %f" % (1.0-score)

def main():
	#inp = raw_input("Reload data? [y/n]: ")
	#if inp is "y":
	for file in os.listdir("../data/"):
		print file
		#try:
		#	load()
		#except Exception as e: #for my OSX notifications
		#	notifier.promptError("Script failed: loading error",str(e))
		#else:
		#	notifier.promptSuccess("Script success!", "No errors")
	try:
		compute()
	except Exception as e:
		notifier.promptError("Script failed: computation error",str(e))
	else:
		notifier.promptSuccess("Script success!", "No errors")


if __name__ == "__main__":
    main()
