import os
import sys
import time
import cPickle

from subprocess import call
from numpy		import *
from numpy.random import *
from tssb		 import *
from alleles	 import *
from util		 import *

import numpy.random

from util2 import *
from params import *
from printo import *

# MCMC settings
burnin		= 1000
num_samples   = 5000
checkpoint	= 1000
dp_alpha	  =25.0
dp_gamma	  = 1.0
init_conc	= 1.0
alpha_decay   = 0.25
top_k = 5

#metropolis-hastings settings
mh_itr = 5000 # no. of iterations in metropolis-hastings
mh_burnin = 0

# number of multiple samples
NTPS = 5


# fin1 and fin2: input data files for ssm and cnv
# fout: output folder for tssb trees/samples
# out2: output file to save top-k trees
# out3: output file to save clonal frequencies
# out4: output file to save log likelihood trace
# num_samples: number of MCMC samples
# mh_itr: number of metropolis-hasting iterations
# dp_alpha: dp alpha
# rand_seed: random seed (initialization)
def run(fin1,fin2,fout='./best',out2='top_k_trees',out3='clonal_frequencies',out4='llh_trace',num_samples=2500,mh_itr=5000,mh_std=100,rand_seed=1):
	if not os.path.exists(fout):
		os.makedirs(fout)
	else:
		call(['rm','-r',fout])
		os.makedirs(fout)

	seed(rand_seed)
	codes, n_ssms, n_cnvs = load_data(fin1,fin2)
	NTPS = len(codes[0].a) # number of samples / time point
	glist = [datum.name for datum in codes if len(datum.name)>0]

	root  = alleles(conc=0.1,ntps=NTPS)
	tssb  = TSSB( dp_alpha=dp_alpha, dp_gamma=dp_gamma, alpha_decay=alpha_decay, root_node=root, data=codes )
	# hack...
	if 1:
		depth=0
		tssb.root['sticks'] = vstack([ tssb.root['sticks'], boundbeta(1, tssb.dp_gamma) if depth!=0 else .999])
		tssb.root['children'].append({ 'node': tssb.root['node'].spawn(),
					'main':boundbeta(1.0, (tssb.alpha_decay**(depth+1))*tssb.dp_alpha) if tssb.min_depth <= (depth+1) else 0.0, 
					'sticks' : empty((0,1)),	
					'children' : [] })
		new_node=tssb.root['children'][0]['node']	
		for n in range(tssb.num_data):	
			tssb.assignments[n].remove_datum(n)
			new_node.add_datum(n)
			tssb.assignments[n] = new_node

		
	####

	
	for datum in codes: datum.tssb=tssb 

	dp_alpha_traces	= zeros((num_samples, 1))
	dp_gamma_traces	= zeros((num_samples, 1))
	alpha_decay_traces = zeros((num_samples, 1))
	conc_traces	   = zeros((num_samples, 1))
	cd_llh_traces	  = zeros((num_samples, 1))

	intervals = zeros((7))
	best_tssb = 0

	# clonal frequencies
	freq = dict([(g,[] )for g in glist])	
	
	print "Starting MCMC run..."
	
	for iter in range(-burnin,num_samples):
		if iter<0: print iter
	
		times = [ time.time() ]
			
		tssb.resample_assignments()
		times.append(time.time())

		tssb.cull_tree()
		times.append(time.time())
		
		# assign node ids
		wts,nodes=tssb.get_mixture()
		for i,node in enumerate(nodes): node.id=i
		
		##################################################
		## some useful info about the tree,
		## used by CNV related computations,
		## to be called only after resampling assignments
		set_node_height(tssb)
		set_path_from_root_to_node(tssb)
		map_datum_to_node(tssb)
		##################################################

		mh_acc = metropolis(tssb,mh_itr,mh_std,mh_burnin,n_ssms,n_cnvs,fin1,fin2,rand_seed,NTPS)
		if float(mh_acc) < 0.08 and mh_std < 10000:
			mh_std = mh_std*2.0
			print "Shrinking MH proposals. Now %f" % mh_std
		if float(mh_acc) > 0.5 and float(mh_acc) < 0.99:
			mh_std = mh_std/2.0
			print "Growing MH proposals. Now %f" % mh_std
		times.append(time.time())
	
		#root.resample_hypers()
		times.append(time.time())
	
		tssb.resample_sticks()
		times.append(time.time())
		
		tssb.resample_stick_orders()
		times.append(time.time())
	
		tssb.resample_hypers(dp_alpha=True, alpha_decay=True, dp_gamma=True)
		times.append(time.time())
 
		intervals = intervals + diff(array(times)) 

		if iter>=0:
			dp_alpha_traces[iter]	= tssb.dp_alpha
			dp_gamma_traces[iter]	= tssb.dp_gamma
			alpha_decay_traces[iter] = tssb.alpha_decay
			conc_traces[iter]	   = root.conc()
			cd_llh_traces[iter]	  = tssb.complete_data_log_likelihood()
	   
		if iter>=0:
			if True or mod(iter, 10) == 0:
				(weights, nodes) = tssb.get_mixture()
				print iter, len(nodes), cd_llh_traces[iter], mh_acc, tssb.dp_alpha, tssb.dp_gamma, tssb.alpha_decay#, " ".join(map(lambda x: "%0.2f" % x, intervals.tolist())) 
				intervals = zeros((7))	  
  
		if iter >= 0 and argmax(cd_llh_traces[:iter+1]) == iter:
			print "\t%f is best per-data complete data likelihood so far." % (cd_llh_traces[iter])

		# save all trees
		if iter >= 0:
			fh = open(fout+'/'+str((cd_llh_traces[iter])[0]), 'w')
			cPickle.dump(tssb, fh)
			fh.close()
		
		wts, nodes = tssb.get_mixture()
		#Save log likelihood:savetxt('loglike',cd_llh_traces)
		savetxt(out4,cd_llh_traces)		

		#log clonal frequencies		
		'''
		if iter >= 0:		
			wts, nodes = tssb.get_mixture()
			for node in nodes:
				data = node.get_data()
				for datum in data:
					if datum.name in freq:
						freq[datum.name].append(float(round(node.params,5)))
		'''
	#save the best tree
	print_top_trees(fout,out2,top_k)

	#save clonal frequencies
	glist = array(freq.keys(),str);glist.shape=(1,len(glist)) 
	savetxt(out3,vstack((glist, array([freq[g] for g in freq.keys()]).T)),fmt='%s',delimiter=', ')


def test():
	tssb=cPickle.load(open('ptree'))
	wts,nodes=tssb.get_mixture()	
	
	for dat in tssb.data:
		print [dat.id, dat.__log_likelihood__(0.5)]

if __name__ == "__main__":
	#run(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4],sys.argv[5],int(sys.argv[6]),int(sys.argv[7]),int(sys.argv[8]))
	run('ssm_data.txt','cnv_data.txt')
	#test()
