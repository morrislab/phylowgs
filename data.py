from numpy import *
import scipy.stats as stat
from scipy.special import gammaln
import util2 as u

class Datum(object):
	def __init__(self, name, id, a, d, mu_r=0, mu_v=0):
		self.name = name # SSM name, blank for CNV
		self.id = id
		self.a = a
		self.d = d
		self.mu_r = mu_r # 1-p_error
		self.mu_v = mu_v
		self._log_bin_norm_const = [u.log_bin_coeff(self.d[tp], self.a[tp]) for tp in arange(len(self.a))]

		## all variables below are used by cnv related computations
		self.nr = 0
		self.nv = 0
		self.node = None # this is the node where the datum resides
		self.cnv = [] # for SSM, this is [(cnv,cp,cm)]
		
		self.tssb = None # this is just a pointer to tssb (tree object), gets initialized in evolve.py
	
	
	# for multiple samples
	def _log_likelihood(self, phi,update_tree=True,new_state=0):
		ntps = len(phi) # multi sample
		return sum([self.__log_likelihood__(phi[tp],tp,update_tree,new_state) for tp in arange(ntps)])	
	
	# new_state is set to 0 or 1 during Metropolis-Hastings updates, defaults to 0 in all other places	
	def __log_likelihood__(self, phi, tp, update_tree=True,new_state=0):	
		if update_tree:
					##################################################
					## some useful info about the tree,
					## used by CNV related computations,
					u.set_node_height(self.tssb)
					u.set_path_from_root_to_node(self.tssb)
					u.map_datum_to_node(self.tssb)
					##################################################
		return self.__log_complete_likelihood__(phi, self.mu_r, self.mu_v, tp, new_state)
	
	
	# for multiple samples
	def _log_complete_likelihood(self, phi, mu_r, mu_v):
		ntps = len(self.a)
		return sum([self.__log_complete_likelihood__(phi, mu_r, mu_v,tp) for tp in arange(ntps)])
	def __log_complete_likelihood__(self, phi, mu_r, mu_v, tp, new_state=0):	
		
		if self.cnv:
			ll = []
			poss_n_genomes = self.compute_n_genomes(tp,new_state)
			poss_n_genomes = [x for x in poss_n_genomes if x[1] > 0]
			for (nr,nv) in poss_n_genomes:
				mu = (nr * mu_r + nv*(1-mu_r) ) / (nr+ nv)
				ll.append(u.log_binomial_likelihood(self.a[tp], self.d[tp], mu) + log(1.0/len(poss_n_genomes)) +  self._log_bin_norm_const[tp])
			if len(poss_n_genomes) == 0:
				ll.append(log(1e-99)) # to handle cases with zero likelihood
			llh = u.logsumexp(ll)
		else: ## CNV datum
			mu = (1 - phi) * mu_r + phi*mu_v # (mu_r=0.999, mu_v=0.5)
			llh = u.log_binomial_likelihood(self.a[tp], self.d[tp], mu) +  self._log_bin_norm_const[tp]
		return 	llh
	
	# computes the binomial parameter
	def compute_n_genomes(self,tp,new_state=0):
	
		def descend(nd,new_state):

			pi = nd.pi1[tp] if new_state else nd.pi[tp] # this is needed for Metropolis-Hastings likelihood computations
			ssm_node = self.node.path[-1]
			mr_cnv = self.find_most_recent_cnv(nd)
			ancestors = nd.get_ancestors()
			if (not ssm_node in ancestors) and (not mr_cnv):
				self.nr1 = self.nr1 + pi * 2
				self.nr2 = self.nr2 + pi * 2
			elif ssm_node in ancestors and (not mr_cnv):
				self.nr1 = self.nr1 + pi
				self.nv1 = self.nv1 + pi
				self.nr2 = self.nr2 + pi
				self.nv2 = self.nv2 + pi
			elif (not ssm_node in ancestors) and mr_cnv:
				self.nr1 = self.nr1 + pi * (mr_cnv[1] + mr_cnv[2])
				self.nr2 = self.nr2 + pi * (mr_cnv[1] + mr_cnv[2])
			elif ssm_node in ancestors and mr_cnv:
				if ssm_node in mr_cnv[0].node.get_ancestors():
					self.nr1 = self.nr1 + pi * mr_cnv[1]
					self.nv1 = self.nv1 + pi * mr_cnv[2]
					self.nr2 = self.nr2 + pi * mr_cnv[2]
					self.nv2 = self.nv2 + pi * mr_cnv[1]
				else:
					self.nr1 = self.nr1 + pi * max(0,(mr_cnv[1]+mr_cnv[2] - 1))
					self.nv1 = self.nv1 + pi * min(1,mr_cnv[1]+mr_cnv[2])
					self.nr2 = self.nr2 + pi * max(0,(mr_cnv[1] + mr_cnv[2] - 1))
					self.nv2 = self.nv2 + pi * min(1,mr_cnv[1]+mr_cnv[2])
			else:
				print "PANIC"
		
		nodes = self.tssb.root['node'].tssb.get_nodes()
		self.nr1 = 0
		self.nv1 = 0
		self.nr2 = 0 
		self.nv2 = 0
		for nd in nodes:
			descend(nd, new_state)
		out = [(self.nr1,self.nv1),(self.nr2,self.nv2)]
 	
		return out
	
	def find_most_recent_cnv(self, nd):
		out = None
		for n in nd.get_ancestors()[::-1]:
			if n in [x[0].node for x in self.cnv]:
				out = [x for x in self.cnv if x[0].node == n][0]
				break
		return out
		
	
	########## old code, not in use, but don't delete #####################
	def __log_complete_likelihood1__(self, phi, mu_r, mu_v, new_state=0):	
		llh = []
		if self.cnv: 
			self.compute_n_genomes(0,new_state) # maternal
			mu = (self.nr * mu_r + self.nv*(1-mu_r) ) / (self.nr+ self.nv)
			llh.append(u.log_binomial_likelihood(self.a, self.d, mu) + log(0.5) +  self._log_bin_norm_const)
			self.compute_n_genomes(1,new_state) # paternal
			mu = (self.nr * mu_r + self.nv*(1-mu_r) ) / (self.nr+ self.nv)
			llh.append(u.log_binomial_likelihood(self.a, self.d, mu) + log(0.5) +  self._log_bin_norm_const)
			llh = u.logsumexp(ll)
		else: ## CNV datum or SSM with no CNV
			mu = (1 - phi) * mu_r + phi*mu_v # (mu_r=0.999, mu_v=0.5)
			llh = u.log_binomial_likelihood(self.a, self.d, mu) +  self._log_bin_norm_const
		return 	llh
	def compute_n_genomes1(self,maternal,new_state=0):
		####### TEMPORARY ONLY ###############
		maternal = True
		self.nr=self.nv=0
	
		wts,nodes = self.tssb.get_mixture()
		ancestors = self.node.get_ancestors() # path from root to ssm node
		
		mr_cnv = self.cnv[0] # CNV corresponding to this SSM
		
		# do this until we encounter the SSM node,
		# i.e., along the path from root to the SSM node
		visited_cnv = False
		for node in ancestors:
			pi = node.pi1 if new_state else node.pi # this is needed for Metropolis-Hastings likelihood computations
		
			if node != mr_cnv[0].node and visited_cnv==False: # until CNV is encountered
				self.nr = self.nr + pi*2
			else:
				visited_cnv = True
				self.nr = self.nr + pi*(mr_cnv[1]+mr_cnv[2])
		
		# do this after the SSM node, i.e, for all nodes in the subtree below the SSM node
		def descend(nd):
			pi = nd.pi1 if new_state else nd.pi # this is needed for Metropolis-Hastings likelihood computations
		
			if nd == mr_cnv[0].node:
				if maternal:
					self.nr = self.nr + pi*mr_cnv[1]
					self.nv = self.nv + pi*mr_cnv[2]
				else:
					self.nr = self.nr + pi*mr_cnv[2]
					self.nv = self.nv + pi*mr_cnv[1]
			else:
				self.nr = self.nr + pi * (mr_cnv[1]+mr_cnv[2] - 1)
				self.nv = self.nv + pi
			
			for child in nd.children():
				descend(child)
		
		# traverse the tree below the ssm node
		for child in node.children(): descend(child)
