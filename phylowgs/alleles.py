from scipy.stats import beta, binom
import scipy.stats as stat
from scipy.misc import comb
from util         import *
from numpy        import *
from node         import *

from util2 import *


class alleles(Node):
	
	init_mean = 0.5
	min_conc = 0.01
	max_conc = 0.1

	def __init__(self, parent=None, tssb=None, conc=0.1, ntps=5):
		super(alleles, self).__init__(parent=parent, tssb=tssb)
		
		if tssb!=None: ntps= len(tssb.data[0].a)

		# pi is a first-class citizen
		self.pi=zeros(ntps); self.params=zeros(ntps) 
		self.params1=zeros(ntps);self.pi1=zeros(ntps) # used in MH	to store old state
		
		self.path = None # set of nodes from root to this node
		self.ht = 0
		
		if parent is None:
			self._conc = conc			
			self.pi = 1.0*ones(ntps)
			self.params = 1.0*ones(ntps)
			
		else:
			self.pi = rand(1)*parent.pi
			parent.pi = parent.pi - self.pi
			self.params = self.pi
			
	def conc(self):
		if self.parent() is None:
			return self._conc
		else:
			return self.parent().conc()			        

	def kill(self):
		if self._parent is not None:
			self._parent._children.remove(self)
		self._parent.pi = self._parent.pi + self.pi
		self._parent   = None
		self._children = None

	def logprob(self, x):
		return x[0]._log_likelihood(self.params)
		
	def complete_logprob(self):
		return sum([self.logprob([data]) for data in self.get_data()])
