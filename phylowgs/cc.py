# implementation of the correlation clustering algorithm
# http://en.wikipedia.org/wiki/Correlation_clustering
# http://www.cs.cornell.edu/People/tj/publications/joachims_hopcroft_05a.pdf

from numpy import *
from numpy.random import *
from scipy import linalg

from cvxopt import matrix, solvers


# W - adjacency matrix of the graph, size mxm
# return the clustering matrix, size mxm
def cc_lp(W):
	m = W.shape[0]
	W.shape = m*m,1	

	#reflexivity constraints
	G1 = zeros((m,m*m))
	r = array([i*m+i for i in arange(m)])
	G1[arange(m),r] = 1
	h1 = ones(m)

	G11 = zeros((m,m*m))
	G11[arange(m),r] = -1
	h11 = -ones(m)

	ids=array([[i,j] for i in arange(m) for j in arange(m)])
	rc=ids[where(ids[:,0]<ids[:,1])]
	r=array([i[0]*m+i[1] for i in rc])

	# symmetry constraints
	nc = len(rc)
	G2 = zeros((nc,m*m))
	G2[arange(nc),r] = 1
	rt = array([i[1]*m+i[0] for i in rc])
	G2[arange(nc),rt] = -1
	h2 = zeros(nc)

	G21 = zeros((nc,m*m))
	G21[arange(nc),r] = -1
	G21[arange(nc),rt] = 1
	h21 = zeros(nc)

	#range [0,1] constraints
	nc = len(r)
	G4 = zeros((nc,m*m))
	G4[arange(nc),r] = -1
	h4 = zeros(nc)
	G41 = zeros((nc,m*m))
	G41[arange(nc),r] = 1
	h41 = ones(nc)

	# transitivity constraints
	ids = array([[i,j,k] for i in arange(m) for j in arange(m) for k in arange(m)])
	ids = ids[where(ids[:,0] != ids[:,1])]
	ids = ids[where(ids[:,1] != ids[:,2])] # there has to be a better way to do this...
	rc = ids[where(ids[:,2] != ids[:,0])]
	nc = len(rc)

	r1=array([i[0]*m+i[1] for i in rc])
	r2=array([i[1]*m+i[2] for i in rc])
	r3=array([i[2]*m+i[0] for i in rc])

	G3 = zeros((nc,m*m))
	G3[arange(nc),r1] = 1
	G3[arange(nc),r2] = 1
	G3[arange(nc),r2] = -1
	h3 = ones(nc)

	G = vstack((G1,G11,G2,G21,G3,G4,G41))
	h = hstack((h1,h11,h2,h21,h3,h4,h41))
	sol = solvers.lp(matrix(-W), matrix(G), matrix(h))
	sol = array(sol['x'])
	sol.shape = m,m
	return sol
	
def knng(X):
	m,n = X.shape
	k = int(ceil(log2(m)))
	XX = dot(X,X.T)
    	D = outer(diag(XX),ones(m))+outer(ones(m),diag(XX))-2*XX
	D = D*1./D.max()
	W = zeros((m,m))
	for i in arange(m):
		ids = argsort(D[i])[1:k+1]
		W[i,ids]=exp(-D[i,ids])
	W[arange(m),arange(m)]=1
	return W