
#include<vector>
#include<cstring>
#include<math.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include "util.hpp"

void sample_cons_params(struct node nodes[],struct config conf,gsl_rng *rand,int tp);
double multi_param_post(struct node nodes[], struct datum data[], int old,struct config conf);
double param_post(struct node nodes[], struct datum data[], int old,struct config conf, int tp);
void update_params(struct node nodes[], struct config conf);
void get_pi(struct node nodes[], double pi[], struct config conf, int old, int tp);

void load_ssm_data(char fname[], struct datum data[], struct config conf);
void load_cnv_data(char fname[], struct datum data[], struct config conf);
void load_data_states(char fname[], struct datum data[], struct node nodes[], struct config conf);

void load_tree(char fname[], struct node nodes[], struct config conf);
void write_params(char fname[], struct node nodes[], struct config conf);

void mh_loop(struct node nodes[], struct datum data[], char* fname, struct config conf);

struct config{
	int MH_ITR;
	float MH_STD;
	
	int N_SSM_DATA; // no. of data points
	int N_CNV_DATA; // no. of data points
	
	int NNODES; // no. of nodes in the tree
	int TREE_HEIGHT; 
	int NTPS; // no. of samples / time points	
};

struct state{
	struct node* nd;
	int nr,nv;
};

struct node{
	int id;
	vector<double> param,pi;
	vector<double> param1,pi1; // dummy	
	int ndata;
	vector<int> dids;	
	int nchild;
	vector<int> cids; // children ids	
	int ht;	
};


struct datum{	
	int id;	
	vector<int> a,d;
	
	double mu_r,mu_v;
	
	vector<double> log_bin_norm_const;//log_bin_coeff(d,a);	

	struct datum* cnv; // for SSM datum, this is a pointer to its CNV datum
	//int cnv;// just an indicator for cnv or ssm datum
	
	// this is used to compute the binomial parameter
	vector <struct state> states1, states2; // maternal and paternal state
	
	double log_ll1111(vector<double> phi, int old){
		double llh = 0.0;
		for(int tp=0; tp<phi.size();tp++)
			llh+=log_complete_ll(phi[tp],mu_r,mu_v,old,tp);
		return llh;
	}

	double log_ll(double phi, int old, int tp){
		return log_complete_ll(phi,mu_r,mu_v,old,tp);
	}
	
	double log_complete_ll(double phi, double mu_r, double mu_v, int old, int tp){
		double nr=0;
		double nv=0;
		double mu = 0;
		double llh = 0;
				
		if(cnv==NULL){	// cnv data
			mu = (1 - phi) * mu_r + phi * mu_v;
			llh = log_binomial_likelihood(a[tp], d[tp], mu) + log_bin_norm_const[tp];			
		}	
		else{ // ssm data
			double ll[2]; // maternal and paternal	
			nr=0;
			nv=0;
			for(int i=0;i<states1.size();i++){
				if(old==0){
					nr += (states1[i].nd->pi1[tp])*states1[i].nr;
					nv += (states1[i].nd->pi1[tp])*states1[i].nv;
				}else{
					nr += (states1[i].nd->pi[tp])*states1[i].nr;
					nv += (states1[i].nd->pi[tp])*states1[i].nv;
				}	
			}
			if(nr+nv>0){
				mu = (nv*(1-mu_r) + nr*mu_r)/(nr+nv);
				ll[0] = log_binomial_likelihood(a[tp], d[tp], mu) + log(0.5) + log_bin_norm_const[tp];
			}else{
				ll[0]=log(pow(10,-99));
			}
			// repetitive...
			nr=nv=0;
			for(int i=0;i<states2.size();i++){
				if(old==0){
					nr += (states2[i].nd->pi1[tp])*states2[i].nr;
					nv += (states2[i].nd->pi1[tp])*states2[i].nv;
				}else{
					nr += (states2[i].nd->pi[tp])*states2[i].nr;
					nv += (states2[i].nd->pi[tp])*states2[i].nv;
				}				
			}
			if(nr+nv>0){
				mu = (nv*(1-mu_r) + nr*mu_r)/(nr+nv);
				ll[1] = log_binomial_likelihood(a[tp], d[tp], mu) + log(0.5) + log_bin_norm_const[tp];
			}else{
				ll[1]=log(pow(10,-99));
			}
			llh = logsumexp(ll,2);
		}
		return llh;
	}	
};
