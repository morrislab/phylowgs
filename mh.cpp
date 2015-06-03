#include <iostream>
#include<fstream>
#include<string>
#include<time.h>
#include<math.h>
#include<cstdlib>	
#include<cstring>
#include <sstream>
#include<map>

#include "mh.hpp"
#include "util.hpp"

using namespace std;

//  g++ -o mh.o  mh.cpp  util.cpp `gsl-config --cflags --libs`
// ./mh.o 5000 100 11 1 8 4 ssm_data.txt cnv_data.txt c_tree_ssm_data_1.txt c_data_states_ssm_data_1.txt c_params.txt 5
//https://www.gnu.org/software/gsl/manual/html_node/Shared-Libraries.html


// done for multi-sample
int main(int argc, char* argv[]){

	// parse command line args
	struct config conf;
	conf.MH_ITR=atoi(argv[1]);//5000
	conf.MH_STD=atof(argv[2]);//100
	conf.N_SSM_DATA=atoi(argv[3]);//12; // no. of ssm data points
	conf.N_CNV_DATA=atoi(argv[4]);//1; // no. of cnv data points	
	conf.NNODES=atoi(argv[5]); //17, no. of nodes in the tree
	conf.TREE_HEIGHT=atoi(argv[6]);//6
	
	// file names
	char* FNAME_SSM_DATA = argv[7];
	char* FNAME_CNV_DATA = argv[8];
	char* FNAME_C_TREE = argv[9];
	char* FNAME_C_DATA_STATES = argv[10];
	char* FNAME_C_PARAMS = argv[11];
	char* FNAME_C_MH_AR = argv[12];

	conf.NTPS = atoi(argv[13]); // no. of samples 
	
	struct datum *data = new datum[conf.N_SSM_DATA+conf.N_CNV_DATA];
	load_ssm_data(FNAME_SSM_DATA, data,conf);
	if (conf.N_CNV_DATA>0)
		load_cnv_data(FNAME_CNV_DATA,data,conf);
	
	
	struct node *nodes = new node[conf.NNODES];
	load_tree(FNAME_C_TREE,nodes,conf);		
	
	load_data_states(FNAME_C_DATA_STATES,data, nodes, conf);
	
	//start MH loop
	mh_loop(nodes,data,FNAME_C_MH_AR,conf);
	
	// write updated params to disk
	write_params(FNAME_C_PARAMS,nodes,conf);		
	
	return 0;	
}


// done for multi-sample
void mh_loop(struct node nodes[],struct datum data[],char* fname, struct config conf){
	gsl_rng *rand = gsl_rng_alloc(gsl_rng_mt19937);
	double ratio=0.0;
	for (int itr=0;itr<conf.MH_ITR;itr++){
	
		for(int tp=0;tp<conf.NTPS;tp++)
			sample_cons_params(nodes,conf,rand,tp);		
			
		double a = multi_param_post(nodes,data,0,conf)-multi_param_post(nodes,data,1,conf);			
		
		//cout<<multi_param_post(nodes,data,0,conf)-multi_param_post(nodes,data,1,conf)<<'\n';
		
		// loop over samples, apply dirichlet correction terms, update a
		double theta[conf.NNODES];// dirichlet params
		double pi[conf.NNODES],pi_new[conf.NNODES];
		for(int tp=0; tp<conf.NTPS;tp++){		
		
			get_pi(nodes,pi_new,conf,0,tp);
			get_pi(nodes,pi,conf,1,tp);		
		
			// apply the dirichlet correction terms			
			for(int i=0;i<conf.NNODES;i++)
				theta[i]=conf.MH_STD*pi_new[i];
			a += gsl_ran_dirichlet_lnpdf(conf.NNODES,theta,pi);
			
			for(int i=0;i<conf.NNODES;i++)
				theta[i]=conf.MH_STD*pi[i];
			a -= gsl_ran_dirichlet_lnpdf(conf.NNODES,theta,pi_new);
		}
			
		double r = gsl_rng_uniform_pos(rand);
		//cout<<log(r)<<'\t'<<a<<'\n';
		if (log(r)<a){
			ratio+=1;
			update_params(nodes,conf);			
		}
	}
	gsl_rng_free(rand);
	
	ofstream dfile;
	dfile.open(fname);	
	dfile<<ratio/conf.MH_ITR;	
	dfile.close();	

}


// done for multi-sample
void sample_cons_params(struct node nodes[], struct config conf, gsl_rng *rand, int tp){

	map <int, int> node_id_map;
	for(int i=0;i<conf.NNODES;i++)
		node_id_map[nodes[i].id]=i;
	
	int NNODES=conf.NNODES;
	double pi[NNODES];
	for(int i=0;i<NNODES;i++)
		pi[i]=nodes[i].pi[tp];

	// randomly sample from a dirichlet
	double pi_new[NNODES],alpha[NNODES];
	for(int i=0;i<NNODES;i++)
		alpha[i]=conf.MH_STD*pi[i]+1;
	dirichlet_sample(NNODES,alpha,pi_new,rand);

	// update the nodes pi1 (new pi)
	for(int i=0;i<NNODES;i++)
		nodes[i].pi1[tp]=pi_new[i];				

	// update the nodes param1 (new param)
	for(int i=0;i<NNODES;i++){		
		double param = nodes[i].pi1[tp];			
		for(int c=0;c<nodes[i].nchild;c++){
			param+=nodes[node_id_map[nodes[i].cids.at(c)]].param1[tp];
		}
		nodes[i].param1[tp]=param;
	}
}


// done for multi-sample
// todo: double check log_ll
double multi_param_post(struct node nodes[], struct datum data[], int old,struct config conf){
	double post=0.0;
	for(int tp=0;tp<conf.NTPS;tp++)
		post+=param_post(nodes,data,old,conf,tp);
	return post;
}		
double param_post(struct node nodes[], struct datum data[], int old,struct config conf,int tp){	
	double llh = 0.0;
	for(int i=0;i<conf.NNODES;i++){
		double p=0;
		if(old==0)
			p=nodes[i].param1[tp];
		else
			p=nodes[i].param[tp];
		for(int j=0;j<nodes[i].ndata;j++){
			llh+=data[nodes[i].dids.at(j)].log_ll(p,old,tp);
		}
	}
	return llh;	
}


// done for multi-sample
void update_params(struct node nodes[],struct config conf){	
	for(int i=0;i<conf.NNODES;i++){
		for(int tp=0;tp<conf.NTPS;tp++){
			nodes[i].param[tp]=nodes[i].param1[tp];
			nodes[i].pi[tp]=nodes[i].pi1[tp];
		}
	}
}

// done for multi-sample
void get_pi(struct node nodes[], double pi[], struct config conf, int old, int tp){
	for(int i=0;i<conf.NNODES;i++){
		if (old==0)
			pi[i]=nodes[i].pi1[tp];
		else
			pi[i]=nodes[i].pi[tp];
	}
}


// done for multi-sample
void write_params(char fname[], struct node *nodes, struct config conf){	
	ofstream dfile;
	dfile.open(fname);	
	for(int i=0;i<conf.NNODES;i++){		
		dfile<<nodes[i].id<<'\t';	
		for(int tp=0;tp<conf.NTPS;tp++)
			dfile<<nodes[i].param[tp]<<',';
		dfile<<'\t';
		for(int tp=0;tp<conf.NTPS;tp++)
			dfile<<nodes[i].pi[tp]<<',';
		dfile<<'\n';
	}	
	dfile.close();	
}


// done for multi-sample
void load_ssm_data(char fname[],struct datum *data, struct config conf){
	string line,token,token1;
	ifstream dfile (fname);
	int ctr=0,id=-1;
	while (getline (dfile,line,'\n')){
		if (id==-1){id+=1;continue;}
		istringstream iss(line);
		ctr=0;
		while(getline(iss,token,'\t')){
			if(ctr==0){
				data->id=id;
			}
			else if(ctr==2){
				istringstream iss(token);
				for(int tp=0;tp<conf.NTPS;tp++){
					getline(iss,token1,',');
					data->a.push_back(atoi(token1.c_str()));
				}
			}
			else if(ctr==3){
				istringstream iss(token);
				for(int tp=0;tp<conf.NTPS;tp++){
					getline(iss,token1,',');
					data->d.push_back(atoi(token1.c_str()));
				}
				for(int tp=0;tp<conf.NTPS;tp++)
					data->log_bin_norm_const.push_back(log_bin_coeff(data->d[tp],data->a[tp]));
			}
			else if(ctr==4){
				data->mu_r=atof(token.c_str());
			}
			else if(ctr==5){
				data->mu_v=atof(token.c_str());
				data->cnv=NULL; // this will be set to CNV address if exists when loading cnvs
			}
			ctr+=1;			
		}		
		data++;
		id+=1;
	}	
	dfile.close();	
}


// done for multi-sample
void load_cnv_data(char fname[],struct datum *data, struct config conf){
	string line,token,token1,token2;
	ifstream dfile (fname);
	int ctr=0,id=-1,did=conf.N_SSM_DATA;
	
	map <int, struct datum*> id_datum_map;
	for(int i=0;i<conf.N_SSM_DATA;i++)
		id_datum_map[data[i].id]=&data[i];
		
	while (getline (dfile,line,'\n')){
		if (id==-1){id+=1;continue;}
		istringstream iss(line);
		ctr=0;
		while(getline(iss,token,'\t')){
			if(ctr==0){
				data[did].id=did;
			}
			else if(ctr==1){
				istringstream iss(token);
				for(int tp=0;tp<conf.NTPS;tp++){
					getline(iss,token1,',');
					data[did].a.push_back(atoi(token1.c_str()));
				}				
			}
			else if(ctr==2){
				istringstream iss(token);
				for(int tp=0;tp<conf.NTPS;tp++){
					getline(iss,token1,',');
					data[did].d.push_back(atoi(token1.c_str()));
				}
				for(int tp=0;tp<conf.NTPS;tp++)
					data[did].log_bin_norm_const.push_back(log_bin_coeff(data[did].d[tp],data[did].a[tp]));
				data[did].mu_r=0.999;
				data[did].mu_v=0.5;
				data[did].cnv=NULL; // todo: is this needed?
			}
			else if(ctr==3){
				istringstream iss1(token);
				while(getline(iss1,token1,';')){
					int ssm_id;
					istringstream iss2(token1);
					getline(iss2,token2,',');
					ssm_id = atoi(token2.substr(1).c_str());
					id_datum_map[ssm_id]->cnv=&data[did];
				}
			}
			ctr++;			
		}
		id++;
		did++;
	}	
	dfile.close();
}

// done for multi-sample
// load the current tree state from the disk
void load_tree(char fname[], struct node *nodes, struct config conf){
	string line,token,token1;
	ifstream dfile (fname);
	int ctr=0;
	while (getline (dfile,line,'\n')){
		istringstream iss(line);
		ctr=0;
		while(getline(iss,token,'\t')){
			if(ctr==0){
				nodes->id=atoi(token.c_str());
			}
			else if(ctr==1){
				istringstream iss(token);
				for(int tp=0;tp<conf.NTPS;tp++){
					getline(iss,token1,',');
					nodes->param.push_back(atof(token1.c_str()));
					nodes->param1.push_back(0.0);
				}				
			}
			else if(ctr==2){
				istringstream iss(token);
				for(int tp=0;tp<conf.NTPS;tp++){
					getline(iss,token1,',');
					nodes->pi.push_back(atof(token1.c_str()));
					nodes->pi1.push_back(0.0);
				}				
			}
			else if(ctr==5){
				nodes->ndata=atoi(token.c_str());
			}
			else if(ctr==6){
				istringstream iss(token);
				for(int i=0;i<nodes->ndata;i++){
					getline(iss,token1,',');
					nodes->dids.push_back(atoi(token1.c_str()));
				}
			}
			else if(ctr==3){
				nodes->nchild=atoi(token.c_str());		
			}			
			else if(ctr==4){
				istringstream iss(token);
				for(int i=0;i<nodes->nchild;i++){
					getline(iss,token1,',');
					nodes->cids.push_back(atoi(token1.c_str()));
				}
			}
			else if(ctr==7){
				nodes->ht=atoi(token.c_str());		
			}			
			ctr+=1;
		}
		nodes++;
	}	
	dfile.close();	
}

// no changes for multi-sample
void load_data_states(char fname[], struct datum* data, struct node* nodes, struct config conf){

	map <int, struct datum*> id_datum_map;
	for(int i=0;i<conf.N_SSM_DATA;i++){
		id_datum_map[data[i].id]=&data[i];
	}
		
	map <int, int> node_id_map;
	for(int i=0;i<conf.NNODES;i++)
		node_id_map[nodes[i].id]=i;
	
	string line,token,token1,token2;
	ifstream dfile (fname);
	int ctr=0;
	while (getline (dfile,line,'\n')){
		
			
		struct datum *dat;
		istringstream iss(line);
		ctr=0;
		while(getline(iss,token,'\t')){
			if(ctr==0){
				dat = id_datum_map[atoi(token.c_str())];
			}
			else if(ctr==1){			
				istringstream iss1(token);
				while(getline(iss1,token1,';')){					
					struct state st8;
					istringstream iss2(token1);
					for(int i=0;i<3;i++){
						getline (iss2,token2,',');
						if(i==0)
							st8.nd = &nodes[node_id_map[atoi(token2.c_str())]];
						else if(i==1)
							st8.nr = atoi(token2.c_str());
						else if(i==2)
							st8.nv = atoi(token2.c_str());
					}
					dat->states1.push_back(st8);
				}
			}
			else if(ctr==2){			
				istringstream iss1(token);
				while(getline(iss1,token1,';')){					
					struct state st8;
					istringstream iss2(token1);
					for(int i=0;i<3;i++){
						getline (iss2,token2,',');
						if(i==0)
							st8.nd = &nodes[node_id_map[atoi(token2.c_str())]];
						else if(i==1)
							st8.nr = atoi(token2.c_str());
						else if(i==2)
							st8.nv = atoi(token2.c_str());
					}
					dat->states2.push_back(st8);
				}
			}
			ctr+=1;
		}		
	}	
	dfile.close();
}
