#define Neuron_no 2400 //(need to be 5*N) Excitatory:2000(Neuron_no*4/5)   Inhibitory:500
#define delta_t 0.05 //ms
#define exc1 2000
#define inh1 2200
#define inh2 2400
#define FFI_size 200
#define memory 200

#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>
#include<string.h>
#include<deque>
#include<vector>
#include <algorithm>

using namespace std;

#include "RandNum.h"


inline double current (double C1, double C2, double tau_d,double tau_r);
inline double dV_dt(double V,double G_kE,double G_kI,double tau, double tau_r, double V_l, double E_E, double E_I);
void evolution (double x, double y, int count, int seed, double xx, double yy, double theta, double gamma, int pulse_number, int learn);


 

inline double current (double C1, double C2,double tau_d,double tau_r)
{
	double S;
	S=(C1-C2)/(tau_d-tau_r);
	return S;
}

inline double dV_dt(double V,double G_kE,double G_kI,double tau, double tau_r, double V_L, double E_E, double E_I)
{
	double V_dot;
	
	V_dot=(V_L-V+G_kE*(E_E-V)+G_kI*(E_I-V))/tau;//unit of delta_t
	
	return V_dot;
}


void evolution (double x, double y, int count, int seed, double yy, double theta, double gamma, int pulse_number, int learn)
{
	int theta_counter;
	
	CRandNum randnum;
	double dt;
	
	double tau_d_E;
	double tau_d_E_NMDA; 
	double tau_d_I;
	double tau_r;
	double tau_r_NMDA;
	double tau_l;
	double tau_E;
	double tau_I;
	
	double E_E;//mV
	double E_I;//mV
	double V_L;//mV
	
	double g_L;
	double g_EO;
	double g_IO;
	double g_EE_L1;
	double g_EE_L2;
	double g_EI;
	double g_IE;
	double g_II;
	
	double link_density;
	double link_density_E;
	double f_ex;
	
	double V_th;
	double V_rest;
	
	double refractory_period_E;
	double refractory_period_I;
	
	double V[Neuron_no];
	int layer[Neuron_no];
	int connection[Neuron_no][10];
	double type[Neuron_no];
	double S_E_C_d[Neuron_no][4];
	double S_E_C_r[Neuron_no][4];
	double S_I_C_d[Neuron_no][4];
	double S_I_C_r[Neuron_no][4];
	double S_O_C_d[Neuron_no];
	double S_O_C_r[Neuron_no];
	double S_E_C_d_NMDA[Neuron_no][4];
	double S_E_C_r_NMDA[Neuron_no][4];
	double sum_layer_S_E_C_d;
	double sum_layer_S_E_C_r;
	double sum_layer_S_I_C_d;
	double sum_layer_S_I_C_r;
	double sum_layer_S_E_C_d_NMDA;
	double sum_layer_S_E_C_r_NMDA;
	double sum_layer_V;
	
	double alpha_NMDA,alpha_AMPA;
	
	
	double G_kE;
	double G_kE_delta_t;
	double G_kI;
	double G_kI_delta_t;
	
	
	deque<double> current_arrival_time[Neuron_no];
	deque<int> spikefrom[Neuron_no];
	
	static double link[Neuron_no][Neuron_no];
	static double link_target[Neuron_no][Neuron_no];
	deque<int> need_modify_link_i;
	deque<int> need_modify_link_j;
	
	deque<double> link_modification;
	
	need_modify_link_i.clear();
 	need_modify_link_j.clear();
	link_modification.clear();
	
	
	double external_input_time[Neuron_no];
	double t_refractory_end[Neuron_no];
	double t_spike;
	
	int count_active_E_layer[4];
	int count_active_I_layer[4];
	
	int count_active_E_M_layer[4];
	int count_active_E_NM_layer[4];
	int count_active_E_memory_layer[4][10];

	
	//STP prarmeter
	double Cu_own[Neuron_no];//own u
	double Cx_own[Neuron_no];//own x
	deque<double> Cu_arrival[Neuron_no];
	deque<double> Cx_arrival[Neuron_no];
	double tau_F;
	double tau_D;
	double U;
	
	
	//memory
	int memory_active_no;
	int active_separation;
	int first_memory_active;
	int active_duration;
	int memory_no[Neuron_no];//0-9; inh=-1;
	double As;
	int number_of_memory;
	
	double trigger_cycle;
	
	
	//STDP
	int STDP_activation;
	double neuron_trace[Neuron_no];
	double neuron_trace_slow[Neuron_no];
	double tau_STDP;
	double tau_STDP_slow;
	double tau_mGluR7;
	double learning_rate;
	
	double beta,beta_IE,delta,delta_IE;
	double neuron_count_trace[Neuron_no];
	double tau_count;

	randnum.InitGenRand(seed);
	

	dt=1;//in unit delta_t
	
	tau_d_E=x/delta_t; //in unit delta_t ,e=2-5ms,i=5-15ms
	tau_d_E_NMDA=100/delta_t;//500
	tau_d_I=y/delta_t; //in unit delta_t
	tau_r=0.5/delta_t; //in unit delta_t
	tau_r_NMDA=20/delta_t;
	tau_l=1/delta_t; //in unit delta_t
	tau_E=20/delta_t; //in unit delta_t
	tau_I=10/delta_t; //in unit delta_t
	tau_mGluR7=100/delta_t;
	g_EO=0.05;
	g_IO=0.08;
	g_EE_L1=0.02;
	g_EE_L2=0.05;
	g_EI=0.9;
	g_IE=0.15; 
	g_II=0.48;
	
	alpha_NMDA=0.8;
	alpha_AMPA=1.1;

	

	link_density=yy;
	link_density_E=yy;
	f_ex=2.5*delta_t*1e-3;

	V_th=-50;//mV
	V_rest=-55;//mV 
	E_E=0;//mV 
	E_I=-70;//mV
	V_L=-70;//mV
	
	As=3.5;
	U=0.2;
	tau_F=1500/delta_t;
	tau_D=200/delta_t;
	first_memory_active=20000/delta_t;
	active_duration=350/delta_t;
	active_separation=300/delta_t;
	number_of_memory=10;
	
	//STDP
	STDP_activation=0;
	tau_STDP=20/delta_t;
	tau_STDP_slow=100/delta_t;
	learning_rate=0.006;
	beta=0.000001;
	beta_IE=0.0002; 
	delta=1*1e-6;
	delta_IE=1*1e-6;
	tau_count=200/delta_t;
	
	
	refractory_period_E=2/delta_t; //in unit delta_t
	refractory_period_I=1/delta_t; //in unit delta_t
	
	FILE *fp0;
	FILE *fp1;
	FILE *fp2;
	FILE *fp3;
	FILE *fp4;
	FILE *fp5;
	FILE *fp6;
	
	int i,j,u,s,t,t_record;
	
	double tau,g_kO,g_kE,g_kI,refractory_period;
	
	int temp_input,temp_input2;
	double temp_S_E_C_d;
	double temp_S_E_C_d_NMDA;
	double temp_S_E_C_r_NMDA;
	double temp_S_E_C_r;
	double temp_S_I_C_d;
	double temp_S_I_C_r;
	double temp_S_O_C_d;
	double temp_S_O_C_r;
	
	double temp_S_E_C_d_m;
	double temp_S_E_C_r_m;
	
	double temp_Cu_own;
	double temp_Cx_own;
	
	double k1,k2,temp,temp1,temp2;
	int temp_int;
	double sum_V_E_layer[4];
	double sum_V_E_M_layer[4];
	double sum_V_E_NM_layer[4];
	double sum_V_E_memory_layer[4][10];
	
	double sum_V_I_layer[4];
	double sum_current_E_layer[4];
	double balance_E[4];
	
	double balance_E_M[4][10];
	double sum_current_E_memory_layer[4][10];

	
	double sum_current_I_layer[4];
	double balance_I[4];
	
	double sum_current_I_M_layer[4];
	double sum_current_I_NM_layer[4];
	double sum_current_I_memory_layer[4][10];
	
	double sum_current_E_I_layer[4];
	double sum_current_I_I_layer[4];
	double sum_U_layer[4];
	double sum_x_layer[4];
	double sum_U_M_layer[4];
	double sum_x_M_layer[4];
	double sum_U_NM_layer[4];
	double sum_x_NM_layer[4];
	double sum_U_memory_layer[4][10];
	double sum_x_memory_layer[4][10];
	
	double update_constant_d_E,update_constant_r,update_constant_d_I;
	double update_constant_Cu,update_constant_Cx;
	double update_constant_neuron_trace,update_constant_neuron_trace_slow,update_constant_neuron_trace_10,update_constant_neuron_trace_5;
	double update_constant_neuron_count_trace;
	double update_constant_d_E_NMDA,update_constant_r_NMDA;
	
	double EM_sum_neuron_count_trace_layer[4];
	double E_memory_sum_neuron_count_trace_layer[4][10];
	
	double E_sum_neuron_count_trace_layer[4];
	double I_sum_neuron_count_trace_layer[4];
	double E_sum_neuron_trace_layer[4];
	double EM_sum_neuron_trace_layer[4];
	double E_memory_sum_neuron_trace_layer[4][10];
	
	double I_sum_neuron_trace_layer[4];
	double EM_mGluR7_sum[4];
	double E_mGluR7_sum[4];
	double E_memory_mGluR7_sum[4][10];
	double update_constant_mGluR7;
	
	int memory_active;
	int counter,counter2;
	int time_count;
	
	double sum_NMDA_2_layer[4];
	double sum_NMDA_layer[4][10];
	
	double C_NMDA=1;
	double sum_hytero;
	
	double a_GE_E;
	double a_GI_E;
	double a_GE_I;
	double a_GI_I;
	double average_trace_slow;
	double average_V;
		
	double temp_link; 
	double max;
	
	deque<int> link_record_EM3_EM1_i[10];
	deque<int> link_record_EM3_EM1_j[10];
	
	deque<int> link_record_E3_EM1_i[10];
	deque<int> link_record_E3_EM1_j[10];
	
	deque<int> link_record_EM3_EM3_i[10];
	deque<int> link_record_EM3_EM3_j[10];
	
	deque<int> link_record_E3_E3_i[10];
	deque<int> link_record_E3_E3_j[10];
	
	deque<int> link_record_E3_EM3_i[10];
	deque<int> link_record_E3_EM3_j[10];
	
	deque<int> link_record_EM3_E3_i[10];
	deque<int> link_record_EM3_E3_j[10];
	
	deque<int> link_record_IM3_EM1_i[10];
	deque<int> link_record_IM3_EM1_j[10];
	
	deque<int> link_record_EM3_IM3_i[10];
	deque<int> link_record_EM3_IM3_j[10];
	deque<int> link_record_EM3_I3_i[10];
	deque<int> link_record_EM3_I3_j[10];
	
	int count_number;
	int flag;
	
	
	for(i=0;i<number_of_memory;i++)
	{
		
		link_record_EM3_EM1_i[i].clear();
		link_record_EM3_EM1_j[i].clear();
		link_record_E3_EM1_i[i].clear();
		link_record_E3_EM1_j[i].clear();
		link_record_EM3_EM3_i[i].clear();
		link_record_EM3_EM3_j[i].clear();
		link_record_E3_E3_i[i].clear();
		link_record_E3_E3_j[i].clear();
		link_record_E3_EM3_i[i].clear();
		link_record_E3_EM3_j[i].clear();
		link_record_EM3_E3_i[i].clear();
		link_record_EM3_E3_j[i].clear();
		link_record_IM3_EM1_i[i].clear();
		link_record_IM3_EM1_j[i].clear();
		link_record_EM3_IM3_i[i].clear();
		link_record_EM3_IM3_j[i].clear();
		link_record_EM3_I3_i[i].clear();
		link_record_EM3_I3_j[i].clear();
		
	}
	
	
	update_constant_d_E=exp(-dt/tau_d_E);
	update_constant_r=exp(-dt/tau_r);
	update_constant_d_I=exp(-dt/tau_d_I);
	update_constant_Cu=exp(-dt/tau_F);
	update_constant_Cx=exp(-dt/tau_D);
	update_constant_neuron_trace=exp(-dt/tau_STDP);
	update_constant_neuron_trace_slow=exp(-dt/tau_STDP_slow);
	update_constant_neuron_count_trace=exp(-dt/tau_count);
	
	update_constant_d_E_NMDA=exp(-dt/tau_d_E_NMDA);
	update_constant_r_NMDA=exp(-dt/tau_r_NMDA);
	update_constant_mGluR7=exp(-dt/tau_mGluR7);

	t_record=(int)(1/delta_t);
	char filename_0[100];
	char filename_1[100];
	char filename_2[100];
	char filename_3[100];
	char filename_4[100];
	char filename_5[100];
	char filename_6[100];

	sprintf(filename_2,"E=%.2fI=%.2fP=%.2f_spike_%d_0_%.2f.txt",x,y,link_density,count,yy);
	fp2=fopen(filename_2,"w");
	fclose(fp2);
	
	for(i=0;i<Neuron_no;i++)
	{
		for (j=0;j<number_of_memory;j++)
		{
			connection[i][j]=0;
		} 
		
	}
	
	
	for(i=0;i<Neuron_no;i++)
	{
		neuron_count_trace[i]=0;
		neuron_trace[i]=0;
		neuron_trace_slow[i]=0;
		
		V[i]=V_rest+(V_th-V_rest)*randnum.GenRandReal_10();
		
		current_arrival_time[i].clear();
		spikefrom[i].clear();
		
		t_refractory_end[i]=-1;
		
		for(j=0;j<3;j++)
		{
			S_E_C_d[i][j]=0;
			S_E_C_r[i][j]=0;
			S_I_C_d[i][j]=0;
			S_I_C_r[i][j]=0;
			S_E_C_d_NMDA[i][j]=0;
			S_E_C_r_NMDA[i][j]=0;
		}
		S_O_C_d[i]=0;
		S_O_C_r[i]=0;
		
		
			
		Cu_own[i]=0;
		Cx_own[i]=0;
		
		Cu_arrival[i].clear();
		Cx_arrival[i].clear();
		
		
		if(i<memory*number_of_memory)
		{
			memory_no[i]=(int)(i/memory);
		}
		else
		{
			memory_no[i]=-1;

			
		}	
		
	}
	
	
	G_kE=0;
	G_kI=0;	
	
	
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////	
	int ii,jj;
	int z=0;
	

	
	for(i=0;i<Neuron_no;i++)
	{
		if(i<inh1)
		{
			layer[i]=0;
			if(i<exc1)
			{
				type[i]=0;//exc
			}
			else
			{
				type[i]=1;//inh
			}
		}
		if(i>=inh1 && i<inh2)
		{
			layer[i]=1;
			type[i]=1;
		}
	}
	

///////////////////////////////////////////////////////////////////////////////////////////////////	
	for(i=0;i<Neuron_no;i++)
	{
		for(j=0;j<Neuron_no;j++)
		{
			
			link[i][j]=0;
			link_target[i][j]=0;
			if(layer[i]==layer[j] && layer[i]==0)
			{
				if(type[i]==1)  
				{
					if(type[j]==1 && randnum.GenRandReal_10()<=link_density && i!=j)
					{
						link[i][j]=g_II;
						link_target[i][j]=g_II;
					}
					if(type[j]==0 && randnum.GenRandReal_10()<=link_density_E && i!=j)
					{
						link[i][j]=g_IE;
						link_target[i][j]=g_IE;
					}
				}
				if(type[i]==0)    
				{
					if(type[j]==0 && randnum.GenRandReal_10()<=link_density_E && i!=j)
					{
						link[i][j]=g_EE_L1;
						link_target[i][j]=g_EE_L1;
						
						
							
					}
					if(type[j]==1 && randnum.GenRandReal_10()<=link_density && i!=j)
					{
						link[i][j]=g_EI;
						link_target[i][j]=g_EI;
					}
				}
			}
			if(layer[i]==layer[j] && layer[i]==1) 
			{
				if(randnum.GenRandReal_10()<=link_density && i!=j)
				{
					link[i][j]=g_II;
				}	
			}
			
			//layer0 exc to layer1 exc
			if(layer[i]!=layer[j])
			{
				if(type[i]==1 && type[j]==0 && randnum.GenRandReal_10()<=link_density_E)    //0.05
				{
					link[i][j]=g_IE;//0.5
					link_target[i][j]=g_IE;
					
	
				}  
				if(type[i]==1 && type[j]==1 && randnum.GenRandReal_10()<=link_density)    //0.05
				{
					link[i][j]=g_II;;//0.5
					link_target[i][j]=g_II;
					

				}  
				if(type[i]==0 && type[j]==1 && randnum.GenRandReal_10()<=link_density_E)    //0.05
				{
					link[i][j]=g_EI;//0.5
					link_target[i][j]=g_EI;
					

				}  
				if(type[i]==0 && type[j]==0 && randnum.GenRandReal_10()<=link_density_E)    //0.05
				{
					link[i][j]=g_EE_L1;//0.5
					link_target[i][j]=g_EE_L1;
					

				}  
			}
			
			
		}
	}
	
	


	
	
	if(learn==0)
	{
		for (int k=0;k<number_of_memory;k++)	
		{
			for (i=k*200;i<(k+1)*200;i++)
			{
				for (j=k*200;j<(k+1)*200;j++)
				{
					if(link[i][j]!=0)
					{
						
						
						link[i][j]=0.5;
						link_target[i][j]=0.5;
					}
				}
			}
		}
	}
	if(learn==1)
	{
		for (int k=1;k<number_of_memory;k++)	
		{
			for (i=k*200;i<(k+1)*200;i++)
			{
				for (j=k*200;j<(k+1)*200;j++)
				{
					if(link[i][j]!=0)
					{
						
						
						link[i][j]=0.5;
						link_target[i][j]=0.5;
					}
				}
			}
		}
	}
	
	
	
	
	
	
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////








	//#pragma omp parallel for
	for(i=0;i<Neuron_no;i++)
	{
		external_input_time[i]=(-1*log(1-randnum.GenRandReal_10())/f_ex/(Neuron_no*4/5*link_density));
	}
	
	//printf("%f\t%f\n",S_O_C_d[2200],external_input_time[2200]);
	//system("PAUSE");
	printf("network finish\n");
	





	

	//Evolution
	memory_active=0;
	time_count=-1;
	
	
	int memory_load_in_last_for=10/delta_t;
	FILE *fp20;
	FILE *fp21;
	FILE *fp22;
	FILE *fp23;
	FILE *fp24;
	FILE *fp25;

		
	sprintf(filename_3,"test_%d_AMPA.txt",count);
	fp20=fopen(filename_3,"w");
	
	
	sprintf(filename_3,"test_%d_NMDA.txt",count);
	fp21=fopen(filename_3,"w");	
	
	
	sprintf(filename_3,"test_%d_GABA_L1.txt",count);
	fp22=fopen(filename_3,"w");	
	
	
	sprintf(filename_3,"test_%d_V.txt",count);
	fp23=fopen(filename_3,"w");	
	
	
	sprintf(filename_3,"test_%d_GABA_L2.txt",count);
	fp25=fopen(filename_3,"w");	

	
	
	
	int retrieve_cue_last_for=140000/delta_t;
	int memory_cue=0;
	int memory_cue_no=-2;
	double pulse_active=0;
	int start;
	if(learn==1)
		start=0;
	else
		start=(int)(100000/delta_t);
	
	for(t=start;t<(100000+2*number_of_memory*15000)/delta_t;t++)
	{

		
		
		

		if(t%(int)(1/delta_t)==0 && t>=(int)((1000)/delta_t))
		{
			sprintf(filename_3,"test_%d_AMPA.txt",count);
			fp20=fopen(filename_3,"a");
			sprintf(filename_3,"test_%d_NMDA.txt",count);
			fp21=fopen(filename_3,"a");
			sprintf(filename_3,"test_%d_GABA_L1.txt",count);
			fp22=fopen(filename_3,"a");
			sprintf(filename_3,"test_%d_V.txt",count);
			fp23=fopen(filename_3,"a");
			sprintf(filename_3,"test_%d_GABA_L2.txt",count);
			fp25=fopen(filename_3,"a");
			
			fprintf(fp20,"%f\t",t*delta_t);
			fprintf(fp21,"%f\t",t*delta_t);
			fprintf(fp22,"%f\t",t*delta_t);
			fprintf(fp23,"%f\t",t*delta_t);
			fprintf(fp25,"%f\t",t*delta_t);
			
		
			for (j=0;j<number_of_memory;j++)
			{
				sum_layer_V=0;
				for (i=j*memory;i<(j+1)*memory;i++)
				{
					sum_layer_V+=V[i];
				}
				fprintf(fp23,"%f\t",sum_layer_V);
			}
			
			sum_layer_V=0;
			for (i=exc1;i<inh1;i++)
			{
				sum_layer_V+=V[i];
			}
			fprintf(fp23,"%f\t",sum_layer_V);
			
			sum_layer_V=0;
			for (i=inh1;i<inh2;i++)
			{
				sum_layer_V+=V[i];
			}
			fprintf(fp23,"%f\t",sum_layer_V);
			
			//AMPA
			for (j=0;j<number_of_memory;j++)
			{
				sum_layer_S_E_C_d=0;
				sum_layer_S_E_C_r=0;

				for (i=j*memory;i<(j+1)*memory;i++)
				{
					sum_layer_S_E_C_d+=S_E_C_d[i][0];
					sum_layer_S_E_C_r+=S_E_C_r[i][0];
				}
				fprintf(fp20,"%f\t",tau_E*alpha_AMPA*current(sum_layer_S_E_C_d,sum_layer_S_E_C_r,tau_d_E,tau_r));
			}
			
			
			sum_layer_S_E_C_d=0;
			sum_layer_S_E_C_r=0;
			for (i=exc1;i<inh1;i++)
			{
				sum_layer_S_E_C_d+=S_E_C_d[i][0];
				sum_layer_S_E_C_r+=S_E_C_r[i][0];
			}
			fprintf(fp20,"%f\t",tau_E*alpha_AMPA*current(sum_layer_S_E_C_d,sum_layer_S_E_C_r,tau_d_E,tau_r));
			
			sum_layer_S_E_C_d=0;
			sum_layer_S_E_C_r=0;
			for (i=inh1;i<inh2;i++)
			{
				sum_layer_S_E_C_d+=S_E_C_d[i][0];
				sum_layer_S_E_C_r+=S_E_C_r[i][0];
			}
			fprintf(fp20,"%f\t",tau_E*alpha_AMPA*current(sum_layer_S_E_C_d,sum_layer_S_E_C_r,tau_d_E,tau_r));
			
			
			//GABA
			for (j=0;j<number_of_memory;j++)
			{
				sum_layer_S_I_C_d=0;
				sum_layer_S_I_C_r=0;
				for (i=j*memory;i<(j+1)*memory;i++)
				{
					sum_layer_S_I_C_d+=S_I_C_d[i][0];
					sum_layer_S_I_C_r+=S_I_C_r[i][0];
				}
				fprintf(fp22,"%f\t",tau_I*current(sum_layer_S_I_C_d,sum_layer_S_I_C_r,tau_d_I,tau_r));
			}
			for (j=0;j<number_of_memory;j++)
			{
				
				sum_layer_S_I_C_d=0;
				sum_layer_S_I_C_r=0;
				for (i=j*memory;i<(j+1)*memory;i++)
				{
					sum_layer_S_I_C_d+=S_I_C_d[i][1];
					sum_layer_S_I_C_r+=S_I_C_r[i][1];
				}
				fprintf(fp25,"%f\t",tau_I*current(sum_layer_S_I_C_d,sum_layer_S_I_C_r,tau_d_I,tau_r));
			}

			sum_layer_S_I_C_d=0;
			sum_layer_S_I_C_r=0;
			for (i=exc1;i<inh1;i++)
			{
				sum_layer_S_I_C_d+=S_I_C_d[i][0];
				sum_layer_S_I_C_r+=S_I_C_r[i][0];
			}
			fprintf(fp22,"%f\t",tau_I*current(sum_layer_S_I_C_d,sum_layer_S_I_C_r,tau_d_I,tau_r));
			sum_layer_S_I_C_d=0;
			sum_layer_S_I_C_r=0;
			for (i=exc1;i<inh1;i++)
			{
				sum_layer_S_I_C_d+=S_I_C_d[i][1];
				sum_layer_S_I_C_r+=S_I_C_r[i][1];
			}
			fprintf(fp25,"%f\t",tau_I*current(sum_layer_S_I_C_d,sum_layer_S_I_C_r,tau_d_I,tau_r));
			
			sum_layer_S_I_C_d=0;
			sum_layer_S_I_C_r=0;
			for (i=inh1;i<inh2;i++)
			{
				sum_layer_S_I_C_d+=S_I_C_d[i][0];
				sum_layer_S_I_C_r+=S_I_C_r[i][0];
			}
			fprintf(fp22,"%f\t",tau_I*current(sum_layer_S_I_C_d,sum_layer_S_I_C_r,tau_d_I,tau_r));
			sum_layer_S_I_C_d=0;
			sum_layer_S_I_C_r=0;
			for (i=inh1;i<inh2;i++)
			{
				sum_layer_S_I_C_d+=S_I_C_d[i][1];
				sum_layer_S_I_C_r+=S_I_C_r[i][1];
			}
			fprintf(fp25,"%f\t",tau_I*current(sum_layer_S_I_C_d,sum_layer_S_I_C_r,tau_d_I,tau_r));
			
			//NMDA
			for (j=0;j<number_of_memory;j++)
			{
				sum_layer_S_E_C_d_NMDA=0;
				sum_layer_S_E_C_r_NMDA=0;
				for (i=j*memory;i<(j+1)*memory;i++)
				{
					sum_layer_S_E_C_d_NMDA+=S_E_C_d_NMDA[i][0];
					sum_layer_S_E_C_r_NMDA+=S_E_C_r_NMDA[i][0];
				}
				fprintf(fp21,"%f\t",tau_E*alpha_NMDA*current(sum_layer_S_E_C_d_NMDA,sum_layer_S_E_C_r_NMDA,tau_d_E_NMDA,tau_r_NMDA));
			}

			sum_layer_S_E_C_d_NMDA=0;
			sum_layer_S_E_C_r_NMDA=0;
			for (i=exc1;i<inh1;i++)
			{
				sum_layer_S_E_C_d_NMDA+=S_E_C_d_NMDA[i][0];
				sum_layer_S_E_C_r_NMDA+=S_E_C_r_NMDA[i][0];
			}
			fprintf(fp21,"%f\t",tau_E*alpha_NMDA*current(sum_layer_S_E_C_d_NMDA,sum_layer_S_E_C_r_NMDA,tau_d_E_NMDA,tau_r_NMDA));
			
			sum_layer_S_E_C_d_NMDA=0;
			sum_layer_S_E_C_r_NMDA=0;
			for (i=inh1;i<inh2;i++)
			{
				sum_layer_S_E_C_d_NMDA+=S_E_C_d_NMDA[i][0];
				sum_layer_S_E_C_r_NMDA+=S_E_C_r_NMDA[i][0];
			}
			fprintf(fp21,"%f\t",tau_E*alpha_NMDA*current(sum_layer_S_E_C_d_NMDA,sum_layer_S_E_C_r_NMDA,tau_d_E_NMDA,tau_r_NMDA));
			
				
			fprintf(fp20,"\n");
			fclose(fp20);
			fprintf(fp21,"\n");
			fclose(fp21);
			fprintf(fp22,"\n");
			fclose(fp22);
			fprintf(fp23,"\n");
			fclose(fp23);
			fprintf(fp25,"\n");
			fclose(fp25);
		}
		
		fp2=fopen(filename_2,"a");
		
		
		if(t>=(int)(1000/delta_t) && t<(int)((100000)/delta_t)) 
		{
			pulse_active=0.5;
			STDP_activation=1;
			if(t==(int)(1000/delta_t))
			{
				
				counter=0;
			}
			if(counter==(int)(0/delta_t))
			{
				memory_active=1;
				memory_cue_no=0;
				memory_cue=2;
			}
			
			
			counter++;
		}
		if(t>=(int)((100000)/delta_t))
		{
			
			if(t==(int)((100000)/delta_t))
			{
				memory_active=1;
				counter=0;
				memory_cue_no=10;
				pulse_active=0.5;
			}
			
			
				
			if(t>=(int)((100000)/delta_t) && t<(int)((100000+number_of_memory*15000)/delta_t))
			{
				STDP_activation=0;
				if(counter==((int)(15000/delta_t)))
				{
					counter=0;
				}
				if(counter==((int)(0/delta_t)))
				{
					memory_active=1;
					memory_cue_no--;
					if(memory_cue_no<0)
						memory_cue_no=number_of_memory;
						
					memory_cue=1;
				}
				if(counter==((int)(5000/delta_t)))
				{
					memory_active=1;
					memory_cue=0;
					
				}
				if(counter==((int)(13000/delta_t)))
				{
					memory_active=2;
				}
				
			}
			else
			{
				memory_active=0;
			}
			counter++;
			 
		}
		
		for(i=0;i<3;i++)
		{
			sum_V_E_layer[i]=0;
			sum_V_I_layer[i]=0;
			count_active_E_layer[i]=0;
			count_active_I_layer[i]=0;
			sum_current_E_layer[i]=0;
			sum_current_I_layer[i]=0;
			sum_current_E_I_layer[i]=0;
			sum_current_I_I_layer[i]=0;
			E_sum_neuron_count_trace_layer[i]=0;
			E_sum_neuron_trace_layer[i]=0;
			I_sum_neuron_count_trace_layer[i]=0;
			I_sum_neuron_trace_layer[i]=0;
			
			sum_U_layer[i]=0;
			sum_x_layer[i]=0;
			
			balance_E[i]=0;
			balance_I[i]=0;;
			
			for(j=0;j<number_of_memory;j++)
			{
				sum_V_E_memory_layer[i][j]=0;
				count_active_E_memory_layer[i][j]=0;
				sum_current_E_memory_layer[i][j]=0;
				sum_current_I_memory_layer[i][j]=0;
				sum_U_memory_layer[i][j]=0;
			    sum_x_memory_layer[i][j]=0;
			    E_memory_sum_neuron_count_trace_layer[i][j]=0;
			    E_memory_sum_neuron_trace_layer[i][j]=0;
			    E_memory_mGluR7_sum[i][j]=0;
			    sum_NMDA_layer[i][j]=0;
			    balance_E_M[i][j]=0;
			}
		}

		a_GE_E=0;
		a_GI_E=0;
		a_GE_I=0;
		a_GI_I=0;
		average_trace_slow=0;
		average_V=0;
		
		
		if((t-(int)(1000/delta_t))%(int)(1000/theta/delta_t)==0 && (memory_active==1 || memory_active==2))
		{
			theta_counter=0;
		}
		
		for(i=0;i<Neuron_no;i++)
		{
			
			if(type[i]==0)
			{
				tau=tau_E;
				g_kO=g_EO;
				refractory_period=refractory_period_E;
			}
			else
			{
				tau=tau_I;
				g_kO=g_IO;
				refractory_period=refractory_period_I;
			}
			
			
			
			
			//voltage update v(t) to v(t+dt)
			sum_layer_S_E_C_d=0;
			sum_layer_S_E_C_r=0;
			sum_layer_S_I_C_d=0;
			sum_layer_S_I_C_r=0;
			sum_layer_S_E_C_d_NMDA=0;
			sum_layer_S_E_C_r_NMDA=0;
			for(j=0;j<3;j++)
			{
				sum_layer_S_E_C_d+=S_E_C_d[i][j];
				sum_layer_S_E_C_r+=S_E_C_r[i][j];
				sum_layer_S_I_C_d+=S_I_C_d[i][j];
				sum_layer_S_I_C_r+=S_I_C_r[i][j];
				sum_layer_S_E_C_d_NMDA+=S_E_C_d_NMDA[i][j];
				sum_layer_S_E_C_r_NMDA+=S_E_C_r_NMDA[i][j];
			}
			if(t>=t_refractory_end[i])
			{
					
				//RK2 V_update
				if(type[i]==0)
				{
					G_kE_delta_t=tau*(current(S_O_C_d[i],S_O_C_r[i],tau_d_E,tau_r)+alpha_AMPA*current(sum_layer_S_E_C_d,sum_layer_S_E_C_r,tau_d_E,tau_r)+alpha_NMDA*current(sum_layer_S_E_C_d_NMDA,sum_layer_S_E_C_r_NMDA,tau_d_E_NMDA,tau_r_NMDA));
				}
				else
					G_kE_delta_t=tau*(current(S_O_C_d[i],S_O_C_r[i],tau_d_E,tau_r)+current(sum_layer_S_E_C_d,sum_layer_S_E_C_r,tau_d_E,tau_r));

				G_kI_delta_t=tau*((current(sum_layer_S_I_C_d,sum_layer_S_I_C_r,tau_d_I,tau_r)));
																				 
				k1=dV_dt(V[i],G_kE_delta_t,G_kI_delta_t,tau, tau_r,V_L, E_E,E_I);
				
				temp_S_E_C_d=sum_layer_S_E_C_d*update_constant_d_E;
				temp_S_E_C_r=sum_layer_S_E_C_r*update_constant_r;
				temp_S_E_C_d_NMDA=sum_layer_S_E_C_d_NMDA*update_constant_d_E_NMDA;
				temp_S_E_C_r_NMDA=sum_layer_S_E_C_r_NMDA*update_constant_r_NMDA;
				temp_S_I_C_d=sum_layer_S_I_C_d*update_constant_d_I;
				temp_S_I_C_r=sum_layer_S_I_C_r*update_constant_r;
				temp_S_O_C_d=S_O_C_d[i]*update_constant_d_E;
				temp_S_O_C_r=S_O_C_r[i]*update_constant_r;
				
				if(type[i]==0)
				{
					G_kE_delta_t=tau*(current(temp_S_O_C_d,temp_S_O_C_r,tau_d_E,tau_r)+alpha_AMPA*current(temp_S_E_C_d,temp_S_E_C_r,tau_d_E,tau_r)+alpha_NMDA*current(temp_S_E_C_d_NMDA,temp_S_E_C_r_NMDA,tau_d_E_NMDA,tau_r_NMDA));
				}		
				else
					G_kE_delta_t=tau*(current(temp_S_O_C_d,temp_S_O_C_r,tau_d_E,tau_r)+current(temp_S_E_C_d,temp_S_E_C_r,tau_d_E,tau_r));
			
				G_kI_delta_t=tau*(current(temp_S_I_C_d,temp_S_I_C_r,tau_d_I,tau_r));
				
					
				k2=dV_dt(V[i]+k1*dt,G_kE_delta_t,G_kI_delta_t,tau, tau_r,V_L, E_E,E_I);
				
				
				if(k1/k2<0 && k2!=0)
						temp=k1*dt;
				else
					temp=((k1+k2)/2)*dt;
			
					
				if(theta_counter%(int)(floor(1000/gamma/delta_t))>0 && theta_counter%(int)(floor(1000/gamma/delta_t))<=(int)(floor(1000/theta/2/delta_t)) && (int)(theta_counter/(floor(1000/gamma/delta_t)))<pulse_number && i>=inh2-FFI_size && (memory_active==1 || memory_active==2) && t>=t_refractory_end[i] && randnum.GenRandReal_10()<pulse_active)
				{
					V[i]=-70;
				}
				
				
				
				//predict firing in next step
				if(V[i]+temp>=V_th)
				{	
					
					//spike calibration
					t_spike=(V_th-V[i])/temp*dt;
					
					
					
					
					
					if(type[i]==0)
					{
						temp_Cu_own=Cu_own[i]*exp(-(t_spike)/tau_F);
						temp_Cx_own=Cx_own[i]*exp(-(t_spike)/tau_D);
					}
					//current_wait_list_update(i,t,linkout,spikefrom,current_arrival_time);
					for(j=0;j<Neuron_no;j++)
					{
						if(link[j][i]!=0)
						{
							spikefrom[j].push_back(i);
							current_arrival_time[j].push_back(t+t_spike+tau_l);
							if(type[i]==0)
							{
								Cu_arrival[j].push_back(temp_Cu_own);
								Cx_arrival[j].push_back(temp_Cx_own);
							}
						}
						//STDP
						if(link[i][j]!=0 && t>1000/delta_t && STDP_activation==1)
						{
							//E->E STDP				
							if(type[i]==0 && type[j]==0 && layer[j]==layer[i] && layer[i]==0)
							{
									need_modify_link_i.push_back(i);
									need_modify_link_j.push_back(j);
									link_modification.push_back((learning_rate)*neuron_trace_slow[i]*neuron_trace[j]*exp(-t_spike*(1/tau_STDP_slow+1/tau_STDP))-beta*(link[i][j]-link_target[i][j])*(link[i][j]-link_target[i][j])*(link[i][j]-link_target[i][j])/link_target[i][j]/link_target[i][j]/link_target[i][j]*neuron_trace[i]*neuron_trace[i]*neuron_trace[i]*exp(-(3*t_spike)/tau_STDP));//add hetero 0.5
								
							}
						}
						
						
					}
					
					V[i]=V_rest;	
					//STDP
					neuron_trace[i]+=exp(t_spike/tau_STDP);
					neuron_trace_slow[i]+=exp(t_spike/tau_STDP_slow);
					neuron_count_trace[i]+=(1/tau_count)*exp(t_spike/tau_count);
					
					//STP
					if(type[i]==0)
					{
						Cx_own[i]+=(U+temp_Cu_own)*(1-temp_Cx_own);
						Cu_own[i]+=U*(1-U-temp_Cu_own);

					}
					fprintf(fp2,"%f\t%d\n",delta_t*(t+t_spike),i);	
					
					t_refractory_end[i]=t+refractory_period+t_spike;
					
					
				}
				else
				{
					if(V[i]+temp<-70)
						V[i]=-70; 
					else
						V[i]+=temp;
 
					
				}
			}
			//refractary time correction v(t) to v(t+dt)
			if(t_refractory_end[i]-t<=dt && t_refractory_end[i]-t>0) 
			{
				//RK2 V_update
				temp_S_E_C_d=sum_layer_S_E_C_d*exp(-(t_refractory_end[i]-t)/tau_d_E);
				temp_S_E_C_r=sum_layer_S_E_C_r*exp(-(t_refractory_end[i]-t)/tau_r);
				temp_S_E_C_d_NMDA=sum_layer_S_E_C_d_NMDA*exp(-(t_refractory_end[i]-t)/tau_d_E_NMDA);
				temp_S_E_C_r_NMDA=sum_layer_S_E_C_r_NMDA*exp(-(t_refractory_end[i]-t)/tau_r_NMDA);
				temp_S_I_C_d=sum_layer_S_I_C_d*exp(-(t_refractory_end[i]-t)/tau_d_I);
				temp_S_I_C_r=sum_layer_S_I_C_r*exp(-(t_refractory_end[i]-t)/tau_r);
				temp_S_O_C_d=S_O_C_d[i]*exp(-(t_refractory_end[i]-t)/tau_d_E);
				temp_S_O_C_r=S_O_C_r[i]*exp(-(t_refractory_end[i]-t)/tau_r);
				

				if(type[i]==0)
				{
					G_kE_delta_t=tau*(current(temp_S_O_C_d,temp_S_O_C_r,tau_d_E,tau_r)+alpha_AMPA*current(temp_S_E_C_d,temp_S_E_C_r,tau_d_E,tau_r)+alpha_NMDA*current(temp_S_E_C_d_NMDA,temp_S_E_C_r_NMDA,tau_d_E_NMDA,tau_r_NMDA));
				}
					
				else
					G_kE_delta_t=tau*(current(temp_S_O_C_d,temp_S_O_C_r,tau_d_E,tau_r)+current(temp_S_E_C_d,temp_S_E_C_r,tau_d_E,tau_r));
				
				G_kI_delta_t=tau*(current(temp_S_I_C_d,temp_S_I_C_r,tau_d_I,tau_r));
				
				k1=dV_dt(V[i],G_kE_delta_t,G_kI_delta_t,tau, tau_r,V_L, E_E,E_I);
				
				temp_S_E_C_d=sum_layer_S_E_C_d*update_constant_d_E;
				temp_S_E_C_r=sum_layer_S_E_C_r*update_constant_r;
				temp_S_E_C_d_NMDA=sum_layer_S_E_C_d_NMDA*update_constant_d_E_NMDA;
				temp_S_E_C_r_NMDA=sum_layer_S_E_C_r_NMDA*update_constant_r;
				temp_S_I_C_d=sum_layer_S_I_C_d*update_constant_d_I;
				temp_S_I_C_r=sum_layer_S_I_C_r*update_constant_r;
				temp_S_O_C_d=S_O_C_d[i]*update_constant_d_E;
				temp_S_O_C_r=S_O_C_r[i]*update_constant_r;
				
				if(type[i]==0)
				{
					G_kE_delta_t=tau*(current(temp_S_O_C_d,temp_S_O_C_r,tau_d_E,tau_r)+alpha_AMPA*current(temp_S_E_C_d,temp_S_E_C_r,tau_d_E,tau_r)+alpha_NMDA*current(temp_S_E_C_d_NMDA,temp_S_E_C_r_NMDA,tau_d_E_NMDA,tau_r_NMDA));
				}	
				else
					G_kE_delta_t=tau*(current(temp_S_O_C_d,temp_S_O_C_r,tau_d_E,tau_r)+current(temp_S_E_C_d,temp_S_E_C_r,tau_d_E,tau_r));
			
				G_kI_delta_t=tau*(current(temp_S_I_C_d,temp_S_I_C_r,tau_d_I,tau_r));
				k2=dV_dt(V[i]+k1*(dt-t_refractory_end[i]+t),G_kE_delta_t,G_kI_delta_t,tau, tau_r,V_L, E_E,E_I);
				
				if(k1/k2<0 && k2!=0)
					temp=k1*(dt-t_refractory_end[i]+t);
				else
					temp=((k1+k2)/2)*(dt-t_refractory_end[i]+t);
				
				
				if(theta_counter%(int)(floor(1000/gamma/delta_t))>0 && theta_counter%(int)(floor(1000/gamma/delta_t))<=(int)(floor(1000/theta/2/delta_t)) && (int)(theta_counter/(floor(1000/gamma/delta_t)))<pulse_number && i>=inh2-FFI_size && (memory_active==1 || memory_active==2) && t>=t_refractory_end[i] && randnum.GenRandReal_10()<pulse_active)
				{
					
					V[i]=-70;
				}
				
				
				
				if(V[i]+temp>=V_th)
				{
					
					t_spike=(V_th-V[i])/temp*(dt-t_refractory_end[i]+t);
					
					
					
					
					if(type[i]==0)
					{
						temp_Cu_own=Cu_own[i]*exp(-(t_refractory_end[i]-t+t_spike)/tau_F);
						temp_Cx_own=Cx_own[i]*exp(-(t_refractory_end[i]-t+t_spike)/tau_D);
					}
					
					for(j=0;j<Neuron_no;j++)
					{
						if(link[j][i]!=0)
						{
							spikefrom[j].push_back(i);
							current_arrival_time[j].push_back(t_refractory_end[i]+t_spike+tau_l);
							
							if(type[i]==0)
							{
								Cu_arrival[j].push_back(temp_Cu_own);
								Cx_arrival[j].push_back(temp_Cx_own);
							}
						}
						//STDP
						if(link[i][j]!=0  && t>1000/delta_t && STDP_activation==1)
						{
							//E->E STDP
							if(type[i]==0 && type[j]==0 && layer[j]==layer[i] && layer[i]==0)
							{
								need_modify_link_i.push_back(i);
								need_modify_link_j.push_back(j);
								link_modification.push_back((learning_rate)*neuron_trace_slow[i]*neuron_trace[j]*exp(-(t_refractory_end[i]-t+t_spike)*(1/tau_STDP_slow+1/tau_STDP))-beta*(link[i][j]-link_target[i][j])*(link[i][j]-link_target[i][j])*(link[i][j]-link_target[i][j])/link_target[i][j]/link_target[i][j]/link_target[i][j]*neuron_trace[i]*neuron_trace[i]*neuron_trace[i]*exp(-3*(t_refractory_end[i]-t+t_spike)/tau_STDP));//add hetor 0.4 
								
							}
							
							
						}
					}
					
					V[i]=V_rest;
					//STDP
					neuron_trace[i]+=exp((t_refractory_end[i]-t+t_spike)/tau_STDP);
					neuron_trace_slow[i]+=exp((t_refractory_end[i]-t+t_spike)/tau_STDP_slow);
					neuron_count_trace[i]+=(1/tau_count)*exp((t_refractory_end[i]-t+t_spike)/tau_count);
					
					//STP
					if(type[i]==0)
					{
						Cx_own[i]+=(U+temp_Cu_own)*(1-temp_Cx_own);
						Cu_own[i]+=U*(1-U-temp_Cu_own);
					}
					fprintf(fp2,"%f\t%d\n",delta_t*(t_refractory_end[i]+t_spike),i);	
					
					t_refractory_end[i]+=refractory_period+t_spike;
					
				}
				else
				{
					if(V[i]+temp<-70)
						V[i]=-70; 
					else
						V[i]+=temp;
					
				}
			}	
		
			
			
			
			
			//update s(t-dt) to s(t)
			for(j=0;j<3;j++)
			{
				S_E_C_d[i][j]*=update_constant_d_E;
				S_E_C_r[i][j]*=update_constant_r;
				S_E_C_d_NMDA[i][j]*=update_constant_d_E_NMDA;
				S_E_C_r_NMDA[i][j]*=update_constant_r_NMDA;
				S_I_C_d[i][j]*=update_constant_d_I;
				S_I_C_r[i][j]*=update_constant_r;
			}
			S_O_C_d[i]*=update_constant_d_E;
			S_O_C_r[i]*=update_constant_r;
			
	
			if(type[i]==0)
			{
				Cu_own[i]*=update_constant_Cu;
				Cx_own[i]*=update_constant_Cx;
			}
			
			//STDP
			neuron_trace[i]*=update_constant_neuron_trace;
			neuron_trace_slow[i]*=update_constant_neuron_trace_slow;
			neuron_count_trace[i]*=update_constant_neuron_count_trace;
			
			
			//check current arrival between t and t+dt	
			if(current_arrival_time[i].size()!=0)
			{	
				while(current_arrival_time[i][0]-t-dt<=0 && current_arrival_time[i].size() != 0)
				{
					
					if(type[spikefrom[i][0]]==0)// from E
					{
						if(type[i]==0)
						{
							S_E_C_d[i][layer[spikefrom[i][0]]]+=link[i][spikefrom[i][0]]*(U+Cu_arrival[i][0])*(1-Cx_arrival[i][0])*exp(-(dt-current_arrival_time[i][0]+t)/tau_d_E);
							S_E_C_r[i][layer[spikefrom[i][0]]]+=link[i][spikefrom[i][0]]*(U+Cu_arrival[i][0])*(1-Cx_arrival[i][0])*exp(-(dt-current_arrival_time[i][0]+t)/tau_r);
							S_E_C_d_NMDA[i][layer[spikefrom[i][0]]]+=link[i][spikefrom[i][0]]*(U+Cu_arrival[i][0])*(1-Cx_arrival[i][0])*exp(-(dt-current_arrival_time[i][0]+t)/tau_d_E_NMDA);
							S_E_C_r_NMDA[i][layer[spikefrom[i][0]]]+=link[i][spikefrom[i][0]]*(U+Cu_arrival[i][0])*(1-Cx_arrival[i][0])*exp(-(dt-current_arrival_time[i][0]+t)/tau_r_NMDA);
							
							
							//STDP
							//E->E
							if(link[i][spikefrom[i][0]]!=0  && t>1000/delta_t && layer[spikefrom[i][0]]==layer[i] && layer[i]==0 && STDP_activation==1)
							{
									need_modify_link_i.push_back(i);
									need_modify_link_j.push_back(spikefrom[i][0]);
									link_modification.push_back(-(learning_rate)*(neuron_trace[i]*exp(-(dt-current_arrival_time[i][0]+t-tau_l)/tau_STDP))+delta);
								
							
							}
						}
						else
						{
							S_E_C_d[i][layer[spikefrom[i][0]]]+=link[i][spikefrom[i][0]]*(U+Cu_arrival[i][0])*(1-Cx_arrival[i][0])*exp(-(dt-current_arrival_time[i][0]+t)/tau_d_E);
							S_E_C_r[i][layer[spikefrom[i][0]]]+=link[i][spikefrom[i][0]]*(U+Cu_arrival[i][0])*(1-Cx_arrival[i][0])*exp(-(dt-current_arrival_time[i][0]+t)/tau_r);
							S_E_C_d_NMDA[i][layer[spikefrom[i][0]]]+=link[i][spikefrom[i][0]]*(U+Cu_arrival[i][0])*(1-Cx_arrival[i][0])*exp(-(dt-current_arrival_time[i][0]+t)/tau_d_E_NMDA);
							S_E_C_r_NMDA[i][layer[spikefrom[i][0]]]+=link[i][spikefrom[i][0]]*(U+Cu_arrival[i][0])*(1-Cx_arrival[i][0])*exp(-(dt-current_arrival_time[i][0]+t)/tau_r_NMDA);
						}
						Cu_arrival[i].pop_front();
						Cx_arrival[i].pop_front();
						
						
						
					}
					if(type[spikefrom[i][0]]==1)// from I
					{
					
						S_I_C_d[i][layer[spikefrom[i][0]]]+=link[i][spikefrom[i][0]]*exp(-(dt-current_arrival_time[i][0]+t)/tau_d_I);
						S_I_C_r[i][layer[spikefrom[i][0]]]+=link[i][spikefrom[i][0]]*exp(-(dt-current_arrival_time[i][0]+t)/tau_r);
					}
					current_arrival_time[i].pop_front();
					spikefrom[i].pop_front();
				}
			}
			
				
			while(external_input_time[i]-t-dt<=0)
			{
				if(layer[i]==0 || layer[i]==1)
				{
					if(memory_active==1 || memory_active==2 || memory_active==3)
					{
						if(type[i]==1)
						{
							S_O_C_d[i]+=g_kO*exp(-(dt-external_input_time[i]+t)/tau_d_E);
							S_O_C_r[i]+=g_kO*exp(-(dt-external_input_time[i]+t)/tau_r);
							
						}
						else
						{
							S_O_C_d[i]+=g_kO*exp(-(dt-external_input_time[i]+t)/tau_d_E);
							S_O_C_r[i]+=g_kO*exp(-(dt-external_input_time[i]+t)/tau_r);
						}
						if(memory_active==2 && type[i]==0)
						{
							external_input_time[i]+=(-1*log(1-randnum.GenRandReal_10())/(f_ex)/(400));
						}
						else
						{
							if(memory_cue==1 && i<(memory_cue_no+1)*memory && i>=memory_cue_no*memory  && type[i]==0)
							{
								external_input_time[i]+=(-1*log(1-randnum.GenRandReal_10())/(5*f_ex)/(400));
							}
							else if(memory_cue==2 && i<(memory_cue_no+1)*memory && i>=memory_cue_no*memory  && type[i]==0)
							{
								external_input_time[i]+=(-1*log(1-randnum.GenRandReal_10())/(3.25*f_ex)/(400));
							}
							else
							{
								external_input_time[i]+=(-1*log(1-randnum.GenRandReal_10())/(f_ex)/(400));
							}
						}
						
						
						
					}
					if(memory_active==0)
					{
						S_O_C_d[i]+=g_kO*exp(-(dt-external_input_time[i]+t)/tau_d_E);
						S_O_C_r[i]+=g_kO*exp(-(dt-external_input_time[i]+t)/tau_r);
						external_input_time[i]+=(-1*log(1-randnum.GenRandReal_10())/(f_ex)/(400));
					}
				}
				else
				{
					S_O_C_d[i]+=g_kO*exp(-(dt-external_input_time[i]+t)/tau_d_E);
					S_O_C_r[i]+=g_kO*exp(-(dt-external_input_time[i]+t)/tau_r);
					external_input_time[i]+=(-1*log(1-randnum.GenRandReal_10())/f_ex/(400));
				}
			}	
		}
		fclose(fp2);
		

		
		if(need_modify_link_i.size()!=0)
		{
			while(need_modify_link_i.size()!=0)
			{
				temp_link=link[need_modify_link_i[0]][need_modify_link_j[0]];
				max=10;
				link[need_modify_link_i[0]][need_modify_link_j[0]]+=link_modification[0];
				if(link[need_modify_link_i[0]][need_modify_link_j[0]]>max)
					link[need_modify_link_i[0]][need_modify_link_j[0]]=max;
				if(link[need_modify_link_i[0]][need_modify_link_j[0]]<0.001)
					link[need_modify_link_i[0]][need_modify_link_j[0]]=0.001;
				
				need_modify_link_i.pop_front();
				need_modify_link_j.pop_front();
				link_modification.pop_front();	
			}
		}

		
		
		theta_counter++;
		
	}



	sprintf(filename_5,"E=%.2fI=%.2fP=%.2f_W_%d_%.0f_%.2f.txt",x,y,link_density,count,t*delta_t,yy);
	fp5=fopen(filename_5,"w");
	for(i=0;i<Neuron_no;i++)
	{
		for(j=0;j<Neuron_no;j++)
		{
			fprintf(fp5,"%f\t",link[i][j]);
				
		}
		fprintf(fp5,"\n");
	}
	fclose(fp5);
}



 
