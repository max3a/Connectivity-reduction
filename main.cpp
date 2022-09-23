#include "run.h"


int main(int argc,char*argv[])
{
	int x,y;
	double xx,yy;
	double link_density;
	double tau_d_e,tau_d_i;
	double gamma,theta,pulse_number;	
	int count;
	int learn;
	char str[100];

	//trial number
	count=80000;
	//allow dual oscillation
	theta=10;//slow oscillation frequency
	gamma=10;//fast oscillation frequency
	
	
	do{
		printf("Rescue (Y\\N)?");
		gets(str);
	}while(strcmp (str,"Y")!=0 && strcmp (str,"N")!=0);
	if(strcmp (str,"Y")==0)
		pulse_number=1;//0: no rescue, 1:rescue
	else
		pulse_number=0;
		
	do{
		printf("Connectivity (0 to 1)?");
		gets(str);
	}while(atof(str)<0 || atof(str)>1);
	yy=atof(str);
	
	do{
		printf("learn (Y\\N)?");
		gets(str);
	}while(strcmp (str,"Y")!=0 && strcmp (str,"N")!=0);
	if(strcmp (str,"Y")==0)
		learn=1;//0: no learn, 1:with learn
	else
		learn=0;

	tau_d_e=3;
	tau_d_i=8;
	evolution(tau_d_e,tau_d_i,count,time(NULL),yy,theta,gamma,pulse_number,learn);

}
 
