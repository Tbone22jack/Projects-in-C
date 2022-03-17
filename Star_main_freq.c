#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
// this allows me to append numbers as chars onto a string
char *
itoa (int value, char *result, int base)
{
    // check that the base if valid
    if (base < 2 || base > 36) { *result = '\0'; return result; }

    char* ptr = result, *ptr1 = result, tmp_char;
    int tmp_value;

    do {
        tmp_value = value;
        value /= base;
        *ptr++ = "zyxwvutsrqponmlkjihgfedcba9876543210123456789abcdefghijklmnopqrstuvwxyz" [35 + (tmp_value - value * base)];
    } while ( value );

    // Apply negative sign
    if (tmp_value < 0) *ptr++ = '-';
    *ptr-- = '\0';
    while (ptr1 < ptr) {
        tmp_char = *ptr;
        *ptr--= *ptr1;
        *ptr1++ = tmp_char;
    }
    return result;
}
//n in the functions
double N = 1.0;
double S_function(double S, double A, double M);
double M_function(double S, double A, double M);
double A_function(double S, double A, double M);
double s_function(double S, double A, double M);
double m_function(double S, double A, double M);
double a_function(double S, double A, double M);
void fourth_order_RungeKutta_1(double I_A, double I_M, double I_S,double delta_t);
void fourth_order_RungeKutta_2(double I_A, double I_M, double I_S,double delta_t,int iteration, FILE *d);
void Print_result(double S, double A, double M,FILE *p,double time);
//These variables express interactions between the masses in the system
//K1 is the constant replenishing of the atomic gas
//K2 is the loss of stellar mass to the solar wind, also replenishes atomic cloud
//K3 is the transformation of the atomic cloud into molecular clouds
//K4 is the star formation rate
double K1, K2, K3, K4;

//Second set of constants for second set of equations
double k1,k2;
//Max number of time steps
#define MAX 3000

int main() {

	//These variables represent the masses within the system.
	//S_active being the active star mass.
	//M_cloud being the mass of the molecular gas clouds
	//A_cloud being the mass of the atomic gas clouds
	double S_active=.4, M_cloud=.4, A_cloud=.2;
	//This variable is the time step
	double delta_t=.001;
	//Constants to be set
	K1=1;K2=.1;K3=4;K4=5;

	//fourth_order_RungeKutta_1(A_cloud, M_cloud, S_active, delta_t);
	//Constants to be set using previous constants
	k1= (K3)/(K1+K2);k2=K4/(K1+K2);
	//Override on constants
	k1=8,k2=15;
	//Reset for next set
	//This would be sO,mO,aO in the paper
	S_active=.3, M_cloud=.3, A_cloud=.4;
	delta_t=.0001;
	FILE *fp = fopen("Frequency_and_Amplitude_vs_n.csv","w");
	fprintf(fp,"Frequency,Amplitude, n \n");
	for(int iteration=0;iteration<80;iteration++)
	fourth_order_RungeKutta_2(A_cloud, M_cloud, S_active, delta_t,iteration,fp);
	fclose(fp);
	return 1;
}
//This is the method I will be using to approximate the system with 4th order delta t accuracy
void fourth_order_RungeKutta_1(double I_A, double I_M, double I_S,double delta_t) {
	double S_active = I_S, M_cloud = I_M, A_cloud = I_A, SS_active, SM_cloud,
			SA_cloud, SSS_active, SSM_cloud, SSA_cloud, SSSS_active, SSSM_cloud,
			SSSA_cloud;
FILE *f =fopen("STAR_RESULT.csv","w");
fprintf(f,"S_active,A_cloud,M_cloud,time\n");
for(double loop=delta_t;loop<MAX;loop+=delta_t)
{
	SA_cloud = A_cloud + (delta_t / 2) * A_function(A_cloud, M_cloud,S_active);
	SM_cloud = M_cloud + (delta_t / 2) * M_function(A_cloud, M_cloud,S_active);
	SS_active = S_active + (delta_t / 2) * S_function(A_cloud, M_cloud,S_active);
	SSA_cloud = A_cloud + (delta_t / 2) * A_function(SA_cloud, SM_cloud,SS_active);
	SSM_cloud = M_cloud + (delta_t / 2) * M_function(SA_cloud, SM_cloud,SS_active);
	SSS_active = S_active + (delta_t / 2) * S_function(SA_cloud, SM_cloud, SS_active);
	SSSA_cloud = A_cloud + delta_t*A_function(SSA_cloud, SSM_cloud, SSS_active);
	SSSM_cloud = M_cloud + delta_t*M_function(SSA_cloud, SSM_cloud,SSS_active);
	SSSS_active = S_active + delta_t*S_function(SSA_cloud, SSM_cloud,SSS_active);
	A_cloud=A_cloud+((delta_t/6)*(A_function(A_cloud, M_cloud,S_active)+ 2*A_function(SA_cloud, SM_cloud,SS_active)+ 2*A_function(SSA_cloud, SSM_cloud,SSS_active)+A_function(SSSA_cloud, SSSM_cloud,SSSS_active)));
	M_cloud=M_cloud+((delta_t/6)*(M_function(A_cloud, M_cloud,S_active)+ 2*M_function(SA_cloud, SM_cloud, SS_active)+ 2*M_function(SSA_cloud, SSM_cloud,SSS_active)+M_function(SSSA_cloud, SSSM_cloud, SSSS_active)));
	S_active=S_active+((delta_t/6)*(S_function(A_cloud, M_cloud,S_active)+ 2*S_function(SA_cloud, SM_cloud, SS_active)+ 2*S_function(SSA_cloud, SSM_cloud,SSS_active)+S_function(SSSA_cloud, SSSM_cloud, SSSS_active)));

	Print_result(S_active,A_cloud,M_cloud,f,N);
}

fclose(f);
}
void fourth_order_RungeKutta_2(double I_A, double I_M, double I_S,double delta_t,int num,FILE *d) {
	double S_active = I_S, M_cloud = I_M, A_cloud = I_A, SS_active, SM_cloud,
			SA_cloud, AP_cloud, SSM_cloud, SSA_cloud, SSSS_active, SSSM_cloud,
			SSSA_cloud;

	char name[40];
	char number [50];
	int nk = num;
	itoa(nk, number, 10);
	strcpy(name,"Val_n/Trial_2/Star_Result_Rep_");
	strncat(name,number,40);
	strncat(name,".csv",40);
	//remove(name);
	//FILE *f =fopen(name,"w");
	//fprintf(f,"S_active,A_cloud,M_cloud,time\n");

int isboring=0,bool=1;
double first=0,last;
double freq,amplitude;

double max,min;
printf("\n%lf\n",N);
for(double loop=delta_t;loop<MAX;loop+=delta_t)
{
	SA_cloud = A_cloud + (delta_t / 2) * a_function(A_cloud, M_cloud,S_active);
	SM_cloud = M_cloud + (delta_t / 2) * m_function(A_cloud, M_cloud,S_active);
	SSA_cloud = A_cloud + (delta_t / 2) * a_function(SA_cloud, SM_cloud,S_active);
	SSM_cloud = M_cloud + (delta_t / 2) * m_function(SA_cloud, SM_cloud,S_active);
	SSSA_cloud = A_cloud + delta_t*a_function(SSA_cloud, SSM_cloud, S_active);
	SSSM_cloud = M_cloud + delta_t*m_function(SSA_cloud, SSM_cloud,S_active);
	AP_cloud=A_cloud+((delta_t/6)*(a_function(A_cloud, M_cloud,S_active)+ 2*a_function(SA_cloud, SM_cloud,S_active)+ 2*a_function(SSA_cloud, SSM_cloud,S_active)+a_function(SSSA_cloud, SSSM_cloud,S_active)));
	M_cloud=M_cloud+((delta_t/6)*(m_function(A_cloud, M_cloud,S_active)+ 2*m_function(SA_cloud, SM_cloud, S_active)+ 2*m_function(SSA_cloud, SSM_cloud,S_active)+m_function(SSSA_cloud, SSSM_cloud, S_active)));
	if(fabs((AP_cloud-A_cloud))/A_cloud < .00001)
		isboring++;
	else
		isboring=0;
	if( isboring>3000000 || A_cloud != A_cloud)
	{
		remove(name);
		bool=3;
		break;
	}
	if(S_active>1-AP_cloud-M_cloud && bool== 0)
	{
					max=1-AP_cloud-M_cloud;
					last=loop;
					bool=1;
					printf("%lf %lf ",first, last);
					if(first!=0)
						freq=1/(last-first);
					first=loop;
					amplitude=(max-min)/2;
					printf("%lf ",freq);
	}

	if(S_active<1-AP_cloud-M_cloud && bool==1)
	{
		min=1-AP_cloud-M_cloud;
		bool=0;
	}
	A_cloud=A_cloud+((delta_t/6)*(a_function(A_cloud, M_cloud,S_active)+ 2*a_function(SA_cloud, SM_cloud,S_active)+ 2*a_function(SSA_cloud, SSM_cloud,S_active)+a_function(SSSA_cloud, SSSM_cloud,S_active)));
	S_active=1-A_cloud-M_cloud;
	//Print_result(S_active,A_cloud,M_cloud,f,loop);
}
if(bool!=3){

fprintf(d,"%lf, %lf, %lf\n",freq,amplitude,N);
}
//fclose(f);
N+=0.01;
}
//Functions that describe the formation behavior
double A_function(double A, double M, double S) {
	return K1 * S + K2 * S - K3 * A * (M * M);
}
double M_function(double A, double M, double S) {
	return K3 * A * (M * M) - K4 * S * pow(M,N);
}
double S_function(double A, double M, double S) {
	return K4 * S * pow(M,N)- K1 * S - K2 * S;
}

//Functions for the second set of equations
double a_function(double A, double M, double S) {
	return 1 - A - M - k1*(M*M)*A;
}
double m_function(double A, double M, double S) {
	return k1*(M*M)*A + (k2*(pow(M,N))) * (A- 1 + M);
}
double s_function(double A, double M, double S) {
	return 1-A-M;
}
void Print_result(double S, double A, double M,FILE *p,double time)
{
	fprintf(p,"%lf,%lf,%lf,%lf\n",S,A,M,time);
	//printf("%lf,%lf,%lf,%lf\n",S,A,M,time);
}
