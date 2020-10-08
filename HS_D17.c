/*CHANGE THE FUNCTION "read_input" if you change the dimensions D. Here Rmin, rd correspond to the minimum and maximum diameter of the spheres respectively.*/
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include <time.h>
#include "ran2.c"

#define eps 0.0000000001
#define D 17
#define Ravg 1.0
#define f 0.008
#define dsigmatol 0.005
#define tr (f*Ravg)
#define swapeq 0.00
#define TD 1000
#define swapprob 0.00
#define IFLAG 1
#define pi (22.0/7.0)

int N,Teq,nsteps,GAP,INIT,MNNEI,tstart,stepstart;
double boxsize,cellsize,boxl,boxl2,vf,drneimax,drneimax2,Rskin,Rskin_in;
char outfile4[200],readfile[200],outfile5[200],outfile6[200];
long int seed=485620;
double *RD,*dispsum[D],dCOM[D];
int **NNList,*Ncross[D];

FILE *fpw;

int find_next_min_image(double x)
{
        int d,fxp,gxp;
        double xp;
        xp=x/boxl2;
        fxp=round(xp);
        if(xp < 0)
                d=-1.0;
        else
                d=1.0;


        if(fabs(fxp) >= fabs(xp))
                gxp=fxp-d;
        else
                gxp=fxp+d;
        if(fxp == 0)
        {
                if(xp < 0)
                        gxp=-1;
                else
                        gxp=1;
        }

        return gxp;
}

int find_max_dr(double dr[D])
{
        int k,ind;
        double max=0;
        for(k=0;k<D;k++)
        {
                if(dr[k] > max)
                {
                        max=dr[k];
                        ind=k;
                }
        }
        return ind;
}

int * image_operation(double array[D])
{
        int index,k,sum,tmprr;
        static int tmprp[D];
        double dtmpr[D];
        sum=0;
        for(k=0;k<D;k++)
        {
                tmprp[k]=round(array[k]/boxl2);
                sum=sum+tmprp[k];
        }
        sum=fabs(sum);
        if(sum % 2 == 1)
        {
                for(k=0;k<D;k++)
                        dtmpr[k]=fabs(array[k]/boxl2-tmprp[k]);
                index=find_max_dr(dtmpr);
                tmprr=find_next_min_image(array[index]);
                tmprp[index]=tmprr;
        }
        return tmprp;
}

void read_input(double r[D][N], double VF)
{
	FILE *fpr;
	int i,k,time;
	time=INIT;
	int c1,c2,c3,c4,c5,c6,c7,c8,c8_,c88,c888,c9_,c10_,c11_,c12_,c13_,c14_;
	double c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,c21,c22,c23,c24,c25;
	i=0;

	FILE *file;
        char outfile7[200];
        sprintf(outfile7,"restart_file/restartfile%dD_N%dBS%1.6lfvf%1.6lfRskin%1.2lfS%1.2lf.txt",D,N,boxsize,vf,Rskin,swapprob);
        if((file = fopen(outfile7,"r"))!=NULL)
        {
                fscanf(file,"%d\n",&tstart);
                fscanf(file,"%d\n",&stepstart);
		fscanf(file,"%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf\n",&dCOM[0],&dCOM[1],&dCOM[2],&dCOM[3],&dCOM[4],&dCOM[5],&dCOM[6],&dCOM[7],&dCOM[8],&dCOM[9],&dCOM[10],&dCOM[11],&dCOM[12],&dCOM[13],&dCOM[14],&dCOM[15],&dCOM[16]);
                for(i=0;i<N;i++)
                        fscanf(file,"%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%lf\n",&r[0][i],&r[1][i],&r[2][i],&r[3][i],&r[4][i],&r[5][i],&r[6][i],&r[7][i],&r[8][i],&r[9][i],&r[10][i],&r[11][i],&r[12][i],&r[13][i],&r[14][i],&r[15][i],&r[16][i],&Ncross[0][i],&Ncross[1][i],&Ncross[2][i],&Ncross[3][i],&Ncross[4][i],&Ncross[5][i],&Ncross[6][i],&Ncross[7][i],&Ncross[8][i],&Ncross[9][i],&Ncross[10][i],&Ncross[11][i],&Ncross[12][i],&Ncross[13][i],&Ncross[14][i],&Ncross[15][i],&Ncross[16][i],&RD[i]);
                fclose(file);
                return;
        }

	if(IFLAG == 1)
	{
		sprintf(readfile,"inflate_out_mono_/configs/fconfig%dD_N%d_BS%1.6lfvf%1.6lf.dat",D,N,boxsize,VF);
		fpr=fopen(readfile,"r");
		while(fscanf(fpr,"%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",&r[0][i],&r[1][i],&r[2][i],&r[3][i],&r[4][i],&r[5][i],&r[6][i],&r[7][i],&r[8][i],&r[9][i],&r[10][i],&r[11][i],&r[12][i],&r[13][i],&r[14][i],&r[15][i],&r[16][i],&RD[i])!=EOF)
			i++;
	}
	else if(IFLAG == 2)
	{
                sprintf(readfile,"./rho_%1.6lfSWP_0.20eq_D17/config%dD_N%d_BS%1.6lf_vf%1.6lfRskin_%1.2lft_%d.dat",vf,D,N,boxsize,vf,Rskin_in,time);
		fpr=fopen(readfile,"r");
		while(fscanf(fpr,"%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",&r[0][i],&r[1][i],&r[2][i],&r[3][i],&r[4][i],&r[5][i],&r[6][i],&r[7][i],&r[8][i],&r[9][i],&r[10][i],&r[11][i],&r[12][i],&r[13][i],&r[14][i],&r[15][i],&r[16][i],&c1,&c2,&c3,&c4,&c5,&c6,&c7,&c8,&c8_,&c88,&c888,&c9_,&c10_,&c11_,&c12_,&c13_,&c14_,&c9,&c10,&c11,&c12,&c13,&c14,&c15,&c16,&c17,&c18,&c19,&c20,&c21,&c22,&c23,&c24,&c25,&RD[i])!=EOF)
			i++;
	}
	fclose(fpr);
}

void next(int v[D],int m)
{
        int ip;
        ip=0;
        while(v[ip] < m)
        {
                v[ip]=v[ip]+1;
                if(v[ip] < m)
                        break;
                else if(v[ip] > (m-1))
                {
                        v[ip]=0;
                        ip=ip+1;
                }
        }
}

void create_neighbour_list(double r[D][N])
{
        int k,j,i,nnei,nnlistbeg;
        double rdiff[D],rdiffsq,rdiffsqrt,sigma;

        int *clnn;
        for(i=0;i<N;i++)
        {
                nnei=0;
                sigma=RD[i];
                for(j=0;j<N;j++)
                {
                        if(j != i)
                        {
                                rdiffsq=0.0;
                                for(k=0;k<D;k++)
                                        rdiff[k]=r[k][i]-r[k][j];
                                clnn = image_operation(rdiff);
                                for(k=0;k<D;k++)
                                {
                                        rdiff[k]=rdiff[k]-*(clnn+k)*boxl2;
                                        rdiffsq=rdiffsq+rdiff[k]*rdiff[k];
                                }
                                rdiffsqrt=sqrt(rdiffsq);
                                if(rdiffsqrt < cellsize)
                                {
                                        nnei++;
                                        NNList[i][nnei]=j;

                                }
                        }
                }
                NNList[i][0]=nnei;
        }
}

int ipow(int base, int exp)
{
	int result = 1;
	while (exp)
	{
		if (exp & 1)
			result *= base;
		exp >>= 1;
		base *= base;
	}
	return result;
}

int check_overlap(int j, double posi[D], double r[D][N], double sigma)
{
        int k,tmpj,nb;
        double rdiff[D],rdiffsq,rdiffsqrt;
        int *clnn;
        for(nb=1;nb<=NNList[j][0];nb++)
        {
                tmpj=NNList[j][nb];
                rdiffsq=0.0;
                for(k=0;k<D;k++)
                        rdiff[k]=posi[k]-r[k][tmpj];
                clnn = image_operation(rdiff);
                for(k=0;k<D;k++)
                {
                        rdiff[k]=rdiff[k]-*(clnn+k)*boxl2;
                        rdiffsq=rdiffsq+rdiff[k]*rdiff[k];
                }
                rdiffsqrt=sqrt(rdiffsq);
                if(rdiffsqrt <((sigma+RD[tmpj])/2.0+eps))
                        return 1;
        }
        return 0;
}

int check_overlap_swap(int j, int j2, double posi[D], double r[D][N], double sigma)
{
        int k,tmpj,nb;
        double rdiff[D],rdiffsq,rdiffsqrt;
        int *clnn;
        for(nb=1;nb<=NNList[j][0];nb++)
        {
                tmpj=NNList[j][nb];
                if(tmpj != j && tmpj != j2)
                {
                        rdiffsq=0.0;
                        for(k=0;k<D;k++)
                                rdiff[k]=posi[k]-r[k][tmpj];
                        clnn = image_operation(rdiff);
                        for(k=0;k<D;k++)
                        {
                                rdiff[k]=rdiff[k]-*(clnn+k)*boxl2;
                                rdiffsq=rdiffsq+rdiff[k]*rdiff[k];
                        }
                        rdiffsqrt=sqrt(rdiffsq);
                        if(rdiffsqrt <((sigma+RD[tmpj])/2.0+eps))
                                return 1;
                }
        }
        return 0;
}


int check_overlap_N(int j, double posi[D], double r[D][N], double sigma)
{
        int k,nb;
        double rdiff[D],rdiffsq,rdiffsqrt;
        int *clnn;
        for(nb=0;nb<N;nb++)
        {
                if(nb != j)
                {
                        rdiffsq=0.0;
                        for(k=0;k<D;k++)
                                rdiff[k]=posi[k]-r[k][nb];
                        clnn = image_operation(rdiff);
                        for(k=0;k<D;k++)
                        {
                                rdiff[k]=rdiff[k]-*(clnn+k)*boxl2;
                                rdiffsq=rdiffsq+rdiff[k]*rdiff[k];
                        }
                        rdiffsqrt=sqrt(rdiffsq);
                        if(rdiffsqrt <((sigma+RD[nb])/2.0+eps))
                                return 1;
                }
        }
        return 0;
}

void find_max_displacement(double dr)
{
	if(dr > drneimax)
	{
		drneimax2=drneimax;
		drneimax=dr;
	}
	else
		if(dr > drneimax2)
			drneimax2=dr;
}

void clear_displacement()
{
	int i,k;
	for(i=0;i<N;i++)
		for(k=0;k<D;k++)
			dispsum[k][i]=0.0;
	drneimax=0.0; drneimax2=0.0;
} 

void restore_dis(int i, double disp[D], double x1, double x2)
{
	int k;
	for(k=0;k<D;k++)
		dispsum[k][i]=dispsum[k][i]-disp[k];
	drneimax2=x2; drneimax=x1;
}

void translate(double r[D][N])
{
        int k,param,param1,j,i,sumclnn,cltmp[D];
        double disp[D],tmpr[D],input_r,sigma,drneisq,drnei,tmp1,tmp2,ur,true_r,cross[D],sum_cross;
        j=(int)(ran2(&seed)*N);

        int *clnn;

        for(k=0;k<D;k++)
        {
                disp[k]=2.0*tr*ran2(&seed)-tr;
                input_r=r[k][j]+disp[k];
                cross[k]=input_r;
        }
        clnn = image_operation(cross);
        for(k=0;k<D;k++)
        {
                tmpr[k]=cross[k]-*(clnn+k)*boxl2;
                cltmp[k]=*(clnn+k);
        }
        sigma=RD[j];
        param=check_overlap(j,tmpr,r,sigma);
        if(param != 0)
                return;
        drneisq=0.0;
        for(k=0;k<D;k++)
        {
                dispsum[k][j]=dispsum[k][j]+disp[k];
                drneisq=drneisq+pow(dispsum[k][j],2.0);
        }
        drnei=sqrt(drneisq);
        tmp1=drneimax; tmp2=drneimax2;
        find_max_displacement(drnei);
        if((drneimax+drneimax2) > (Rskin-eps))
        {
                param1=check_overlap_N(j,tmpr,r,sigma);
                if(param1 == 0)
                {
                        sumclnn=0;
                        for(k=0;k<D;k++)
                        {
                                r[k][j]=tmpr[k];
                                dCOM[k]=dCOM[k]+disp[k];
                                sumclnn=sumclnn+fabs(cltmp[k]);
			}

                        if((sumclnn) > 0)
                        {
                                for(k=0;k<D;k++)
                                        Ncross[k][j]=Ncross[k][j]+cltmp[k];
                        }
                        create_neighbour_list(r);
                        clear_displacement();
                }
                else
                        restore_dis(j,disp,tmp1,tmp2);
        }
        else
        {
                sumclnn=0;
                for(k=0;k<D;k++)
                {
                        r[k][j]=tmpr[k];
                        dCOM[k]=dCOM[k]+disp[k];
                        sumclnn=sumclnn+fabs(cltmp[k]);
                }
                if((sumclnn) > 0)
                        for(k=0;k<D;k++)
                                Ncross[k][j]=Ncross[k][j]+cltmp[k];
        }
}

void swap(double r[D][N])
{
        int i1,i2,param1,param2,k;
        double sigma1,sigma2,dsigma,tmpr1[D],tmpr2[D];

        dsigma=1.0;
        while(dsigma > dsigmatol)
        {
                i1=(int)(ran2(&seed)*N);
                i2=(int)(ran2(&seed)*N);
                if(i1 != i2)
                        dsigma=fabs(RD[i1]-RD[i2]);
        }

        for(k=0;k<D;k++)
        {
                tmpr1[k]=r[k][i1];
                tmpr2[k]=r[k][i2];
        }
        sigma1=RD[i1]; sigma2=RD[i2];
        param1=check_overlap_swap(i1,i2,tmpr1,r,sigma2);
        param2=check_overlap_swap(i2,i1,tmpr2,r,sigma1);
        if((param1+param2) == 0)
        {
                RD[i1]=sigma2; RD[i2]=sigma1;
        }
}

void evolve(double r[D][N], double pswp)
{
	int i;
	//printf("I AM HERE\n");
	for(i=0;i<N;i++)
	{
		if(ran2(&seed) < pswp)
			swap(r);
		else
			translate(r);
	}
}

double vol_measure(double diam)
{
        double vol,g,dd,ag;
        dd=1.0*D;
        ag=D/2.0+1.0;
        vol=pow(pi,dd/2.0)*pow(diam,dd)/tgamma(ag);
        return vol;
}

double measure_packing_frac()
{
        int i;
        double vol,est;
        vol=0.0;
        for(i=0;i<N;i++)
                vol=vol+1.0*vol_measure(RD[i]);
        vol=vol/pow(boxl2,1.0*D)/ipow(2,D)/2.0;
        printf("packing fraction = %lf\n",vol);
        return vol;
}

void print_configuration(double r[D][N], int time)
{
        int i,k;
        double input_r,tmpr[D],rr[D];
        sprintf(outfile4,"./rho_%1.6lfSWP_%1.2lfeq_D17/config%dD_N%d_BS%1.6lf_vf%1.6lfRskin_%1.2lft_%d.dat",vf,swapprob,D,N,boxsize,vf,Rskin,time);
        fpw=fopen(outfile4,"w");
        int *clnn;
        for(i=0;i<N;i++)
        {
                for(k=0;k<D;k++)
                {
                        input_r=r[k][i]-1.0*dCOM[k]/N;
                        rr[k]=input_r;
                }
                clnn = image_operation(rr);
                for(k=0;k<D;k++)
                {
                        tmpr[k]=rr[k]-*(clnn+k)*boxl2;
                        fprintf(fpw,"%.16e\t",tmpr[k]);
                }
                for(k=0;k<D;k++)
                        fprintf(fpw,"%d\t",Ncross[k][i]);
                for(k=0;k<D;k++)
                        fprintf(fpw,"%.16e\t",rr[k]+boxl2*Ncross[k][i]);
                fprintf(fpw,"%.16e\n",RD[i]);
        }
        fclose(fpw);
}

void print_restartfile(double r[D][N], int te, int ts)
{
        int i,k;
        FILE *fpwr;
        sprintf(outfile4,"restart_file/restartfile%dD_N%dBS%1.6lfvf%1.6lfRskin%1.2lfS%1.2lf.txt",D,N,boxsize,vf,Rskin,swapprob);
        fpwr=fopen(outfile4,"w");
        fprintf(fpwr,"%d\n",te);
        fprintf(fpwr,"%d\n",ts);
        for(k=0;k<D;k++)
                fprintf(fpwr,"%.16e\t",dCOM[k]);
        fprintf(fpwr,"\n");
        for(i=0;i<N;i++)
        {
                for(k=0;k<D;k++)
                        fprintf(fpwr,"%.16e\t",r[k][i]);
                for(k=0;k<D;k++)
                        fprintf(fpwr,"%d\t",Ncross[k][i]);
                fprintf(fpwr,"%.16e",RD[i]);
                fprintf(fpwr,"\n");
        }
        fclose(fpwr);
}

void wrap_coordinates(double r[D][N])
{
        int *clnn;
        int i,k;
        double tmp[D];
        for(i=0;i<N;i++)
        {
                for(k=0;k<D;k++)
                        tmp[k]=r[k][i];
                clnn = image_operation(tmp);
                for(k=0;k<D;k++)
                        r[k][i]=tmp[k]-*(clnn+k)*boxl2;
        }
}

int main(void)
{
	FILE *fp,*fpt;
	int i,it,gapres;

	scanf("%d",&N);
	scanf("%lf",&boxsize);
	scanf("%d",&nsteps);
	scanf("%d",&Teq);
	scanf("%d",&INIT);
	scanf("%lf",&vf);
	scanf("%lf",&Rskin);

	Rskin_in=Rskin;

//	N=4000; boxsize=2.6590858290740664; nsteps=0;
//	Teq=50; vf=0.000325675727635; INIT=0; Rskin_in=0.00; 
//	Rskin=0.28;
//	printf("Rskin=%lf\n",Rskin);

	boxl=boxsize; boxl2=((boxl)/2.0); GAP=(nsteps/TD);
	cellsize=(Ravg+Rskin);
	MNNEI=N;
	//gapres=nsteps/TD;
	gapres=1000;

	RD=(double *) malloc (N * sizeof(double));
	NNList=(int **) malloc (N * sizeof(int *));
	for(it=0;it<D;it++)
	{
		dispsum[it]=(double *) malloc (N * sizeof(double));
		Ncross[it]=(int *) malloc (N * sizeof(int));
	}
	for(it=0;it<N;it++)
		NNList[it]=(int *) malloc (MNNEI * sizeof(int));
	double r[D][N];
	sprintf(outfile5,"time_tracker_%dD_N%dBS%1.6lf_vf%1.6lfSWP%1.2lfRskin%1.2lf.dat",D,N,boxsize,vf,swapprob,Rskin);
	read_input(r,vf);
	wrap_coordinates(r);
	
	measure_packing_frac();
	create_neighbour_list(r);

	for(it=tstart+1;it<=Teq;it++)
	{
		evolve(r,swapeq);
		if(it % GAP == 0)
		{
			fp=fopen(outfile5,"w");
			fprintf(fp,"# it = %d\n",it);
			fclose(fp);	
		}
		if(it % gapres == 0)
                        print_restartfile(r,it,0);
	}	
	
	for(it=stepstart+1;it<=nsteps;it++)
	{
		evolve(r,swapprob);
		if(it % GAP == 0)
			print_configuration(r,it);
		if(it % gapres == 0)
                        print_restartfile(r,Teq,it);
	}
}
