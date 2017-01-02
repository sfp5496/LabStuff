#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "double3.h"
#include "double2.h"
#include "mysyslib.h"
//#include "d3func.c"

unsigned int npart;
double ALPHA,R,H;
unsigned int timeskip;

double* sigma;
double* l;
double3* r;
double3* u;
double3* r1;
double3* r2;

double3 d3add(double3 a, double3 b) {
	double3 ret;
	ret.x=a.x+b.x;
	ret.y=a.y+b.y;
	ret.z=a.z+b.z;
	return ret;
}

double3 d3sub(double3 a, double3 b) {
	double3 ret;
	ret.x=a.x-b.x;
	ret.y=a.y-b.y;
	ret.z=a.z-b.z;
	return ret;
}

double3 d3multscal(double3 a, double b) {
	double3 ret;
	ret.x=b*a.x;
	ret.y=b*a.y;
	ret.z=b*a.z;
	return ret;
}

double3 d3multscal2(double3 a, double b) {
	double3 ret;
	ret.x=b*b*a.x;
	ret.y=b*b*a.y;
	ret.z=b*b*a.z;
}

double3 d3divscal(double3 a, double b) {
	double3 ret;
	ret.x=a.x/b;
	ret.y=a.y/b;
	ret.z=a.z/b;
	return ret;
}

double d3dotp(double3 a, double3 b) {
	return (a.x*b.x)+(a.y*b.y)+(a.z*b.z);
}

double3 d3crossp(double3 a, double3 b) {
	double3 ret;
	ret.x=(a.y*b.z)-(a.z*b.y);
	ret.y=-(a.x*b.z)+(a.z*b.x);
	ret.z=(a.x*b.y)-(a.y*b.x);
	return ret;
}

double d3mag(double3 a) {
	return sqrt(a.x*a.x+a.y*a.y+a.z*a.z);
}

double d3dist(double3 a, double3 b) {
	return d3mag(d3sub(a,b));
}

double3 d3unit(double3 a) {
	return d3divscal(a,d3mag(a));
}

//THESE FUNCTIONS ARE FOR THE SPHEROCYLINDER CODE//

double d3SCdist(double3 ri, double3 rj, double3 ui, double3 uj, double li, double lj) {
	double3 dij;
	dij=d3sub(d3sub(d3add(rj,d3multscal(uj,lj)),ri),d3multscal(ui,li));
	return d3mag(dij);
}

double L_i(double3 ri, double3 rj, double3 ui, double3 uj, double li) {
	double retn,retd;
	retn=d3dotp(ui,d3sub(rj,ri))-d3dotp(ui,uj)*d3dotp(uj,d3sub(rj,ri));
	retd=1.0-(d3dotp(ui,uj)*d3dotp(ui,uj));
	/*if(retd=0.0) {
		printf("L_i: Divide by 0 Error\n");
	}*/
	if((retn/retd)>(li/2.0)) {
		return (li/2.0);
	}
	else if((retn/retd)<(-li/2.0)) {
		return (-li/2.0);
	}
	else {
		return (retn/retd);
	}
}

double L_j(double3 ri, double3 rj, double3 ui, double3 uj, double lj) {
	double retn,retd;
	retn=d3dotp(ui,uj)*d3dotp(ui,d3sub(rj,ri))-d3dotp(uj,d3sub(rj,ri));
	retd=1.0-(d3dotp(ui,uj)*d3dotp(ui,uj));
	/*if(retd=0.0) {
		printf("L_j: Divide by 0 Error\n");
	}*/
	if((retn/retd)>(lj/2.0)) {
		return (lj/2.0);
	}
	else if((retn/retd)<(-lj/2.0)) {
		return (-lj/2.0);
	}
	else {
		return (retn/retd);
	}
}
void ends(int i) {
	//r1[i]=d3add(r[i],d3multscal(u[i],0.5*l[i]));
	r1[i].x=r[i].x+0.5*u[i].x*l[i];
	r1[i].y=r[i].y+0.5*u[i].y*l[i];
	r1[i].z=r[i].z+0.5*u[i].z*l[i];
	//r2[i]=d3sub(r[i],d3multscal(u[i],0.5*l[i]));
	r2[i].x=r[i].x+0.5*u[i].x*l[i];
	r2[i].y=r[i].y+0.5*u[i].y*l[i];
	r2[i].z=r[i].z+0.5*u[i].z*l[i];
	r1[i].x=hypot(r1[i].x,r1[i].y);
	r2[i].x=hypot(r2[i].x,r2[i].y);
}

double packfrac() {
	double Vp=0.0;
	double Vb=M_PI*R*R*H;
	int i=0;
	Vp=npart*(M_PI*(sigma[i]/2.0)*(sigma[i]/2.0)*l[i]+(4.0/3.0)*M_PI*(sigma[i]/2.0)*(sigma[i]/2.0)*(sigma[i]/2.0));
	return Vp/Vb;
	}

double PotEnergy() {
	int m;
	for (m=0;m<npart;m++) {
		double ret=0.0;
		ends(m);
		int k;
		for(k=0;k<npart;k++) {
			if((m!=k) && (d3dist(r[m],r[k])<(l[m]+sigma[m]+l[k]+sigma[k])/2.0)) {
				double lambda_i, lambda_j;
				lambda_i=L_i(r[m],r[k],u[m],u[k],l[m]);
				lambda_j=L_j(r[m],r[k],u[m],u[k],l[m]);
				double d;
				d=d3SCdist(r[m],r[k],u[m],u[k],lambda_i,lambda_j);
				if(d<(sigma[m]+sigma[k])/2.0) {
					ret+=0.25*((sigma[m]+sigma[k])/2.0-d)*((sigma[m]+sigma[k])/2.0-d);
				}
			}	
		}	
		if(r1[m].x>R-(sigma[m]/2.0)) {
			ret+=0.5*(r1[m].x-(R-sigma[m]/2.0))*(r1[m].x-(R-sigma[m]/2.0));
		}
		if(r2[m].x>R-(sigma[m]/2.0)) {
			ret+=0.5*(r2[m].x-(R-sigma[m]/2.0))*(r2[m].x-(R-sigma[m]/2.0));
		}
		if(r1[m].z-sigma[m]/2.0<0.0) {
			ret+=0.5*(r1[m].z-sigma[m]/2.0)*(r1[m].z-sigma[m]/2.0);
		}
		if(r1[m].z+sigma[m]/2.0>H) {
			ret+=0.5*(r1[m].z+sigma[m]/2.0-H)*(r1[m].z+sigma[m]/2.0-H);
		}
		if(r2[m].z-sigma[m]/2.0<0.0) {
			ret+=0.5*(r2[m].z-sigma[m]/2.0)*(r2[m].z-sigma[m]/2.0);
		}
		if(r2[m].z+sigma[m]/2.0>H) {
			ret+=0.5*(r2[m].z+sigma[m]/2.0-H)*(r2[m].z+sigma[m]/2.0-H);
		}
		return ret;
	}
}	
double Pressure() {
	int m;
	double3 F; //F.x is the force on the radial wall, F.y is on the top wall, F.z -> bottom
	double3 P;
	F.x=0.0;
	F.y=0.0;
	F.z=0.0;
	P.x=0.0;
	P.y=0.0;
	P.z=0.0;
	for(m=0;m<npart;m++) {
		ends(m);
		if(r1[m].x>R-(sigma[m]/2.0)) {
			F.x+=r1[m].x-(R-sigma[m]/2.0);
		}
		if(r2[m].x>R-(sigma[m]/2.0)) {
			F.x+=r2[m].x-(R-sigma[m]/2.0);
		}
		if(r1[m].z-sigma[m]/2.0<0.0) {
			F.z+=r1[m].z-sigma[m]/2.0;
		}
		if(r1[m].z+sigma[m]/2.0>H) {
			F.y+=r1[m].z+sigma[m]/2.0-H;
		}
		if(r2[m].z-sigma[m]/2.0<0.0) {
			F.z+=r2[m].z-sigma[m]/2.0;
		}
		if(r2[m].z+sigma[m]/2.0>H) {
			F.y+=r2[m].z+sigma[m]/2.0-H;
		}
	}
	P.x=F.x/(2.0*M_PI*R*H);
	P.y=F.y/(M_PI*R*R);
	P.z=F.z/(M_PI*R*R);
	return sqrt(P.x*P.x+P.y*P.y+P.z*P.z);
}


//unsigned int npart;
//double ALPHA,R,H;
//unsigned int timeskip;

//double* sigma;
//double* l;
//double3* r;
//double3* u;
//double3* r1;
//double3* r2;

void start() {
	sigma=malloc(npart*sizeof(double));
	l=malloc(npart*sizeof(double));
	r=malloc(npart*sizeof(double3));
	u=malloc(npart*sizeof(double3));
	r1=malloc(npart*sizeof(double3));
	r2=malloc(npart*sizeof(double3));
}

/*
 * It seems like the argv stuff that I'm going to use here is going to get out of hand real quick
 * So I might as well find a place to keep what each one means as I'm writing them so I don't forget
 *
 * argv[1] : the input file, probably generated by an SCsim.cu file
 * argv[2] : timeskip, number of frames you skip, any int > 0
 * argv[3] : if you have something here you do the analysis, else you just turn the data into something people can read, just put like a '1' there for simplicity sake
 */

int main(int argc, char *argv[]) {
	if(argc>1) {
		doOpen(argv[1],"r");
	}
	else {
		fprintf(stderr,"No File Name Given\n");
		return 1;
	}

	/*Probably won't need timeskip stuff yet, but it might come in handy when I do animations*/
	if(argc>2) {
		timeskip=(unsigned int)argv[2]; //skip as many frames as the timeskip tells you to
	}
	else {
		timeskip=0; //Look at all frames
	}
	
	/*Check to see if it's the right file type*/
	const char header[4];
	const char thisn[4]="SP2";
	doRead(&header,4,1);
	if(strcmp(header,thisn)==0) {
		doRead(&npart,sizeof(unsigned int),1);
		doRead(&ALPHA,sizeof(double),1);
		doRead(&R,sizeof(double),1);
		doRead(&H,sizeof(double),1);
		printf("\n========================================\n");
		printf("Your parameters for this run were:\n");
		printf("Number of Particles: %i\n",npart);
		printf("Particle Aspect Ratio: %lf\n",ALPHA);
		printf("Container Radius: %lf\n",R);
		printf("Container Height: %lf\n",H);
		printf("========================================\n\n");
	}
	else {
		fprintf(stderr,"wrong file type, make sure that this is an 'SP2' file\n");
		return 1;
	}
	start();
	int i;
	int j = 0;
	while(1) {
		for(i=0;i<npart;i++) {
			doRead((sigma+i),sizeof(double),1);
			doRead((l+i),sizeof(double),1);
			doRead((r+i),sizeof(double3),1);
			doRead((u+i),sizeof(double3),1);
		}
		if(argc>3) {
			//Analysis code here...
			//Basically just steal it from the program that generated the data
				
			double U = 0.0;
			double P = 0.0;
			double p = 0.0;

			U=PotEnergy();
			P=Pressure();
			p=packfrac();
				
			fprintf(stdout,"%i	%lf	%lf	%lf\n",j,p,P,U);
		}
		else {
			for(i=0;i<npart;i++) {
				fprintf(stdout,"%i	%i	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf\n",j,i,sigma[i],l[i],r[i].x,r[i].y,r[i].z,u[i].x,u[i].y,u[i].z);
			}
		}
		j++;
	}
}
