#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "SC.h"
#include <time.h>
#include <string.h>
#include "myhelpers.h"
//#include "mysyslib.h"

double3* r;
double3* rhost;
double3* rchost;
double3* r1;
double3* r1host;
double3* r2;
double3* r2host;
double3* u;
double3* uhost;
double2* theta;
double2* thetahost;
double* l;
double* lhost;
double* sigma;
double* sigmahost;
double3* GU1;
double3* GU1host;
double2* GU1A;
double2* GU1Ahost;
double3* GU0;
double2* GU0A;
double3* dr;
double2* dtheta;
double3* h;
double2* hA;
double* GUI;
double3* GUII;
double2* GUIIA;
double2* UW;
double* GS;
double* GSA;
double* U;
double* p1;
double* GShost;
double* GSAhost;
double* Uhost;
double* p1host;
double* Vars; //Vars[0] == GS; Vars[1] == GSA; Vars[2] == U; Vars[3] == p1;
double* Varshost;

size_t pitch;

dim3 blocks, threads;
dim3 grid, block;

FILE *datfile;

void doOpen(const char *fname, const char *mode) {
	fprintf(stderr, "Opening %s with mode %s\n", fname, mode);
	datfile = fopen(fname, mode);
	if (!datfile) {
		perror("fopen");
		exit(2);
	}
}

void doClose() {
	if (fclose(datfile))
		perror("fclose");
}

int safeRead(void *buf, size_t size, size_t nr) {
	size_t num = fread(buf, size, nr, datfile);
	if (num < nr) {
		if (feof(datfile)) {
			fprintf(stderr, "Reached EOF\n");
			return 1;
		}
		else {
			perror("fread");
			exit(-1);
		}
	}
	return 0;
}

int doRead(void *buf, size_t size, size_t nr) {
	size_t num = fread(buf, size, nr, datfile);
	if (num < nr) {
		if (feof(datfile)) {
			fprintf(stderr, "Reached EOF, exitting...\n");
			exit(0);
		}
		else {
			perror("fread");
			exit(-1);
		}
	}
	return 0;
}

void doWrite(const void *buf, size_t size, size_t nr) {
	if (!fwrite(buf, size, nr, datfile)) {
		perror("fwrite");
		exit(-1);
	}
}


__attribute__ ((malloc)) void * xmalloc(size_t size) {
	void *p;
	p=malloc(size);
	if (!p) {
		perror("xmalloc");
		exit(-1);
	}
	return p;
}

void xmemcpy(void *out, void *in, size_t size) {
        if (!memcpy(out, in, size)) {
		perror("xmemcpy");
		exit(1);
        }
}


__device__ __host__ double2 d2add(double2 a, double2 b) {
	double2 ret;
	ret.x=a.x+b.x;
	ret.y=a.y+b.y;
	return ret;
}

__device__ __host__ double2 d2sub(double2 a, double2 b) {
	double2 ret;
	ret.x=a.x-b.x;
	ret.y=a.y-b.y;
	return ret;
}

__device__ __host__ double2 d2multscal(double2 a, double b) {
	double2 ret;
	ret.x=b*a.x;
	ret.y=b*a.y;
	return ret;
}

__device__ __host__ double2 d2multscal2(double2 a, double b) {
	double2 ret;
	ret.x=b*b*a.x;
	ret.y=b*b*a.y;
	return ret;
}

__device__ __host__ double2 d2divscal(double2 a, double b) {
	double2 ret;
	ret.x=a.x/b;
	ret.y=a.y/b;
	return ret;
}

__device__ __host__ double d2dotp(double2 a, double2 b) {
	return (a.x*b.x)+(a.y*b.y);
}

__device__ __host__ double d2mag(double2 a) {
	return sqrt((a.x)*(a.x)+(a.y)*(a.y));
}

__device__ __host__ double d2dist(double2 a, double2 b) {
	return d2mag(d2sub(a,b));
}

__device__ __host__ double2 d2unit(double2 a) {
	return d2divscal(a, d2mag(a));
}

__device__ __host__ double3 d3add(double3 a, double3 b) {
	double3 ret;
	ret.x=a.x+b.x;
	ret.y=a.y+b.y;
	ret.z=a.z+b.z;
	return ret;
}

__device__ __host__ double3 d3sub(double3 a, double3 b) {
	double3 ret;
	ret.x=a.x-b.x;
	ret.y=a.y-b.y;
	ret.z=a.z-b.z;
	return ret;
}

__device__ __host__ double3 d3multscal(double3 a, double b) {
	double3 ret;
	ret.x=b*a.x;
	ret.y=b*a.y;
	ret.z=b*a.z;
	return ret;
}

__device__ __host__ double3 d3multscal2(double3 a, double b) {
	double3 ret;
	ret.x=b*b*a.x;
	ret.y=b*b*a.y;
	ret.z=b*b*a.z;
	return ret;
}

__device__ __host__ double3 d3divscal(double3 a, double b) {
	double3 ret;
	ret.x=a.x/b;
	ret.y=a.y/b;
	ret.z=a.z/b;
	return ret;
}

__device__ __host__ double d3dotp(double3 a, double3 b) {
	return (a.x*b.x)+(a.y*b.y)+(a.z*b.z);
}

__device__ __host__ double3 d3crossp(double3 a, double3 b) {
	double3 ret;
	ret.x=(a.y*b.z)-(a.z*b.y);
	ret.y=-(a.x*b.z)+(a.z*b.x);
	ret.z=(a.x*b.y)-(a.y*b.x);
	return ret;
}

__device__ __host__ double d3mag(double3 a) {
	return sqrt(a.x*a.x+a.y*a.y+a.z*a.z);
}

__device__ __host__ double d3dist(double3 a, double3 b) {
	return d3mag(d3sub(a,b));
}

__device__ __host__ double3 d3unit(double3 a) {
	return d3divscal(a,d3mag(a));
}

//THESE FUNCTIONS ARE FOR THE SPHEROCYLINDER CODE//

__device__ __host__ double d3SCdist(double3 ri, double3 rj, double3 ui, double3 uj, double li, double lj) {
	double3 dij;
	dij=d3sub(d3sub(d3add(rj,d3multscal(uj,lj)),ri),d3multscal(ui,li));
	return d3mag(dij);
}

__device__ __host__ double L_i(double3 ri, double3 rj, double3 ui, double3 uj, double li) {
	double retn,retd;
	retn=d3dotp(ui,d3sub(rj,ri))-d3dotp(ui,uj)*d3dotp(uj,d3sub(rj,ri));
	retd=1.0-(d3dotp(ui,uj)*d3dotp(ui,uj));
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

__device__ __host__ double L_j(double3 ri, double3 rj, double3 ui, double3 uj, double lj) {
	double retn,retd;
	retn=d3dotp(ui,uj)*d3dotp(ui,d3sub(rj,ri))-d3dotp(uj,d3sub(rj,ri));
	retd=1.0-(d3dotp(ui,uj)*d3dotp(ui,uj));
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

double unitrand() {
	return (((double)rand())/((double)(RAND_MAX)));
}

void initcond(double3 *x, double2 *ang, double *len, double *diam) {
	int i;
	srand(time(NULL));
	for(i=0;i<npart;i++) {
		x[i].x=R*unitrand();
		x[i].y=2*M_PI*unitrand();
		x[i].z=H*unitrand();
		len[i]=L;
		diam[i]=SIGMA;
		ang[i].x=unitrand()*2*M_PI;
		ang[i].y=unitrand()*2*M_PI;
	}
}

__device__ __host__ void sptoca(int i, double2* theta, double3* u) {
	u[i].x=sin(theta[i].y)*cos(theta[i].x);
	u[i].y=sin(theta[i].y)*sin(theta[i].x);
	u[i].z=cos(theta[i].y);
}

__device__ __host__ void ends(int i,double3* r,double3* u,double* l,double3* r1,double3* r2) {
	r1[i]=d3add(r[i],d3multscal(u[i],0.5*l[i]));
	r2[i]=d3sub(r[i],d3multscal(u[i],0.5*l[i]));
	r1[i].x=hypot(r1[i].x,r1[i].y);
	r2[i].x=hypot(r2[i].x,r2[i].y);
}

void start() {
	cudaSetDevice(1);

	if(npart<1024) {
		blocks.x=1;
		threads.x=npart;

		grid.x=npart;
		grid.y=npart;
		grid.z=1;
		block.x=1;
		block.y=1;
		block.z=1;
	}
	else {
		blocks.x=npart/512;
		threads.x=512;

		grid.x=512;
		grid.y=512;
		grid.z=1;
		block.x=npart/512;
		block.y=npart/512;
		block.z=1;
	}
	HANDLE_ERROR(cudaMalloc(&r,npart*sizeof(double3)));
	HANDLE_ERROR(cudaMalloc(&r1,npart*sizeof(double3)));
	HANDLE_ERROR(cudaMalloc(&r2,npart*sizeof(double3)));
	HANDLE_ERROR(cudaMalloc(&u,npart*sizeof(double3)));
	HANDLE_ERROR(cudaMalloc(&theta,npart*sizeof(double2)));
	HANDLE_ERROR(cudaMalloc(&l,npart*sizeof(double)));
	HANDLE_ERROR(cudaMalloc(&sigma,npart*sizeof(double)));
	HANDLE_ERROR(cudaMalloc(&GU1,npart*sizeof(double3)));
	HANDLE_ERROR(cudaMalloc(&GU1A,npart*sizeof(double2)));
	HANDLE_ERROR(cudaMalloc(&GU0,npart*sizeof(double3)));
	HANDLE_ERROR(cudaMalloc(&GU0A,npart*sizeof(double3)));
	HANDLE_ERROR(cudaMalloc(&dr,npart*sizeof(double3)));
	HANDLE_ERROR(cudaMalloc(&dtheta,npart*sizeof(double2)));
	HANDLE_ERROR(cudaMalloc(&h,npart*sizeof(double3)));
	HANDLE_ERROR(cudaMalloc(&hA,npart*sizeof(double2)));
	HANDLE_ERROR(cudaMalloc(&GUI,npart*sizeof(double)));
	HANDLE_ERROR(cudaMalloc(&GUII,npart*sizeof(double3)));
	HANDLE_ERROR(cudaMalloc(&GUIIA,npart*sizeof(double2)));
	HANDLE_ERROR(cudaMalloc(&UW,npart*sizeof(double2)));
	HANDLE_ERROR(cudaMalloc(&GS,sizeof(double)));
	HANDLE_ERROR(cudaMalloc(&GSA,sizeof(double)));
	HANDLE_ERROR(cudaMalloc(&U,sizeof(double)));
	HANDLE_ERROR(cudaMalloc(&p1,sizeof(double)));
	rhost=(double3*)malloc(npart*sizeof(double3));
	r1host=(double3*)malloc(npart*sizeof(double3));
	r2host=(double3*)malloc(npart*sizeof(double3));
	uhost=(double3*)malloc(npart*sizeof(double3));
	thetahost=(double2*)malloc(npart*sizeof(double2));
	lhost=(double*)malloc(npart*sizeof(double));
	sigmahost=(double*)malloc(npart*sizeof(double));
	GU1host=(double3*)malloc(npart*sizeof(double3));
	GU1Ahost=(double2*)malloc(npart*sizeof(double2));
	HANDLE_ERROR(cudaMallocHost(&Varshost,5*sizeof(double)));
	HANDLE_ERROR(cudaMalloc(&Vars,5*sizeof(double)));
	HANDLE_ERROR(cudaMallocHost(&GShost,sizeof(double)));
	HANDLE_ERROR(cudaMallocHost(&GSAhost,sizeof(double)));
	HANDLE_ERROR(cudaMallocHost(&Uhost,sizeof(double)));
	HANDLE_ERROR(cudaMallocHost(&p1host,sizeof(double)));
	rchost=(double3*)malloc(npart*sizeof(double3));
	initcond(rchost,thetahost,lhost,sigmahost);
	int i;
	for(i=0;i<npart;i++) {
		rhost[i].x=rchost[i].x*cos(rchost[i].y);
		rhost[i].y=rchost[i].x*sin(rchost[i].y);
		rhost[i].z=rchost[i].z;
		sptoca(i,thetahost,uhost);
	}
	HANDLE_ERROR(cudaMemcpy(r,rhost,npart*sizeof(double3),cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpy(theta,thetahost,npart*sizeof(double2),cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpy(u,uhost,npart*sizeof(double3),cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpy(l,lhost,npart*sizeof(double),cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpy(sigma,sigmahost,npart*sizeof(double),cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemset(dr,0.0,npart*sizeof(double3)));
	HANDLE_ERROR(cudaMemset(GU1,0.0,npart*sizeof(double3)));
	HANDLE_ERROR(cudaMemset(GU1A,0.0,npart*sizeof(double2)));
}

__device__ __host__ double PotEnergyI(int i, int j, double3* r, double3* u, double2* theta, double* l, double* sigma, double3* r1, double3* r2) {
	double ret=0.0;
	sptoca(i,theta,u);
	sptoca(j,theta,u);
	ends(i,r,u,l,r1,r2);
	ends(j,r,u,l,r1,r2);
	if((i!=j) && (d3dist(r[i],r[j])<(l[i]+sigma[i]+l[j]+sigma[j])/2.0)) {
		double lambda_i, lambda_j;
		lambda_i=L_i(r[i],r[j],u[i],u[j],l[i]);
		lambda_j=L_j(r[i],r[j],u[i],u[j],l[j]);
		double d;
		d=d3SCdist(r[i],r[j],u[i],u[j],lambda_i,lambda_j);
		if(d<(sigma[i]+sigma[j])/2.0) {
			ret+=0.5*((sigma[i]+sigma[j])/2.0-d)*((sigma[i]+sigma[j])/2.0-d);
		}
	}
	return ret;
}

__device__ __host__ double WallEnergy(int i, double3* r, double3* u, double2* theta,  double* l, double* sigma, double3* r1, double3* r2) {
	double ret=0.0;
	sptoca(i,theta,u);
	ends(i,r,u,l,r1,r2);
	
	if(r1[i].x>R-(sigma[i]/2.0)) {
		ret+=0.5*(r1[i].x-(R-sigma[i]/2.0))*(r1[i].x-(R-sigma[i]/2.0));
	}
	if(r2[i].x>R-(sigma[i]/2.0)) {
		ret+=0.5*(r2[i].x-(R-sigma[i]/2.0))*(r2[i].x-(R-sigma[i]/2.0));
	}
	if(r1[i].z-sigma[i]/2.0<0.0) {
		ret+=0.5*(r1[i].z-sigma[i]/2.0)*(r1[i].z-sigma[i]/2.0);
	}
	if(r1[i].z+sigma[i]/2.0>H) {
		ret+=0.5*(r1[i].z+sigma[i]/2.0-H)*(r1[i].z+sigma[i]/2.0-H);
	}
	if(r2[i].z-sigma[i]/2.0<0.0) {
		ret+=0.5*(r2[i].z-sigma[i]/2.0)*(r2[i].z-sigma[i]/2.0);
	}
	if(r2[i].z+sigma[i]/2.0>H) {
		ret+=0.5*(r2[i].z+sigma[i]/2.0-H)*(r2[i].z+sigma[i]/2.0-H);
	}
	return ret;
}

__global__ void Zero(double* GUI, double3* GUII, double2* GUIIA) {
	int tx=threadIdx.x+blockIdx.x*blockDim.x;
	*(GUI+tx)=0.0;
	((GUII+tx)->x)=0.0;
	((GUII+tx)->y)=0.0;
	((GUII+tx)->z)=0.0;
	((GUIIA+tx)->x)=0.0;
	((GUIIA+tx)->y)=0.0;
}

__global__ void GradI(double3* r, double3* u, double2* theta, double* l, double* sigma, double3* r1, double3* r2, double* GUI, double3* GUII, double2* GUIIA) {
	int tx=threadIdx.x+blockIdx.x*blockDim.x;
	int ty=threadIdx.y+blockIdx.y*blockDim.y;

	*(GUI+tx)+=PotEnergyI(tx,ty,r,u,theta,l,sigma,r1,r2);

	double dx=.0001;
	((r+tx)->x)+=dx;
	((GUII+tx)->x)+=PotEnergyI(tx,ty,r,u,theta,l,sigma,r1,r2);
	((r+tx)->x)-=dx;

	double dy=.0001;
	((r+tx)->y)+=dy;
	((GUII+tx)->y)+=PotEnergyI(tx,ty,r,u,theta,l,sigma,r1,r2);
	((r+tx)->y)-=dy;

	double dz=.0001;
	((r+tx)->z)+=dz;
	((GUII+tx)->z)+=PotEnergyI(tx,ty,r,u,theta,l,sigma,r1,r2);
	((r+tx)->z)-=dz;

	double dt=.0001;
	((theta+tx)->x)+=dt;
	((GUIIA+tx)->x)+=PotEnergyI(tx,ty,r,u,theta,l,sigma,r1,r2);
	((theta+tx)->x)-=dt;

	double dp=.0001;
	((theta+tx)->y)+=dp;
	((GUIIA+tx)->y)+=PotEnergyI(tx,ty,r,u,theta,l,sigma,r1,r2);
	((theta+tx)->y)-=dp;
}

__global__ void GradWall(double3* r, double3* u, double2* theta, double* l, double* sigma, double3* r1, double3* r2, double* GUI, double3* GUII, double2* GUIIA) {
	int tx=threadIdx.x+blockIdx.x*blockDim.x;

	*(GUI+tx)+=WallEnergy(tx,r,u,theta,l,sigma,r1,r2);

	double dx=.0001;
	((r+tx)->x)+=dx;
	((GUII+tx)->x)+=WallEnergy(tx,r,u,theta,l,sigma,r1,r2);
	((r+tx)->x)-=dx;

	double dy=.0001;
	((r+tx)->y)+=dy;
	((GUII+tx)->y)+=WallEnergy(tx,r,u,theta,l,sigma,r1,r2);
	((r+tx)->y)-=dy;

	double dz=.0001;
	((r+tx)->z)+=dz;
	((GUII+tx)->z)+=WallEnergy(tx,r,u,theta,l,sigma,r1,r2);
	((r+tx)->z)-=dz;

	double dt=.0001;
	((theta+tx)->x)+=dt;
	((GUIIA+tx)->x)+=WallEnergy(tx,r,u,theta,l,sigma,r1,r2);
	((theta+tx)->x)-=dt;

	double dp=.0001;
	((theta+tx)->y)+=dp;
	((GUIIA+tx)->y)+=WallEnergy(tx,r,u,theta,l,sigma,r1,r2);
	((theta+tx)->y)-=dp;
}

__global__ void GradII(double* GUI, double3* GUII, double2* GUIIA, double3* GU1, double2* GU1A,double* sigma) {
	int tx=threadIdx.x+blockIdx.x*blockDim.x;

	double dx=.0001;
	((GU1+tx)->x)=(((GUII+tx)->x)-(*(GUI+tx)))/dx;

	double dy=.0001;
	((GU1+tx)->y)=(((GUII+tx)->y)-(*(GUI+tx)))/dy;

	double dz=.0001;
	((GU1+tx)->z)=(((GUII+tx)->z)-(*(GUI+tx)))/dz;

	double dt=.0001;
	((GU1A+tx)->x)=(((GUIIA+tx)->x)-(*(GUI+tx)))/dt;
	
	double dp=.0001;
	((GU1A+tx)->y)=(((GUIIA+tx)->y)-(*(GUI+tx)))/dp;

	if(sigma[tx]==0.0) {
		((GU1+tx)->x)=0.0;
		((GU1+tx)->y)=0.0;
		((GU1+tx)->z)=0.0;
		((GU1A+tx)->x)=0.0;
		((GU1A+tx)->y)=0.0;
	}
}

__global__ void ConGradI(double3* dr, double3* h, double3* GU1) {
	int tx=threadIdx.x+blockIdx.x*blockDim.x;
	
	double eta=-.001;
	h[tx]=GU1[tx];
	dr[tx]=d3multscal(h[tx],eta);
}

__global__ void ConGradII(double3* GU1, double3* GU0, double3* dr, double3* h) {
	int tx=threadIdx.x+blockIdx.x*blockDim.x;

	((dr+tx)->x)=0.0;
	((dr+tx)->y)=0.0;
	((dr+tx)->z)=0.0;
	if(d3mag(GU1[tx])>0.0) {
		double gamma;
		double eta=-.001;
		gamma=d3dotp(GU1[tx],GU1[tx])/d3dotp(GU0[tx],GU0[tx]);
		if(gamma>1.0) {
			gamma=1.0;
		}
		if(gamma<-1.0) {
			gamma=-1.0;
		}
		h[tx]=d3add(GU1[tx],d3multscal(h[tx],gamma));
		dr[tx]=d3multscal(h[tx],eta);
	}
}

__global__ void ConGradAI(double2* dtheta, double2* hA, double2* GU1A) {
	int tx=threadIdx.x+blockIdx.x*blockDim.x;

	double eta=-.001;
	hA[tx]=GU1A[tx];
	dtheta[tx]=d2multscal(hA[tx],eta);
}

__global__ void ConGradAII(double2* GU1A, double2* GU0A, double2* dtheta, double2* hA) {
	int tx=threadIdx.x+blockIdx.x*blockDim.x;

	((dtheta+tx)->x)=0.0;
	((dtheta+tx)->y)=0.0;
	if(d2mag(GU1A[tx])>0.0) {
		double gamma;
		double eta=-.001;
		gamma=d2dotp(GU1A[tx],GU1A[tx])/d2dotp(GU0A[tx],GU0A[tx]);
		if(gamma>1.0) {
			gamma=1.0;
		}
		if(gamma<-1.0) {
			gamma=-1.0;
		}
		hA[tx]=d2add(GU1A[tx],d2multscal(hA[tx],gamma));
		dtheta[tx]=d2multscal(hA[tx],eta);
	}
}

__global__ void ConGradMove(double3* r, double3* dr, double2* theta, double2* dtheta, double3* GU1, double3* GU0, double2* GU1A, double2* GU0A) {
	int tx=threadIdx.x+blockIdx.x*blockDim.x;
	
	if(tx>npart) {
		return;
	}

	r[tx]=d3add(r[tx],dr[tx]);
	GU0[tx]=GU1[tx];

	if(dtheta[tx].x>1.0) {
		dtheta[tx].x=1.0;
	}
	else if(dtheta[tx].x<-1.0) {
		dtheta[tx].x=-1.0;
	}
	if(dtheta[tx].y>1.0) {
		dtheta[tx].y=1.0;
	}
	else if(dtheta[tx].y<-1.0) {
		dtheta[tx].y=-1.0;
	}
	theta[tx]=d2add(theta[tx],dtheta[tx]);
	GU0A[tx].x=GU1A[tx].x;
	GU0A[tx].y=GU1A[tx].y;

	if((theta[tx].x>2.0*M_PI) || (theta[tx].x<0.0)) {
		double t,td;
		td=theta[tx].x/(2.0*M_PI);
		t=floor(td);
		theta[tx].x=(2.0*M_PI)*(td-t);
	}
	if((theta[tx].y>2.0*M_PI) || (theta[tx].y<0.0)) {
		double t,td;
		td=theta[tx].y/(2.0*M_PI);
		t=floor(td);
		theta[tx].y=(2.0*M_PI)*(td-t);
	}
}

__global__ void PotEnergy(double3* r, double3* u, double2* theta, double* l, double* sigma, double3* r1, double3* r2, double* U) {
	int tx=threadIdx.x+blockIdx.x*blockDim.x;
	
	int ty;
	for(ty=0;ty<npart;ty++) {
		U[2]+=PotEnergyI(tx,ty,r,u,theta,l,sigma,r1,r2);
	}
	U[2]+=WallEnergy(tx,r,u,theta,l,sigma,r1,r2);
}

double contacts() {
	int i, j;
	double ret=0.0;
	for(i=0;i<npart;i++) {
		for(j=0;j<npart;j++) {
			sptoca(i,thetahost,uhost);
			ends(i,rhost,uhost,lhost,r1host,r2host);
			if((i!=j) && (d3dist(rhost[i],rhost[j])<(lhost[i]+sigmahost[i]+lhost[j]+sigmahost[j])/2.0)) {
				double lambda_i, lambda_j;
				lambda_i=L_i(rhost[i],rhost[j],uhost[i],uhost[j],lhost[i]);
				lambda_j=L_j(rhost[i],rhost[j],uhost[i],uhost[j],lhost[j]);
				double d;
				d=d3SCdist(rhost[i],rhost[j],uhost[i],uhost[j],lambda_i,lambda_j);
				if(d<(sigmahost[i]+sigmahost[j])*(.9)/2.0) {
					ret+=1.0;
				}
			}
		}
		sptoca(i,thetahost,uhost);
		ends(i,rhost,uhost,lhost,r1host,r2host);
		if(r1host[i].x>R-(sigmahost[i]/2.0)*(.95)) {
			ret+=1.0;
		}
		if(r2host[i].x>R-(sigmahost[i]/2.0)*(.95)) {
			ret+=1.0;
		}
		if(r1host[i].z-(sigmahost[i]*.95)/2.0<0.0) {
			ret+=1.0;
		}
		if(r1host[i].z+(sigmahost[i]*.95)/2.0>H) {
			ret+=1.0;
		}
		if(r2host[i].z-(sigmahost[i]*.95)/2.0<0.0) {
			ret+=1.0;
		}
		if(r2host[i].z+(sigmahost[i]*.95)/2.0>H) {
			ret+=1.0;
		}
	}
	return ret;	
}

double packfrac() {
	double Vp=0.0;
	double Vb=M_PI*R*R*H;
	int i=0;
	Vp=npart*(M_PI*(sigmahost[i]/2.0)*(sigmahost[i]/2.0)*lhost[i]+(4.0/3.0)*M_PI*(sigmahost[i]/2.0)*(sigmahost[i]/2.0)*(sigmahost[i]/2.0));
	return Vp/Vb;
}

__global__ void Pressure(double3* r, double3* u, double2* theta, double* l, double* sigma, double3* r1, double3* r2, double* p) {
	int tx=threadIdx.x+blockIdx.x*blockDim.x;

	double3 F;
	double3 P;
	F.x=0.0;
	F.y=0.0;
	F.z=0.0;
	P.x=0.0;
	P.y=0.0;
	P.z=0.0;

	sptoca(tx,theta,u);
	ends(tx,r,u,l,r1,r2);
	if(r1[tx].x>R-(sigma[tx]/2.0)) {
		F.x+=r1[tx].x-(R-sigma[tx]/2.0);
	}
	if(r2[tx].x>R-(sigma[tx]/2.0)) {
		F.x+=r2[tx].x-(R-sigma[tx]/2.0);
	}
	if(r1[tx].z-sigma[tx]/2.0<0.0) {
		F.z+=r1[tx].z-sigma[tx]/2.0;
	}
	if(r1[tx].z+sigma[tx]/2.0>H) {
		F.y+=r1[tx].z+sigma[tx]/2.0-H;
	}
	if(r2[tx].z-sigma[tx]/2.0<0.0) {
		F.z+=r2[tx].z-sigma[tx]/2.0;
	}
	if(r2[tx].z+sigma[tx]/2.0>H) {
		F.y+=r2[tx].z+sigma[tx]/2.0-H;
	}
	P.x=F.x/(2.0*M_PI*R*H);
	P.y=F.y/(M_PI*R*R);
	P.z=F.z/(M_PI*R*R);
	p[3]+=fabs(P.x)+fabs(P.y)+fabs(P.z);
}

__global__ void GradSum(double3* GU1, double* GS) {
	int tx=threadIdx.x+blockIdx.x*blockDim.x;
	
	GS[0]+=d3mag(GU1[tx]);
}

__global__ void GradSumA(double2* GU1A, double* GSA) {
	int tx=threadIdx.x+blockIdx.x*blockDim.x;

	GSA[1]+=d2mag(GU1A[tx]);
}

__global__ void ZeroVars(double* Vars) {
	int i;
	for(i=0;i<4;i++) {
		Vars[i]=0.0;
	}
}

void doAnIter() {
	//clock_t t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15;
	//t0=clock();
	Zero<<<blocks, threads>>>(GUI,GUII,GUIIA);
	//t1=clock();
	GradI<<<grid, block>>>(r,u,theta,l,sigma,r1,r2,GUI,GUII,GUIIA);
	//t2=clock();
	GradWall<<<blocks, threads>>>(r,u,theta,l,sigma,r1,r2,GUI,GUII,GUIIA);
	//t3=clock();
	GradII<<<blocks, threads>>>(GUI,GUII,GUIIA,GU1,GU1A,sigma);
	//t4=clock();
	ConGradI<<<blocks, threads>>>(dr,h,GU1);
	//t5=clock();
	ConGradAI<<<blocks, threads>>>(dtheta,hA,GU1A);
	//t6=clock();
	ConGradMove<<<blocks, threads>>>(r,dr,theta,dtheta,GU1,GU0,GU1A,GU0A);
	//t7=clock();
	ZeroVars<<<1,1>>>(Vars);
	//t8=clock();
	GradSum<<<blocks, threads>>>(GU1,Vars);
	//t9=clock();
	GradSumA<<<blocks, threads>>>(GU1A,Vars);
	//t10=clock();
	Pressure<<<threads, blocks>>>(r,u,theta,l,sigma,r1,r2,Vars);
	PotEnergy<<<threads, blocks>>>(r,u,theta,l,sigma,r1,r2,Vars);
	HANDLE_ERROR(cudaMemcpy(Varshost,Vars,5*sizeof(double),cudaMemcpyDeviceToHost));
	HANDLE_ERROR(cudaMemcpy(rhost,r,npart*sizeof(double3),cudaMemcpyDeviceToHost));
	HANDLE_ERROR(cudaMemcpy(uhost,u,npart*sizeof(double3),cudaMemcpyDeviceToHost));
	//t11=clock();
	int i;
	Varshost[4]=1.0;
	while((fabs(Varshost[0])/(npart*lhost[1])>0.000000001) && (fabs(Varshost[1])/(npart*lhost[1])>0.00000001)) {
		Zero<<<blocks, threads>>>(GUI,GUII,GUIIA);
		GradI<<<grid, block>>>(r,u,theta,l,sigma,r1,r2,GUI,GUII,GUIIA);
		GradWall<<<blocks, threads>>>(r,u,theta,l,sigma,r1,r2,GUI,GUII,GUIIA);
		GradII<<<blocks, threads>>>(GUI,GUII,GUIIA,GU1,GU1A,sigma);
		//t12=clock();
		ConGradII<<<blocks, threads>>>(GU1,GU0,dr,h);
		//t13=clock();
		ConGradAII<<<blocks,threads>>>(GU1A,GU0A,dtheta,hA);
		//t14=clock();
		ConGradMove<<<blocks,threads>>>(r,dr,theta,dtheta,GU1,GU0,GU1A,GU0A);
		ZeroVars<<<1,1>>>(Vars);
		GradSum<<<blocks, threads>>>(GU1,Vars);
		GradSumA<<<blocks, threads>>>(GU1A,Vars);
		PotEnergy<<<threads, blocks>>>(r,u,theta,l,sigma,r1,r2,Vars);
		Pressure<<<threads, blocks>>>(r,u,theta,l,sigma,r1,r2,Vars);
		Varshost[4]=Varshost[4]+1.0;
		HANDLE_ERROR(cudaMemcpy(Varshost,Vars,5*sizeof(double),cudaMemcpyDeviceToHost));
		HANDLE_ERROR(cudaMemcpy(rhost,r,npart*sizeof(double3),cudaMemcpyDeviceToHost));
		HANDLE_ERROR(cudaMemcpy(uhost,u,npart*sizeof(double3),cudaMemcpyDeviceToHost));
		fprintf(stderr,"Grad Sum: %lf\n",Varshost[0]);
	}
	//t15=clock();
	double c=contacts();
	double p=packfrac();
	Uhost[0]=Varshost[2];
	p1host[0]=Varshost[3];
	fprintf(stderr,"Rod Radius: %lf\n", sigmahost[0]/2.0);
	fprintf(stderr,"Rod Length: %lf\n", lhost[0]);
	fprintf(stderr,"Packing Fraction: %lf\n",p);
	fprintf(stderr,"Average Number of Contacts: %lf\n",c/npart);
	fprintf(stderr,"Potential Energy: %lf\n",Uhost[0]);
	fprintf(stderr,"Average Overlap: %lf\n",sqrt((Uhost[0]/npart)*2.0)/(sigmahost[0]/2.0));
	fprintf(stderr,"Pressure: %lf\n",p1host[0]);
	fprintf(stderr,"GradSum: %lf\n",Varshost[0]);
	fprintf(stderr,"GradSumA: %lf\n",Varshost[1]);
	fprintf(stdout,"%lf	%lf	%lf\n",p,p1host[0],Uhost[0]);
	/*for(i=0;i<npart;i++) {
		fprintf(stdout,"%i	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf\n",i,rhost[i].x,rhost[i].y,rhost[i].z,uhost[i].x,uhost[i].y,uhost[i].z,sigmahost[i],lhost[i]);
	}*/
	if(Varshost[0]+Varshost[1]==0.0) {
		for(i=0;i<npart;i++) {
			lhost[i]+=0.001;
			sigmahost[i]=(1.0/ALPHA)*lhost[i];
		}
	}
	else {
		for(i=0;i<npart;i++) {
			lhost[i]+=0.0001;
			sigmahost[i]=(1.0/ALPHA)*lhost[i];
		}
	}
	HANDLE_ERROR(cudaMemcpy(l,lhost,npart*sizeof(double),cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpy(sigma,sigmahost,npart*sizeof(double),cudaMemcpyHostToDevice));
}

int main(int argc, char *argv[]) {
	if(argc>1) {
		doOpen(argv[1],"w");
	}
	else {
		fprintf(stderr,"No File Name Given\n");
	}
	const char header[4]="SP2";

	unsigned int n=npart;
	double alpha=ALPHA;
	double radius=R;
	double height=H;

	doWrite(header,4,1);
	doWrite(&n,sizeof(unsigned int),1);
	doWrite(&alpha,sizeof(double),1);
	doWrite(&radius,sizeof(double),1);
	doWrite(&height,sizeof(double),1);

	clock_t t2, t3;
	t2=clock();
	start();
	HANDLE_ERROR(cudaMemcpy(sigmahost,sigma,npart*sizeof(double),cudaMemcpyDeviceToHost));
	HANDLE_ERROR(cudaMemcpy(lhost,l,npart*sizeof(double),cudaMemcpyDeviceToHost));
	int i;
	clock_t t0, t1;
	while(packfrac()<PHI) {
		t0=clock();
		doAnIter();
		t1=clock();
		int q;
		for(q=0;q<npart;q++) {
			doWrite((sigmahost+q),sizeof(double),1);
			doWrite((lhost+q),sizeof(double),1);
			doWrite((rhost+q),sizeof(double3),1);
			doWrite((uhost+q),sizeof(double3),1);
		}
		fprintf(stderr,"Iteration Time: %lf\n",(t1-t0)/((double)CLOCKS_PER_SEC));
	}
	for(i=0;i<npart;i++) {
		fprintf(stdout,"%lf	%lf	%lf	%lf	%lf	%lf\n",rhost[i].x,rhost[i].y,rhost[i].z,uhost[i].x,uhost[i].y,uhost[i].z);
	}
	t3=clock();
	fprintf(stderr,"Number of Particles: %i\n",npart);
	fprintf(stderr,"Total Time Taken: %lf\n",(t3-t2)/((double)CLOCKS_PER_SEC));
	doClose();
	return 0;
}
