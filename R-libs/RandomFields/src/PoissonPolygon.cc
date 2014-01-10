/*******************************************************************************
 * PoissonPolygon.c
 *
 * written by Felix Ballani, 2012
 ****/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "PoissonPolygon.h"
#include "convhull2D.h"

#include "RF.h"


/*******************************************************************************
 * double scProd(const double *x, const double *y)
 *
 * returns the scalar product of two two-dimensional vectors x and y
 ****/
double scProd(const double *x, const double *y)
{
	return x[0]*y[0]+x[1]*y[1];
}

/******************************************************************************* 
 * int compareAngles(const void * a, const void * b)
 *
 * returns -1 if a<b, 1 if a>b, or 0 if a=b
 ****/
int compareAngles(const void * a, const void * b)
{
	double diff = ( *(double*)a - *(double*)b );
	if(diff<0.0) return -1;
	if(diff>0.0) return 1;
	return 0;
}

/*******************************************************************************
 * void rTriangle(double *phi)
 *
 * returns three random angles from (0,2*PI) according to a joint distribution
 * which has a density w.r.t. the i.i.d. case proportional to the area of the 
 * triangle spanned by the corresponding points on the unit circle (multiplied 
 * by the indicator that this triangle contains the origin)
 *
 * needs a function 
 *   double runif(void)
 * which generates a uniform random number from (0,1)
 ****/
void rTriangle(double *phi)
{
	int onceagain=1, ok=0;
	double a, amax=1.2990381056766579701455847561294; 
	// amax = maximum area of a triangle with vertices on the unit circle 
	//      = 3*sqrt(3)/4
	while(onceagain)
	{
		ok = 0;
		while(!ok)
		{
			ok = 1;
			phi[0] = 2*PI*UNIFORM_RANDOM;
			phi[1] = 2*PI*UNIFORM_RANDOM;
			phi[2] = 2*PI*UNIFORM_RANDOM;
			qsort(phi, 3, sizeof(double), compareAngles);
			if(phi[2]-phi[0]<PI) ok=0;
			else if(phi[1]<phi[2]-PI || phi[1]>phi[0]+PI) ok=0;

		}
		a = 0.5*(fabs(sin(phi[2]-phi[1]))+fabs(sin(phi[0]-phi[2]))+fabs(sin(phi[1]-phi[0])));
		if(amax*UNIFORM_RANDOM<a) onceagain=0;
	}
}

/*******************************************************************************
 * int rPoissonPolygon(struct polygon *P, double lambda)
 * 
 * generates the typical cell of a stationary and isotropic Poisson line 
 * tessellation with intensity lambda>0, see
 * - Schneider, R. and Weil, W. (2008). Stochastic and Integral Geometry.
 *   Springer-Verlag, Berlin, Heidelberg, pp. 497f.
 * - Calka, P. (2009). Tessellations. Contributed chapter to New Perspectives in
 *   Stochastic Geometry. Edited by W. Kendall and I. Molchanov. Oxford 
 *   University Press
 *
 * returns an error code err
 * 
 * needs 
 * - a function 
 *     double runif(void)
 *   which generates a uniform random number from (0,1)
 * - a function 
 *     double rexp(double mu)
 *   which generates a random number from the exponential distribution with
 *   mean 1/mu
 * - a function 
 *     double rpois(double mu)
 *   which generates a random number from the Poisson distribution with mean mu
 *
 * ! no errors are catched so far
 ****/
int rPoissonPolygon(struct polygon *P, double lambda)
{
	double R,T,RRight,RLeft,RMax,uprim[2],pprim,udual[2],pdual;
	double sqrlen,ulen,maxlen;
	int n=0,nold=0,nneu,ok=0,i,j,nV=0;
	int err=0;
	// arrays
	double phi[3], psi, **vdual, **vdtemp, **vsave;
	struct vertex *vprim;

	// generate initial triangle
	T = rexp(2.0);
	n = 3;
	rTriangle(phi);
	vdual = (double **) CALLOC(n, sizeof(double *));
	vsave = (double **) CALLOC(n, sizeof(double *));
	for(i=0;i<n;i++){
		uprim[0] = cos(phi[i]);
		uprim[1] = sin(phi[i]);
		vdual[i] = (double *) CALLOC(2, sizeof(double));
		for(j=0;j<2;j++) vdual[i][j] = uprim[j]/T;
		vsave[i] = vdual[i];
	}
	nold = n;
	nV = 3;
	vprim = (struct vertex *) CALLOC(nV, sizeof(struct vertex));
	// determine vertices
	for(i=0;i<nV;i++){
		j = (i+1) % nV;
		udual[0] = vdual[j][1]-vdual[i][1];
		udual[1] = vdual[i][0]-vdual[j][0];
		ulen = sqrt(scProd(udual,udual));
		for(j=0;j<2;j++) udual[j] /= ulen;
		pdual = scProd(udual,vdual[i]);
		for(j=0;j<2;j++) vprim[i].x[j] = udual[j]/pdual;
	}

	// determine distance of the farthest vertex
	RMax = T/cos(0.5*(phi[1]-phi[0]));
	R = T/cos(0.5*(phi[2]-phi[1]));
	if(R>RMax) RMax = R;
	R = T/cos(PI+0.5*(phi[0]-phi[2]));
	if(R>RMax) RMax = R;

	// rest of the line process
	R = 10.0; // seems to be a good initial value
	RLeft = T-R;
	RRight = T;
	while(!ok){
		RLeft += R;
		if(RRight+R>RMax){
			ok = 1;
			R = RMax-RRight;
		}
		RRight += R;
		
		nneu = rpois(2.0*R);		
		if(nneu>0){
			if (vprim != NULL) free(vprim); 
			n += nneu;
			vdtemp = (double **) realloc(vdual, (n+1)*sizeof(double *));
			if(vdtemp != NULL) vdual = vdtemp; // else 'error handling'
			vdtemp = (double **) realloc(vsave, n*sizeof(double *));
			if(vdtemp != NULL) vsave = vdtemp; // else 'error handling'
			for(i=nold;i<n;i++){
				psi = 2*PI*UNIFORM_RANDOM;
				uprim[0] = cos(psi);
				uprim[1] = sin(psi);
				pprim = RLeft+R*UNIFORM_RANDOM;
				vdual[i] = (double *) CALLOC(2, sizeof(double));
				for(j=0;j<2;j++) vdual[i][j] = uprim[j]/pprim;
				vsave[i] = vdual[i];
			}
			nV = ch2d(vdual,n);
			nold = n;
			vprim = (struct vertex *) CALLOC(nV, sizeof(struct vertex));
			// determine vertices
			for(i=0;i<nV;i++){
				j = (i+1) % nV;
				udual[0] = vdual[j][1]-vdual[i][1];
				udual[1] = vdual[i][0]-vdual[j][0];
				ulen = sqrt(scProd(udual,udual));
				for(j=0;j<2;j++) udual[j] /= ulen;
				pdual = scProd(udual,vdual[i]);
				for(j=0;j<2;j++) vprim[i].x[j] = udual[j]/pdual;
			}
			// check whether the farthest vertex is closer than RRight
			maxlen = 0.0;
			for(i=0;i<nV;i++){
				sqrlen = scProd(vprim[i].x,vprim[i].x);
				if(sqrlen>maxlen) maxlen = sqrlen;
			}
			if(maxlen<RRight*RRight) ok = 1;
		}
	}
	// free vdual and vsave
	for(i=0;i<n;i++) if(vsave[i] != NULL) free(vsave[i]);
	if (vdual != NULL) free (vdual);
	if (vsave != NULL) free(vsave);

	// determine the final polygon
	P->n = nV; 
	P->box0[0] = 0.0;
	P->box0[1] = 0.0;
	P->box1[0] = 0.0;
	P->box1[1] = 0.0;
	P->v = (struct vertex *) CALLOC(P->n, sizeof(struct vertex));
	for(i=0;i<P->n;i++){
		for(j=0;j<2;j++){
			P->v[i].x[j] = vprim[i].x[j]/lambda;
			if(P->v[i].x[j]<P->box0[j]) P->box0[j] = P->v[i].x[j];
			if(P->v[i].x[j]>P->box1[j]) P->box1[j] = P->v[i].x[j];
		}
	}
	P->e = (struct edge *) CALLOC(P->n, sizeof(struct edge));
	for(i=0;i<P->n;i++) {
		j = (i+1) % nV;
		P->e[i].u[0] = P->v[j].x[1]-P->v[i].x[1];
		P->e[i].u[1] = P->v[i].x[0]-P->v[j].x[0];
		ulen = sqrt(scProd(P->e[i].u,P->e[i].u));
		P->e[i].u[0] /= ulen;
		P->e[i].u[1] /= ulen;
		P->e[i].p = scProd(P->e[i].u,P->v[i].x);
	}
	
	// free vprim
	if (vprim != NULL) free(vprim);
	return err;
}

/*******************************************************************************
 * void freePolygon(struct polygon *P)
 *
 * frees allocated memory
 ****/
void freePolygon(struct polygon *P)
{
	if (P->e != NULL) free (P->e); 
	if (P->v != NULL) free (P->v); 
}

/*******************************************************************************
 * int isInside(struct polygon *P, double *x)
 *
 * returns 1 if the two-dimensional vectors x is inside the polygon P, else 0
 ****/
int isInside(struct polygon *P, double *x)
{
	int i, ok=1;
	for(i=0;i<P->n;i++)
		if(scProd(x,P->e[i].u) > P->e[i].p){
			ok = 0;
			break;
		}
	return ok;
}

/*******************************************************************************
 * double getArea(struct polygon *P)
 *
 * returns the area of the polygon P
 ****/
double getArea(struct polygon *P)
{
	int i, j;
	double a=0.0, dx, dy;
	for(i=0;i<P->n;i++){
		j = (i+1) % P->n;
		dx = P->v[i].x[0]-P->v[j].x[0];
		dy = P->v[i].x[1]-P->v[j].x[1];
		a += 0.5*P->e[i].p*sqrt(dx*dx+dy*dy);
	}
	return a;
}
