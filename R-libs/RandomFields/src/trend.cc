
/*
 Authors 
 Martin Schlather, schlather@math.uni-mannheim.de

 Handling of the different possibilities to pass the trend

Note:
 * Never use the below functions directly, but only by the functions indicated 
   in RFsimu.h, since there is no error check (e.g. initialization of RANDOM)
 * VARIANCE, SCALE are not used here 
 * definitions for the random coin method can be found in MPPFcts.cc
 * definitions for genuinely anisotropic or nondomain models are in
   SophisticatedModel.cc; hyper models also in Hypermodel.cc

 Copyright (C) 2011 - 2013 Marco Oesting & Martin Schlather

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 3
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.  
*/

/*Note:
the parameter 'polycoeff' is hidden in R; it is a vector consisting of the
coefficients of the corresponding trend polynomials:
the first choose(polydeg[1]+d,d) components belong to the first polynomial,
the next choose(polydeg[2]+d,d) to the second one,...
the corresponding monomial functions are of the following order:
1, x, x^2, ..., x^k, y, x y, ... x^(k-1) y, y^2, x y^2..., y^k,
z, x z, x^2 z, ...., x^(k-1) z, y z, x y z, x^2 y z, ..., z^k
*/


#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <math.h>
 
#include "RF.h"
#include "primitive.h"
#include <R_ext/Lapack.h>
#include <R_ext/Linpack.h>
#include "Covariance.h"



//////////////////////////////////////////////////////////////////////
//    mixed
//////////////////////////////////////////////////////////////////////

#define MIXED_CONSTANT 0

void mixed(double *x, cov_model *cov, double *v) {
  // bis auf den letzten Schluss alles schon auf bel vdim geschrieben.
  cov_model  *next = cov->sub[0];
  location_type *loc = Loc(cov);
  mixed_storage *s = cov->Smixed;
  listoftype *X = (listoftype*) (cov->p[MIXED_X]);
  double //distTF, 
    var = 1.0, 
    *covmatrix=NULL;
  int i, err, // Xnrow, 
    //    err=NOERROR,
    element = ((int*) cov->p[MIXED_ELMNT])[0],
    vdim = cov->vdim,
    vdimsq = vdim * vdim;
  
  error("mixed not programmed yet");

 if (cov->nsub == 0) {
    for (i=0; i<vdimsq; i++) v[i] = 0.0;
    return;
  }
  if (X == NULL) { // no X given, but a submodel	
    COV(x, next, v);
    return;
  }
    
  if (element < 0 || element >= cov->nrow[MIXED_X]) { 
    // realex element <0 ??
    GERR1("Illegal call of the internal function '%s'", NICK(cov)); 
  }

  // PMI(cov->calling);
  // crash(cov);
  assert(false); // hier irgendwo stehen geblieben

  if (cov->q[MIXED_CONSTANT]) {
    // Xu , u ~ N(0, C)
    // classical random effect
    cov_model *sub = next;
    while (isDollar(sub)) {
      var *= sub->p[DVAR][0];
      sub=sub->sub[0];
    }
    listoftype *list= (listoftype*) (next->p[CONSTANT_M]);
    covmatrix = list->p[element];
    if (X->ncol[element] != list->nrow[element]) {
      GERR("mixed model: X and covariance matrix M do not match");
    }
  } else {    
    // Xu , u ~ random field
    // random effect, geostatistical covariance model

    CovarianceMatrix(cov->key, s->mixedcov);
    covmatrix = s->mixedcov;
  } // sub->nr != CONSTANT

  // ab hier vdim =1  
  *v = XkCXtl(X->p[element], covmatrix, 
	      X->nrow[element], X->ncol[element], 
	      loc->i_row, loc->i_col);
  *v *= var;
  return;

 ErrorHandling:
  XERR(err);
  
}

void mixed_nonstat(double *x, double *y, cov_model *cov, double *v){
  int 
    element = ((int*) cov->p[MIXED_ELMNT])[0];

  error("mixed_constant not programmed yet");
  
  if (element < 0 || element >= cov->nrow[MIXED_X]) { 
    // relax element <0 ??
    char msg[100];
    sprintf(msg, "Illegal call of the internal function '%s'", 
	    NICK(cov)); 
    error(msg);
  }


   if (cov->nsub != 0 && cov->p[MIXED_X] == NULL) { // no X given, but a submodel	
    cov_model *next = cov->sub[0];
    NONSTATCOV(x, y, next, v);
  } else mixed(x, cov, v);
}

void covmatrix_mixed(cov_model *cov, double *v) {
  cov_model *sub = cov->sub[0];
  int element = ((int*) cov->p[MIXED_ELMNT])[0];
  
  if (cov->ncol[MIXED_X] == 0 || element < 0) {
    CovList[sub->nr].covmatrix(sub, v);
    return;
  }
  if (element >= cov->nrow[MIXED_X]) ERR("value of 'element' is too large");

  listoftype *X = (listoftype*) (cov->p[MIXED_X]);
  double *C=NULL;
  int nrow = X->nrow[element],
    ncol=X->ncol[element]; 
  C = (double*) MALLOC(sizeof(double) * ncol * ncol);
  if (C==NULL) {
    StandardCovMatrix(cov, v);
    return;
  }
  
  CovList[sub->nr].covmatrix(sub, C);
  XCXt(X->p[element], C, v, nrow, ncol);
  Loc(cov)->totalpoints = nrow;

  // PMI(cov);
  
  /*
  int i,j,k, tot=Loc(cov)->totalpoints;
  printf("\nStart mixed C %d  ---- check within comment\n", element); //
  for (k=i=0; i<ncol * ncol; i+=tot) {
     for (j=0; j<ncol; j++) printf("%f ", C[k++]); //
     printf("\n");//
  }

  printf("\nStart mixed t(X) %d\n", element); //
  for (k=i=0; i<tot*tot; i+=tot) {
  for (j=0; j<tot; j++) printf("%f ", X->p[element][k++]);//
   printf("\n");//
  }

  printf("\nStart mixed t(v) %d\n", element); //
    for (k=i=0; i<tot*tot; i+=tot) {
    for (j=0; j<tot; j++) printf("%f ", v[k++]);//
    printf("--- end check within comment \n");//
  }
  */
  
  free(C);
}

int set_mixed_constant(cov_model *cov) {
  cov_model 
    *next = cov->sub[0],
    *sub = next;   
  bool simple = true;
  location_type *loc=Loc(cov);
  int i,
    totalpoints = loc->totalpoints,
    *ncol = cov->ncol,
    *nrow = cov->nrow; 
  listoftype *X = (listoftype*) (cov->p[MIXED_X]);

  cov->q = (double*) MALLOC(sizeof(double) * 1);
  cov->qlen = 1;

  // PMI(cov);

  while (sub != NULL && isDollar(sub) &&
	 ((simple = sub->p[DPROJ] == NULL 
	   && (sub->p[DSCALE] == NULL || sub->p[DSCALE][0] == 1.0)
	   && sub->p[DANISO] == NULL))) sub=sub->sub[0];
  
  if ((cov->q[MIXED_CONSTANT] = sub != NULL && sub->nr == CONSTANT)) {
    //next->delflag = DEL_COV - 6;
    if (isDollar(next) && next->nrow[DVAR]==0) {
      //next->delflag = DEL_COV - 6;
      if (!simple) 
	SERR1("'%s' not allowed together with an anisotropic structrue",
	     NICK(cov));
    }
    
    int constant_size;
    
    for (i=0; i<nrow[MIXED_X]; i++) {
      constant_size = ((listoftype*) (sub->p[CONSTANT_M]))->nrow[i];
      
      if (ncol[MIXED_X] > 0 && X->ncol[i] != constant_size) {
	SERR5("%dth matrix 'X' (%d x %d) and (%d x %d) constant matrix 'M' do not match", i, X->nrow[i], X->ncol[i], constant_size, constant_size);
      }
    }
  } else {
    for (i=0; i<nrow[MIXED_X]; i++) {
      if (X->nrow[i] != X->ncol[i])
	SERR3("%dth  matrix is not symmetric (%d x %d)",
	      i+1, X->nrow[i], X->ncol[i]);
    }
  }
  
  if (false) // test ist wichtig, kollidiert z.Zt. aber mit SetAndGetModelInfo
    for (i=0; i<nrow[MIXED_X]; i++) {
      if (X->nrow[i] != totalpoints) {
	//PMI(cov);
	SERR3("number of rows of entry %d of 'X' (%d) are different from the number of locations (%d)", i+1, X->nrow[i], totalpoints);
      }
    }
  return NOERROR;
}

char iscovmatrix_mixed(cov_model *cov) {
  int err;
  // printf("iscov %d %ld\n", cov->qlen, cov->q);
  if (cov->q == NULL && (err = set_mixed_constant(cov)) != NOERROR) XERR(err);
  // printf("iscov %d %ld\n", cov->qlen, cov->q);
  return 2 * (int) (cov->nsub > 0) * (int) (cov->q[MIXED_CONSTANT] || 
					    cov->ncol[MIXED_X] > 0);
}

void kappamixed(int i, cov_model  VARIABLE_IS_NOT_USED *cov, int *nr, int *nc){
  if (i==MIXED_DIM || i==MIXED_ELMNT) *nc = *nr = 1;
  else if (i==MIXED_BETA || i==MIXED_DIST) {
    *nc = 1;
    *nr = SIZE_NOT_DETERMINED; 
  } else if (i==MIXED_X || i==MIXED_COORD) *nc = *nr = SIZE_NOT_DETERMINED; 
  else *nc = *nr = -1;
}
int checkmixed(cov_model *cov) {
    // cov_model *sub;
  //location_type *loc=Loc(cov);
  int i, err, 
    //totalpoints = loc->totalpoints,
    nsub = cov->nsub,
    *ncol = cov->ncol,
    *nrow = cov->nrow; // taken[MAX DIM],
  char msg[300];
  listoftype *X = (listoftype*) (cov->p[MIXED_X]);
  
  ROLE_ASSERT(ROLE_GAUSS || cov->role == ROLE_COV);

  cov->vdim=1; //falls kein submodel vorhanden (Marco)
  cov->maxdim=INFDIM;
  cov->matrix_indep_of_x = true;

  kdefault(cov, MIXED_ELMNT, 0);

  if (ncol[MIXED_BETA] > 0) { // b is given  
    if (nsub != 0) 
      SERR("b and a covariance model may not be given at the same time");
    if (ncol[MIXED_X] == 0) SERR("if 'b' is given 'X' must be given");
    for (i=0; i<nrow[MIXED_X]; i++) {
      if (X->ncol[i] != nrow[MIXED_BETA]) {
	sprintf(msg,
       	"%dth set: (%d x %d) matrix X and (%d x %d) vector b do not match",
	i, X->nrow[0], X->ncol[i], nrow[MIXED_BETA], ncol[MIXED_BETA]);
	SERR(msg);
      }
    }
  } else if (nsub == 0) { // only X is given -- then only a deterministic 
	//                                 constant is given
    if (ncol[MIXED_BETA] == 0) 
      SERR("if no covariance model is given then 'b' must be given");
    if (ncol[MIXED_X] != 1) // deterministic effect
      SERR("X must have one column");
    kdefault(cov, MIXED_BETA, 1);
  } else { // submodel is given
    cov_model 
      *next = cov->sub[0];   
    //    double var = 1.0;

    if (cov->tsdim != cov->xdimprev || cov->tsdim != cov->xdimown) 
      return ERRORDIM;
    if ((err = CHECK(next, cov->tsdim, cov->xdimown, PosDefType,
		       cov->domown, 
                       cov->isoown, SUBMODEL_DEP, ROLE_COV)) != NOERROR) {
        // print("error\n");
      return err;
    }

    if (cov->q == NULL && (err = set_mixed_constant(cov)) != NOERROR) {
      return err;
    }
    
    // warning("some checks in model `mixed' are missing");
    // ob X mit C zusammengeht.
    
    setbackward(cov, next);
  }

  if (cov->vdim > 1) 
    SERR("multivariate version of mixed not programmed yet");

  if ((cov->p[MIXED_DIST] == NULL) xor (cov->p[MIXED_DIM] == NULL))
    SERR("if 'dim' and 'dist' must be given at the same time");
  if ((cov->p[MIXED_DIST] != NULL) xor (cov->p[MIXED_COORD] != NULL))
    SERR("'dist' and 'coord' may not be given together");

  // incorrect. but save !!
  cov->semiseparatelast = false; // taken[tsxdim - 1] <= 1;
  cov->separatelast = false;     // taken[tsxdim - 1] <= 1; ?? 
  return NOERROR;
}


  void rangemixed(cov_model  VARIABLE_IS_NOT_USED *cov, range_type *range){
  int i;

  range->min[MIXED_ELMNT] = 0;
  range->max[MIXED_ELMNT] = MAXELEMENTS;
  range->pmin[MIXED_ELMNT] = 0;
  range->pmax[MIXED_ELMNT] = MAXELEMENTS;
  range->openmin[MIXED_ELMNT] = false;
  range->openmax[MIXED_ELMNT] = false;


  for (i=MIXED_X; i<=MIXED_COORD; i++) {
    range->min[i] = RF_NEGINF;
    range->max[i] = RF_INF;
    range->pmin[i] = -1e10;
    range->pmax[i] = 1e10;
    range->openmin[i] = true;
    range->openmax[i] = true;
  }

  i=MIXED_DIST;
  range->min[i] = 0;
  range->max[i] = RF_INF;
  range->pmin[i] = 1e-10;
  range->pmax[i] = 1e10;
  range->openmin[i] = false;
  range->openmax[i] = true;

}


int initmixed(cov_model *cov, storage  VARIABLE_IS_NOT_USED *S) {
  location_type *loc = Loc(cov);
  mixed_storage *s;
  char errorloc_save[nErrorLoc];
  double 
    *b = cov->p[MIXED_BETA],
    *coord = cov->p[MIXED_COORD],
    *dist = cov->p[MIXED_DIST];
  int 
    Cn = -1, 
    Cdim = -1, 
    list_element=0,
    err = NOERROR,
    vdim = cov->vdim,
    dim = coord!=NULL ? cov->ncol[MIXED_COORD] : ((int*) cov->p[MIXED_DIM])[0],
    totalpoints = loc->totalpoints,
    total = vdim * totalpoints;

  listoftype 
    *X = (listoftype*) (cov->p[MIXED_X]);
  bool distTF;


  return ERRORFAILED; // muss in zerlegt werden in init und struct
  // und unten struct aufrufen richtig praeparieren.


  ROLE_ASSERT_GAUSS;
  
  // cholesky zerlegung + bereitstellung von b

  strcpy(errorloc_save, ERROR_LOC);
  sprintf(ERROR_LOC, "%s%s: ", errorloc_save, "init mixed model");
  
  assert(cov->nr = MIXEDEFFECT);
 

  // submodel exists:
  if ((cov->Smixed = (mixed_storage*) MALLOC(sizeof(mixed_storage)))==NULL){
      err=ERRORMEMORYALLOCATION; goto ErrorHandling;
  }  
  s = cov->Smixed;
  MIXED_NULL(s);
  assert(false);

  if (s->Xb != NULL) free(s->Xb);
  if (cov->ncol[MIXED_BETA] > 0) { // b is given
    // X is given, but no covariance model
    if (cov->nrow[MIXED_X] > 1) {
      warning("using first element of list X only");
    }
    
    if ((s->Xb = (double*) MALLOC(sizeof(double) * X->nrow[list_element])) == 
	NULL) { err=ERRORMEMORYALLOCATION; goto ErrorHandling; }  
   
    Ax(X->p[list_element], b, X->nrow[list_element], X->ncol[list_element], 
       s->Xb);

  } else { // submodel is given

    assert(false) ; //if (cov->q[MIXED_CONSTANT]) error("not 'constant'").

    cov_model *sub, *next = cov->sub[0];
    
    if ((s->Xb = (res_type*) MALLOC(sizeof(res_type) * total)) == 0) {
      err=ERRORMEMORYALLOCATION; goto ErrorHandling;
    }  
 
    if ((distTF = dist!=NULL)) {
      // Xu , u ~ geostatmodel, u i.a. viel laenger als X (tierzucht)
      // stehende Vektoren, aneinandergereiht
      Cn = cov->ncol[MIXED_DIST];
      Cn = (int) (0.5 + 0.5 * (1.0 + sqrt(1 + 8.0 * Cn)));
      Cdim = cov->nrow[MIXED_DIST];
      assert(Cdim==1);
      if ((err = covcpy(&(cov->key), next, dist, NULL, dim, dim, Cn, 
			false /*Time */,
			false, distTF)) != NOERROR) goto ErrorHandling;
    } else if (coord != NULL) {
      // Xu , u ~ geostatmodel, u i.a. viel laenger als X (tierzucht)
      // stehende Vektoren, aneinandergereiht
      Cn = cov->ncol[MIXED_COORD];
      Cdim = cov->nrow[MIXED_COORD];
      
      if (Cdim > Cn)
	warning("The dimension of the coordinates is higher than the number of points");
      if ((err = covcpy(&(cov->key), next)) != NOERROR) goto ErrorHandling;
    } else {
      // Xu , u ~ geostatmodel, X quadratisch	
      Cn = loc->totalpoints;
      if ((err = covcpy(&(cov->key), next, coord, NULL, dim, dim, Cn, false, 
			false, distTF)) != NOERROR) goto ErrorHandling;
    }
    if (cov->key->nr != GAUSSPROC) addModel(&(cov->key), GAUSSPROC);   
    cov->key->calling = cov->key;
    cov->stor = (storage *) MALLOC(sizeof(storage));
    STORAGE_NULL(cov->stor);
    if ((err = INIT(cov->key, 0, cov->stor)) != NOERROR) goto ErrorHandling;
   
    int Xnrow, Xncol;
    Xnrow = Xncol = Cn * sub->vdim;
    if (loc->i_row==0 && loc->i_col==0) {
      if (s->mixedcov != NULL) free(s->mixedcov);   
      s->mixedcov = (double*) MALLOC(sizeof(double) * Xnrow * Xnrow);
      if (s->mixedcov == NULL) {
	err = ERRORMEMORYALLOCATION;
	goto ErrorHandling;
      }
    }

  } // end of submodel
 
  cov->initialised = true;
  FieldReturn(cov);

 ErrorHandling: 
  return err;
}


static int keeprandomeffect = false;
void domixed(cov_model *cov, storage  VARIABLE_IS_NOT_USED *S){
  location_type *loc = Loc(cov);
  mixed_storage *s = cov->Smixed;
  double *res  = cov->rf;
  int i,
    list_element=0,
    vdim = cov->vdim,
    totalpoints = loc->totalpoints,
    total = vdim * totalpoints;
 
  listoftype 
    *X = (listoftype*) (cov->p[MIXED_X]);

  if (cov->ncol[MIXED_BETA] > 0) { // b is given
    // X is given, but no covariance model
    if (total == X->nrow[list_element]) {
      for (i=0; i<total; i++) res[i] = s->Xb[i];
    } else {
      assert(X->nrow[list_element] == 1);
      for (i=0; i<total; i++) res[i] = s->Xb[0];
    }
     //     print("%f\n", res[0]); assert(false);
    
  } else { // submodel is given
    cov_model *key = cov->key;
    if (!keeprandomeffect || !s->initialized) {
      do_gaussprocess(key, cov->stor);
    }
     
    if (X != NULL) {
      AxResType(X->p[list_element], cov->key->rf, X->nrow[list_element], 
		X->ncol[list_element], res); 
    } else {
      res_type *rf = cov->key->rf;
      for (i=0; i<total; i++) res[i] = rf[i];
    } 
  }
}
    






//////////////////////////////////////////////////////////////////////
//    trend
//////////////////////////////////////////////////////////////////////

void trend(double  VARIABLE_IS_NOT_USED *x, cov_model *cov, double *v){
  // todo : ueberlegen ob hier nicht der Trend-Wert berechnet werden sollte

  ///BUG; // sollte durch plusStat & plusNonStat abgeblock sein ?
  // to do: stattdessen sollte hier der Trendwert zurueckgeliefert werden ?
  // nein: model kann auch nur aus Trend bestehen!

  int i,
    vSq = cov->vdim * cov->vdim;
  for (i=0; i<vSq; i++) v[i]=0.0;
  //print("trend %d %f\n", vSq, *v);
}

void trend_nonstat(double  VARIABLE_IS_NOT_USED *x, double  VARIABLE_IS_NOT_USED *y, cov_model *cov, double *v){
  int i,
    vSq = cov->vdim * cov->vdim;
 
  if (cov->role == ROLE_COV)  for (i=0; i<vSq; i++) v[i]=0.0;

  else ERR("trend is called unexpectately.");
}


void kappatrend(int i, cov_model *cov, int *nr, int *nc){
  // i nummer des parameters
  int k;
  
  switch(i) {
    case TREND_MEAN: //mu
      *nr = SIZE_NOT_DETERMINED; 
      *nc = 1;
    break;
  
    case TREND_LINEAR: //plane
      *nr = cov->tsdim;
      *nc = SIZE_NOT_DETERMINED;
    break;
    
    case TREND_POLY: //polydeg
      *nr = SIZE_NOT_DETERMINED;
      *nc = 1;
     break;
    
    case TREND_PARAM_POLY: //polycoeff
      if(cov->p[TREND_POLY] == NULL) {
        *nr = -1; 
      } else {
        *nr = 0;
	
	//APMI(cov);
        for(k=0; k < cov->nrow[TREND_POLY]; k++) {
	  //printf("%d\n", *nr);
	  *nr += binomialcoeff(((int*) cov->p[TREND_POLY])[k] + 
			      cov->tsdim, cov->tsdim);
	}
      }
      *nc = 1;
    break;
    
    case TREND_FCT: //arbitraryfct
      *nr = 1;
      *nc = 1;
    break;
    
    case TREND_PARAM_FCT: //fctcoeff
      *nr = 1;
      *nc = 1;
    break;

    default:
      *nr = -1; 
      *nc = -1;
  }
}


int checktrend(cov_model *cov){
  // if (cov->ncol[TREND_LINEAR] > 0 || cov->ncol[TREND_FCT]>0
  //    || cov->ncol[TREND_PARAM_FCT] > 0)
  //  return(ERRORNOTPROGRAMMED);
  double *mu = cov->p[TREND_MEAN],
         *plane = cov->p[TREND_LINEAR],
         *polydeg = cov->p[TREND_POLY],
         *arbitraryfct = cov->p[TREND_FCT];
  int i, 
    vdim = 0, 
    tsdim= cov->tsdim,
    basislen = 0;
  //SEXP Rx, fctbody, envir;

  ROLE_ASSERT(ROLE_GAUSS || cov->role == ROLE_COV);
 
  if (mu != NULL) {
    vdim = cov->nrow[TREND_MEAN];
    cov->matrix_indep_of_x = true;
  }
  
  if (plane != NULL) {
    if(vdim>0 && (vdim != cov->ncol[TREND_LINEAR])) {
      SERR("trend parameters have different multivariate dimensions");
    } else vdim = cov->ncol[TREND_LINEAR];
    cov->matrix_indep_of_x = false;
  }
  
  if (polydeg != NULL) {
    if (vdim>0) {
      if (vdim != cov->nrow[TREND_POLY])
	SERR("trend parameters have different multivariate dimensions");
      SERR("polynomials and free functions in trend may not be mixed with other trend definitions. Please use a sum of trends.");
    }
    vdim = cov->nrow[TREND_POLY];
    if (cov->p[TREND_PARAM_POLY] == NULL) {
      for (i=0; i<vdim; i++) 
	basislen += binomialcoeff(((int*) polydeg)[i] + tsdim, tsdim);
      if ((cov->p[TREND_PARAM_POLY] = 
	   (double *) MALLOC(sizeof(double) * basislen))==NULL) {
	return ERRORMEMORYALLOCATION;
      }
      for(i=0; i<basislen; i++) (cov->p[TREND_PARAM_POLY])[i] = NA_REAL;
      cov->ncol[TREND_PARAM_POLY] = 1;
      cov->nrow[TREND_PARAM_POLY] = basislen;
    }
    cov->matrix_indep_of_x = false;
  }
  
  if (arbitraryfct != NULL) { //brauche hier: construct.fct$vdim
    cov->matrix_indep_of_x = false;

    error("arbitrary function not programmed yet");

    assert(false);
//     kdefault(cov, TREND_PARAM_FCT, 1.0);
//     PROTECT(envir = allocSExp(ENVSXP));
//     SET_ENCLOS(envir, R_GlobalEnv);
//     
//     PROTECT(fctbody = allocSExp(CLOSXP));
//     fctbody = BODY(*((SEXP *) arbitraryfct));
//     
//     PROTECT(Rx = allocVector(REALSXP,xlen));
//     for(j=0; j<xlen; j++) REAL(Rx)[j] = 0;
//     defineVar(install("x"), Rx, envir);
//     defineVar(install("y"), ScalarReal(0), envir);
//     defineVar(install("z"), ScalarReal(0), envir);
//     defineVar(install("T"), ScalarReal(0), envir);
//     
//     int vdimproposal = length(eval(fctbody, envir));       
//     UNPROTECT(3);
//     
//     if (vdim > 0) {
//       if (vdim != vdimproposal)
// 	GERR("trend parameters have different multivariate dimensions");
//       GERR("polynomials and free functions in trend may not be mixed with other trend definitions. Please use a sum of trends.");
//     }
//     vdim = vdimproposal;
  }
  
  if (vdim <= 0) {
    vdim = cov->calling->vdim;
    if (vdim <= 0) 
      SERR("multivariate dimension for trend cannot be determined.");
    if ((cov->p[TREND_MEAN] =
	 (double *) MALLOC(vdim*sizeof(double)))==NULL) {
      return ERRORMEMORYALLOCATION;
    }
    for(i=0; i<vdim; i++) (cov->p[TREND_MEAN])[i] = 0.0;
    cov->nrow[TREND_MEAN] = vdim;
    cov->ncol[TREND_MEAN] = 1;
  }
  cov->vdim = vdim;
  cov->isoown = cov->matrix_indep_of_x ? ISOTROPIC : NO_ROTAT_INV;
  
 return NOERROR;

}


void rangetrend(cov_model  VARIABLE_IS_NOT_USED *cov, range_type *range){
  //cov->p[TREND_MEAN]: mu / mean
  
  range->min[TREND_MEAN] = RF_NEGINF;
  range->max[TREND_MEAN] = RF_INF;
  range->pmin[TREND_MEAN] = -10^10;
  range->pmax[TREND_MEAN] = 10^10;
  range->openmin[TREND_MEAN] = true;
  range->openmax[TREND_MEAN] = true;
  
  //cov->p[TREND_LINEAR]: plane
  range->min[TREND_LINEAR] = RF_NEGINF;
  range->max[TREND_LINEAR] = RF_INF;
  range->pmin[TREND_LINEAR] = -10^10;
  range->pmax[TREND_LINEAR] = 10^10;
  range->openmin[TREND_LINEAR] = true;
  range->openmax[TREND_LINEAR] = true;
  
  //cov->p[TREND_POLY]: polydeg / polynomial degree
  range->min[TREND_POLY] = 0;
  range->max[TREND_POLY] = RF_INF;
  range->pmin[TREND_POLY] = 0;
  range->pmax[TREND_POLY] = 10;
  range->openmin[TREND_POLY] = false;
  range->openmax[TREND_POLY] = false;
  
  //cov->p[TREND_PARAM_POLY]: polycoeff / coefficients of polynomial
  range->min[TREND_PARAM_POLY] = RF_NEGINF;
  range->max[TREND_PARAM_POLY] = RF_INF;
  range->pmin[TREND_PARAM_POLY] = -10^10;
  range->pmax[TREND_PARAM_POLY] = 10^10;
  range->openmin[TREND_PARAM_POLY] = true;
  range->openmax[TREND_PARAM_POLY] = true;
 
  //cov->p[TREND_FCT]: arbitraryfct / arbitrary function
  range->min[TREND_FCT] = RF_NEGINF;
  range->max[TREND_FCT] = RF_INF;
  range->pmin[TREND_FCT] = -10^10;
  range->pmax[TREND_FCT] = 10^10;
  range->openmin[TREND_FCT] = true;
  range->openmax[TREND_FCT] = true;
  
  //cov->p[TREND_PARAM_FCT]: fctcoeff / coefficient of arbitrary function
  range->min[TREND_PARAM_FCT] = RF_NEGINF;
  range->max[TREND_PARAM_FCT] = RF_INF;
  range->pmin[TREND_PARAM_FCT] = -10^10;
  range->pmax[TREND_PARAM_FCT] = 10^10;
  range->openmin[TREND_PARAM_FCT] = true;
  range->openmax[TREND_PARAM_FCT] = true;
}

				
double GetInternalMean(cov_model *cov){
  if (cov->nr == TREND) {
    if (cov->ncol[TREND_MEAN]==1) {
      if (cov->nrow[TREND_MEAN] != 1) {
	return(RF_NAN); // only scalar allowed !
      }
      return(cov->p[TREND_MEAN][0]);
    } 
  }
  double sum=0.0;
  if (cov->nr == PLUS || cov->nr == TREND) {
    int i;
    for (i=0; i<cov->nsub; i++)
      sum += GetInternalMean(cov->sub[i]);
  }
  return sum;
}



int init_trend(cov_model *cov, storage *S) {
  
  long err = NOERROR;
  trend_storage *s;
  double *polydeg =  cov->p[TREND_POLY],
         *arbitraryfct = cov->p[TREND_FCT];
  int i, basislen=0,
      tsdim = cov->tsdim,
      vdim = cov->vdim;
  //SEXP fctformals, argnames;

  //assert(false);

  ROLE_ASSERT_GAUSS;
   
  if(polydeg != NULL) {
    for(i=0; i<vdim; i++) 
      basislen += binomialcoeff(((int*) polydeg)[i]+tsdim,tsdim);
  }

  if (cov->Strend != NULL) free(cov->Strend);
  if ((cov->Strend = (trend_storage *) MALLOC(sizeof(trend_storage)))==NULL) {
    err=ERRORMEMORYALLOCATION;   
    goto ErrorHandling;
  }
  s = cov->Strend;
  TREND_NULL(s);

  if ((s->x = (double *) MALLOC(sizeof(double) * tsdim))==NULL ||
      (s->xi = (int *) MALLOC(sizeof(int) * tsdim))==NULL ||
      (s->evalplane = (double *) MALLOC(sizeof(double) * vdim))==NULL || 
      (basislen > 0 && 
       (s->powmatrix = (int *) MALLOC(sizeof(int) * basislen * tsdim))==NULL)){
    err=ERRORMEMORYALLOCATION;   
    goto ErrorHandling;
  }
  
  if (basislen > 0) poly_basis(cov, S); //generates basis of monomials
  //each row consists of one basis element
  //the j-th column consists of the power of the j-th space-time-dimension
  
  
  if (arbitraryfct != NULL) { //hier werden Argumente von arbitraryfct ueberprueft
    assert(false);
//       fctformals = getAttrib(FORMALS(*((SEXP *) arbitraryfct)), R_NamesSymbol);
//       nargs = length(fctformals);
//       PROTECT(argnames = allocVector(STRSXP,1));
//       
//       SET_STRING_ELT(argnames,0,mkChar("y"));
//       matchind = 0;
//       for(j=0; j<nargs; j++) 
// 	 matchind += (STRING_ELT(argnames,0) == STRING_ELT(fctformals,j));
//       if (matchind>0) {
//         if((meth->loc->spatialdim < 2) || (meth->loc->xvectorvalued))
// 	  GERR("The variable y does not match to the locations.\n");
//       }
//       
//       SET_STRING_ELT(argnames,0,mkChar("z"));
//       matchind = 0;
//       for(j=0; j<nargs; j++) 
// 	 matchind += (STRING_ELT(argnames,0) == STRING_ELT(fctformals,j));
//       if (matchind>0) {
// 	if((meth->loc->spatialdim < 3) || (meth->loc->xvectorvalued))
// 	  GERR("The variable z does not match to the locations.\n")	
//       }
//       
//       SET_STRING_ELT(argnames,0,mkChar("T"));
//       matchind = 0;
//       for(j=0; j<nargs; j++) 
// 	 matchind += (STRING_ELT(argnames,0) == STRING_ELT(fctformals,j));
//       if (matchind>0) {
//         if (meth->loc->Time == false) 
// 	  GERR("The variable T may be used for time only.\n")		  
//       }
//       UNPROTECT(1);
  }

  err = FieldReturn(cov);

  return NOERROR;
  
  ErrorHandling:
   return err;
}


void do_trend(cov_model *cov, storage  VARIABLE_IS_NOT_USED *s){
  location_type *loc = Loc(cov);
  char errorloc_save[nErrorLoc];
  trend_storage *S = cov->Strend;
  //PMI(cov->calling->calling);
  assert(S != NULL);
  double t,
    *mu = cov->p[TREND_MEAN],
    *plane   = cov->p[TREND_LINEAR],
    *polydeg = cov->p[TREND_POLY],
    *polycoeff = cov->p[TREND_PARAM_POLY],
    *arbitraryfct = cov->p[TREND_FCT],
    //    *fctcoeff = cov->p[TREND_PARAM_FCT],
    **xgr = loc->xgr,
    *x = S->x,
    *evalplane = S->evalplane;
  int i, j, k, v, w,
    basislen, startindex,
    totalpoints = loc->totalpoints,
    vdim = cov->vdim,
    tsdim = cov->tsdim,
    spatialdim = loc->spatialdim,
    *len = loc->length,
    total = totalpoints * vdim,
    *xi = S->xi,
    *powmatrix = S->powmatrix;
  //SEXP fctbody, tempres, envir, Rx;
  res_type *res = cov->rf;
  
 
  strcpy(errorloc_save, ERROR_LOC);
  sprintf(ERROR_LOC, "%s%s: ", errorloc_save, "add trend model");
 
  // print("%s\n", ERROR_LOC);

  if (mu != NULL) {  
    for (k=0; k<total; ) {
      for (v=0; v<vdim; res[k++] = mu[v++]); 
    } 
  } else {
    for (k=0; k<total; res[k++]=0.0);
  }
  
  if (plane != NULL) {  
    if (loc->grid) {
      for (w=0; w<tsdim; w++)  {
	x[w]=xgr[w][XSTART];
	xi[w]=0;
      }
      int endfor = totalpoints * vdim;
      for (k=0; k<endfor; ) {
	xA(x, plane, cov->nrow[TREND_LINEAR], cov->ncol[TREND_LINEAR],
	   evalplane);
	for(v=0; v<vdim; v++) res[k++] += evalplane[v];

	i = 0;
	(xi[i])++;
	x[i] += xgr[i][XSTEP];
	while(xi[i]>=len[i]) {
	  xi[i] = 0;
	  x[i] = xgr[i][XSTART];
	  if (i<tsdim-1) {
	    i++;
	    (xi[i])++;
	    x[i] += xgr[i][XSTEP];
	  } else {
	    assert(k==endfor);
	  }
	}
      }	
    } else if (loc->Time) {
      int m, 
	endfor= loc->length[loc->timespacedim-1];
      for(k=m=0, t=loc->T[XSTART]; m<endfor; m++, t+=loc->T[XSTEP]) {
	// for(t=loc->T[XSTART], i=0; t < loc->T[XEND]; t += loc->T[XSTEP]) {
	for(m=0, j=0; m < loc->spatialtotalpoints; m++) {
	  for(w=0; w<spatialdim; w++) x[w] = (loc->x)[j++];
	  x[spatialdim] = t;
	  xA(x, plane, cov->nrow[TREND_LINEAR], cov->ncol[TREND_LINEAR],
	     evalplane);
	  for(v=0;v<vdim;v++) res[k++] += evalplane[v];
	}
      }
    } else {
      int m,
	endfor = totalpoints * tsdim;
      for (k=m=0; m<endfor; m+=tsdim) {
	xA(loc->x + m, plane, cov->nrow[TREND_LINEAR], 
	   cov->ncol[TREND_LINEAR], evalplane);
	for(v=0;v<vdim;v++) res[k++] += evalplane[v];
      }
    }
  }
  
  if(polydeg != NULL) {
    startindex=0;
    int end_k = totalpoints * vdim;
    for(v=0; v<vdim; v++) { 
      basislen = binomialcoeff(((int *) polydeg)[v] + tsdim, tsdim);
      if(isnan(polycoeff[0])) {
	ERR("Error: cannot evaluate polynomial without coefficients.\n");
      }
      if (loc->grid) {
	for (w=0; w<tsdim; w++)  {
	  x[w]=xgr[w][XSTART];
	  xi[w]=0;
	}
	for(k=v; k<end_k; k+=vdim) {
	  //evaluation of trend polynomial
	  res[k] += evalpoly(x, powmatrix + startindex*tsdim,
			     polycoeff + startindex, basislen, tsdim);
	  i = 0;
	  (xi[i])++;
	  x[i] += xgr[i][XSTEP];
	  while(xi[i]>=len[i]) {
	    xi[i] = 0;
	    x[i] = xgr[i][XSTART];
	    if (i<tsdim-1) {
	      i++;
	      (xi[i])++;
	      x[i] += xgr[i][XSTEP];
	    } else {
	      assert(k==v + end_k - vdim);
	    }
	  }
	}	
      } else if(loc->Time) {
	 for(t=loc->T[XSTART], k=v; k<end_k; t+=loc->T[XSTEP]) {
	   for(i=0, j=0; i<loc->spatialtotalpoints; i++, k+=vdim) {
             for(w=0; w<spatialdim; w++)  x[w] = (loc->x)[j++];
	     x[spatialdim] = t;
	     //evaluation of trend polynomial
             res[k] += evalpoly(x, powmatrix + startindex*tsdim,
				    polycoeff + startindex, basislen, tsdim);
           }
	 }
      } else {
	for(k=v; k<end_k; k+=vdim) {
	  //evaluation of trend polynomial
	  res[k] += evalpoly(x, powmatrix + startindex*tsdim,
			     polycoeff + startindex, basislen, tsdim);
	}
      }
      startindex += basislen;
    }
  }
  
  if (arbitraryfct != NULL) { //muss hier arbitraryfct auswerten
    BUG;
//      if (isnan(fctcoeff[0])) {
// 	ERR("Error: cannot evaluate function without coefficient.\n");
//      }
//      
//      PROTECT(envir = allocSExp(ENVSXP));
//      SET_ENCLOS(envir, R_GlobalEnv);
//      
//      PROTECT(fctbody = allocSExp(CLOSXP));
//      fctbody = BODY(*((SEXP *) arbitraryfct));
//      
//      PROTECT(tempres = allocVector(REALSXP, vdim));
//      PROTECT(Rx = allocVector(REALSXP,xlen));     
//      
//      if (loc->grid) {
//          for (w=0; w<tsdim; w++)  {
// 	    x[w]=xgr[w][XSTART];
// 	    xi[w]=0;
//          }
//          for(k=0; k<totalpoints; k++) {
// 	   //evaluation of trend polynomial
//            for(w=0; w<xlen; w++) REAL(Rx)[w] = x[w];
//            defineVar(install("x"), Rx, envir);
//            if (spatialdim>1)
//               defineVar(install("y"), ScalarReal(x[1]), envir);
//            if (spatialdim>2)
//               defineVar(install("z"), ScalarReal(x[2]), envir);
//            defineVar(install("T"), ScalarReal(x[tsdim-1]), envir);
// 	   tempres = eval(fctbody, envir);
// 	   for (v=0; v<vdim; v++)
// 	     res[k*vdim+v] += fctcoeff[0]*REAL(tempres)[v];
//            i = 0;
//            (xi[i])++;
//            x[i] += xgr[i][XSTEP];
//            while(xi[i]>=len[i]) {
//              xi[i] = 0;
//              x[i] = xgr[i][XSTART];
//              if (i<tsdim-1) {
// 	       i++;
// 	       (xi[i])++;
// 	       x[i] += xgr[i][XSTEP];
//              } else {
// 	       assert(k==totalpoints-1);
//              }
//            }
//          }	
//      } else if(loc->Time) {
// 	 int k, 
// 	   endfor= loc->length[loc->timespacedim-1];
// 	 for(k=0, t=loc->T[XSTART], i=0; k<endfor; k++, t += loc->T[XSTEP]) {
// 	   //	 for(t=loc->T[XSTART], i=0; t < loc->T[XEND]; t+=loc->T[XSTEP]) {
// 	   for(k=0, j=0; k < loc->spatialtotalpoints; k++, i++) {
//              for(w=0; w<spatialdim; w++) x[w] = (loc->x)[j++];
// 	     x[spatialdim] = t;
// 	     //evaluation of trend polynomial
//              for(w=0; w<xlen; w++) REAL(Rx)[w] = x[w];
//              defineVar(install("x"), Rx, envir);
//              if (spatialdim > 1)
//                defineVar(install("y"), ScalarReal(x[1]), envir);
//              if (spatialdim > 2)
//                defineVar(install("z"), ScalarReal(x[2]), envir);
//              if (loc->Time)
//                defineVar(install("T"), ScalarReal(x[tsdim-1]), envir);
//              tempres = eval(fctbody, envir);
// 	     for (v=0; v<vdim; v++)
// 	       res[i*vdim+v] += fctcoeff[0]*REAL(tempres)[v];
//            }
// 	 }
//        } else {
// 	 for(k=0; k<totalpoints; k++) {
// 	   //evaluation of trend polynomial
// 	   for(w=0; w<spatialdim; w++) x[w] = (loc->x)[k*spatialdim+w];
//            for(w=0; w<xlen; w++) REAL(Rx)[w] = x[w];
//            defineVar(install("x"), Rx, envir);
//            if (spatialdim>1)
//               defineVar(install("y"), ScalarReal(x[1]), envir);
//            if (spatialdim>2)
//               defineVar(install("z"), ScalarReal(x[2]), envir);
//            defineVar(install("T"), ScalarReal(x[tsdim-1]), envir);
// 	   tempres = eval(fctbody, envir);
// 	   for (v=0; v<vdim; v++)
// 	     res[k*vdim+v] += fctcoeff[0]*REAL(tempres)[v];
// 	   }
//        }
//        UNPROTECT(4);
  }  
  
  return;

}


int binomialcoeff(int n, int k) {
  //programmed as in wikipedia
  int i, res=1;
  
  if((k < 0) || (k > n)) return 0;
  if(k > n-k) k = n-k; //symmetry
  for(i=0; i<k; i++) {
     res *= n-i;
     res /= i+1; //two steps because of integer division
  }
  return res;
}

void poly_basis(cov_model *cov, storage  VARIABLE_IS_NOT_USED *s) {
  
  trend_storage *S = cov->Strend;
  int basislen=0, powsum, d, i, j, k, v,
    dim = cov->tsdim,
    vdim = cov->vdim,
    *powmatrix = S->powmatrix,
      *dimi=NULL,
      err=NOERROR;
  double *polydeg = cov->p[TREND_POLY];
  
  dimi = (int *) MALLOC(dim * sizeof(int));
  if (dimi == NULL) {
     err = ERRORMEMORYALLOCATION;
     goto ErrorHandling;
  }
  
  for(v=0, j=0; v<vdim; v++) {
    basislen = binomialcoeff(((int *) polydeg)[v]+dim,dim);
    //   print("v: %d\n", v);
    for(d=0; d<dim; d++) dimi[d] = 0;
    for(k=0; k<basislen; k++) {
      for(d=0; d<dim; d++) powmatrix[j++] = dimi[d];
      // for(d=0;d<dim;d++) print("%d ", powmatrix[j*dim+d]);
      // print("\n");
      i=0;
      (dimi[i])++;
      powsum = 0;
      for(d=0; d<dim; d++) powsum += dimi[d];
      while(powsum > ((int*) polydeg)[v]) {
        dimi[i] = 0;
        if (i < dim-1) {
	  i++;
	  (dimi[i])++;
        }
        powsum = 0;
        for(d=0; d<dim; d++) powsum += dimi[d];
      }
    }
    // print("\n");
  }
  //print("\n");
  
  ErrorHandling:
  if (dimi != NULL) free(dimi);   
  dimi = NULL;
  if (err != NOERROR) XERR(err);
   
  return;
  
}

double evalpoly(double *x, int *powmatrix, double *polycoeff, int basislen,
	      int dim) {
  int i, d, j;
  double res = 0, tempres;
  for(j=i=0; i<basislen; i++) {
     tempres=1;
     for(d=0; d<dim; d++) tempres *= pow(x[d], powmatrix[j++]);
     res += polycoeff[i] * tempres;
  }
  return(res);
}

