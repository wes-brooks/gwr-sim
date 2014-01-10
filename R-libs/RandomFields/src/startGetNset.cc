/* 
 Authors
 Martin Schlather, schlather@math.uni-mannheim.de

 (library for simulation of random fields)

 Copyright (C) 2001 -- 2013 Martin Schlather, 

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 3
of the License, or (at your option) any later version.
RO
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/

#include <R.h>
#include <Rdefines.h>
#include <math.h>  
#include <stdio.h>  
#include <stdlib.h>
 
#include <string.h>
#include "RF.h"
#include "primitive.h"
#include <R_ext/Linpack.h>
#include "Covariance.h"


void addkappa(int i, const char *n, SEXPTYPE t, Types ParamType) {
  cov_fct *C = CovList + currentNrCov - 1;

  // printf("%s\n", C->name);

  assert(n[0] != '\0' && 
	 (n[0] != 'k' || n[1] != '\0') // reserved for standard parameter names
	 && (n[0] != 'C' || (n[1] != '\0' && (n[1] < '1' ||  n[1] > '9')))
                                        // reserved for standard submodel names
	 );
  assert(i < C->kappas);
  strcopyN(C->kappanames[i], n, PARAMMAXCHAR);
  C->kappatype[i] = t;
  C->kappaParamType[i] = ParamType;
  assert(strcmp(n, FREEVARIABLE) || C->internal);
  if (t >= LISTOF) assert(strcmp(C->kappanames[0], "element") == 0 
			  ||  strcmp(C->name, "select") == 0
			  );
}


void addkappa(int i, const char *n, SEXPTYPE t) {
  addkappa(i, n, t, RandomType);
}

void kappanames(const char* n1, SEXPTYPE t1) {
  assert({cov_fct *C = CovList + currentNrCov - 1; C->kappas == 1;});
  addkappa(0, n1, t1);
}
void kappanames(const char* n1, SEXPTYPE t1, const char* n2, SEXPTYPE t2) {
  assert({cov_fct *C = CovList + currentNrCov - 1; C->kappas == 2;});
  addkappa(0, n1, t1);
  addkappa(1, n2, t2);
}
void kappanames(const char* n1, SEXPTYPE t1, const char* n2, SEXPTYPE t2,
		const char* n3, SEXPTYPE t3) {
  assert({cov_fct *C = CovList + currentNrCov - 1; C->kappas == 3;});
  addkappa(0, n1, t1);
  addkappa(1, n2, t2);
  addkappa(2, n3, t3);
}
void kappanames(const char* n1, SEXPTYPE t1, const char* n2, SEXPTYPE t2,
		const char* n3, SEXPTYPE t3, const char* n4, SEXPTYPE t4) {
  assert({ cov_fct *C = CovList + currentNrCov - 1;C->kappas == 4;});
  addkappa(0, n1, t1);
  addkappa(1, n2, t2);
  addkappa(2, n3, t3);
  addkappa(3, n4, t4);
}
void kappanamesX(const char* n1, SEXPTYPE t1, const char* n2, SEXPTYPE t2,
		const char* n3, SEXPTYPE t3, const char* n4, SEXPTYPE t4,
		const char* n5, SEXPTYPE t5) {
  addkappa(0, n1, t1);
  addkappa(1, n2, t2);
  addkappa(2, n3, t3);
  addkappa(3, n4, t4);
  addkappa(4, n5, t5);
}
void kappanames(const char* n1, SEXPTYPE t1, const char* n2, SEXPTYPE t2,
		const char* n3, SEXPTYPE t3, const char* n4, SEXPTYPE t4,
		const char* n5, SEXPTYPE t5) {
  assert({ cov_fct *C = CovList + currentNrCov - 1;C->kappas == 5;});
  kappanamesX(n1, t1, n2, t2, n3, t3, n4, t4, n5, t5);
}

void kappanames(const char* n1, SEXPTYPE t1, const char* n2, SEXPTYPE t2,
		const char* n3, SEXPTYPE t3, const char* n4, SEXPTYPE t4,
		const char* n5, SEXPTYPE t5, const char* n6, SEXPTYPE t6) {
  assert({ cov_fct *C = CovList + currentNrCov - 1;C->kappas == 6;});
  kappanamesX(n1, t1, n2, t2, n3, t3, n4, t4, n5, t5);
  addkappa(5, n6, t6);
}
void kappanames(const char* n1, SEXPTYPE t1, const char* n2, SEXPTYPE t2,
		const char* n3, SEXPTYPE t3, const char* n4, SEXPTYPE t4,
		const char* n5, SEXPTYPE t5, const char* n6, SEXPTYPE t6,
		const char* n7, SEXPTYPE t7) {
  assert({ cov_fct *C = CovList + currentNrCov - 1; C->kappas == 7;});
  kappanamesX(n1, t1, n2, t2, n3, t3, n4, t4, n5, t5);
  addkappa(5, n6, t6);
  addkappa(6, n7, t7);
}
void kappanames(const char* n1, SEXPTYPE t1, const char* n2, SEXPTYPE t2,
		const char* n3, SEXPTYPE t3, const char* n4, SEXPTYPE t4,
		const char* n5, SEXPTYPE t5, const char* n6, SEXPTYPE t6,
		const char* n7, SEXPTYPE t7, const char* n8, SEXPTYPE t8) {
  assert({ cov_fct *C = CovList + currentNrCov - 1;C->kappas == 8;});
  kappanamesX(n1, t1, n2, t2, n3, t3, n4, t4, n5, t5);
  addkappa(5, n6, t6);
  addkappa(6, n7, t7);
  addkappa(7, n8, t8);
}
void kappanames(const char* n1, SEXPTYPE t1, const char* n2, SEXPTYPE t2,
		const char* n3, SEXPTYPE t3, const char* n4, SEXPTYPE t4,
		const char* n5, SEXPTYPE t5, const char* n6, SEXPTYPE t6,
	       	const char* n7, SEXPTYPE t7, const char* n8, SEXPTYPE t8,
		const char* n9, SEXPTYPE t9) {
  assert({ cov_fct *C = CovList + currentNrCov - 1;C->kappas == 9;});
 kappanamesX(n1, t1, n2, t2, n3, t3, n4, t4, n5, t5);
  addkappa(5, n6, t6);
  addkappa(6, n7, t7);
  addkappa(7, n8, t8);
  addkappa(8, n9, t9);
}
void kappanamesX(const char* n1, SEXPTYPE t1, const char* n2, SEXPTYPE t2,
		const char* n3, SEXPTYPE t3, const char* n4, SEXPTYPE t4,
		const char* n5, SEXPTYPE t5, const char* n6, SEXPTYPE t6,
	       	const char* n7, SEXPTYPE t7, const char* n8, SEXPTYPE t8,
		const char* n9, SEXPTYPE t9, const char* n10, SEXPTYPE t10) {
  kappanamesX(n1, t1, n2, t2, n3, t3, n4, t4, n5, t5);
  addkappa(5, n6, t6);
  addkappa(6, n7, t7);
  addkappa(7, n8, t8);
  addkappa(8, n9, t9);
  addkappa(9, n10, t10);
}

void kappanames(const char* n1, SEXPTYPE t1, const char* n2, SEXPTYPE t2,
		const char* n3, SEXPTYPE t3, const char* n4, SEXPTYPE t4,
		const char* n5, SEXPTYPE t5, const char* n6, SEXPTYPE t6,
	       	const char* n7, SEXPTYPE t7, const char* n8, SEXPTYPE t8,
		const char* n9, SEXPTYPE t9, const char* n10, SEXPTYPE t10) {
  assert({ cov_fct *C = CovList + currentNrCov - 1;C->kappas == 10;});
  kappanamesX(n1, t1, n2, t2, n3, t3, n4, t4, n5, t5, n6, t6, n7, t7, n8, t8,
	     n9, t9, n10, t10);
}

void kappanames(const char* n1, SEXPTYPE t1, const char* n2, SEXPTYPE t2,
		const char* n3, SEXPTYPE t3, const char* n4, SEXPTYPE t4,
		const char* n5, SEXPTYPE t5, const char* n6, SEXPTYPE t6,
	       	const char* n7, SEXPTYPE t7, const char* n8, SEXPTYPE t8,
		const char* n9, SEXPTYPE t9, const char* n10, SEXPTYPE t10,
		const char* n11, SEXPTYPE t11) {
  assert({ cov_fct *C = CovList + currentNrCov - 1; C->kappas == 11;});
  kappanamesX(n1, t1, n2, t2, n3, t3, n4, t4, n5, t5, n6, t6, n7, t7, n8, t8,
	     n9, t9, n10, t10);
  addkappa(10, n11, t11);
}

void kappanames(const char* n1, SEXPTYPE t1, const char* n2, SEXPTYPE t2,
		const char* n3, SEXPTYPE t3, const char* n4, SEXPTYPE t4,
		const char* n5, SEXPTYPE t5, const char* n6, SEXPTYPE t6,
	       	const char* n7, SEXPTYPE t7, const char* n8, SEXPTYPE t8,
		const char* n9, SEXPTYPE t9, const char* n10, SEXPTYPE t10,
		const char* n11, SEXPTYPE t11, const char* n12, SEXPTYPE t12) {
  assert({ cov_fct *C = CovList + currentNrCov - 1;C->kappas == 12;});
  kappanamesX(n1, t1, n2, t2, n3, t3, n4, t4, n5, t5, n6, t6, n7, t7, n8, t8,
	     n9, t9, n10, t10);
  addkappa(10, n11, t11);
  addkappa(11, n12, t12);
}

void kappanames(const char* n1, SEXPTYPE t1, const char* n2, SEXPTYPE t2,
		const char* n3, SEXPTYPE t3, const char* n4, SEXPTYPE t4,
		const char* n5, SEXPTYPE t5, const char* n6, SEXPTYPE t6,
	       	const char* n7, SEXPTYPE t7, const char* n8, SEXPTYPE t8,
		const char* n9, SEXPTYPE t9, const char* n10, SEXPTYPE t10,
		const char* n11, SEXPTYPE t11, const char* n12, SEXPTYPE t12,
		const char* n13, SEXPTYPE t13) {
   assert({ cov_fct *C = CovList + currentNrCov - 1;C->kappas == 13;});
  kappanamesX(n1, t1, n2, t2, n3, t3, n4, t4, n5, t5, n6, t6, n7, t7, n8, t8,
	     n9, t9, n10, t10);
  addkappa(10, n11, t11);
  addkappa(11, n12, t12);
  addkappa(12, n13, t13);
}

void kappanames(const char* n1, SEXPTYPE t1, const char* n2, SEXPTYPE t2,
		const char* n3, SEXPTYPE t3, const char* n4, SEXPTYPE t4,
		const char* n5, SEXPTYPE t5, const char* n6, SEXPTYPE t6,
	       	const char* n7, SEXPTYPE t7, const char* n8, SEXPTYPE t8,
		const char* n9, SEXPTYPE t9, const char* n10, SEXPTYPE t10,
		const char* n11, SEXPTYPE t11, const char* n12, SEXPTYPE t12,
		const char* n13, SEXPTYPE t13, const char* n14, SEXPTYPE t14) {
  assert({ cov_fct *C = CovList + currentNrCov - 1;C->kappas == 14;});
  kappanamesX(n1, t1, n2, t2, n3, t3, n4, t4, n5, t5, n6, t6, n7, t7, n8, t8,
	     n9, t9, n10, t10);
  addkappa(10, n11, t11);
  addkappa(11, n12, t12);
  addkappa(12, n13, t13);
  addkappa(13, n14, t14);
}

void kappanamesX(const char* n1, SEXPTYPE t1, const char* n2, SEXPTYPE t2,
		const char* n3, SEXPTYPE t3, const char* n4, SEXPTYPE t4,
		const char* n5, SEXPTYPE t5, const char* n6, SEXPTYPE t6,
	       	const char* n7, SEXPTYPE t7, const char* n8, SEXPTYPE t8,
		 const char* n9, SEXPTYPE t9, const char* n10, SEXPTYPE t10,
		const char* n11, SEXPTYPE t11, const char* n12, SEXPTYPE t12,
		const char* n13, SEXPTYPE t13, const char* n14, SEXPTYPE t14, 
		const char* n15, SEXPTYPE t15) {
  kappanamesX(n1, t1, n2, t2, n3, t3, n4, t4, n5, t5, n6, t6, n7, t7, n8, t8,
	     n9, t9, n10, t10);
  addkappa(10, n11, t11);
  addkappa(11, n12, t12);
  addkappa(12, n13, t13);
  addkappa(13, n14, t14);
  addkappa(14, n15, t15);
}

void kappanames(const char* n1, SEXPTYPE t1, const char* n2, SEXPTYPE t2,
		const char* n3, SEXPTYPE t3, const char* n4, SEXPTYPE t4,
		const char* n5, SEXPTYPE t5, const char* n6, SEXPTYPE t6,
	       	const char* n7, SEXPTYPE t7, const char* n8, SEXPTYPE t8,
		const char* n9, SEXPTYPE t9, const char* n10, SEXPTYPE t10,
		const char* n11, SEXPTYPE t11, const char* n12, SEXPTYPE t12,
		const char* n13, SEXPTYPE t13, const char* n14, SEXPTYPE t14, 
		const char* n15, SEXPTYPE t15) {
  assert({ cov_fct *C = CovList + currentNrCov - 1;C->kappas == 15;});
  kappanamesX(n1, t1, n2, t2, n3, t3, n4, t4, n5, t5, n6, t6, n7, t7, n8, t8,
	      n9, t9, n10, t10, n11, t11, n12, t12, n13, t13, n14, t14, n15, 
	      t15);  
}

void kappanames(const char* n1, SEXPTYPE t1, const char* n2, SEXPTYPE t2,
		const char* n3, SEXPTYPE t3, const char* n4, SEXPTYPE t4,
		const char* n5, SEXPTYPE t5, const char* n6, SEXPTYPE t6,
	       	const char* n7, SEXPTYPE t7, const char* n8, SEXPTYPE t8,
		const char* n9, SEXPTYPE t9, const char* n10, SEXPTYPE t10,
		const char* n11, SEXPTYPE t11, const char* n12, SEXPTYPE t12,
		const char* n13, SEXPTYPE t13, const char* n14, SEXPTYPE t14, 
		const char* n15, SEXPTYPE t15, const char* n16, SEXPTYPE t16) {
  assert({ cov_fct *C = CovList + currentNrCov - 1;C->kappas == 16;});
  kappanamesX(n1, t1, n2, t2, n3, t3, n4, t4, n5, t5, n6, t6, n7, t7, n8, t8,
	      n9, t9, n10, t10, n11, t11, n12, t12, n13, t13, n14, t14, n15, 
	      t15);  
  addkappa(15, n16, t16);
}

sortsofparam standard_paramtype(int k, int VARIABLE_IS_NOT_USED row, int VARIABLE_IS_NOT_USED col) {
  // for MLE
// k non-negative: k-th parameter
// k negative: (i-1)th submodel
  return k<0 ? VARPARAM : CRITICALPARAM;
// k<0: varparam used to indicate that variance is allowed for submodel,
// see recursive call getnaposition; E.g. not allowed for second submodel
// of nsst: the variance parameter is something between scale and variance
}

sortsofparam ignoreall_paramtype(int VARIABLE_IS_NOT_USED k, int VARIABLE_IS_NOT_USED row, int VARIABLE_IS_NOT_USED col) {
  // for MLE
// k non-negative: k-th parameter
// k negative: (i-1)th submodel
  return IGNOREPARAM;
// k<0: varparam used to indicate that variance is allowed for submodel,
// see recursive call getnaposition; E.g. not allowed for second submodel
// of nsst: the variance parameter is something between scale and variance
}

sortsofparam uncritical_paramtype(int k, int VARIABLE_IS_NOT_USED row, int VARIABLE_IS_NOT_USED col) {
  // for MLE
// k non-negative: k-th parameter
// k negative: (i-1)th submodel
  return k<0 ? VARPARAM : ANYPARAM;
// k<0: varparam used to indicate that variance is allowed for submodel,
// see recursive call getnaposition; E.g. not allowed for second submodel
// of nsst: the variance parameter is something between scale and variance
}



void add_paramtype(paramtype_fct paramtype) {
  cov_fct *C = CovList + currentNrCov - 1;
  C->paramtype = paramtype;
}

void addsub(int i, const char *n) {
  cov_fct *C = CovList + currentNrCov - 1;
  int j;
  assert(n[0] != 'k' && n[0] != 'C');
  assert(i < MAXSUB);

  strcopyN(C->subnames[i], n, PARAMMAXCHAR);
  C->subintern[i] = false;
  for (j=0; j<C->kappas; j++) {
    if (C->kappanames[j] == NULL)
      error("give all names for the parameters first");
    if ((C->subintern[i] = strcmp(C->kappanames[j], C->subnames[i]) == 0))
      break;
  }
}

void subnames(const char* n1) {
  addsub(0, n1);
}
void subnames(const char* n1, const char* n2) {
  addsub(0, n1);
  addsub(1, n2);
}
void subnames(const char* n1, const char* n2, const char* n3) {
  addsub(0, n1);
  addsub(1, n2);
  addsub(2, n3);
}
void subnames(const char* n1, const char* n2, const char* n3, const char* n4) {
  addsub(0, n1);
  addsub(1, n2);
  addsub(2, n3);
  addsub(3, n4);
}
void subnames(const char* n1, const char* n2, const char* n3, const char* n4,
	      const char* n5) {
  addsub(0, n1);
  addsub(1, n2);
  addsub(2, n3);
  addsub(3, n4);
  addsub(4, n5);
}
void subnames(const char* n1, const char* n2, const char* n3, const char* n4,
	      const char* n5, const char* n6) {
  addsub(0, n1);
  addsub(1, n2);
  addsub(2, n3);
  addsub(3, n4);
  addsub(4, n5);
  addsub(5, n6);
}


int xxx(int x) {return (int) pow(10, (double) x);}
void crash(cov_model *cov) {
    PRINTF("crash!!!!\n");
    char m[1]; MEMCOPY(m, cov, xxx(3));  PRINTF("%s\n", m);
}
void crash() {cov_model *cov=NULL; crash(cov);}

void ErrCov(double VARIABLE_IS_NOT_USED *x, cov_model *cov,
	    double VARIABLE_IS_NOT_USED *v) {
    // 
  PRINTF("\nErrCov %s [%d] gatter=%d:\n", 
	 NICK(cov), cov->nr, cov->gatternr); // ok
  if (PL >=  PL_ERRORS){
    PMI(cov, "ErrCov");//
    crash(cov);
  }
  ERR("unallowed or undefined call of function");
}
void ErrLogCov(double VARIABLE_IS_NOT_USED *x, cov_model *cov, 
	       double VARIABLE_IS_NOT_USED *v, 
	       double VARIABLE_IS_NOT_USED *sign) {
  // 
  PRINTF("\nErrlogCov %s:\n", NICK(cov));
  if (PL >=  PL_ERRORS) {
    PMI(cov, "ErrlogCov");//
    crash(cov); 
  }
  ERR("unallowed or undefined call of function (log)");
}
void ErrCovNonstat(double VARIABLE_IS_NOT_USED *x, 
		   double VARIABLE_IS_NOT_USED *y, cov_model *cov, 
		   double VARIABLE_IS_NOT_USED *v) {
  PRINTF("\nErrCovNonstat %s: (%d)\n", NICK(cov), cov->nr);
  //  APMI(cov->calling);
  if (PL >= PL_ERRORS) {
    PMI(cov->calling, "ErrCovNonstat"); //
    crash(cov);
  }
  ERR("unallowed or undefined call of non-domain function");
}
void ErrLogCovNonstat(double VARIABLE_IS_NOT_USED *x, 
		      double VARIABLE_IS_NOT_USED *y, cov_model *cov,
		      double VARIABLE_IS_NOT_USED *v, 
		      double VARIABLE_IS_NOT_USED *sign) {
  PRINTF("\nErrlogCovNonstat %s: (%d)\n", NICK(cov), cov->nr);
  if (PL >= PL_ERRORS) {
    PMI(cov->calling, "ErrLogCovNonstat"); //
    crash(cov);
  }
  ERR("unallowed or undefined call of non-domain function (log)");
}
void Errspectral(cov_model *cov,
		 storage VARIABLE_IS_NOT_USED *s, 
		 double VARIABLE_IS_NOT_USED *e) {
  PRINTF("\nErrlogCovNonstat %s: (%d)\n", NICK(cov), cov->nr);
 if (PL >= PL_ERRORS) {
   PMI(cov->calling, "Errspectral"); //
    crash(cov);
  }
  ERR("unallowed or undefined call of spectral function");
} 

void ErrInverseNonstat(double VARIABLE_IS_NOT_USED *v, cov_model *cov,
		       double *x, double *y) {
  x[0] = -INVERSE_UNKNOWN;
  y[0] = INVERSE_UNKNOWN;
  return;

  PRINTF("\nErrCovNonstat %s: (%d)\n", NICK(cov), cov->nr);
  //  APMI(cov->calling);
  if (PL >= PL_ERRORS) {
    PMI(cov->calling, "ErrCovNonstat"); //
    crash(cov);
  }
  ERR("unallowed or undefined call of non-domain function");
}

void ErrInverse(double VARIABLE_IS_NOT_USED *v, 
		       cov_model VARIABLE_IS_NOT_USED *cov,
		       double *x) {
  x[0] = INVERSE_UNKNOWN;
  return;
}

char InternalName[]="-";
void kappasize1(int VARIABLE_IS_NOT_USED i,
		cov_model VARIABLE_IS_NOT_USED *cov, int *nrow, int *ncol) {
  *nrow = *ncol = 1;
}


void rangeOK(cov_model VARIABLE_IS_NOT_USED *cov, 
	     range_type VARIABLE_IS_NOT_USED *range) { 
  assert(CovList[cov->nr].kappas == 0);
}

int structOK(cov_model VARIABLE_IS_NOT_USED *cov, cov_model VARIABLE_IS_NOT_USED **newmodel){ return NOERROR;}

int struct_statiso(cov_model *cov, cov_model **newmodel) {
  cov_fct *C = CovList + cov->nr;

  //   printf("struct_statiso %s role=%s\n", NICK(cov), ROLENAMES[cov->role]);
   
  ASSERT_NEWMODEL_NOT_NULL;
  

  switch (cov->role) {
  case ROLE_POISSON_GAUSS :
    SERR1("Unexpected call of'struct' by '%s'", C->nick);
    break;      
  case ROLE_POISSON : case ROLE_MAXSTABLE :
    cov->mpp.maxheight = 1.0;
    if (C->finiterange == true) {
      return addUnifModel(cov, 1.0, newmodel);
    } else SERR2("The function '%s' has inifinite support use '%s' to truncate the support.", C->nick, CovList[TRUNCSUPPORT].nick);
  default : ILLEGAL_ROLE_STRUCT;
  }    
  return NOERROR;
}

int struct_failed(cov_model *cov, cov_model VARIABLE_IS_NOT_USED **newmodel) {
  //  PMI(cov);
  //  printf("%ld %ld %ld\n", CovList[cov->nr].Struct, CovList[DISTRIBUTION].Struct, struct_failed);  crash();
  SERR3("initialization failed --  model '%s' (%d) does not fit (yet) the properties required by '%s'.",
	NICK(cov), cov->nr, cov->calling == NULL ? "<null>" : NICK(cov->calling)); 
}

int initOK(cov_model *cov, storage *s) {
  cov_fct *C = CovList + cov->nr; // nicht gatternr  
  int i, err = NOERROR,
    nk = C->kappas;
  bool random = false;
  for (i=0; i<nk; i++) {
    cov_model *ks = cov->kappasub[i];
    if (ks != NULL) {
      if (isRandom(C->kappaParamType[i])) {
	random = true;
	if ((err = INIT(ks, cov->mpp.moments, s)) != NOERROR) return err;
      } else {
	SERR2("%s : parameter %s is not of random type", 
	      C->nick, C->kappanames[i]);
      }
    }
  }
  if (random) SERR("'initOK' not programmed yet for 'random'");
  return err; 
}

int init_failed(cov_model *cov, storage VARIABLE_IS_NOT_USED *s) {
  if (PL >= PL_ERRORS) PRINTF("init failed cov=%s:\n", NICK(cov));
  SERR("Init failed. Compound Poisson fields are essentially only programmed for simple, domain and isotropic shape functions");
  return ERRORFAILED;
}

int init_statiso(cov_model *cov, storage *s) {
  // only domain and isotropic models
  cov_fct *C = CovList + cov->nr; // nicht gatternr
  int err;

  if (cov->finiterange == true) {
    if  (C->finiterange == true) {
      // ?? to do ??
      // cov->mpp.refradius = 1.0; // jetzt ueber Inverse abgefangen
    }
  }
  
  if ((err = initOK(cov, s)) == NOERROR) return err;
 
    if (cov->role == ROLE_POISSON) { 
      return NOERROR;
  }


  if (PL >= PL_ERRORS) PRINTF("init failed cov=%s:\n", NICK(cov));

  //  PMI(cov);

  SERR("Call of init: Compound Poisson fields are essentially only programmed for domain and isotropic functions");

 
  return NOERROR;
}


void doOK(cov_model VARIABLE_IS_NOT_USED *cov, storage VARIABLE_IS_NOT_USED *s) { }

void do_failed(cov_model *cov, storage VARIABLE_IS_NOT_USED *s) {
  //PMI(cov->calling);
  if (PL >= PL_ERRORS) PRINTF("do failed %s:\n", NICK(cov));
  ERR("Call of do: Compound Poisson fields are essentially only programmed for domain and isotropic shape functions");
}

void do_random_failed(cov_model *cov, double VARIABLE_IS_NOT_USED *v) {
  if (PL >= PL_ERRORS) PRINTF("do_random failed %s:\n", NICK(cov));
  ERR("Call of do: Compound Poisson fields are essentially only programmed for domain and isotropic shape functions");
}

void do_statiso(cov_model *cov, storage VARIABLE_IS_NOT_USED *s) {
 if (cov->role == ROLE_POISSON) return;
  if (PL >= PL_ERRORS) ERR("Call of do: Compound Poisson fields are essentially only programmed for domain and isotropic shape functions");
 
}

static int badname = -1;
void nickname(const char *name, int nr) {
  char dummy[MAXCHAR];
  cov_fct *C = CovList + nr; // nicht gatternr 
  
  int sl = strlen(CAT_TYPENAMES[C->Type]);  
  strcopyN(dummy, name, MAXCHAR-sl);
  sprintf(C->nick, "%s%s", CAT_TYPENAMES[C->Type], dummy);
  strcpy(CovNickNames[nr], C->nick);

  if ((int) strlen(name) >= (int) MAXCHAR - sl) {
    badname = nr;
  } else {
    if (badname >= 0 && badname != nr) 
      PRINTF("Warning! Nick name is truncated to `%s'.\n", 
	     CovList[badname].nick);
    badname = -1;
  }
}

void insert_name(int curNrCov, const char *name) {
  cov_fct *C = CovList + curNrCov;
  char dummy[MAXCHAR];
  strcopyN(dummy, name, MAXCHAR);
  strcpy(CovNames[curNrCov], dummy);
  strcpy(C->name, dummy);
  if (strlen(name)>=MAXCHAR) {
    PRINTF("Warning! Covariance name is truncated to `%s'.\n", C->name);
  }
  assert(strcmp(InternalName, name));
  nickname(name, curNrCov);
}

void StandardCovariance(cov_model *cov, double *v){
  CovVario(cov, true, false, v); 
}
void StandardCovMatrix(cov_model *cov, double *v) { 
  CovarianceMatrix(cov, v); 
}
void StandardInverseCovMatrix(cov_model *cov, double *v) { 
  InverseCovMatrix(cov, v); 
}
void StandardVariogram(cov_model *cov, double *v) {
  //assert(false);
  CovVario(cov, false, false, v); 
}
void StandardPseudoVariogram(cov_model *cov, double *v) {
  CovVario(cov, false, true, v); 
}
void StandardSelectedCovMatrix(cov_model *cov, int *sel, int nsel, double *v) { 
  SelectedCovMatrix(cov, sel, nsel, v); 
}

void StandardInverseNonstat(double *v, cov_model *cov,
			    double *left, double *right) {
  assert(CovList[cov->nr].inverse != NULL);
  double x;
  int d,
    dim = cov->xdimprev;
  
  // PMI(cov); printf("%d %d\n", dim, cov->tsdim);

  assert(dim == cov->tsdim);
  
  INVERSE(v, cov, &x);

  for (d=0; d<dim; d++) {
    //    printf("inverse %d %f %f\n", d, *v, x);

    left[d] = -x;
    right[d] = x;
  }
}

char isTrue(cov_model VARIABLE_IS_NOT_USED *cov) {error("isTrue may not be used anymore"); return true;}
char isFalse(cov_model VARIABLE_IS_NOT_USED *cov) {return false;}

void InverseFiniteRange(double VARIABLE_IS_NOT_USED *x, cov_model VARIABLE_IS_NOT_USED *cov, double *v){ *v = 1.0; }

void InverseIsoMon(double VARIABLE_IS_NOT_USED *x, cov_model *cov, double VARIABLE_IS_NOT_USED *v){
  searchInverse(CovList[cov->nr].cov, cov, 1.0, *x, 0.001);
}


void createmodel(const char *name, Types type, int kappas, size_fct kappasize,	
		 domain_type domain, isotropy_type isotropy,
		 checkfct check, rangefct range, int vdim, pref_shorttype pref,
		 int maxdim, ext_bool finiterange, int monotone) {
  int i;
  cov_fct *C = CovList + currentNrCov;
  bool 
    stat_iso = domain == XONLY && isotropy == ISOTROPIC;

  assert(domain >= 0 && isotropy >= 0);

  if (currentNrCov==-1) InitModelList(); 
  assert(CovList!=NULL);
  assert(currentNrCov>=0);
  if (currentNrCov>=MAXNRCOVFCTS) {
    char msg[150];
    sprintf(msg, "maximum number of covariance models reached. Last included  model is '%s'.", CovList[MAXNRCOVFCTS-1].name);
    warning(msg);
  }
  
  if (PL >= PL_DETAILS) PRINTF("%d %s vdim=%d\n", currentNrCov, name, vdim); 
  C->Type = type; // ganz vorne 
  assert(type >=0 && type <= OtherType);
  C->TypeFct = NULL;

  insert_name(currentNrCov, name);

  //if (domain == PARAM_DEP) printf("%s\n", C->nick);

  // nachfolgendes unbedingt lassen, da sonst vor Aufruf von Interface
  // reduzierende Gatter gesetzt werden koennten
  assert(type != InterfaceType || (domain == XONLY && isotropy==NO_ROTAT_INV));

  C->kappas = kappas;
  assert(kappas >= 0 && kappas <= MAXPARAM);
  C->minsub = C->maxsub = 0;
  C->domain = domain; 
  C->isotropy = isotropy;
  C->vdim = vdim; assert(vdim != MISMATCH);
  C->maxdim = maxdim;
  C->maxmoments = 0;
  assert(maxdim != MISMATCH);
  for (i=0; i<kappas; i++) {
    sprintf(C->kappanames[i], "k%d", i); // default (repeated twice)
    C->kappatype[i] = REALSXP;
  }
  C->kappasize = (kappasize == NULL) ? kappasize1 : kappasize;
  C->paramtype=standard_paramtype;
    
  if (kappas==0) {
    assert(range == NULL);
    C->range=rangeOK; 
  } else {
    assert(range != NULL);
    C->range= range;
  }
  C->check = check;
  assert(check != NULL);
  for (i=0; i<Forbidden; i++)
    C->implemented[i] = NOT_IMPLEMENTED;
  C->internal = false;
  assert(finiterange != MISMATCH);
  assert(finiterange != PARAM_DEP || check==checktbmop);
  C->finiterange = finiterange;
  assert(monotone != MISMATCH || 
	 (monotone == MISMATCH && (type==RandomType || type== PointShapeType
				   || type== InterfaceType 
				   || type== GaussMethodType
				   || type== BrMethodType 
				   || type==ProcessType)));
  assert(monotone != PARAM_DEP);
  C->Monotone = monotone;


  MEMCOPY(C->pref, pref, sizeof(pref_shorttype));
 
  C->cov = C->D = C->D2 = C->D3 = C->D4 =C->tbm2 = C->nabla = C->hess = 
    C->random = ErrCov;
  C->inverse = finiterange == true ? InverseFiniteRange : ErrInverse;
  C->nonstat_inverse =  C->nonstat_inverse_D = ErrInverseNonstat;
  C->density = MISMATCH;
  C->log = C->logD = ErrLogCov;
  C->F_derivs = C->RS_derivs = isProcess(type) ? 0 : MISMATCH;
  C->nonstat_cov = C->nonstat_D = C->nonstat_random = ErrCovNonstat;
  C->nonstatlog = ErrLogCovNonstat;
  C->aux_cov = NULL;
  C->coinit = C->ieinit = NULL;
  C->alternative = NULL;

  C->spectral=Errspectral;

  C->drawmix = NULL;
  C->logmixdens = NULL;

  C->Struct = stat_iso ? struct_statiso : struct_failed;
  C->Init = stat_iso ? init_statiso : init_failed;
  //  printf("%s failed=%d\n", C->nick, !stat_iso);
  C->Do = stat_iso ? do_statiso : do_failed;
  C->minmaxeigenvalue = NULL;

  C->hyperplane=NULL;  
  C->primitive = true;
  C->covariance = StandardCovariance;
  C->covmatrix = StandardCovMatrix;
  C->inversecovmatrix = StandardInverseCovMatrix;
  C->variogram = StandardVariogram;
  C->pseudovariogram = StandardPseudoVariogram;
  C->selectedcovmatrix = StandardSelectedCovMatrix;
  C->is_covariance=C->is_covmatrix=C->is_inversecovmatrix=C->is_variogram
    = C->is_pseudovariogram = isFalse;

  C->TaylorN = C->TailN = MISMATCH;
 
  currentNrCov++;
}

int CopyModel(const char *name, int which) {
  memcpy(CovList + currentNrCov, CovList + which, sizeof(cov_fct));  
  insert_name(currentNrCov, name);
  currentNrCov++;
  return currentNrCov - 1;
}


int CopyModel(const char *name, int which, Types type) {
  CopyModel(name, which);
  int nr = currentNrCov - 1;
  CovList[nr].Type = type;
  return nr;
}


int CopyModel(const char *name, int which, checkfct check) {  
  CopyModel(name, which);
  int nr = currentNrCov - 1;  
  CovList[nr].check = check;
  return nr;
}

/*
  int CopyModel(const char *name, int which, checkfct check) {
  CopyModel(name, which);
  int nr = currentNrCov - 1;
  CovList[nr].check = check;
  }
*/

void nickname(const char *name) {
  int nr = currentNrCov - 1;
  nickname(name, nr);
}

int IncludePrim(const char *name, Types type,  int kappas, 
		domain_type domain, isotropy_type isotropy,	
		checkfct check, rangefct range, int maxdim, 
		ext_bool finiterange, int monotonicity) {  
  createmodel(name, type, kappas, NULL, domain, isotropy, check, range,
	      SCALAR, PREF_ALL, maxdim, finiterange, monotonicity);
  return currentNrCov - 1;
}
int IncludePrim(const char *name,Types type,  int kappas, size_fct kappasize,
		domain_type domain, isotropy_type isotropy,	
		checkfct check, rangefct range, int maxdim, 
		ext_bool finiterange, int monotonicity) {  
  createmodel(name, type, kappas, kappasize, domain, isotropy, check, range,
	      SCALAR, PREF_ALL, maxdim, finiterange, monotonicity);
  return currentNrCov - 1;
}

int IncludePrim(const char *name, Types type,  int kappas, 
		domain_type domain, isotropy_type isotropy,	
		checkfct check, rangefct range, int vdim, 
		int maxdim, ext_bool finiterange, int monotonicity) {  
  createmodel(name, type, kappas, NULL, domain, isotropy, check, range,
	      vdim, PREF_ALL, maxdim, finiterange, monotonicity);
  return currentNrCov - 1;
}
int IncludePrim(const char *name,Types type,  int kappas, size_fct kappasize,
		domain_type domain, isotropy_type isotropy,	
		checkfct check, rangefct range, int vdim, 
		int maxdim, ext_bool finiterange, int monotonicity) {  
  createmodel(name, type, kappas, kappasize, domain, isotropy, check, range,
	      vdim, PREF_ALL, maxdim, finiterange, monotonicity);
  return currentNrCov - 1;
}

int IncludePrim(const char *name,Types type,  int kappas, 
		domain_type domain, isotropy_type isotropy,	
		checkfct check, rangefct range, pref_type pref,
		int vdim, int maxdim, ext_bool finiterange, 
		int monotonicity) {  
  createmodel(name, type, kappas, NULL, domain, isotropy, check, range,
	      vdim, pref, maxdim, finiterange, monotonicity);
  return currentNrCov - 1;
}
int IncludePrim(const char *name,Types type,  int kappas, size_fct kappasize,
		domain_type domain, isotropy_type isotropy,	
		checkfct check, rangefct range, pref_type pref,
		int vdim, int maxdim, ext_bool finiterange, int monotonicity) { 
  createmodel(name, type, kappas, kappasize, domain, isotropy, check, range,
	      vdim, pref, maxdim, finiterange, monotonicity);
  return currentNrCov - 1;
}


void make_internal() {  
  int nr = currentNrCov - 1;  
  cov_fct *C = CovList + nr; // nicht gatternr
  C->internal = true; 
}






// extern ?!
int IncludeModel(const char *name, Types type, 
		 char minsub, char maxsub, int kappas,
		 size_fct kappasize,
		 domain_type domain, isotropy_type isotropy,
		 checkfct check, rangefct range, pref_type pref, 
		 bool internal, int vdim, int maxdim, ext_bool finiterange,
		 int monotonicity) {  
  createmodel(name, type, kappas, kappasize, domain, isotropy, check, range,
	      vdim, pref, maxdim, finiterange, monotonicity);
  //    assert(maxsub > 0); // check deleted 25. nov 2008 due to nugget 
  assert(maxsub >= minsub && maxsub <= MAXSUB);
  assert(check != checkOK || maxsub==0);
  int i, 
    nr = currentNrCov - 1;  
  cov_fct *C = CovList + nr; // nicht gatternr
  C->minsub = minsub;
  C->maxsub = maxsub;  

  /*
    if (maxsub == 0) printf("not primitive: %s \n", name); //
    not primitive: biWM 
    not primitive: constant 
    not primitive: epsC 
    not primitive: matern 
    not primitive: nugget 
    not primitive: trend 
    not primitive: whittle 
    (and others?)
  */

  assert(minsub <= maxsub);
  if (PL>=PL_SUBDETAILS && maxsub == 0) 
    PRINTF("note: %s has no submodels\n", C->name);
  C->primitive = false; // muss falls sein, sonst wird kappacheck 
  // aufgerufen
  C->internal = internal;

  //if (internal) printf("internal: %s\n", C->name);

  if (maxsub <= 2) {
    if (maxsub >= 1) {
      addsub(0, "phi");
      if (maxsub >= 2) {
	addsub(1, "psi");
      }
    }
  } else {
    for (i=0; i<maxsub; i++) {      
      sprintf(C->subnames[i], "C%d", i); // default (repeated twice)
      C->subintern[i] = false;
    }
  }
  return nr;
}

int IncludeModel(const char *name, Types type, 
		 char minsub, char maxsub, int kappas, 
		 domain_type domain, isotropy_type isotropy,
		 checkfct check, rangefct range, pref_type pref, 
		 int maxdim, ext_bool finiterange, int monotonicity) {
  return
    IncludeModel(name, type, minsub, maxsub, kappas,
		 NULL, domain, isotropy, check, 
		 range, pref, false,
		 SCALAR, maxdim, finiterange, monotonicity);
}

#define IMPLEMENTED_SEQUENTIAL C->vdim <= 1 &&				\
     (C->domain == PREVMODELD ||	\
      ((isPosDef(C->Type) || isUndefinedType(C->Type)) && C->domain == XONLY))

#define IMPLEMENTED_CE							\
  (((isPosDef(C->Type) || isUndefinedType(C->Type)) && C->domain == XONLY) \
   || C->domain == PREVMODELD)

void addCov(int F_derivs, covfct cf, covfct D, covfct D2, covfct inverse,
	    nonstat_inv nonstat_inverse) {
  int nr = currentNrCov - 1;

  assert(nr>=0 && nr<currentNrCov && cf!=NULL);
  cov_fct *C = CovList + nr; // nicht gatternr
  bool stat_iso = C->domain == XONLY && C->isotropy == ISOTROPIC;

  C->cov = cf;  
  if (C->RS_derivs < 0) C->RS_derivs = 0;
  assert(C->nonstat_cov == ErrCovNonstat);

  if (D != NULL) {
    assert(cf != NULL);
    if (C->cov!=NULL && C->RS_derivs < 1) C->RS_derivs = 1;
    C->D=D;    
    assert(C->isotropy == ISOTROPIC || C->isotropy == SPACEISOTROPIC ||
	   C->isotropy == PREVMODELI || C->isotropy == PARAM_DEP);
    // C->implemented[TBM2] = NUM_APPROX;
    C->implemented[TBM] = true; 
  }
  if (D2 != NULL) {
    assert(D != NULL);
    C->D2 = D2;
   if (C->cov!=NULL && C->D != NULL && C->RS_derivs < 2)  C->RS_derivs = 2;
  }
  if (inverse != NULL) C->inverse = inverse;
  else if (isMonotone(C->Monotone) && C->isotropy == ISOTROPIC && 
	   C->inverse==ErrCov)
    C->inverse = InverseIsoMon;

  C->nonstat_inverse = nonstat_inverse!=NULL ? nonstat_inverse :
    stat_iso && inverse != NULL ? StandardInverseNonstat : ErrInverseNonstat;
  C->implemented[Direct] = cf != NULL;
  C->implemented[CircEmbed] = cf != NULL && IMPLEMENTED_CE;

  // printf("%s %d %d\n", C->nick, C->implemented[CircEmbed], C->pref[CircEmbed]);

  C->implemented[Sequential] = IMPLEMENTED_SEQUENTIAL;
  C->F_derivs = F_derivs >= 0 ?  F_derivs : C->RS_derivs;
  assert(C->F_derivs <= C->RS_derivs);
}

void addCov(covfct cf, covfct D, covfct D2, covfct inverse) {
  addCov(MISMATCH, cf, D, D2, inverse, NULL);
}

void addCov(covfct cf, covfct D, covfct D2, nonstat_inv inverse) {
  addCov(MISMATCH, cf, D, D2, NULL, inverse);
}

void addCov(covfct cf, covfct D, covfct D2, 
	    covfct inverse, nonstat_inv nonstat_inverse) {
  addCov(MISMATCH, cf, D, D2, inverse, nonstat_inverse);
}

void addCov(int F_derivs, covfct cf, covfct D, covfct inverse) {
  addCov(F_derivs, cf, D, NULL, inverse, NULL);
}


void addCov(covfct cf, covfct D, covfct inverse) {
  addCov(MISMATCH, cf, D, inverse);
}

void addCov(int F_derivs, covfct cf, covfct D, covfct D2, covfct D3, covfct D4,
	    covfct inverse) {
  addCov(cf, D, D2, inverse, NULL);
  cov_fct *C = CovList + currentNrCov - 1;
  C->D3 = D3;
  C->D4 = D4;
  assert(C->RS_derivs == 2 && D3!=NULL && D4!=NULL);
  C->RS_derivs = 4;
  C->F_derivs = F_derivs >= 0 ?  F_derivs : C->RS_derivs;
  assert(C->F_derivs <= C->RS_derivs);
}

void addCov(covfct cf, covfct D, covfct D2, covfct D3, covfct D4,
	    covfct inverse) {
  addCov(-1, cf, D, D2, D3, D4, inverse);
}

void addCov(int F_derivs, nonstat_covfct cf) {
  int nr = currentNrCov - 1;
  assert((nr>=0) && (nr<currentNrCov) && cf!=NULL);
  cov_fct *C = CovList + nr; // nicht gatternr

  C->nonstat_cov = cf;
  C->implemented[Direct] = true;
  C->implemented[CircEmbed] = IMPLEMENTED_CE;
  C->implemented[Sequential] = IMPLEMENTED_SEQUENTIAL;

  if (C->RS_derivs < 0) {
    C->RS_derivs = 0;
    C->F_derivs = F_derivs >= 0 ?  F_derivs : C->RS_derivs;
    assert(C->F_derivs <= C->RS_derivs);
  }
}
void addCov(nonstat_covfct cf) {
  addCov(-1, cf);
}

void addCov(aux_covfct auxcf){
  int nr = currentNrCov - 1;
  assert((nr>=0) && (nr<currentNrCov) && auxcf!=NULL);
  cov_fct *C = CovList + nr; // nicht gatternr
  assert(C->cov == ErrCov && C->nonstat_cov==ErrCovNonstat);
  C->aux_cov = auxcf;
}

void addCov(covfct distrD, logfct logdistrD, nonstat_inv Dinverse,
	    covfct distrP, nonstat_covfct distrP2sided,
	    covfct distrQ, covfct distrR, nonstat_covfct distrR2sided) {
  int nr = currentNrCov - 1;
  cov_fct *C = CovList + nr; // nicht gatternr
  assert(nr>=0 && nr<currentNrCov);
  assert(CovList[nr].domain == STAT_MISMATCH);
  assert(CovList[nr].isotropy == ISO_MISMATCH);
  assert(distrD != NULL && logdistrD!=NULL && Dinverse != NULL &&
	 distrP != NULL && distrQ != NULL &&
	 distrR != NULL);

  C->RS_derivs = 1;
  C->F_derivs = 0;

  C->cov = distrP;
  C->D = distrD; 
  C->logD = logdistrD;
  C->nonstat_inverse_D = Dinverse;
  C->inverse = distrQ;
  C->random = distrR;
  if (distrP2sided != NULL || distrR2sided !=NULL) {
    assert(distrP2sided != NULL && distrR2sided !=NULL);
    C->nonstat_cov = distrP2sided;
    C->nonstat_random = distrR2sided;
    C->domain = PREVMODELD;
  } else {
    assert(distrP2sided == NULL && distrR2sided ==NULL);
    C->domain = XONLY;
  }
  
  C->isotropy = NO_ROTAT_INV;
}

void addlogCov(logfct log, nonstat_logfct nonstatlog) {
  int nr = currentNrCov - 1;
  assert((nr>=0) && (nr<currentNrCov));
  cov_fct *C = CovList + nr; // nicht gatternr
  if (log != NULL) C->log = log;
  else assert(nonstatlog != NULL);
  if (nonstatlog != NULL) C->nonstatlog = nonstatlog;
}
 
void addlogCov(logfct log) {
  addlogCov(log, NULL);
}

void addlogCov(nonstat_logfct nonstatlog) {
  addlogCov(NULL, nonstatlog);
}


int addFurtherCov(int F_derivs, covfct cf, covfct D, covfct D2) {
  assert(currentNrCov > 0);
  cov_fct *C = CovList + currentNrCov;
  MEMCOPY(C, C - 1, sizeof(cov_fct));
  assert(C->vdim == SCALAR || C->vdim == SUBMODEL_DEP ||
	 C->vdim == PARAM_DEP || D == NULL);

  strcopyN(CovNames[currentNrCov], InternalName, MAXCHAR);
  C->name[0] = InternalName[0];
  strcopyN(C->name + 1, CovList[currentNrCov-1].name, MAXCHAR - 1);
  if (cf != NULL) {
    C->cov = cf;
    C->RS_derivs = 0;
  }
  if (D != NULL) {
    assert(cf != NULL);
    C->D = D;
    C->RS_derivs = 1;

    assert(C->isotropy == ISOTROPIC || C->isotropy == SPACEISOTROPIC ||
	   C->isotropy == ZEROSPACEISO ||
	   C->isotropy == PREVMODELI);
 
    //C->implemented[TBM2] = NUM_APPROX;
    C->implemented[TBM] = IMPLEMENTED; 
  }
  if (D2 != NULL) {
    assert(D != NULL);
    C->D2 = D2;
    C->RS_derivs = 2;
  }
  C->F_derivs = F_derivs >= 0 ?  F_derivs : C->RS_derivs;
  assert(C->F_derivs <= C->RS_derivs);
  C->internal = true; // addCov is used without previous call of IncludeModel
  currentNrCov++;
  return currentNrCov - 1;
}
int addFurtherCov(covfct cf, covfct D, covfct D2) {
  return addFurtherCov(-1, cf, D, D2);
}

int addFurtherCov(covfct cf, covfct D) {
  return addFurtherCov(cf, D, NULL);
}
int addFurtherCov(int F_derivs, covfct cf, covfct D) {
  return addFurtherCov(F_derivs, cf, D, NULL);
}
int addFurtherCov(int F_derivs, nonstat_covfct cf) {
  assert(currentNrCov > 0);
  cov_fct *C = CovList + currentNrCov;
  MEMCOPY(C, C - 1, sizeof(cov_fct));
  strcopyN(CovNames[currentNrCov], InternalName, MAXCHAR);
  C->name[0] = InternalName[0];
  strcopyN(C->name + 1, CovList[currentNrCov-1].name, MAXCHAR - 1);
  C->RS_derivs = -1;
  if (cf != NULL) {
    C->nonstat_cov = cf;
    C->RS_derivs = 0;
  }
  C->F_derivs = F_derivs >= 0 ?  F_derivs : C->RS_derivs;
  assert(C->F_derivs <= C->RS_derivs);
  C->D = ErrCov;
  C->internal = true; // addCov is used without previous call of IncludeModel
  currentNrCov++;
  return currentNrCov - 1;
}

int addFurtherCov(nonstat_covfct cf) {
  return addFurtherCov(-1, cf);
}


void addTypeFct(type_fct TypeFct) {
  int nr = currentNrCov - 1;
  cov_fct *C = CovList + nr; // nicht gatternr
  C->TypeFct = TypeFct;
}


//int addFurtherCov(nonstat_covfct cf, covfct D, covfct D2) {
//  return addFurtherCov(-1, cf, D, D2);
//}

void nablahess(covfct nabla, covfct hess) {
  int nr = currentNrCov - 1;
  cov_fct *C = CovList + nr; // nicht gatternr

  assert((nr>=0) && (nr<currentNrCov));
  assert(C->vdim == SCALAR || C->vdim == SUBMODEL_DEP || 
	 C->vdim == PARAM_DEP || nabla==NULL);
  assert(C->cov != NULL && nabla!=NULL && hess != NULL);
  
  C->nabla=nabla;    
  C->hess = hess;

}

void addLocal(getlocalparam coinit, getlocalparam ieinit) {
  int nr = currentNrCov - 1;
  assert(nr>=0 && nr < currentNrCov) ;
  cov_fct *C = CovList + nr; // nicht gatternr
  assert(C->D!=ErrCov);
  if ((C->implemented[CircEmbedIntrinsic] = ieinit != NULL)) {
    assert(C->D2 != NULL);
    C->ieinit = ieinit;
  } 
  if ((C->implemented[CircEmbedCutoff] = coinit != NULL)) {
    C->coinit = coinit;
  }
}
void addCallLocal(altlocalparam alt) {
  int nr = currentNrCov - 1;
  assert(nr>=0 && nr < currentNrCov) ;
  cov_fct *C = CovList + nr; // nicht gatternr
  C->alternative = alt;
}

#define ASSERT_TBM  \
  assert(C->domain == PREVMODELD ||					\
     ((isPosDef(C->Type)|| isUndefinedType(C->Type)) && C->domain == XONLY)); \
  assert(C->vdim == SCALAR || C->vdim == SUBMODEL_DEP ||		\
	 C->vdim == PARAM_DEP)

int addTBM(covfct tbm2) {
  // must be called always AFTER addCov !!!!
  int nr = currentNrCov - 1;
  assert((nr>=0) && (nr<currentNrCov));
  cov_fct *C = CovList + nr; // nicht gatternr
  ASSERT_TBM;
  C->tbm2=tbm2;
  if (tbm2 != NULL) {
    // addTBM is called from the other addTBM's -- so tbm2 might
    // be NULL
    assert(C->isotropy==ISOTROPIC || C->isotropy==SPACEISOTROPIC || 
	   C->isotropy==PREVMODELI);
    assert(C->D != ErrCov);
    C->implemented[TBM] = IMPLEMENTED;
  }
  // IMPLEMENTED must imply the NUM_APPROX to simplify the choice
  // between TBM2 and Tbm2Num
  return nr;
}

void addTBM(covfct tbm2, initfct Init, spectral_do spectral) {
  int nr = addTBM(tbm2);
  cov_fct *C = CovList + nr; // nicht gatternr
  ASSERT_TBM;
  C->spectral=spectral;
  C->Init=Init;
  C->implemented[SpectralTBM] = true;
}

void addTBM(initfct Init, spectral_do spectral) {
  addTBM((covfct) NULL, Init, spectral);
}
	
	
void addHyper(hyper_pp_fct hyper_pp) {
  int nr = currentNrCov - 1;
  cov_fct *C = CovList + nr; // nicht gatternr
  assert((nr>=0) && (nr<currentNrCov));
  C->hyperplane=hyper_pp;
  C->implemented[Hyperplane] = hyper_pp!=NULL;
}
		   
//void addSpecialMeth(initstandard initspecial, dometh special)  {
///  int nr = currentNrCov - 1;
//  cov_fct *C = CovList + nr; // nicht gatternr
//  C->initspecial=initspecial;
//  C->special=special;
//  if ((special!=NULL) || (initspecial!=NULL)) 
//    assert((special!=NULL) && (initspecial!=NULL));
//  C->implemented[Special] = true;
//}

void addMarkov(int *variable) {
  int nr = currentNrCov - 1;
  cov_fct *C = CovList + nr; // nicht gatternr
  C->implemented[Markov] = true; 
  *variable = nr;
}

void RandomShape(int maxmoments, structfct Struct, initfct Init, 
		 dofct Do, do_random_fct DoRandom, 
		 bool average, bool randomcoin, bool specific){
 // rein auf shape-Function bezogen, ohne Kenntnis von irgendeinem
  // Kovarianzmodell
  //
  // init and do der elementaren shape-Funktion
  // insbesondere notwendig wenn die elementare shape-Function zufaellig ist.
  // Erst dann koennen werte abgefragt werden
  int nr = currentNrCov - 1;
  cov_fct *C = CovList + nr; // nicht gatternr

  assert(Init != NULL && Do != NULL && DoRandom != NULL && Struct!=NULL);

  C->maxmoments = maxmoments;
  C->Struct = Struct;
  C->Do = Do;
  C->DoRandom = DoRandom;
  C->Init = Init;
  C->implemented[Average] = average; 
  C->implemented[RandomCoin] = randomcoin; 
  C->implemented[Specific] = specific; 
}

void RandomShape(int maxmoments, structfct Struct, initfct Init, dofct Do,
		 bool average, bool randomcoin, bool specific){
 RandomShape(maxmoments, Struct, Init, Do, do_random_failed,
	     average, randomcoin, specific);
}

void RandomShape(int maxmoments, structfct Struct, initfct Init,
		 dofct Do){
  RandomShape(maxmoments, Struct, Init, Do, do_random_failed, 
	      false, false, false);
}

void RandomShape(int maxmoments, structfct Struct, initfct Init, 
		 do_random_fct DoRandom){
  RandomShape(maxmoments, Struct, Init, do_failed, DoRandom, 
	      false, false, false);
}

void RandomShape(int maxmoments, initfct Init, dofct Do, 
		 bool average){
  RandomShape(maxmoments, structOK, Init, Do, do_random_failed,
	      average, !average, false);
}

void RandomShape(int maxmoments, initfct init, dofct Do){
   RandomShape(maxmoments, structOK, init, Do, do_random_failed, 
 	      false, true, false); 
}

void RandomShape(structfct Struct, bool average){
   RandomShape(-1, Struct, initOK, doOK,  do_random_failed, 
	      average, !average, false);
}

void RandomShape(structfct Struct){
  RandomShape(-1, Struct, initOK, doOK, do_random_failed, false, true, false);
}


void addSpecial(minmaxfct minmaxeigen){
  int nr = currentNrCov - 1;
  cov_fct *C = CovList + nr; // nicht gatternr
  C->minmaxeigenvalue = minmaxeigen;
}

void addGaussMixture(draw_random drawmix,
		     log_mixdens logmixdens) {
  int nr = currentNrCov - 1;
  cov_fct *C = CovList + nr; // nicht gatternr
  assert(drawmix != NULL && logmixdens != NULL);
  C->drawmix = drawmix;
  C->logmixdens = logmixdens;
}

void addReturns(return_fct Covariance, ext_bool_ret_fct isCovariance, 
		return_fct CovMatrix, ext_bool_ret_fct isCovMatrix,
		return_fct InverseCovMatrix,ext_bool_ret_fct isInverseCovMatrix,
		return_fct Variogram, ext_bool_ret_fct isVariogram,
		return_fct PseudoVariogram, ext_bool_ret_fct isPseudoVariogram,
		returnX_fct SelectedCovMatrix,
		ext_bool_ret_fct isSelectedCovMatrix) {
  int nr = currentNrCov - 1;
  cov_fct *C = CovList + nr; // nicht gatternr
  if (Covariance!=NULL) {
    C->covariance=Covariance;
    C->is_covariance = isCovariance==NULL ? isTrue : isCovariance;
  } else assert(isCovariance==NULL);
  if (CovMatrix!=NULL) {
    C->covmatrix=CovMatrix;
    C->is_covmatrix = isCovMatrix==NULL ? isTrue : isCovMatrix;
  } else assert(isCovMatrix==NULL);
  if (InverseCovMatrix!=NULL) {
    C->inversecovmatrix=InverseCovMatrix;
    C->is_inversecovmatrix =
      isInverseCovMatrix==NULL ? isTrue : isInverseCovMatrix;
  } else assert(isInverseCovMatrix==NULL);
  if (Variogram!=NULL) {
    C->variogram=Variogram;
    C->is_variogram = isVariogram==NULL ? isTrue : isVariogram;
  } else assert(isVariogram==NULL);
  if (PseudoVariogram!=NULL) {
    C->pseudovariogram=PseudoVariogram;
    C->is_pseudovariogram =
      isPseudoVariogram==NULL ? isTrue : isPseudoVariogram;
  } else assert(isPseudoVariogram==NULL);
  if (SelectedCovMatrix!=NULL) {
    C->selectedcovmatrix=SelectedCovMatrix;
    C->is_selectedcovmatrix = 
      isSelectedCovMatrix==NULL ? isTrue : isSelectedCovMatrix;
  } else assert(isSelectedCovMatrix==NULL);
}


void Taylor(double c, double pow) {
  int 
    nr = currentNrCov - 1;
  cov_fct *C = CovList + nr; // nicht gatternr

  assert(pow != 0.0);

  C->TaylorN = 0;
  if (isPosDef(C->Type) || isUndefinedType(C->Type)) {
    C->Taylor[C->TaylorN][TaylorConst] = 1.0;
    C->Taylor[C->TaylorN][TaylorPow] = 0.0;    
    C->TaylorN++;
  } 

  C->Taylor[C->TaylorN][TaylorConst] = c;
  C->Taylor[C->TaylorN][TaylorPow] = pow;
  C->TaylorN++;
  assert(C->TaylorN <= MAXTAYLOR);

  if (C->finiterange == true) TailTaylor(0, 0, 0, 0);

}


void Taylor(double c, double pow, double c1, double pow1) {
  int  
    nr = currentNrCov - 1;
  cov_fct *C = CovList + nr; // nicht gatternr

  Taylor(c, pow);
  C->Taylor[C->TaylorN][TaylorConst] = c1;
  C->Taylor[C->TaylorN][TaylorPow] = pow1;
  C->TaylorN++;
  assert(C->TaylorN <= MAXTAYLOR);
}

void Taylor(double c, double pow, double c1, double pow1, 
	    double c2, double pow2) {
  int 
    nr = currentNrCov - 1;
  cov_fct *C = CovList + nr; // nicht gatternr

  Taylor(c, pow, c1, pow1);
  C->Taylor[C->TaylorN][TaylorConst] = c2;
  C->Taylor[C->TaylorN][TaylorPow] = pow2;
  C->TaylorN++;

  assert(C->TaylorN <= MAXTAYLOR);

}

void TailTaylor(double t, double tpow, double texpc, double texppow) {
  int
    nr = currentNrCov - 1;
  cov_fct *C = CovList + nr; // nicht gatternr

  if (C->finiterange == true) {
    assert(t == 0 && tpow==0 && texpc==0 && texppow==0);
  }


  C->TailN = 0;
  C->Tail[C->TailN][TaylorConst] = t;
  C->Tail[C->TailN][TaylorPow] = tpow;
  C->Tail[C->TailN][TaylorExpConst] = texpc;
  C->Tail[C->TailN][TaylorExpPow] = texppow;
  C->TailN++;

  assert(C->TailN <= MAXTAYLOR);
}

