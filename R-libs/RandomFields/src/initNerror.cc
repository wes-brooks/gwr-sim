/* 
 Authors
 Martin Schlather, schlather@math.uni-mannheim.de

 library for simulation of random fields -- init part and error messages

 Copyright (C) 2001 -- 2013 Martin Schlather

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


#include <math.h>  
#include <stdio.h>  
#include <stdlib.h>
//#include <sys/timeb.h>
 
#include <string.h>
#include "RF.h"
#include "Covariance.h"
#include <unistd.h>


int gaussmethod[Forbidden+1];

cov_model *KEY[MODEL_MAX+1];
double  ZERO[MAXSIMUDIM], 
  ONE = 1,
//    *userdefinedCovMatrix[MAXDEFMATRIX][MAXMAKEEXPLICITE],
    *OutX=NULL, 
    *OutY=NULL;
int GENERALISEDCAUCHY, STABLE,  BROWNIAN, CAUCHY, 
  GAUSS, NUGGET, PLUS, TBM2NR, BALL, ECF, MULT,
  DISTRIBUTION, DETERM_DISTR, GAUSS_DISTR, SETPARAM, COVFCTN,
  COVMATRIX, RFGET, STROKORB_MONO, STROKORB_BALL_INNER, RECTANGULAR,
  POLYGON,
  MULT_INVERSE,
  TRUNCSUPPORT, SHAPESTP, SHAPEAVE, BROWNRESNICK, UNIF, MPPPLUS, CUTOFF, STEIN,
  BRSHIFTED_USER, BRMIXED_USER, BRORIGINAL_USER, 
  BRSHIFTED_INTERN, BRMIXED_INTERN, BRORIGINAL_INTERN,   
  EXTREMALGAUSSIAN, RANDOMSIGN,  
  ARCSQRT_DISTR,
  PTS_GIVEN_SHAPE, STATIONARY_SHAPE, STANDARD_SHAPE,
  LOC, SCALESPHERICAL, TBM_OP, USER,
  MIXEDEFFECT, // MLEMIXEDEFFECT,
  VARIOGRAM_CALL, 
  MISSING_COV,
  DOLLAR_PROC, PLUS_PROC,
  BINARYPROC, BROWNRESNICKPROC,
  GAUSSPROC, POISSONPROC,  SCHLATHERPROC, SMITHPROC, CHI2PROC,
  NUGGET_USER, NUGGET_INTERN,
  CIRCEMBED,  SPECTRAL_PROC_USER, SPECTRAL_PROC_INTERN,
  DIRECT, SEQUENTIAL, SPECIFIC, SELECT,
  AVERAGE_USER, AVERAGE_INTERN, HYPERPLANE_USER, HYPERPLANE_INTERN,
  RANDOMCOIN_USER, CE_CUTOFFPROC_USER, CE_CUTOFFPROC_INTERN, 
  CE_INTRINPROC_USER, CE_INTRINPROC_INTERN, TBM_PROC_USER, TBM_PROC_INTERN, 
  VECTOR,
  ISO2ISO, SP2SP, SP2ISO, S2ISO, S2SP, S2S, SId, FIRST_GATTER, LAST_GATTER;
// userdefinedCM_RC[MAXDEFMATRIX][MAXMAKEEXPLICITE], 
bool
    NAOK_RANGE=false;
char CovNames[MAXNRCOVFCTS][MAXCHAR], CovNickNames[MAXNRCOVFCTS][MAXCHAR];
char MSG[1000], NEWMSG[1000], BUG_MSG[250];

cov_fct *CovList=NULL;
int currentNrCov=-1,
  CONSTANT = -1,
  OUT = -1,
  DOLLAR = -1,
  TREND = -1,
  LASTDOLLAR = -1,
  NATSC = -1;

int True=1; // never change
int False=0; // never change
char *FREEVARIABLE= (char*) "...";
#define MAX_CE_MEM 16777216
#define R_PRINTLEVEL 1
#define C_PRINTLEVEL 1
#define NAT_SCALE 0
int PL=C_PRINTLEVEL, 
    NS=NAT_SCALE;


globalparam GLOBAL = {
  {'.', false, false, false, false, false, true,
   normal, R_PRINTLEVEL, C_PRINTLEVEL, NAT_SCALE, 1, 0, 
   0, MODEL_KRIGE, MODEL_COND, MODEL_ERR, MODEL_GUI,
   {0, 2} /* chol, SVD */,
   NA_REAL,  1e-6, 1e-6, NA_REAL
  },// general;
  {NA_REAL, 0.05, false, false, 800},  // gauss
  {'A', false, false, true,
   8000, {5000, 200, 1000}, 2}, // krige
  {false, true, false, TRIVIALSTRATEGY, 3, MAX_CE_MEM, MAX_CE_MEM,
   -1e-7, 1e-3, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
   1.0}, //ce_param ce, 13 NULLEN
  {false, -2, 3, 0, {1, 60, 500}, RF_NAN, 2.0, 0.0, 
   {RF_NAN, RF_NAN, RF_NAN, RF_NAN}},  // TBM 
  {true, false, 50, 0.0, {2500, 2500, 2500, 2500}},  // spectral_param
  {Cholesky, 1e-12, 8192 },// direct_param direct;
  {5000, 5, -10}, //sequ_param sequ;
  {2, false, 250000, 1.0}, //markov_param markov; 
  {}, // ave
  {0.0}, // nugget_param nugget;
  {50000, // mpp,   n_estim_E 
   {1.0, 1.0, 1.0, 1.0}, // intens; 
    1e-4 /* about zero */
  },  //mpp_param mpp;
  {700, 1000, HYPER_UNIFORM, RF_NAN}, // hyper_param hyper;
  {}, //special_param special;
  {MAXINT, 30, FLAT_UNDETERMINED, // int
   4.0,  // oder 5 !; jedoch ist 3 zu wenig
   1.0, 0.0},//extremes_storage extremes; maxstable
  {10000000, 7, 2000000000, 0, 0.1, -2.0, 0.01, 5.0}, // br (BrownResnick)
  {0.08, 0.1, 1e-20, 1e5, 100, 8, 20, 15, 1000}, // distr (rectangular) // todo should be 500 and better algorithm for approximation!
  {0.5, 3.0, 5.0, 3.0, 10.0, 1e4, // 6
   1E-10, 1000.0, 0.001, 0.02, 1/1000, 1000,  // 12
   false, -10.0, 10.0, NA_REAL, NA_REAL, 0.1,  // 18 
   1e-7, 
   50, 20, 1, 1, 20, 1,            // 6
   10, 5000, {3000, 200, 1000}, 2, 2000,  // 12
   0, -1, 1, // besser 3, 0, 2
   true, true, false, true, true,
   4}, // fit
  {0.0, 0.0, 1e-13, false, true}, // empvario (automatically chosen as -pi/n/2))
  {true, Nothing, {1024, 64}},  // gui
  {true, 1}, // graphics
  {true, true, true, true, true, 
   true, true, false}, // warnings
};
 
//globalorig GLOBALORIG = {false, {}};
int PrInL=-1;				


pref_type PREF_ALL = {PREF_BEST, PREF_BEST, PREF_BEST, PREF_BEST, PREF_BEST, 
		      PREF_BEST, PREF_BEST, PREF_BEST, PREF_BEST, PREF_BEST,
		      PREF_BEST, PREF_BEST, PREF_NONE, // specific
		                                       PREF_BEST, PREF_BEST},
  PREF_NOTHING = {PREF_NONE, PREF_NONE, PREF_NONE, PREF_NONE, PREF_NONE, 
		  PREF_NONE, PREF_NONE, PREF_NONE, PREF_NONE, PREF_NONE,
		  PREF_NONE, PREF_NONE, PREF_NONE, PREF_BEST, // nothing
		                                              PREF_NONE},
  PREF_AUX = {PREF_NONE, PREF_NONE, PREF_NONE, PREF_NONE, PREF_NONE, 
	      PREF_NONE, PREF_NONE, PREF_NONE, PREF_NONE, PREF_NONE,
	      PREF_NONE, PREF_NONE, PREF_NONE, PREF_NONE, PREF_NONE};



/*
  record containing all user specified directions for choosing
  a simulation method. Currently only one element is contained.
*/

double GENERAL_PRECISION = 5e-15;
double EIGENVALUE_EPS = 1e-15;

char ERRORSTRING[MAXERRORSTRING], ERRORSTRING_OK[MAXERRORSTRING], 
  ERRORSTRING_WRONG[MAXERRORSTRING],
    ERROR_LOC[nErrorLoc]="";
int ERRORMODELNUMBER = -1;

char PREF_FAILURE[100 * Nothing];


const char *METHODNAMES[Forbidden+1]={"circulant", //0
				      "cutoff",
				      "intrinsic",
				      "tbm", 
				      "spectral", //4
				      "direct",
				      "sequential",
				      "markov",
				      "average",
				      "nugget", //9
				      "coins",
				      "hyperplane",
				      "specific",
				      "any method", // nothing
				      "forbidden"},
  *STATNAMES[LAST_STAT + 1] = {
    "single variable", "kernel", "calling model", "mismatch"},
  *ISONAMES[LAST_ISO + 1] =  { 
    "isotropic", "space-isotropic", "zero-space-isotropic", "vector-isotropic",
    "symmetric", "no rotation invariance", "parameter dependent", "<mismatch>"},
  *ROLENAMES[ROLE_LAST + 1] = {
    "<none>",                                                     // 0
    "covariance model", "Gauss", "max-stable", "BrownResnick", "Smith",  // 5
    "Schlather", "Poisson", "PoissonGauss", "Bernoulli", "distribution", // 10
    "<rotten>", "<undefined>"},
  *TYPENAMES[OtherType + 1] = {
    "tail correlation function", "positive definite", "negative definite", 
    "process", 
    "method for Gauss processes", "method for Brown-Resnick processes",
    "point-shape function",
    "distribution family", "shape function", "trend", "interface",
    "undefined", "other type"},
  *MONOTONE_NAMES[BERNSTEIN + 1 - MISMATCH] = {
    "mismatch in monotonicity", "submodel dependent monotonicity",
    "previvous model dependent monotonicity",
    "parameter dependent monotonicity",
    "not monotone", "monotone", "Gneiting-Schaback class", 
    "completely monotone",  "normal mixture", 
    "Bernstein"     
  },
  *CAT_TYPENAMES[OtherType + 1] = {
    // TcfType, PosDefType, NegDefType, ProcessType, GaussMethodType, 
    // BrMethodType, PointShapeType, RandomType, ShapeType, TrendType, 
    // InterfaceType,
    // UnDefinedType, OtherType
    "RM", "RM", "RM", "RP", "RP",
    "RP", "RM", "RR", "RM", "RM", 
    "RF", 
    "RM", "RO"},
  *REGNAMES[MODEL_MAX+1] = {"reg0", "reg1", "reg2", "reg3", "reg4", 
			    "reg5", "reg6", "reg7", "reg8", "reg9",
			    "user", "unused", "intern", "split", "gui",
			    "mle", "mlesplit", "mletrend", "mlebounds",
			    "kriging", "conditional", "error model"},
  *MODENAMES[nr_modes] = {"careless", "sloppy", "easygoing", "normal", 
			  "precise", "pedantic", "neurotic"};  

char
  STANDARDPARAM[MAXPARAM][MAXCHAR],
  STANDARDSUB[MAXSUB][MAXCHAR];

bool RELAX_UNKNOWN_RFOPTION=false; // auf keinen Fall aendern!

void errorMSG(int err, char* m, int len) {
  if (err >= ERRORM && err <= ERRORMEND) err = ERRORM;

  switch (err) {
  case NOERROR : strcpy(m,"none"); break;
  case NOERROR_REPEAT : strcpy(m,"none; looking for further covariances applicable to the same method");break;
  case NOERROR_ENDOFLIST : strcpy(m,"none; end of list");break;
  case ERRORDUMMY : strcpy(m,"none (dummy)"); break;
  case ERRORNOTDEFINED :       
    strcpy(m,"specified method undefined for the given model or no simulation method found for the given model");break;
  case ERRORNOTPROGRAMMED :    
    strcpy(m,"not programmed yet. Sorry."); break;
  case ERRORVDIMNOTPROGRAMMED :    
    strcpy(m,"multivariate version not programmed yet. Sorry."); break;
  case ERRORTYPECONSISTENCY :
     strcpy(m,"incorrect choice of submodel (type inconsistency)"); break;
  case ERRORFAILED: 
   strcpy(m,"algorithm failed (partially)");break;
  case ERRORMEMORYALLOCATION: 
    strcpy(m,"memory allocation error"); break;
  case ERRORNOTINITIALIZED: 
    strcpy(m,"not initialized or storing=FALSE");break;
  case ERRORDECOMPOSITION:
    strcpy(m,"covariance function does not seem to be (strictly) positive definite");break;
  case ERRORCOVFAILED: 
    sprintf(m, "model and method only valid for %s. Got %s",
	    ERRORSTRING_OK,ERRORSTRING_WRONG);
    break;
  case ERRORNOMULTIVARIATE :
    strcpy(m, "multivariate models not allowed (yet)"); 
    break;
  case ERROR_MATRIX_SQUARE :
    strcpy(m, "square matrix expected"); break;
  case ERROR_MATRIX_VDIM :
    strcpy(m, "size of matrix is not a multiple of the multivariate dimension"); break;
  case ERROR_MATRIX_POSDEF :
    strcpy(m, "matrix does not seem to be strictly positive definite"); break;
    //  case ERROR_MATRIX_ :   strcpy(m, ""); break;
  case ERRORDIM: 
    //
    //    { printf("error dimension\n"); cov_model *cov; crash(cov); }
    sprintf(m,"dimension specification not in [1,%d] or dimension of coordinates larger than that the supposed spatio-temporal process",
	    MAXSIMUDIM);break;
  case ERRORWAVING :
    strcpy(m,"Rescaling not possible (waving or large nugget effect?)");break;
  case ERRORRESCALING:
    strcpy(m,"practical range not defined");
    break;
  case ERRORNOSTATMATCH : 
    strcpy(m,"no matching assumption found for the domains");
    break;
  case ERRORANISO:
    strcpy(m,"anisotropic call not allowed"); break; 
  case ERRORUNKNOWNMETHOD:
    strcpy(m,"Unknown method in this context or unallowed mixture of methods"); 
    break;
  case ERRORWRONGDIM:
    strcpy(m,"wrong dimension"); break;
   case ERRORUNKOWNSXPTYPE:
    strcpy(m, "parameter value of unknown SXP type");
    break;
  case ERROROUTOFMETHODLIST:
    sprintf(m, 
	    "Running out of list of methods. %s%s",
	    GLOBAL.general.skipchecks
	    ? "Did you try an invalid parameter combination?"
	    : "Are the RFoptions() too restrictive?",
	    PL <= 2
	    ? "\n You get (more) internal information if you set RFoptions(cPrintlevel=3) before running your code."
	    : ""
	    );
    break;
  case ERRORREGISTER: 
    strcpy(m, "register number out of range");
    break;
  case ERRORINCORRECTSTATISO: 
    strcpy(m, "distances only allow domain and isotropic frame  or  Gaussian fields need 'covariance' and 'anisotropic' as input");
    strcpy(m, "wrong domown or wrong isoown");
    strcpy(m, "wrong domain/isotropy or wrong/missing parameters");
    break;
  case ERRORM: 
    strcpy(m, ERRORSTRING);
    break;
  case ERRORWRONG:
    sprintf(m, "%s", ERRORSTRING_WRONG);
    break;
    
    // extremes:
  case ERRORSUBMETHODFAILED:
    sprintf(m, "no good submethods exist");
  case ERRORONLYISOTROPIC :
    strcpy(m, "only isotropic fields are allowed");
    break;
  case  ERRORSTATVARIO:
    strcpy(m, 
	   "negative definite function expected depending on 1 variable only");
    break;
   case ERRORNOVARIOGRAM:
    strcpy(m, "Variogram model not allowed in this context");
    break;
  case ERRORMARKOVPARAMETER :
    sprintf(m, "GMRF method not available for the given parameters (%s)",
	    ERRORSTRING_WRONG);
    break;
  case ERRORNORMALMIXTURE:
    strcpy(m, "only normal mixtures as first submodel allowed (Gneiting, 2002)");
    break;
  case ERRORMAXDIMMETH:
    strcpy(m, "maximal dimension of variables for the method exceeded");
    break;
  case ERRORPREVDOLLAR:
    strcpy(m, "method may not be initialised by preceding initS");
    break;
  case ERRORSPECTRAL: 
    strcpy(m, "submodel does not have spectral representation");
    break;    
  case ERRORTBMCOMBI: 
    strcpy(m, "the given combination of 'fulldim' and 'reduceddim' is not possible yet.");
    break;    

  case ERRORINVALIDMODEL : // gauss distribution, no method
    strcpy(m, "Invalid covariance model: did you wrongly use an auxiliary function to construct the model?");
    break;    
  case ERRORODDMODEL : // gauss distribution, no method
    strcpy(m, "Odd covariance model: the use of auxiliary functions and/or your choice of the parameters lead to a covariance model for which no simulation methods exist.");
    break;    
  case ERRORANISO_T :
    strcpy(m, "'anisoT' may not be given at the same time with 'Aniso' or 'proj'");
    break;
  case ERRORDIAMETERNOTGIVEN:
    strcpy(m, "Diameter must always be given");
    break;
  case ERRORPREFNONE:
    strcpy(m, "the simulation method does not allow for the given model.");
    break;
    
    //    case : strcpy(m,"");break;
    //
    // Poisson:
  case ERRORUNKNOWNMAXTYPE :
    strcpy(m, "unknown type of max-stable process");
    break;
 
  case ERRORATOMP :
    strcpy(m, "p must be given everywhere or nowhere");
    break;
   
  case ERRORKRIGETOL :
    strcpy(m,"sigma must be at most KRIGE_TOLERANCE");
    break;


  case MSGLOCAL_OK :
    strcpy(m,"fine");
    break;
  case MSGLOCAL_JUSTTRY :
    strcpy(m,
	   "unclear whether algorithm will work for specified parameters");
    break;
  case MSGLOCAL_NUMOK :
    strcpy(m,"fine. Algorithm should work for specified parameters");
    break;
  case MSGLOCAL_ENDOFLIST :
    strcpy(m,"end of list for variants of the algorithm");
    break;
  case MSGLOCAL_SIGNPHI :
    strcpy(m,"wrong sign of covariance function at break point");
    break;
  case MSGLOCAL_SIGNPHIFST :
    strcpy(m, "wrong sign of 1st derivative of the covariance function at the break point");
    break;
  case MSGLOCAL_SIGNPHISND :
    strcpy(m, "wrong sign of 2nd derivative of the covariance function at the break point");
    break;
  case MSGLOCAL_INITINTRINSIC :
    strcpy(m,"one of a2, b or a0+phi(0) has wrong sign");
    break;
  case ERRORUNSPECIFIED :
    strcpy(m,"(unspecified)");
    break;
   default : 
     PRINTF(" error=%d\n", err); 
     // crash();
     BUG;
  }
  
  if (strlen(m) > (unsigned int) len && len > 6) {    
    //  printf("%s %d %d\n", m, strlen(m), len);
    m[len-2] = m[len-3] = m[len-4] = '.';
    m[len-5] = ' ';
    m[len-1] ='\0';
    // printf("out %s %d %d\n", m, strlen(m), len);
   }
  if (PL >= PL_ERRORS) {
    PRINTF("error code %d [%s]\n", err, m);
  }
}

void errorMSG(int err, char* m) {
  errorMSG(err, m, 100000);
}

void ErrorStop(int err) {
  char m[1000];
  errorMSG(err, m);
  error(m);
}

int checkOK(cov_model VARIABLE_IS_NOT_USED *cov){
   return NOERROR;
}

int checkMissing(cov_model *cov){
  if (cov->calling == NULL) ERR("missing may not be called by the user");
  char S[100];
  cov_model *prev=cov->calling;
  sprintf(S, "'%s' does have not enough submodels", NICK(prev));
  ERR(S);
  return ERRORFAILED; // damit compiler keine Warnung bringt
}

int checkNotOK(cov_model VARIABLE_IS_NOT_USED *cov){
   return ERRORFAILED;
}

void ScaleOne(double *x, cov_model VARIABLE_IS_NOT_USED *cov, double *v){ 
  *v = *x <= 0.05 ? 1.0 : RF_NAN;
} 


sortsofparam paramtypeAny(int VARIABLE_IS_NOT_USED k, int VARIABLE_IS_NOT_USED row, int VARIABLE_IS_NOT_USED col) { return ANYPARAM; }

double *EinheitsMatrix(int dim) {
  // Einheitsmatrizen
  double *mem;
  if ((mem = (double*) CALLOC(dim * dim, sizeof(double))) != NULL) {
    int d;
    for (d=0; d<dim; d+=dim+1) mem[d] = 1.0;
  }
  return mem;
}


bool SimpleChecks() {
  ///////////////////////////////////////////////////////////////////////
  // SIMPLE CHECKS 
  ///////////////////////////////////////////////////////////////////////
  
  // must be the very last one!
  int oldcurrentNrCov = currentNrCov;
  assert(currentNrCov < MAXNRCOVFCTS); // otherwise we have already reached
  //                the maximum number
  //                of models; note that it is planned that further
  //                 functions might be added by users
  while (currentNrCov < MAXNRCOVFCTS) addFurtherCov(NULL, NULL);
  currentNrCov = oldcurrentNrCov;

 
  if (PL > PL_STRUCTURE) PRINTF("init: end of checking model definitions\n");
  
  cov_fct *C;
  cov_model *cov = NULL;
  char biWM[] = "biWM";
  int nr, i, err,
    biwm = getmodelnr(biWM);
  bool skipchecks = GLOBAL.general.skipchecks;
  if (PL >= PL_STRUCTURE) PRINTF("init: checking model definitions\n");

    
  for (nr=0; nr<currentNrCov; nr++) { 
    if (nr == MISSING_COV) continue;
     if (nr == biwm || nr == CONSTANT || nr == MIXEDEFFECT || nr == MPPPLUS)
	continue;
   
    C = CovList + nr; // nicht gatternr
    if (PL >= PL_DETAILS) 
      PRINTF("\n nr =%d (%d) %s  ",nr, currentNrCov, C->name); 
    // general check

    int correct;
    for (correct=false; correct<true; correct++) {
      
      GLOBAL.general.skipchecks = !correct;

      cov = (cov_model *) MALLOC(sizeof(cov_model));
      COV_NULL(cov);
      cov->ownloc = (location_type *)  MALLOC(sizeof(location_type));
      LOC_NULL(cov->ownloc);
      cov->nr = nr;
      cov->maxdim = 1;
      
      for (i=0; i<Nothing; i++) cov->pref[i] = PREF_NONE;
      
      if (correct) {
	for (i=0; i<C->kappas; i++) {
	  int ncol, nrow, tot;
	  C->kappasize(i, cov, &nrow, &ncol);
	  tot = nrow * ncol;
	  cov->p[i] = (double*) MALLOC(sizeof(double) * tot);
	  int j;
	  for (j=0; j<tot; j++) cov->p[i][j] =0.1;
	  cov->ncol[i] = ncol;
	  cov->nrow[i] = nrow;
	}
      } else {
	for (i=0; i<C->kappas; i++) {
	  int j, ncol=3, nrow=5, tot=ncol*nrow;
	  cov->p[i] = (double*) MALLOC(sizeof(double) * tot);
	  for (j=0; j<tot; j++) cov->p[i][j] =0.1;
	  cov->ncol[i] = ncol;
	  cov->nrow[i] = nrow;
	}
      }

      for (i=0; i<C->minsub; i++) {
	addModel(cov->sub + i, MISSING_COV);
	cov->sub[i]->calling = cov;
      }

      int oldPL = PL;      
      PL = 0;

      //  PMI(cov);

      err = C->check(cov); // no check on error; !!
      PL = oldPL;
      if (PL >= PL_DETAILS) {
	if (err == NOERROR) PRINTF("  OK");
	else {
	  errorMSG(err, MSG); 
	  if (correct) PRINTF("!!! %s%s", ERROR_LOC, MSG);
	  else PRINTF(" %s%s", ERROR_LOC, MSG);
	}
      }
      for (i=0; i<Nothing; i++) {
	if(!C->implemented[i] && cov->pref[i]>0) {
	  PRINTF("%d: %s  pos=%d, value %d > 0\n", nr, C->name, i,
		 cov->pref[i]);
	  assert(false);
	}
      }
      COV_DELETE(&cov);
    }
  }
    
  if (PL > PL_STRUCTURE) PRINTF("init: end of checking model definitions\n");
  
  GLOBAL.general.skipchecks = skipchecks;
  return true;
}


void InitModelList() {
  assert(currentNrCov=-1);  // otherwise something went wrong with the call
  assert(MODEL_MAX == 21); // otherwise change REGNAMES
  assert(ROLE_LAST == 12); // otherwise change ROLENAMES
  assert(OtherType == 12); // otherwise change TYPENAMES, 
  //                                           CAT_TYPENAMES[OtherType + 1]
  //                                           TypeConsistency in getNset
  //                                           RF_GLOBAL on ../R/
  assert(MAXMPPDIM <= MAXSIMUDIM); // ZERO
  // assert(CUTOFF_THEOR == 4);/* do not change this value as used in RFmethods.Rd */

  int i;

  for (i=0; i<MAXSIMUDIM; i++) ZERO[i] = 0.0;

  for (i=0; i<MAXPARAM; i++) sprintf(STANDARDPARAM[i], "k%d", i+1);
  for (i=0; i<MAXSUB; i++) sprintf(STANDARDSUB[i], "u%d", i+1);

  //  assert (KEY == NULL);
  //  KEY = (cov_model*) MALLOC(sizeof(cov_model **) * (MODEL_MAX+1));

   // init models
  for (i=0; i<=MODEL_MAX; i++) {
    KEY[i] = NULL; 
    MEM_NAS[i] = -1;
  }

  if (CovList!=NULL) {
    PRINTF("List of covariance functions looks already initiated.\n"); 
    return;
  }
  CovList = (cov_fct*) MALLOC(sizeof(cov_fct) * (MAXNRCOVFCTS+1));
  // + 1 is necessary because of COVINFO_NULL that uses the last + 
  currentNrCov = 0;

  // *******************
  // **** trend-models ****
  // *******************
 
  MIXEDEFFECT = 
    IncludeModel("mixed", TrendType, 0, 1, 6, kappamixed,
		 XONLY, PREVMODELI,  // todo !!
		 checkmixed, rangemixed, PREF_NOTHING,
		 false, SUBMODEL_DEP, SUBMODEL_DEP,
		 SUBMODEL_DEP, NOT_MONOTONE);
  make_internal();
  // if element is negative and SpaceEffect then only the covariance of the subsequent model is returned (as if 'X' were not given)
  kappanames("element", INTSXP, "X", LISTOF+REALSXP, "beta", REALSXP,
	     "coord", REALSXP, "dist", REALSXP, "dim", INTSXP);
  subnames("cov");
  addCov(mixed, NULL, NULL);
  addCov(mixed_nonstat);
  RandomShape(0, initmixed, domixed);
  addReturns(NULL, NULL, covmatrix_mixed, iscovmatrix_mixed,
	     NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
  //  MLEMIXEDEFFECT = addFurtherCov(MLEmixed, ErrCov); 
  // addCov(MLEmixed_nonstat);


  pref_type ptrend = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10, 5};
        //           CE CO CI TBM Sp di sq Ma av n mpp Hy spf any
  TREND = IncludeModel("trend", TrendType,  0, 0, 6, kappatrend, 
		       XONLY, PREVMODELI,
		       checktrend, 
		       rangetrend,
		       ptrend,
		       false, PARAM_DEP, INFDIM, false, NOT_MONOTONE);
  kappanames("mean", REALSXP, "plane", REALSXP, "polydeg", INTSXP, "polycoeff",
	     REALSXP, "arbitraryfct", CLOSXP, "fctcoeff", REALSXP); 
  addCov(trend, NULL, NULL);
  addCov(trend_nonstat);
  RandomShape(0, structOK, init_trend, do_trend, false, false, true);



  // *******************
  // **** RO-models ****
  // *******************

  pref_type pGatter=  {5, 0, 0,  5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5};
  //                  CE CO CI TBM Sp di sq Ma av n mpp Hy spf any
  FIRST_GATTER = ISO2ISO = // 2
      IncludeModel("#",  OtherType, 1, 1, 0, NULL, PREVMODELD, PREVMODELI,
		   checkNotOK, NULL, pGatter, true, SUBMODEL_DEP,
		   SUBMODEL_DEP, SUBMODEL_DEP, SUBMODEL_DEP);
  addCov(iso2iso, D_2, DD_2, inverse2, nonstatinverse2);
  addlogCov(logiso2iso);
  RandomShape(INFTY, struct2, init2, do2, dorandom2, true, true, false); 
  addReturns(NULL, NULL, covmatrixS,  iscovmatrixS, 
  	     NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
    
  SP2SP = addFurtherCov(spiso2spiso, D_2, DD_2); // 3
  addlogCov(logspiso2spiso);
 
  SP2ISO = addFurtherCov(spacetime2iso, D_2, DD_2); // 4
  addlogCov(logspacetime2iso);

  S2ISO = addFurtherCov(Stat2iso, ErrCov); // 5
  addCov(Nonstat2iso);// 
  addlogCov(logStat2iso, logNonstat2iso);

  assert(CovList[S2ISO].Init != NULL);
  assert(S2ISO == 5);

  S2SP = addFurtherCov(Stat2spacetime, ErrCov);// 6
  addCov(Nonstat2spacetime);// 
  addlogCov(logStat2spacetime, logNonstat2spacetime);
  
  S2S = addFurtherCov(Stat2Stat, ErrCov);// 7
  addCov(Nonstat2Stat);// 
  addlogCov(logStat2Stat, logNonstat2Stat);
  // printf("# %ld %ld %ld\n", Stat2Stat, CovList[currentNrCov-1].cov, Stat2iso);
  LAST_GATTER = SId = addFurtherCov(Stat2Stat, ErrCov);// 8
  addCov(Nonstat2Nonstat);// 
  addlogCov(logStat2Stat, logNonstat2Nonstat);
  assert(SId == 8);

  //  addFurtherCov(iso2iso_MLE, D_2, DD_2); // 12
  //  addFurtherCov(spiso2spiso_MLE, D_2, DD_2); // 13
  //  addFurtherCov(spacetime2iso_MLE, D_2, DD_2); // 14
  //  addFurtherCov(Stat2iso_MLE, ErrCov); // 15
  //  addCov(Nonstat2iso_MLE);// 
  //  addFurtherCov(Stat2spacetime_MLE, ErrCov);// 16
  //  int mleGatter = addFurtherCov(Stat2Stat_MLE, ErrCov);// 17
  //  addCov(Nonstat2Stat_MLE);// 
  //  if (mleGatter != S2S + (S2S - ISO2ISO + 1)) {
  //    error("mleGatter has the wrong number");
  //  }

  
  MISSING_COV =
    IncludePrim("missing", OtherType, 0, XONLY, SYMMETRIC,
		checkMissing,  NULL, INFDIM, true, NOT_MONOTONE);
  make_internal(); 


  // *******************
  // **** definite functions  ****
  // *******************

  SELECT =  // to do: replace by parameter in '+', selecting the 'type' or
    // 'models'
    IncludeModel("select", TcfType, 1, MAXSUB, 1, NULL,
		 PREVMODELD, PREVMODELI,
		 checkselect, rangeconstant, PREF_ALL,
		 true, PARAM_DEP, INFDIM, SUBMODEL_DEP, NOT_MONOTONE);
  kappanames("subnr", INTSXP);
  addCov(select, NULL, NULL); 
  addReturns(NULL, NULL, covmatrix_select, iscovmatrix_select, 
	     NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
   
  pref_type pplus =  {5, 0, 0,  5, 0, 5, 5, 0, 0, 0, 0, 0, 5, 5};
  //                  CE CO CI TBM Sp di sq Ma av n mpp Hy spf any
  PLUS = 
    IncludeModel("+", UndefinedType, 1, MAXSUB, 0, NULL, PREVMODELD, PREVMODELI,
		 checkplus, NULL, pplus, 
		 false, SUBMODEL_DEP, SUBMODEL_DEP, SUBMODEL_DEP, SUBMODEL_DEP);
  // Achtung in covmatrix_plus wird SELECT_SUBNR verwendet!
  nickname("plus");
  addCov(plusStat, Dplus, DDplus, NULL, NULL);
  addCov(plusNonStat);
  addTBM(NULL, spectralplus);
  RandomShape(0, structplus, initplus, doplus, false, false, true);
  addReturns(NULL, NULL, covmatrix_plus, iscovmatrix_plus, 
	     NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
  addTypeFct(Typeplus);



 pref_type pmal =  {5, 0, 0,  5, 0, 5, 5, 0, 0, 0, 0, 0, 0, 5};
  //                 CE CO CI TBM Sp di sq Ma av n mpp Hy spf any
  MULT = IncludeModel("*", TcfType,  1, MAXSUB, 0, NULL, PREVMODELD, PREVMODELI,
		      checkmal, NULL, pmal, false, SUBMODEL_DEP,
		      SUBMODEL_DEP, SUBMODEL_DEP, SUBMODEL_DEP);
  nickname("mult");
  addCov(malStat, Dmal, NULL);
  addCov(malNonStat);
  addlogCov(logmalStat, logmalNonStat);
  //  RandomShape(structplusmal, initmal, domal, NULL);
  addTypeFct(Typemal);

  pref_type pS=  {5, 0, 0,  5, 5, 5, 5, 0, 0, 5, 0, 0, 1, 5};
  //        CE CO CI TBM Sp di sq Ma av n mpp Hy spf any
  DOLLAR = IncludeModel("$",  UndefinedType, // to do: tcftype durch einen allgemeinen Type ersetzen, da auch Trend dem "$" folgen kann. Z.Z. nicht moeglich.

			1, 1, 5, kappaS, // kappadollar,
			PREVMODELD, PREVMODELI, checkS, rangeS, pS,
			false, SUBMODEL_DEP, SUBMODEL_DEP, SUBMODEL_DEP,
			SUBMODEL_DEP);
  // do not change Order!!
  nickname("S");
  kappanames("var", REALSXP, "scale", REALSXP, "anisoT", REALSXP,
	     "Aniso", REALSXP, "proj", INTSXP);
  addkappa(3, "Aniso", REALSXP, ShapeType);
  subnames("phi");
  addTypeFct(TypeS);
  addCov(Siso, DS, DDS, inverseS); // unterscheidung nur wegen der 
  //  geschwindigkeit, also Siso ist ein sehr haeufiger Spezialfall von Sstat
  addCov(Snonstat);
  addlogCov(logSiso);
  // addLocal(coinitS, ieinitS);  
  addTBM(tbm2S, NULL, spectralS);
  nablahess(nablaS, hessS);
  RandomShape(INFTY, structS, initS, doS, true, true, true);
  addReturns(NULL, NULL, covmatrixS, iscovmatrixS, 
  	     NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
  Taylor(RF_NAN, RF_NAN, RF_NAN, RF_NAN);
  TailTaylor(RF_NAN, RF_NAN, RF_NAN, RF_NAN);
   
  LASTDOLLAR = addFurtherCov(Sstat, ErrCov); // 3
  addCov(Snonstat);
  addlogCov(logSstat, logSnonstat);
  RandomShape(INFTY, structS, initS, doS, true, true, false);
  


  // at the origin: prespecified distribution
  // old RandomFields ave1,ave2
  IncludeModel("ave",  PosDefType, 1, 1, 3, kappa_ave, XONLY, SYMMETRIC,
	       checkave, rangeave, PREF_ALL, 
	       false, SCALAR, AveMaxDim, false, NOT_MONOTONE);
  kappanames("A", REALSXP, "z", REALSXP, "spacetime", INTSXP);
  addCov(ave, NULL, NULL);
  RandomShape(structAve, true);



  pref_type
    pbessel = {2, 0, 0,  0, 5, 4, 5, 0, 5, 0, 5, 0, 0, 5};
  //            CE CO CI TBM Sp di sq Ma av n mpp Hy spf any
  IncludePrim("bessel",  PosDefType, 1, XONLY, ISOTROPIC, 
	      checkBessel, rangeBessel,
	      pbessel, SCALAR, INFDIM, false, NOT_MONOTONE);
  kappanames("nu", REALSXP);
  addCov(Bessel, NULL, NULL);
  addTBM(initBessel, spectralBessel);	       
  add_paramtype(uncritical_paramtype);
  
  IncludeModel("bigneiting", PosDefType, 0, 0, 8, kappa_biGneiting, XONLY,
	       ISOTROPIC, checkbiGneiting, rangebiGneiting, PREF_ALL, 
	       false, 2, PARAM_DEP, true, NOT_MONOTONE);
  addCov(biGneiting, DbiGneiting, DDbiGneiting, NULL, NULL);
  kappanames("kappa", INTSXP,
	     "mu", REALSXP,
	     "s", REALSXP, "sred12", REALSXP,
	     "gamma", REALSXP,
	     "cdiag", REALSXP, "rhored", REALSXP, "c", REALSXP);
  add_paramtype(paramtype_biGneiting);
 

  pref_type
    pbernoulli = {5, 0, 0,  0, 0, 5, 5, 0, 0, 0, 0, 0, 0, 5};
  //             CE CO CI TBM Sp di sq Ma av n mpp Hy spf any
  IncludeModel("bernoulli", TcfType, 1, 1, 1, NULL, PREVMODELD, PREVMODELI,
	       checkbinary, rangebinary, pbernoulli,
	       false, SCALAR, SUBMODEL_DEP, SUBMODEL_DEP, SUBMODEL_DEP);
  kappanames("threshold", REALSXP);
  addCov(binary, NULL, NULL);
  

  IncludeModel("biWM",  PosDefType, 0, 0, 8, kappa_biWM, XONLY, ISOTROPIC,
	       checkbiWM2, rangebiWM2, PREF_ALL,
	       false, 2, INFDIM, false, NOT_MONOTONE);
  nickname("biwm");
  addCov(biWM2, biWM2D, NULL);
  kappanames("nudiag", REALSXP, "nured12", REALSXP, 
	     "nu", REALSXP, // or lower triangle
	     "s", REALSXP,  // lower triangle definition
	     "cdiag", REALSXP, "rhored", REALSXP,
	     "c", REALSXP,  // or lower triangle
	     "notinvnu", INTSXP);
  add_paramtype(paramtype_biWM);

  pref_type
    pbrownresnick = {5, 0, 0,  5, 0, 5, 5, 0, 0, 0, 0, 0, 0, 5};
  //                CE CO CI TBM Sp di sq Ma av n mpp Hy spf any
  IncludeModel("brownresnick", TcfType, 1, 1, 0, NULL, XONLY, PREVMODELI,
	       checkbrownresnick, NULL , pbrownresnick, false,
	       SCALAR, SUBMODEL_DEP, false, SUBMODEL_DEP);
  addCov(brownresnick, Dbrownresnick, DDbrownresnick, NULL, NULL);
  RandomShape(0, struct_brownresnick, init_brownresnick, do_brownresnick);
  //  Taylor(0, 0, 0, 0, 0, 0, 0, 0);

  
  IncludeModel("br2bg",  PosDefType, 1, 1, 0, XONLY, PREVMODELI, 
	       check_BR2BG, NULL, PREF_ALL, SUBMODEL_DEP, false, SUBMODEL_DEP);
  addCov(BR2BG, NULL, NULL);
  add_paramtype(paramtypeAny);

  IncludeModel("br2eg", PosDefType, 1, 1, 0,  XONLY, PREVMODELI, 
	       check_BR2EG, NULL, PREF_ALL, SUBMODEL_DEP, false, SUBMODEL_DEP);
  addCov(BR2EG, NULL, NULL);
  add_paramtype(paramtypeAny);

 
  pref_type pcauchy=  {2, 0, 0,  3, 0, 4, 0, 0, 0, 0, 0, 0, 0, 5};
  //        CE CO CI TBM Sp di sq Ma av n mpp Hy spf any
  IncludePrim("cauchy", TcfType, 1, XONLY, ISOTROPIC, 
	      checkCauchy, rangeCauchy, pcauchy, 
	      SCALAR, INFDIM, false, NORMAL_MIXTURE);
  kappanames("gamma", REALSXP);
  addCov(Cauchy, DCauchy, DDCauchy, InverseCauchy);
  addlogCov(logCauchy);
  addTBM(TBM2Cauchy);
  addLocal(coinitCauchy, NULL);
  addGaussMixture(DrawMixCauchy, LogMixDensCauchy);
	       
  //  pref_type pctbm={2, 0, 0,  5, 0, 5, 5, 0, 0, 0, 0, 0, 0, 5};
  //     //           CE CO CI TBM Sp di sq Ma av n mpp Hy spf any
  // IncludePrim("cauchytbm", PosDefType,  3, XONLY, ISOTROPIC, checkOK,
  //	      rangeCauchytbm, pctbm, SCALAR, INFDIM, false);
  // kappanames("alpha", REALSXP, "beta", REALSXP, "gamma", REALSXP);
  // addCov(Cauchytbm, DCauchytbm, InverseCauchy); // scale not correct, but
  // should be an approximation that is good enough
	      
  IncludePrim("circular",  TcfType, 0, XONLY, ISOTROPIC,
	      checkOK, NULL, 2, false, GNEITING_MON);
  addCov(circular, Dcircular, ScaleOne);
  RandomShape(structcircular);


  // IncludePrim("cone", PosDefType,  3, XONLY, ISOTROPIC, checkcone, rangecone);
  //  kappanames("r", REALSXP, "socle", REALSXP, "height", REALSXP);
  // RandomShape(init_cone, mppget_cone, sd_standard, MPP_POISS);

  IncludePrim("CDeWijsian",  NegDefType, 2, NULL, XONLY, ISOTROPIC, 
	      checkdewijsian,  rangeDeWijsian, PREF_NOTHING, 
	      1, INFDIM, false, MONOTONE); 
  nickname("cdewijs");
  make_internal();
  kappanames("alpha", REALSXP, "range", REALSXP);
  addCov(DeWijsian, NULL, NULL, InverseDeWijsian); 

  CONSTANT = 
    IncludeModel("constant", TcfType, 0, 0, 3, NULL, XONLY, ISOTROPIC,
		 //  PREVMODELD, PREVMODELI, 
		 //wegen Variogramm berechnung in stat. Fall
		 checkconstant, rangeconstant, PREF_ALL,
		 false, PARAM_DEP, INFDIM, false, COMPLETELY_MON);
  kappanames("element", INTSXP, "M", LISTOF+REALSXP, "vdim", INTSXP);  
  addCov(constant, NULL, NULL); 
  addCov(constant_nonstat);
  add_paramtype(uncritical_paramtype);
  addReturns(NULL, NULL, covmatrix_constant, iscovmatrix_constant, 
	     NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);

  pref_type pcox={2, 0, 0,  5, 0, 5, 5, 0, 0, 0, 0, 0, 0, 5};
  //           CE CO CI TBM Sp di sq Ma av n mpp Hy spf any
  IncludeModel("coxisham",  PosDefType, 1, 1, 3, kappa_cox, 
	       XONLY, ZEROSPACEISO, 
	       checkcox, rangecox, pcox,
	       false, SCALAR, CoxMaxDim, false, NOT_MONOTONE);
  kappanames("mu", REALSXP, "D", REALSXP, "beta", REALSXP);  
  addCov(cox, NULL, NULL);
  addTBM(initcox, spectralcox);
  nablahess(coxnabla, coxhess);

  IncludePrim("cubic",  TcfType, 0, XONLY, ISOTROPIC, 
	      checkOK, NULL, 3, false, MONOTONE);
  addCov(cubic, Dcubic, ScaleOne);
	       
  pref_type pcurl= {2, 0, 0, 0, 0, 5, 5, 0, 0, 0, 0, 0, 0, 5};
  //           CE CO CI TBM Sp di sq Ma av n mpp Hy spf any
  IncludeModel("curlfree",  PosDefType, 1, 1, 0, NULL, XONLY, SYMMETRIC,
	       checkdivcurl, NULL, pcurl,
	       false, PARAM_DEP, SUBMODEL_DEP, SUBMODEL_DEP, NOT_MONOTONE);
  addCov(curl, NULL, NULL);
 
  pref_type plocal={5, 0, 0, 0, 0, 5, 5, 0, 0, 0, 0, 0, 0, 5};
  //            CE CO CI TBM Sp di sq Ma av n mpp Hy spf any
  CUTOFF =  
    IncludeModel("cutoff",  PosDefType, 1, 1,2, NULL, XONLY, ISOTROPIC,
		 check_co, range_co, plocal,
		 false, SCALAR, MAXCEDIM,  true, MONOTONE);
  kappanames("diameter", REALSXP, "a", REALSXP);  
  addCov(co, NULL, NULL);
  addCallLocal(alternativeparam_co);
 
  //  warning("siehe Primitive.cc/biWM: cutoff funktioniert nicht bei MLE, vereinheitlichung mit natsc und verbesserung von biWM notwendig");
 

  IncludeModel("dagum",  PosDefType, 0, 0, 2, NULL, XONLY, ISOTROPIC,
	       checkOK, rangedagum, PREF_ALL, false, 1, INFDIM, false, MONOTONE);
  kappanames("beta", REALSXP, "gamma", REALSXP);
  addCov(dagum, Ddagum, Inversedagum);



  IncludePrim("dampedcosine",  PosDefType, 1, XONLY, ISOTROPIC,
	      checkdampedcosine, rangedampedcosine, PARAM_DEP,
	      false, NOT_MONOTONE);
  nickname("dampedcos");
  kappanames("lambda", REALSXP);
  addCov(dampedcosine, Ddampedcosine, Inversedampedcosine);
  // addlogCov(logdampedcosine);

  IncludePrim("DeWijsian", NegDefType,  1, XONLY, ISOTROPIC,
	      checkOK, rangedewijsian, INFDIM, false, MONOTONE);
  nickname("dewijsian");
  kappanames("alpha", REALSXP);
  addCov(dewijsian, Ddewijsian, NULL, Inversedewijsian); 

 

  pref_type pdiv= {2, 0, 0, 0, 0, 5, 5, 0, 0, 0, 0, 0, 0, 5};
  //           CE CO CI TBM Sp di sq Ma av n mpp Hy spf any
  IncludeModel("divfree", PosDefType, 1, 1, 0, NULL, XONLY, SYMMETRIC, 
	       checkdivcurl, NULL, pdiv, 
	       false, PARAM_DEP, SUBMODEL_DEP, SUBMODEL_DEP, NOT_MONOTONE);
  addCov(div, NULL, NULL);


 
  // epsC has been for internal reasons only ; essentially
  // the gencauchy model, except that 1 in the denominator 
  // is replaced by epsilon



  pref_type pepsC = {2, 5, 5, 5, 0, 5, 5, 0, 0, 0, 0, 0, 0, 5};
  //           CE CO CI TBM Sp di sq Ma av n mpp Hy spf any
  IncludeModel("epsC",  PosDefType, 0, 0, 3, NULL, XONLY, ISOTROPIC,
	       checkepsC, rangeepsC, pepsC, 
	       false, SCALAR, INFDIM, false, NORMAL_MIXTURE);
  nickname("epscauchy");
  kappanames("alpha", REALSXP, "beta", REALSXP, "eps", REALSXP);
  addCov(epsC, DepsC, DDepsC, NULL, NULL);
  addlogCov(logepsC);

  IncludePrim("exponential", TcfType, 0, XONLY, ISOTROPIC,
	      checkexponential, NULL, INFDIM, false, COMPLETELY_MON);
  nickname("exp");
  addCov(0, exponential, Dexponential, DDexponential, Inverseexponential, NULL);
  addlogCov(logexponential);
  addLocal(coinitExp, ieinitExp);
  //addMarkov(&EXPONENTIAL);
  addHyper(hyperexponential);
  // numerisches Ergebnis passt nicht !
  addGaussMixture(DrawMixExp, LogMixDensExp);
  addTBM(TBM2exponential, NULL, spectralexponential);
  RandomShape(0, initexponential, do_exp);
  Taylor(-1, 1.0, 0.5, 2.0);
  TailTaylor(1, 0, 1, 1);
  
  // operator, replacing any covariance fct C by exp(C) (componentwise)
  IncludeModel("Exp", 
	       PosDefType, 1, 1, 2, PREVMODELD, PREVMODELI, checkExp,
	       rangeExp, PREF_ALL, SUBMODEL_DEP, false, NOT_MONOTONE);
  nickname("exponential");
  kappanames("n", INTSXP, "standardised", INTSXP);
  addCov(Exp, DExp, DDExp, NULL, NULL);
  add_paramtype(paramtypeAny);
  
  pref_type
    pextrgauss = {5, 0, 0,  0, 0, 5, 5, 0, 0, 0, 0, 0, 0, 5};
  //             CE CO CI TBM Sp di sq Ma av n mpp Hy spf any
  IncludeModel("extremalgauss", TcfType, 1, 1, 0, NULL,
	       XONLY, PREVMODELI,
	       check_extrgauss, NULL, pextrgauss, false,
	       SCALAR, SUBMODEL_DEP, SUBMODEL_DEP, NOT_MONOTONE);
  nickname("schlather");
  addCov(extrgauss, NULL, NULL);
  

  IncludePrim("FD", PosDefType,  1, XONLY, ISOTROPIC,
	      checkOK, rangeFD, 1, false, NOT_MONOTONE); 
  nickname("fractdiff");
  kappanames("a", REALSXP);
  addCov(FD, NULL, NULL);


  BROWNIAN = 
    IncludePrim("fractalB", NegDefType, 1, XONLY, ISOTROPIC, 
		checkfractalBrownian, rangefractalBrownian, INFDIM, false,
		BERNSTEIN); // todo BERNSTEIN
  nickname("fbm");
  kappanames("alpha", REALSXP);
  addCov(fractalBrownian, DfractalBrownian, DDfractalBrownian, 
	 D3fractalBrownian, D4fractalBrownian, 
	 InversefractalBrownian);
  addlogCov(logfractalBrownian);
  addLocal(NULL, ieinitBrownian);
  add_paramtype(uncritical_paramtype);
  Taylor(-1, RF_NAN, 0, 0);
  TailTaylor(-1, RF_NAN, 0, 0);
  
  IncludePrim("fractgauss", PosDefType, 1, XONLY, ISOTROPIC,
	      checkOK, rangefractGauss, 1, false, NOT_MONOTONE);
  kappanames("alpha", REALSXP);
  addCov(fractGauss, NULL, NULL);

  pref_type pgauss= {2, 0, 0, 5, 5, 5, 5, 5, 0, 0, 5, 0, 0, 5};
  //                CE CO CI TBM Sp di sq Ma av n mpp Hy spf any
  GAUSS = 
    IncludePrim("gauss",  PosDefType, 0, XONLY, ISOTROPIC,
		checkOK, NULL, pgauss,
		SCALAR, INFDIM, false, NORMAL_MIXTURE);
  addCov(Gauss, DGauss, DDGauss, D3Gauss, D4Gauss, InverseGauss);
  addlogCov(logGauss);
  //addMarkov(&GAUSS);
  addTBM(NULL, spectralGauss);
  RandomShape(INFTY, struct_Gauss, initGauss, do_Gauss, false, true, false);
  addGaussMixture(DrawMixGauss, LogMixDensGauss);
  Taylor(-1.0, 2.0);
  TailTaylor(1, 0, 1.0, 2.0);

  
  IncludePrim("genB", NegDefType, 2, XONLY, ISOTROPIC, 
	      checkOK, rangegenBrownian, INFDIM, false, MONOTONE);
  nickname("genfbm");
  kappanames("alpha", REALSXP, "delta", REALSXP);
  addCov(genBrownian, NULL, NULL, InversegenBrownian); 
  addlogCov(loggenBrownian);
  
  pref_type pgenc = {2, 5, 5, 5, 0, 5, 0, 0, 0, 0, 0, 0, 0, 5};
  //           CE CO CI TBM Sp di sq Ma av n mpp Hy spf any
  IncludePrim("gencauchy", UndefinedType, 2, XONLY, ISOTROPIC,
	      checkgeneralisedCauchy, rangegeneralisedCauchy, pgenc,
	      SCALAR, INFDIM, false, NORMAL_MIXTURE); // todo part is even
  // LAPLACE
  kappanames("alpha", REALSXP, "beta", REALSXP);
  addCov(generalisedCauchy, DgeneralisedCauchy, DDgeneralisedCauchy,
	 InversegeneralisedCauchy);
  addlogCov(loggeneralisedCauchy);
  addLocal(coinitgenCauchy, ieinitgenCauchy);
  addTypeFct(TypegeneralisedCauchy);

  IncludePrim("gengneiting",  PosDefType, 2, XONLY, ISOTROPIC, 
	      checkgenGneiting, rangegenGneiting, INFDIM-1, true,
	      MONOTONE); // GNEITING_MON ??
  // not INFDIM, also not normalscale mixture and alpha will be void
  kappanames("kappa", INTSXP, "mu", REALSXP);
  addCov(genGneiting, DgenGneiting, DDgenGneiting, ScaleOne);

  IncludePrim("gneiting", PosDefType, 0, XONLY, ISOTROPIC,
	      checkOK, NULL, 3, true, MONOTONE);  // GNEITING_MON ??
  addCov(Gneiting, DGneiting, ScaleOne);
 
  pref_type phyper= {2, 0, 0, 3, 0, 4, 5, 0, 5, 0, 5, 0, 0, 5};
  //           CE CO CI TBM Sp di sq Ma av n mpp Hy spf any
  IncludePrim("hyperbolic",  PosDefType, 3, XONLY, ISOTROPIC,
	      checkhyperbolic, rangehyperbolic, phyper,
	      SCALAR, INFDIM, false, NORMAL_MIXTURE);
  kappanames("nu", REALSXP, "lambda", REALSXP, "delta", REALSXP);
  addCov(hyperbolic, Dhyperbolic, NULL); // InversehyperbolicSq);
  addlogCov(loghyperbolic);

  IncludePrim("iacocesare",  PosDefType, 3, XONLY, SPACEISOTROPIC, 
	      checkOK, rangeIacoCesare, INFDIM, false, NOT_MONOTONE);
  nickname("iaco");
  kappanames("nu", REALSXP, "lambda", REALSXP, "delta", REALSXP);
  addCov(IacoCesare, NULL, NULL);
  
  IncludeModel("identity", UndefinedType, 1, 1, 1, NULL, PREVMODELD, PREVMODELI,
	       checkId, rangeId, PREF_ALL, 
	       false, SUBMODEL_DEP, SUBMODEL_DEP, SUBMODEL_DEP, SUBMODEL_DEP);
  nickname("id");
  kappanames("vdim", INTSXP);
  addCov(IdStat, DId, DDId, IdInverse);
  addCov(IdNonStat);
  addTBM(TBM2Id, initId, spectralId);
  addLocal(coinitId, ieinitId);
  addTypeFct(Typesetparam);

  IncludePrim("kolmogorov",  NegDefType, 0, XONLY, VECTORISOTROPIC,
	      checkKolmogorov, NULL, 3, 3, false, NOT_MONOTONE);
  addCov(Kolmogorov, NULL, NULL);

  pref_type plgd1= {2, 0, 0, 5, 0, 5, 5, 0, 0, 0, 0, 0, 0, 5};
  //           CE CO CI TBM Sp di sq Ma av n mpp Hy spf any
  IncludePrim("lgd1",  PosDefType, 2, NULL, XONLY, ISOTROPIC, 
	      checklgd1, rangelgd1, plgd1, 
	      SCALAR, PARAM_DEP, false, MONOTONE);
  nickname("lgd");
  kappanames("alpha", REALSXP, "beta", REALSXP);
  addCov(lgd1, Dlgd1, NULL); // Inverselgd1);



  // stimmt so nicht, siehe Gneiting 1998, on a alpha-sym multiv. char. fct:
  //  IncludeModel("lp", PosDefType,  1, 1, 1, XONLY, SYMMETRIC, 
  //	       checklp, rangelp,
  //	       (pref_type) {5, 0, 0, 0, 0, 5, 0, 0, 0, 0, 0, 0, 0, 5}
  //	       //          CE CO CI TBM Sp di sq Ma av n mpp Hy spf any
  //             );
  //  kappanames("p", REALSXP);
  //  addCov(lp, NULL, NULL);
 
  IncludeModel("mastein",  PosDefType, 1, 1, 2, XONLY, SPACEISOTROPIC, 
	       check_MaStein, range_MaStein, PREF_ALL, 
	       SUBMODEL_DEP, false, NOT_MONOTONE);
  kappanames("nu", REALSXP, "delta", REALSXP);
  addCov(MaStein, NULL, NULL);

    
  IncludeModel("ma1", PosDefType,  1, 1, 2, XONLY, SYMMETRIC,
	       checkma1, rangema1, PREF_ALL, 
	       SUBMODEL_DEP, false, NOT_MONOTONE);
  nickname("ma");
  kappanames("alpha", REALSXP, "theta", REALSXP);
  addCov(ma1, NULL, NULL);


  IncludeModel("ma2",  PosDefType, 1, 1, 0, XONLY, SYMMETRIC,
	       checkma2, NULL, PREF_ALL, SUBMODEL_DEP, false, NOT_MONOTONE);
  nickname("intexp");
  addCov(ma2, NULL, NULL);

  IncludeModel("M",  PosDefType, 1, 1, 1, kappaM, PREVMODELD, PREVMODELI,
	       checkM, rangeM, PREF_ALL,
	       false, PARAM_DEP, SUBMODEL_DEP, SUBMODEL_DEP, NOT_MONOTONE);
  nickname("matrix");
  kappanames("M", REALSXP);
  addCov(Mstat, NULL, NULL);
  addCov(Mnonstat);
  add_paramtype(paramtype_M);

	  
  IncludeModel("matern", UndefinedType, 0, 0, 2, XONLY, ISOTROPIC, 
	       checkMatern, rangeWM, PREF_ALL, INFDIM, false, NORMAL_MIXTURE);
  kappanames("nu", REALSXP, "notinvnu", INTSXP);
  addCov(Matern, DMatern, DDMatern, D3Matern, D4Matern, InverseMatern);
  addlogCov(logMatern);
  addTBM(initMatern, spectralMatern);
  addLocal(coinitWM, ieinitWM);
  addTypeFct(TypeWM);

  //  addGaussMixture(DrawMixWM, LogMixDensWM);


  IncludeModel("mqam", PosDefType,
  	       2, 10, 1, kappamqam, XONLY, SYMMETRIC,
  	       checkmqam, rangemqam, PREF_ALL, 
	       false, PARAM_DEP, SUBMODEL_DEP, false, NOT_MONOTONE);
  kappanames("theta", REALSXP);
  subnames("phi");
  addCov(mqam, NULL, NULL);
  add_paramtype(paramtype_qam);


  NATSC = 
    IncludeModel("natsc", TcfType,  1, 1, 0, NULL, XONLY, ISOTROPIC,
		 checknatsc, NULL, PREF_ALL,
		 false, 1, SUBMODEL_DEP, SUBMODEL_DEP, SUBMODEL_DEP);
  // nie einen Parameter !
  addCov(natsc, Dnatsc, DDnatsc, Inversenatsc);
  addLocal(coinitnatsc, ieinitnatsc);
  addTBM(tbm2natsc, initnatsc, spectralnatsc);

  IncludeModel("nonstWM", PosDefType, 0, 0, 1, kappaNonStWM, KERNEL, SYMMETRIC,
  	       checkNonStWM, rangeNonStWM, PREF_ALL, false, SCALAR,
	       INFDIM, false, NOT_MONOTONE);
  nickname("nonstwm");
  addCov(NonStWMQ); // anders herum gibt es fehler in addCov(aux_covfct auxcf),
  addCov(NonStWM);  // da auxiliary ia.. nicht mit cov hand in hand gehen kann
  addkappa(0, "nu", REALSXP, ShapeType);
  // subnames("nu"); // i.e. nu can be a constant or a submodel !!!
  //                 see GetSubNames in userinterface, how it is programmed
  //  addGaussMixture(DrawMixNonStWM, LogMixDensNonStWM);
  add_paramtype(paramtype_nonstWM);


  IncludeModel("nsst",  PosDefType, 2, 2, 1, XONLY, SPACEISOTROPIC,
	       checknsst, rangensst, PREF_ALL,
	       SUBMODEL_DEP, false, NOT_MONOTONE);
  kappanames("delta", REALSXP);
  subnames("phi", "psi");
  add_paramtype(paramtype_nsst);
  addCov(nsst, Dnsst, NULL);
  addTBM(TBM2nsst);
 
  //  IncludePrim("nsst2", 7, checknsst2, SPACEISOTROPIC, 
  //		   rangensst2);
  //  addCov(nsst2, Dnsst2, NULL);
  //  addTBM(NULL, NULL /* TBM3nsst2 */);

  pref_type pnugget= { 4, 0, 0, 0, 0, 4, 4, 0, 0, 5, 0, 0, 0, 5};
  //                  CE CO CI TBM Sp di sq Ma av n mpp Hy spf any
  NUGGET  = 
    IncludeModel("nugget",  TcfType, 0, 0, 2, NULL, XONLY, ISOTROPIC,
		 check_nugget, range_nugget, pnugget, 
		 false, PREVMODEL_DEP, INFDIM, true, MONOTONE);
  kappanames("tol", REALSXP, "vdim", INTSXP);
  add_paramtype(ignoreall_paramtype);
  addCov(nugget, NULL, Inversenugget);
  addReturns(NULL, NULL, covmatrix_nugget, iscovmatrix_nugget, 
	     NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);


  IncludeModel("parsWM", PosDefType, 0, 0, 1, kappa_parsWM, 
	       XONLY, ISOTROPIC,
	       checkparsWM, rangeparsWM, PREF_ALL,
	       false, PARAM_DEP, INFDIM, false, NOT_MONOTONE);
  nickname("parswm");
  addCov(parsWM, parsWMD, NULL);
  kappanames("nudiag", REALSXP);
  add_paramtype(paramtype_parsWM);

  IncludePrim("penta", PosDefType, 0, XONLY, ISOTROPIC,
	      checkOK, NULL, 3, true, MONOTONE);
  addCov(penta, Dpenta, ScaleOne);

  IncludePrim("power", UndefinedType,  1, XONLY, ISOTROPIC, 
	      checkpower, rangepower, INFDIM-1, true, MONOTONE);
  nickname("askey");
  kappanames("alpha", REALSXP);
  addCov(power, Dpower, ScaleOne);	
 addTypeFct(Typepower);
 
  IncludeModel("Pow", PosDefType, 1, 1, 1, PREVMODELD, PREVMODELI,
	       checkPow, rangePow, PREF_ALL, SUBMODEL_DEP, false, NOT_MONOTONE);
  nickname("power");
  addCov(Pow, DPow, DDPow, InversePow); 
  kappanames("alpha", REALSXP);

  
  IncludeModel("qam",  PosDefType, 2, MAXSUB, 1, kappaqam, XONLY, ISOTROPIC,
	       checkqam, rangeqam, PREF_ALL, 
	       false, SCALAR, SUBMODEL_DEP, false, NOT_MONOTONE);
  kappanames("theta", REALSXP);
  subnames("phi");
  addCov(qam, NULL, NULL);
  add_paramtype(paramtype_qam);
  

  IncludePrim("qexponential",  PosDefType, 1, XONLY, ISOTROPIC, 
	      checkOK, rangeqexponential, INFDIM, false, NOT_MONOTONE);
  nickname("qexp");
  kappanames("alpha", REALSXP);
  addCov(qexponential, Dqexponential, Inverseqexponential);

  IncludeModel("schur",  PosDefType, 1, 1, 3, kappaS, PREVMODELD, PREVMODELI, 
	       checkSchur, rangeSchur, PREF_ALL,
	       false, SUBMODEL_DEP, SUBMODEL_DEP, SUBMODEL_DEP, NOT_MONOTONE);
  kappanames("M", REALSXP, "diag", REALSXP, "rhored", REALSXP);
  addCov(Mstat, NULL, NULL);
  addCov(Mnonstat);
  add_paramtype(paramtype_M); 
 

  IncludeModel("shift", PosDefType, 1, 1, 1, kappashift, XONLY, SYMMETRIC,
	       checkshift, rangeshift, PREF_ALL, 
	       false, PARAM_DEP, SUBMODEL_DEP, SUBMODEL_DEP, NOT_MONOTONE);
  nickname("delay"); // delayeffect
  addCov(shift, NULL, NULL);
  kappanames("s", REALSXP);


  IncludePrim("spherical", TcfType, 0, NULL, XONLY, ISOTROPIC, 
	      checkOK, NULL, 3, true, GNEITING_MON);
  nickname("spheric");
  addCov(spherical, Dspherical, DDspherical, ScaleOne);
  addTBM(TBM2spherical);
  RandomShape(0, structspherical, initspherical, dospherical,
	      false, true, false);
  Taylor(-3.0, 1.0, 0.5, 3.0);


  IncludePrim("stable",  UndefinedType, 1, XONLY, ISOTROPIC, 
	      checkstable, rangestable, INFDIM, false, NORMAL_MIXTURE);
  kappanames("alpha", REALSXP);
  addCov(stable, Dstable, DDstable, Inversestable);
  addlogCov(logstable);
  addLocal(coinitstable, ieinitstable);
  addTypeFct(Typestable);

  // SPACEISOTROPIC variant of stable -- used for testing purposes only
  //  IncludePrim("stableX", 1, checkOK, SPACEISOTROPIC, 
  //		  rangestable);
  //  addCov(stableX, DstableX, Inversestable);
  //  addTBM(NULL, NULL)


  STEIN =  
    IncludeModel("Stein", PosDefType,  1, 1, 2, NULL, XONLY, ISOTROPIC,
		 check_Stein, range_Stein, plocal,
		 false, SCALAR, MAXCEDIM, true, NOT_MONOTONE);
  nickname("intrinsic");
  kappanames("diameter", REALSXP, "rawR", REALSXP);  
  add_paramtype(ignoreall_paramtype);
  addCov(Stein, NULL, NULL);
  addCallLocal(alternativeparam_Stein);
  //  RandomShape(struct_ce_approx, init_ce_approx, do_ce_approx);

 
  IncludePrim("steinst1",  PosDefType, 2, kappaSteinST1, XONLY, SYMMETRIC,
	      checkSteinST1, rangeSteinST1, INFDIM, false, NOT_MONOTONE);
  nickname("stein");
  kappanames("nu", REALSXP, "z", REALSXP);
  addCov(SteinST1, NULL, NULL);
  addTBM(initSteinST1, spectralSteinST1);
 

  IncludeModel("stp", PosDefType, 1, 2, 3, kappa_stp, KERNEL, SYMMETRIC,
	       checkstp, rangestp, PREF_ALL,
	       false, SCALAR, StpMaxDim, false, NOT_MONOTONE);
  addCov(stp);
  kappanames("S", REALSXP, "z", REALSXP, "M", REALSXP);
  addkappa(0, "S", REALSXP, ShapeType);
  RandomShape(structStp, true);
  subnames("xi2", "phi"); // H ueberall wo U-x steht. dort U-H(x)
  //                           gedoppelte immer zum Schluss!
  
  TBM_OP = // old RandomFields tbm2, tbm3
    IncludeModel("tbm",  PosDefType, 1, 1, 3, NULL, XONLY, PREVMODELI,
		 checktbmop, rangetbmop, PREF_ALL,
		 false, SUBMODEL_DEP, PARAM_DEP, PARAM_DEP, NOT_MONOTONE);
  kappanames("fulldim", INTSXP, "reduceddim", INTSXP, "layers", REALSXP); 
  addCov(tbm, NULL, NULL); // Dtbm, NULL); todo


  //  { int i; for (i=0; i<=Nothing; i++) printf("%d %d\n ", i, CovList[TREND].pref[i]); assert(false); }

  USER =
    IncludeModel("U", UndefinedType, 0, 0, 16, kappaUser, 
		 PREVMODELD, PREVMODELI,
		 checkUser, rangeUser, PREF_AUX, 
		 true,// FREEVARIABLE vorhanden. Muss extra in SpecialRMmodel.R
		 // definiert und nicht ueber generatemodels.R
		 PARAM_DEP, INFDIM, false, // per default.
		 NOT_MONOTONE);
  nickname("user");
  kappanames("type", INTSXP, "domain", INTSXP,  "isotropy", INTSXP,
	     "vdim", INTSXP, "beta", REALSXP, "variab.names", INTSXP,
	     "fctn", LANGSXP, "fst", LANGSXP, "snd", LANGSXP,
	     "envir", LANGSXP,
	     FREEVARIABLE, REALSXP, FREEVARIABLE, REALSXP, 
	     FREEVARIABLE, REALSXP, FREEVARIABLE, REALSXP,
	     FREEVARIABLE, REALSXP, FREEVARIABLE, REALSXP
	     //, "trd", LANGSXP
	     ); 
  // H ueberall wo U-x steht. dort U-H(x)
  addCov(User, DUser, DDUser, NULL, NULL);
  addCov(UserNonStat);
  addTypeFct(TypeUser);

  VECTOR = 
    IncludeModel("vector",  PosDefType, 1, 1, 2, NULL, XONLY, SYMMETRIC,
		 checkvector, rangevector, PREF_ALL, 
		 false, PARAM_DEP, SUBMODEL_DEP, SUBMODEL_DEP, NOT_MONOTONE);
  addCov(vector, NULL, NULL);
  kappanames("a", REALSXP, "Dspace", INTSXP);
  addFurtherCov(vectorAniso, NULL); 


  pref_type pwave = {2, 0, 0, 0, 5, 4, 5, 0, 0, 0, 0, 0, 0, 5};
  //           CE CO CI TBM Sp di sq Ma av n mpp Hy spf any
  IncludePrim("wave", PosDefType, 0, XONLY, ISOTROPIC, 
	      checkOK, NULL, pwave, SCALAR, 3, false, NOT_MONOTONE);
  addCov(wave, NULL, Inversewave);
  addTBM(initwave, spectralwave);
 
  IncludeModel("whittle",  UndefinedType, 0,0, 2, XONLY, ISOTROPIC, 
	       checkWM, rangeWM, PREF_ALL, INFDIM, false, NORMAL_MIXTURE);
  kappanames("nu", REALSXP, "notinvnu", INTSXP);
  addCov(Whittle, DWhittle, DDWhittle, D3Whittle, D4Whittle, InverseWhittle);
  addlogCov(logWhittle);
  addTBM(initWhittle, spectralWhittle);
  //addMarkov(&WHITTLE);
  addLocal(coinitWM, ieinitWM);
  addGaussMixture(DrawMixWM, LogMixDensW);
  addTypeFct(TypeWM);
 
 
  // *******************
  // **** shape types  ****
  // *******************
  
  BALL= IncludePrim("ball",  ShapeType, 0,  NULL, 
		    XONLY, ISOTROPIC, checkOK, NULL, PREF_NOTHING, 
		    SCALAR, INFDIM-1, true, MONOTONE);
  addCov(ball, NULL, NULL, Inverseball);
  RandomShape(INFTY, struct_ball, init_ball, do_ball);   
  
  IncludePrim("EAxxA", ShapeType,  2, kappa_EAxxA, XONLY, NO_ROTAT_INV,
	      checkEAxxA, rangeEAxxA, PREF_NOTHING, 
	      PARAM_DEP, EaxxaMaxDim, false, NOT_MONOTONE);
  nickname("eaxxa");
  addCov(EAxxA, NULL, NULL);
  kappanames("E", REALSXP, "A", REALSXP);
  addSpecial(minmaxEigenEAxxA);

  IncludePrim("EtAxxA",  ShapeType, 3, kappa_EtAxxA, XONLY, NO_ROTAT_INV,
	      checkEtAxxA, rangeEtAxxA, 3, EaxxaMaxDim, false, NOT_MONOTONE);
  nickname("etaxxa");
  addCov(EtAxxA, NULL, NULL);
  kappanames("E", REALSXP, "A", REALSXP, "alpha", REALSXP);
  addSpecial(minmaxEigenEtAxxA);

  MULT_INVERSE =
    IncludeModel("mult_inverse", ShapeType, 1, 1, 0, NULL,
		 PREVMODELD, PREVMODELI,
		 checkmult_inverse, NULL, PREF_NOTHING,
		 true, SCALAR, SUBMODEL_DEP, SUBMODEL_DEP, SUBMODEL_DEP);
  addCov(mult_inverse, NULL, NULL); 
  addCov(mult_inverseNonstat);
 
  POLYGON =
    IncludeModel("polygon",  ShapeType, 0, 0, 2, NULL, XONLY, NO_ROTAT_INV, 
	       check_polygon, range_polygon, PREF_NOTHING,
	       false, SCALAR, 2, true, MONOTONE);
  kappanames("lambda", REALSXP, "safetyfactor", REALSXP);
  addCov(Polygon, NULL, NULL, Inversepolygon, InversepolygonNonstat);
  RandomShape(INFTY, struct_polygon, init_polygon, do_polygon);   
 
  IncludePrim("rational", ShapeType, 2, kappa_rational,
	      XONLY, NO_ROTAT_INV,
	      checkrational, rangerational, INFDIM, false, NOT_MONOTONE);
  addCov(rational, NULL, NULL);
  kappanames("A", REALSXP, "a", REALSXP);
  addSpecial(minmaxEigenrational);
  
  IncludePrim("rotat",  ShapeType, 2, kappa_rotat, XONLY, NO_ROTAT_INV,
	      checkrotat, rangerotat, PREF_NOTHING, SCALAR, 3, false, NOT_MONOTONE);
  addCov(rotat, NULL, NULL);
  kappanames("speed", REALSXP, "phi", REALSXP);
  addSpecial(minmaxEigenrotat);

  IncludePrim("Rotat",  ShapeType, 1, kappa_Rotat, XONLY, NO_ROTAT_INV,
	      checkRotat, rangeRotat, PARAM_DEP, 3, false, NOT_MONOTONE);
  nickname("rotation");
  addCov(Rotat, NULL, NULL);
  kappanames("phi", REALSXP);

  RANDOMSIGN = 
    IncludeModel("sign",  ShapeType, 1, 1, 1, NULL, XONLY, PREVMODELI,
		 check_randomsign, range_randomsign, PREF_NOTHING,
		 false, SCALAR, SUBMODEL_DEP, SUBMODEL_DEP, NOT_MONOTONE);
  //nickname("");
  kappanames("p", REALSXP);
  addCov(randomsign, NULL, NULL, randomsignInverse, randomsignNonstatInverse);
  addlogCov(lograndomsign);
  RandomShape(1, struct_randomsign, init_randomsign, do_randomsign,
	      true, true, false); 
 
  SETPARAM = 
    IncludeModel("setparam", UndefinedType, 1, 1, 1, NULL, 
		 PREVMODELD, PREVMODELI,
		 checksetparam,  range_setparam, PREF_ALL, 
		 true, SUBMODEL_DEP, SUBMODEL_DEP, SUBMODEL_DEP, SUBMODEL_DEP);
  // Achtung in covmatrix_setparam wird SELECT_SUBNR verwendet!
  nickname("setparam");
  kappanames("performDo", INTSXP);
  addCov(setparamStat, Dsetparam, DDsetparam, D3setparam, D4setparam, 
	 Inverse_setparam);
  addCov(setparamNonStat);
  addTBM(NULL, spectralsetparam);
  RandomShape(0, struct_failed, initOK, dosetparam, false, false, true);
  addReturns(NULL, NULL, covmatrix_setparam, iscovmatrix_setparam, 
	     NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
  addTypeFct(Typesetparam);

  

  SHAPEAVE =
    IncludeModel("shape.ave", ShapeType, 1, 2, 3, kappa_ave,
		 XONLY, NO_ROTAT_INV,  
		 check_shapeave, rangeave, PREF_NOTHING,
		 true, SCALAR, INFDIM-1, true, NOT_MONOTONE);
  kappanames("A", REALSXP, "z", REALSXP, "spacetime", INTSXP);
  subnames("phi", "gauss");
  addlogCov(logshapeave);
  RandomShape(0, init_shapeave, do_shapeave, true);

  SHAPESTP = 
    IncludeModel("shape.stp",  ShapeType, 1, 4, 3, kappa_stp, KERNEL, 
		 NO_ROTAT_INV, check_shapestp, rangestp, PREF_NOTHING,
		 true, SCALAR, StpMaxDim, true, NOT_MONOTONE);
  kappanames("S", REALSXP, "z", REALSXP, "M", REALSXP); 
  addlogCov(logshapestp);
  subnames("xi2", "phi", "S", "gauss"); // hier gedoppeltes S nicht am Schluss 
  //                                       da auf checkstp zugreifend
  RandomShape(0, init_shapestp, do_shapestp);

  STROKORB_MONO =
    IncludeModel("strokorbMono", ShapeType, 1, 1, 0, NULL, XONLY, ISOTROPIC,
		 checkstrokorb, NULL, PREF_NOTHING,
		 false, SCALAR, 3, SUBMODEL_DEP, SUBMODEL_DEP);
  addCov(strokorb, NULL, NULL); 
  RandomShape(1, structOK, init_strokorb, do_strokorb,
	      false, false, false);

  IncludeModel("strokorbBall", ShapeType, 1, 1, 0, NULL, XONLY, ISOTROPIC,
		 checkstrokorbBall, NULL, PREF_AUX,
		 false, SCALAR, 3, true, MONOTONE);
  // addCov(strokorb, NULL, NULL); 
  RandomShape(1, struct_strokorbBall, init_failed, do_failed, do_random_failed,
	      false, false, false);
 
  STROKORB_BALL_INNER = // !! inverse scale gegenueber paper
    IncludeModel("strokorbBscale", ShapeType, 1, 1, 1, NULL,
		 XONLY, NO_ROTAT_INV,
		 check_strokorbBallInner, range_strokorbBallInner, PREF_AUX,
		 true, 1, 1, true, NOT_MONOTONE);
  kappanames("dim", INTSXP);
  addCov(strokorbBallInner, NULL, NULL);


  /*
    fktioniert so nicht !!
    da wiederum gewichtet und zwar mit b^2 falls b die intensitaet.
    kann nicht in dichte function g(b) reingezogen werden, da
    b^2 g(b) nicht integrierbar. Stattdessen darf f (Dichte im Raum)
    nicht die Gleichverteilung sein, sondern bei grossen b um
    die zu simulierenden Punkte zusammenschrumpfen.
    Dabei nimmt man an, dass die Radien ein Vielfaches des mittleren
    Radius nicht ueberschreiten. Dies ist OK, da ungefaehr 
    exponentielles Abfallen der WK.

  IncludeModel("strokorbPoly", ShapeType, 1, 1, 0, NULL, XONLY, ISOTROPIC,
		 checkstrokorbPoly, NULL, PREF_AUX,
		 false, SCALAR, 3, true, MONOTONE);
  // addCov(strokorb, NULL, NULL); 
  RandomShape(1, struct_strokorbPoly, init_failed, do_failed, do_random_failed,
	      false, false, false);
  */

  TRUNCSUPPORT =
    IncludeModel("truncsupport", ShapeType, 
		 1, 1, 1, NULL, XONLY, PREVMODELI, checktruncsupport,
		 rangetruncsupport, PREF_NOTHING, false, SCALAR,
		 SUBMODEL_DEP, SUBMODEL_DEP, SUBMODEL_DEP);
  kappanames("radius", REALSXP); // neg value == auto
  addCov(truncsupport, NULL, NULL, truncsupportInverse, StandardInverseNonstat);
  RandomShape(0, struct_truncsupport, init_truncsupport, do_truncsupport, false,
	      false, false);
  

  //////////////////////////////////////////////////
  // families of multivariate distribution; used in 

  // ACHTUNG!! addCov muss ganz zum Schluss !!

 ARCSQRT_DISTR = 
    IncludeModel("arcsqrt", RandomType, 0, 0, 0, NULL, 
		 STAT_MISMATCH, ISO_MISMATCH,
		 checkOK, NULL, PREF_AUX, 
		 true, 1, 1, false, MISMATCH);
  addCov(arcsqrtD, arcsqrtDlog, arcsqrtDinverse, 
	 arcsqrtP, NULL, arcsqrtQ, arcsqrtR, NULL);
 
  DETERM_DISTR = 
    IncludeModel("determ", RandomType, 0, 0, 1, kappa_determ, 
		 STAT_MISMATCH, ISO_MISMATCH,
		 check_determ, range_determ, PREF_AUX,
		 false, SUBMODEL_DEP, INFDIM-1, SUBMODEL_DEP, MISMATCH);
  kappanames("mean", REALSXP);
  RandomShape(INFTY, structOK, init_determ, do_determ); 
  addCov(determD, determDlog, determDinverse, determP, determP2sided, determQ, 
	 determR, determR2sided);

 
  DISTRIBUTION = // FREEVARIABLE vorhanden. Muss extra in SpecialRMmodel.R
		 // definiert und nicht ueber generatemodels.R 
    IncludeModel("distr", RandomType, 0, 0, 16, kappa_distr, 
		 STAT_MISMATCH, ISO_MISMATCH,
		 check_distr, range_distr, PREF_AUX,
		 true, PARAM_DEP, INFDIM-1, false, MISMATCH);
  kappanames("ddistr", LANGSXP, 
	     "pdistr", LANGSXP,
	     "qdistr", LANGSXP, 
	     "rdistr", LANGSXP,
	     "nrow", INTSXP,
	     "ncol", INTSXP,
	     "envir", LANGSXP, 
	     FREEVARIABLE, REALSXP, // wird nie verwendet -- Puffer fuer 
	     // einfachen Algorithmus
	     FREEVARIABLE, REALSXP, FREEVARIABLE, REALSXP, 
	     FREEVARIABLE, REALSXP, FREEVARIABLE, REALSXP,
	     FREEVARIABLE, REALSXP, FREEVARIABLE, REALSXP,
	     FREEVARIABLE, REALSXP, FREEVARIABLE, REALSXP
	     ); // 7 free ones are remaining !	     
  RandomShape(0, structOK, init_distr, do_distr_do);
  addCov(distrD, distrDlog, distrDinverse, distrP, distrP2sided, distrQ,
	 distrR, distrR2sided);
 
  GAUSS_DISTR = 
    IncludeModel("normal", RandomType, 0, 0, 3, kappa_gauss_distr, 
		 STAT_MISMATCH, ISO_MISMATCH,
		 check_gauss_distr, range_gauss_distr, PREF_AUX, 
		 false, PARAM_DEP, INFDIM-1, false, MISMATCH);
  nickname("gauss");
  kappanames("mu", REALSXP, "sd", REALSXP, "log", REALSXP);
  RandomShape(INFTY, structOK, init_gauss_distr, do_gauss_distr);
  addCov(gaussD, gaussDlog, gaussDinverse, 
	 gaussP, gaussP2sided, gaussQ, gaussR, gaussR2sided);
 
  /*
  SET_DISTR = 
    IncludeModel("set", RandomType, 1, 1, 1, NULL, STAT_MISMATCH, ISO_MISMATCH,
		 check_setParam, range_setParam, PREF_AUX, 
		 true, SUBMODEL_DEP, SUBMODEL_DEP, SUBMODEL_DEP, MISMATCH);
  subnames("to");
  kappanames("variant", INTSXP);
  RandomShape(INFTY, structOK, init_setParam, do_setParam);
  addCov(setParamD, setParamDlog, setParamDinverse, 
	 setParamP, setParamP2sided, setParamQ,
	 setParamR, setParamR2sided);
  */

  LOC =
    IncludeModel("loc", RandomType, 1, 1, 2, kappa_loc, 
		 STAT_MISMATCH, ISO_MISMATCH, 
		 check_loc, range_loc, PREF_AUX, 
		 false, SUBMODEL_DEP, SUBMODEL_DEP, SUBMODEL_DEP, MISMATCH);
  kappanames("mu", REALSXP, "scale", REALSXP);
  RandomShape(INFTY, structOK, init_loc, do_loc);
  addCov(locD, locDlog, locDinverse, locP, locP2sided, locQ, locR, locR2sided);

  RECTANGULAR =
    IncludeModel("rectangular", RandomType, 1, 1, 11, NULL, 
		 // ACHTUNG! Model kann auch ueber cov->q uebergeben werden.
		 //          dies vereinfacht die Verwendung von zufaelligen
		 //          Huetchen, da keine Parameter kopiert werden 
		 //          muesen, sondern direkt auf das Huetchen zugegriffen
		 //         
		 STAT_MISMATCH, ISO_MISMATCH,
		 check_rectangular, range_rectangular, PREF_AUX, 
		 false, PARAM_DEP,  INFDIM-1, true, MISMATCH);
  kappanames("safety", REALSXP, "minsteplen", REALSXP, "maxsteps", INTSXP,
	     "parts", INTSXP, "maxit", INTSXP, "innermin", REALSXP,
	     "outermax", REALSXP, "mcmc_n", INTSXP,
	     "normed", INTSXP, "approx", INTSXP, "onesided", INTSXP
	     );
  RandomShape(INFTY, structOK, init_rectangular, do_rectangular); 
  addCov(rectangularD, rectangularDlog, rectangularDinverse, rectangularP, 
	 rectangularP2sided, rectangularQ, rectangularR, rectangularR2sided);

  SCALESPHERICAL = 
    IncludeModel("spherical", RandomType, 0, 0, 2, NULL,
		 STAT_MISMATCH, ISO_MISMATCH,
		 check_RRspheric, range_RRspheric, PREF_AUX,
		 false, 1, 1, true, MISMATCH);
  kappanames("spacedim", INTSXP, "balldim", INTSXP);
  RandomShape(INFTY, structOK, init_RRspheric, do_RRspheric);
  addCov(sphericD, sphericDlog, sphericDinverse, sphericP, NULL, sphericQ,
	 sphericR, NULL);
 
 
  UNIF = IncludeModel("unif", RandomType, 0, 0, 2, kappa_unif, 
		      STAT_MISMATCH, ISO_MISMATCH,
		      check_unif, range_unif, PREF_AUX, 
		      false, PARAM_DEP,  INFDIM-1, true, MISMATCH);
  kappanames("min", REALSXP, "max", REALSXP);
  RandomShape(INFTY, structOK, init_unif, do_unif); 
  addCov(unifD, unifDlog, unifDinverse, unifP, unifP2sided, unifQ, 
	 unifR, unifR2sided);

  
  // -----------------------------
  // shape + locations 
  // they *take* all very detailed roles like ROLE_SMITH and pass
  // ROLE_MAXSTABLE to the submodel, in general
  // storage always pgs_storage !!
  PTS_GIVEN_SHAPE = 
    IncludeModel("ptsGivenShape", PointShapeType, 2, 2, 2, NULL, 
		 XONLY, NO_ROTAT_INV,
		 check_pts_given_shape, range_pts_given_shape, PREF_AUX, 
		 true, SUBMODEL_DEP, SUBMODEL_DEP, SUBMODEL_DEP, MISMATCH);
  kappanames("density.ratio", REALSXP,  // stemming from gauss
	     "flat", INTSXP	  
	     );
  subnames("shape", "loc");
  addCov(pts_given_shape, NULL, NULL);
  addlogCov(logpts_given_shape);
  RandomShape(SUBMODEL_DEP, struct_pts_given_shape, init_pts_given_shape, 
	      do_pts_given_shape, do_random_failed, true, true, false); 

  STANDARD_SHAPE = 
    IncludeModel("standardShape", PointShapeType, 1, 2, 0, NULL, 
		 XONLY, NO_ROTAT_INV,
		 check_standard_shape, NULL, PREF_AUX, 
		 true, SUBMODEL_DEP, SUBMODEL_DEP, SUBMODEL_DEP, MISMATCH);
  subnames("shape");
  addCov(standard_shape, NULL, NULL);
  addlogCov(logstandard_shape);
  RandomShape(SUBMODEL_DEP, struct_standard_shape, init_standard_shape, 
	      do_standard_shape, do_random_failed, true, true, false); 

  STATIONARY_SHAPE = 
    IncludeModel("poissonShape", PointShapeType, 1, 1, 0, NULL, 
		 XONLY, NO_ROTAT_INV,
		 check_stationary_shape, NULL, PREF_AUX, 
		 true, SUBMODEL_DEP, SUBMODEL_DEP, SUBMODEL_DEP, MISMATCH);
  subnames("shape");
  addCov(stationary_shape, NULL, NULL);
  addlogCov(logstationary_shape);
  RandomShape(SUBMODEL_DEP, struct_stationary_shape, init_stationary_shape, 
	      do_stationary_shape, do_random_failed, true, true, false); 

  pref_type pmppp =  {0, 0, 0, 0, 0, 0, 0, 0, 5, 0, 5, 0, 10, 5};
  //           CE CO CI TBM Sp di sq Ma av n mpp Hy spf any
  MPPPLUS =
    IncludeModel("++",  PointShapeType, 1, MAXSUB, 1, kappamppplus, 
		 PREVMODELD, PREVMODELI, // NO_ROTAT_INV,
		 checkmppplus, rangempplus, pmppp, 
		 false, SUBMODEL_DEP, SUBMODEL_DEP, SUBMODEL_DEP, MISMATCH);
  nickname("mppplus");
  kappanames("p", REALSXP);
  addCov(mppplus, NULL, NULL);
  //  addlogCov(logmppplus);
  RandomShape(SUBMODEL_DEP, struct_mppplus, init_mppplus, do_mppplus,  
	      true, true, false); 


  // -----------------------------
  // Interfaces
  // -----------------------------
  IncludeModel("Simulate", InterfaceType, 1, 1, 0, NULL, 
	       XONLY, NO_ROTAT_INV, 
	       check_simulate, NULL, PREF_AUX, 
	       true, SUBMODEL_DEP, SUBMODEL_DEP, SUBMODEL_DEP, MISMATCH);
  addCov(simulate, NULL, NULL);
  RandomShape(struct_simulate); 
 
  COVFCTN =
    IncludeModel("Cov", InterfaceType, 1, 1, 0, NULL, XONLY, NO_ROTAT_INV, 
		 check_cov, NULL, PREF_AUX, 
		 true, SUBMODEL_DEP, SUBMODEL_DEP, SUBMODEL_DEP, MISMATCH);
  addCov(Cov, NULL, NULL);
  RandomShape(struct_cov); 

  COVMATRIX = 
    IncludeModel("CovMatrix", InterfaceType, 1, 1, 0, NULL, XONLY,
		 NO_ROTAT_INV, //NO_ROTAT_INV,ISOTROPIC dependening on whether
		 // distances are givenXONLY, NO_ROTAT_INV,
		 check_covmatrix, NULL, PREF_AUX, 
		 true, SUBMODEL_DEP, SUBMODEL_DEP, SUBMODEL_DEP, MISMATCH);
  addCov(CovMatrix, NULL, NULL);
  RandomShape(struct_cov); 

  IncludeModel("Dummy", InterfaceType, 1, 1, 0, NULL, XONLY, NO_ROTAT_INV, 
	       check_dummy, NULL, PREF_AUX, 
	       true, SUBMODEL_DEP, SUBMODEL_DEP, SUBMODEL_DEP, MISMATCH);
  addCov(Dummy, NULL, NULL);
  RandomShape(struct_dummy); 

  RFGET = 
    IncludeModel("get", InterfaceType, 1, 1, 2, NULL, XONLY, NO_ROTAT_INV,
		 check_RFget, range_RFget, PREF_AUX, 
		 true, SUBMODEL_DEP, SUBMODEL_DEP, SUBMODEL_DEP, MISMATCH);
  kappanames("up", INTSXP, "register", INTSXP);
  addCov(RFget, NULL, NULL);
  RandomShape(struct_RFget); 

 
  IncludeModel("Fctn", InterfaceType, 1, 1, 0, NULL, XONLY, NO_ROTAT_INV, 
		 check_fctn, NULL, PREF_AUX, 
	       true, SUBMODEL_DEP, SUBMODEL_DEP, SUBMODEL_DEP, MISMATCH);
  addCov(Cov, NULL, NULL);
  RandomShape(structOK); 

  IncludeModel("Distr", InterfaceType, 1, 1, 5, kappa_EvalDistr,
	       XONLY, NO_ROTAT_INV,
	       check_EvalDistr, range_EvalDistr, PREF_AUX, 
	       true, SUBMODEL_DEP, SUBMODEL_DEP, SUBMODEL_DEP, MISMATCH);
  kappanames("x", REALSXP, "q", REALSXP, "p", REALSXP, "n", REALSXP,
	     "dim", INTSXP);
  addCov(EvalDistr, NULL, NULL);
  RandomShape(struct_EvalDistr); 


  IncludeModel("Pseudovariogram", InterfaceType, 1, 1, 0, NULL, 
	       XONLY, NO_ROTAT_INV,
	       check_cov, NULL, PREF_AUX, 
	       true, SUBMODEL_DEP, SUBMODEL_DEP, SUBMODEL_DEP, MISMATCH);
  nickname("Pseudovario");
  addCov(Pseudovariogram, NULL, NULL);
  RandomShape(struct_variogram); 
  
  VARIOGRAM_CALL =
    IncludeModel("Variogram", InterfaceType, 1, 1, 0, NULL, 
		 XONLY, NO_ROTAT_INV,
		 check_vario, NULL, PREF_AUX, 
		 true, SUBMODEL_DEP, SUBMODEL_DEP, SUBMODEL_DEP, MISMATCH);
  addCov(Variogram, NULL, NULL);
  RandomShape(struct_variogram); 


  
  // ----------------------------
  // processes  //        CE CO CI TBM Sp di sq Ma av n mpp Hy spf any

  DOLLAR_PROC 
    = IncludeModel("proc$", ProcessType,		   
		   1, 1, 5, kappaS, // kappadollar,
		   XONLY, NO_ROTAT_INV, checkS, rangeS, PREF_ALL,
		   true, SUBMODEL_DEP, SUBMODEL_DEP, SUBMODEL_DEP,
		   SUBMODEL_DEP);
  // do not change Order!!
  nickname("S");
  kappanames("var", REALSXP, "scale", REALSXP, "anisoT", REALSXP,
	     "Aniso", REALSXP, "proj", INTSXP);
  add_paramtype(ignoreall_paramtype);
  subnames("phi");
  RandomShape(2, structSproc, initSproc, doSproc, true, true, true);

  
  PLUS_PROC = 
    IncludeModel("plusproc", ProcessType, 1, MAXSUB, 0, NULL, 
		  XONLY, NO_ROTAT_INV,
		 checkplusproc, NULL, PREF_ALL, 
		 true, SUBMODEL_DEP, SUBMODEL_DEP, SUBMODEL_DEP, SUBMODEL_DEP);
  // Achtung in covmatrix_plus wird SELECT_SUBNR verwendet!
  nickname("plus");
  RandomShape(2, structplusproc, initplusproc, doplusproc, false, false, true);

   
  AVERAGE_USER = 
    IncludeModel("average", GaussMethodType, 1, 2, 2, NULL, 
		 XONLY, NO_ROTAT_INV,
		 check_randomcoin, range_randomcoin, PREF_NOTHING,  
		 false, SCALAR, MAXMPPDIM, false, MISMATCH);
  kappanames("loggauss",  INTSXP, "intensity", REALSXP);  
  add_paramtype(ignoreall_paramtype);
  subnames("phi", "shape");
  // addCov(coin, NULL, NULL, coinInverse);
  RandomShape(2, struct_extractdollar, init_gaussprocess, do_gaussprocess); 

  AVERAGE_INTERN = 
    CopyModel("averageIntern", AVERAGE_USER);
  make_internal(); 
  RandomShape(2, struct_randomcoin, init_randomcoin, dompp, true, true, false);


  CIRCEMBED = // und die anderen fehlen auch noch !!
    IncludeModel("circulant", GaussMethodType, 1, 1, 12, kappa_ce,
		 XONLY, NO_ROTAT_INV,
		 check_ce, range_ce, PREF_NOTHING,
		 false, SUBMODEL_DEP, MAXCEDIM, false, MISMATCH);
  kappanames("loggauss", INTSXP, 
	     "force", INTSXP, "mmin", REALSXP, "strategy", INTSXP,
	     "maxmem", INTSXP, "tolIm", REALSXP, "tolRe", REALSXP,
	     "trials", INTSXP, "useprimes", INTSXP, "dependent", INTSXP,
	     "approx_step", REALSXP,  "approx_maxgrid", INTSXP);
  add_paramtype(ignoreall_paramtype);
  RandomShape(2, struct_ce_approx, init_ce_approx, do_ce_approx);
 
  CE_CUTOFFPROC_USER =  
    IncludeModel("cutoff", GaussMethodType, 1, 1, 14, kappa_localproc,
		 XONLY, NO_ROTAT_INV,
		 check_local_proc, range_co_proc, PREF_NOTHING,
		 false, SCALAR, MAXCEDIM,  false, MISMATCH);
  kappanames("loggauss", INTSXP, 
	     "force", INTSXP, "mmin", REALSXP, "strategy", INTSXP,
	     "maxmem", INTSXP, "tolIm", REALSXP, "tolRe", REALSXP,
	     "trials", INTSXP, "useprimes", INTSXP, "dependent", INTSXP,
	     "approx_step", REALSXP,  "approx_maxgrid", INTSXP, 
	     "diameter", REALSXP, "a", REALSXP);  
  add_paramtype(ignoreall_paramtype);
  RandomShape(2, struct_extractdollar, init_gaussprocess, do_gaussprocess); 
 
  CE_CUTOFFPROC_INTERN = CopyModel("cutoffIntern", CE_CUTOFFPROC_USER);
  make_internal();
  RandomShape(2, struct_ce_approx, init_ce_approx, do_ce_approx);

  CE_INTRINPROC_USER =  
    IncludeModel("intrinsic", GaussMethodType, 1, 1, 14, kappa_localproc, 
		 XONLY, NO_ROTAT_INV,
		 check_local_proc, range_intrinCE, PREF_NOTHING,
		 false, SCALAR, MAXCEDIM, false, MISMATCH);
  kappanames("loggauss", INTSXP, 
	     "force", INTSXP, "mmin", REALSXP, "strategy", INTSXP,
	     "maxmem", INTSXP, "tolIm", REALSXP, "tolRe", REALSXP,
	     "trials", INTSXP, "useprimes", INTSXP, "dependent", INTSXP,
	     "approx_step",REALSXP, "approx_maxgrid",INTSXP, "diameter",REALSXP,
	     "rawR", REALSXP);  
  add_paramtype(ignoreall_paramtype);
  RandomShape(2, struct_extractdollar, init_gaussprocess, do_gaussprocess); 
 
  CE_INTRINPROC_INTERN = CopyModel("intrinsIntern", CE_INTRINPROC_USER);
  make_internal();
  RandomShape(2, struct_ce_approx, init_ce_approx, do_ce_approx);


  DIRECT = 
    IncludeModel("direct",GaussMethodType,  1, 1, 4, NULL, XONLY, NO_ROTAT_INV,
		 check_directGauss, range_direct, PREF_NOTHING,
		 false,  SUBMODEL_DEP, INFDIM-1, false, MISMATCH);
  kappanames("loggauss", INTSXP, 
	     "root_method", INTSXP, "svdtolerance", REALSXP, "max_variab", INTSXP);
  add_paramtype(ignoreall_paramtype);
  RandomShape(2, init_directGauss, do_directGauss);


  HYPERPLANE_USER  = 
    IncludeModel("hyperplane", GaussMethodType, 1, 1, 5, NULL, 
		 XONLY, NO_ROTAT_INV,
		 check_hyperplane, range_hyperplane, PREF_NOTHING,  
		 false, SCALAR, 2, false, MISMATCH);
  kappanames("loggauss", INTSXP,
	     "superpos", INTSXP, "maxlines", INTSXP, "mar_distr", INTSXP, 
	     "mar_param", REALSXP);
  add_paramtype(ignoreall_paramtype);
  //  addCov(IdStat, NULL, NULL, IdInverse);
  //  addCov(IdNonStat);
  RandomShape(2, struct_extractdollar, init_gaussprocess, do_gaussprocess); 

  HYPERPLANE_INTERN =   
    CopyModel("hyperIntern", HYPERPLANE_USER, check_hyperplane_intern);
  make_internal();
  RandomShape(2, struct_hyperplane, init_hyperplane, do_hyperplane);
  //  printf("%ld %ld\n", CovList[HYPERPLANE_INTERN].nonstat_cov,
  //	 IdNonStat); assert(false);

  NUGGET_USER  = 
    IncludeModel("nugget", GaussMethodType, 1, 1, 3, NULL, XONLY, NO_ROTAT_INV,
		 check_nugget_proc, range_nugget_proc, PREF_NOTHING, false, 
		 PREVMODEL_DEP, INFDIM, true, MISMATCH);
  kappanames("loggauss", INTSXP, 
	     "tol", REALSXP, "vdim", INTSXP);
  add_paramtype(ignoreall_paramtype);
  RandomShape(2, struct_extractdollar, init_gaussprocess, do_gaussprocess); 


  NUGGET_INTERN =   
    CopyModel("NuggetIntern", NUGGET_USER);
  make_internal();
  RandomShape(2, struct_nugget, init_nugget, do_nugget);

  /* see simu.cc, CMbuild for special treatment of nugget when
     users choice is given */
  
  /* cf. convert.R, PrepareModel, near end of function */
  /* deleted 27.12.08
     { 
     int nr = currentNrCov - 1;
     for (i=0; i<MethodsLength; i++) {
     CovList[nr].implemented[i] = GIVEN_METH_IGNORED;
     }
     CovList[nr].implemented[Nugget] = 
     CovList[nr].implemented[Direct] =
     CovList[nr].implemented[Sequential] =
     CovList[nr].implemented[CircEmbed] =IMPLEMENTED;
     }
  */
 
  RANDOMCOIN_USER = CopyModel("coins", AVERAGE_USER);

 
  SEQUENTIAL = 
    IncludeModel("sequential", GaussMethodType, 1, 1, 4, NULL, 
		 XONLY, NO_ROTAT_INV,
		 check_sequential, range_sequential, PREF_NOTHING,  
		 false, SCALAR, INFDIM-1, false, MISMATCH);
  kappanames("loggauss", INTSXP, 
	     "max_variables", INTSXP, "back_steps", INTSXP, "initial", INTSXP);
  add_paramtype(ignoreall_paramtype); 
  RandomShape(2, init_sequential, do_sequential);


  SPECTRAL_PROC_USER = 
    IncludeModel("spectral", GaussMethodType,  1, 1, 5, NULL,
		 XONLY, NO_ROTAT_INV,
		 check_spectral, range_spectral, PREF_NOTHING,  
		 false, SCALAR, MAXTBMSPDIM, false, MISMATCH);
  kappanames("loggauss", INTSXP, 
	     "sp_lines", INTSXP, "sp_grid", INTSXP,
	     "prop_factor", REALSXP, "sigma", REALSXP );
  add_paramtype(ignoreall_paramtype);
  RandomShape(2, struct_extractdollar, init_gaussprocess, do_gaussprocess); 
 
  SPECTRAL_PROC_INTERN = 
    CopyModel("spectralIntern", SPECTRAL_PROC_USER);
  make_internal();
  RandomShape(2, struct_spectral, init_spectral, do_spectral);


  SPECIFIC = 
    IncludeModel("specific", GaussMethodType, 1, 1, 1, NULL, 
		 XONLY, NO_ROTAT_INV,
		 check_specificGauss, range_common_gauss, PREF_NOTHING,  
		 false, SUBMODEL_DEP, MAXTBMSPDIM, false, MISMATCH);
  kappanames("loggauss", INTSXP);
  RandomShape(2, struct_specificGauss, init_specificGauss, do_specificGauss);
  

  TBM_PROC_USER = 
    IncludeModel("tbm", GaussMethodType,  1, 1, 9, tbm_kappasproc, 
		 XONLY, NO_ROTAT_INV, 
		 checktbmproc, rangetbmproc, PREF_NOTHING,
		 false, SCALAR, SUBMODEL_DEP, false, MISMATCH);
  kappanames("loggauss", INTSXP, 
	     "fulldim", INTSXP, "reduceddim", INTSXP, "layers", REALSXP,
	     "lines", INTSXP, "linessimufactor", REALSXP,
	     "linesimustep",  REALSXP,  // "grid", INTSXP,
	     "center", REALSXP, "points", INTSXP); 
  add_paramtype(ignoreall_paramtype);
  // addFurtherCov(tbm2num, NULL);
  RandomShape(2, struct_extractdollar, init_gaussprocess, do_gaussprocess); 
 
  TBM_PROC_INTERN = 
    CopyModel("tbmIntern", TBM_PROC_USER);
  make_internal();
  RandomShape(2, struct_tbmproc, init_tbmproc, do_tbmproc); 
 

  gaussmethod[CircEmbed] = CIRCEMBED;
  gaussmethod[CircEmbedCutoff]= CE_CUTOFFPROC_INTERN;
  gaussmethod[CircEmbedIntrinsic] = CE_INTRINPROC_INTERN;
  gaussmethod[TBM] = TBM_PROC_INTERN;
  gaussmethod[SpectralTBM] = SPECTRAL_PROC_INTERN;
  gaussmethod[Direct] = DIRECT;
  gaussmethod[Sequential] = SEQUENTIAL;
  gaussmethod[Markov] = -1;
  gaussmethod[Average] = AVERAGE_INTERN;
  gaussmethod[Nugget] = NUGGET_INTERN;
  gaussmethod[RandomCoin] = AVERAGE_INTERN;
  gaussmethod[Hyperplane] = HYPERPLANE_INTERN;
  gaussmethod[Specific] = SPECIFIC;
  gaussmethod[Nothing] =  gaussmethod[Forbidden] = -1;
 

  // non sub-gaussian processe
  BRORIGINAL_USER =
    IncludeModel("brorig", BrMethodType, 1, 2, 3, NULL, 
		 XONLY, NO_ROTAT_INV,
		 checkBrownResnick, range_mpp, PREF_NOTHING,
		 false, SUBMODEL_DEP, MAXMPPDIM, false, MISMATCH);
  addkappa(0, "xi", REALSXP, RandomType);
  addkappa(1, "mu", REALSXP, ShapeType);
  addkappa(2, "s",  REALSXP, ShapeType);
  subnames("phi", "tcf");
  add_paramtype(ignoreall_paramtype);
  RandomShape(0, structBRuser, initBRuser, dompp); 

  BRMIXED_USER =
    IncludeModel("brmixed", BrMethodType, 1, 2, 12, kappaBRmixed, 
		 XONLY, NO_ROTAT_INV, 
		 check_BRmixed, range_BRmixed, PREF_NOTHING,
		 false, SUBMODEL_DEP, MAXMPPDIM, false, MISMATCH);
  kappanames("xi", REALSXP, "mu", REALSXP, "s",  REALSXP, 
	     "meshsize", REALSXP,
             "lowerbound_optimarea", REALSXP, "vertnumber", INTSXP,
             "optim_mixed", INTSXP, "optim_mixed_tol", REALSXP, 
	     "optim_mixed_maxpoints", INTSXP, "lambda", REALSXP,
	     "areamat", REALSXP, "variobound", REALSXP);
  addkappa(1, "mu", REALSXP, ShapeType);
  addkappa(2, "s",  REALSXP, ShapeType);
  subnames("phi", "tcf");
  add_paramtype(ignoreall_paramtype);
  RandomShape(0, structBRuser, initBRuser, dompp);
  
  BRSHIFTED_USER = 
    IncludeModel("brshifted", BrMethodType, 
		 1, 2, 3, NULL, 
		 XONLY, NO_ROTAT_INV, 
		 checkBrownResnick, range_mpp, PREF_NOTHING,
		 false, SUBMODEL_DEP, MAXMPPDIM, false, MISMATCH);
  addkappa(0, "xi", REALSXP, RandomType);
  addkappa(1, "mu", REALSXP, ShapeType);
  addkappa(2, "s",  REALSXP, ShapeType);
  subnames("phi", "tcf");
  add_paramtype(ignoreall_paramtype);
  RandomShape(0, structBRuser, initBRuser, dompp);
  
  BRORIGINAL_INTERN =
    CopyModel("brorigIntern", BRORIGINAL_USER, PointShapeType);
  make_internal();
  RandomShape(SUBMODEL_DEP, structBRintern, init_BRorig, do_BRorig);
  
  BRMIXED_INTERN =
    CopyModel("brmixedIntern", BRMIXED_USER, PointShapeType); 
  make_internal();
  RandomShape(SUBMODEL_DEP, structBRintern, init_BRmixed, do_BRmixed);
  
  BRSHIFTED_INTERN =
    CopyModel("brshiftIntern", BRSHIFTED_USER, PointShapeType); 
  make_internal();
  RandomShape(SUBMODEL_DEP, structBRintern, init_BRshifted, do_BRshifted);
  
  // distributions

  
  BINARYPROC = 
    IncludeModel("binaryprocess", ProcessType, 1, 1, 3, NULL, 
		 XONLY, NO_ROTAT_INV, 
		 checkbinaryprocess, rangebinaryprocess, PREF_NOTHING,
		 false, SCALAR, MAXSIMUDIM, false, MISMATCH);
  nickname("bernoulli");
  kappanames("loggauss", INTSXP,  
	     "stationary_only", INTSXP, 
	     "threshold", REALSXP);
  add_paramtype(ignoreall_paramtype);
  RandomShape(INFTY, struct_binaryprocess, init_binaryprocess,
	      do_binaryprocess);
  
  BROWNRESNICKPROC =
    IncludeModel("brownresnick",ProcessType,  1, 2, 3, NULL, 
		 XONLY, NO_ROTAT_INV, 
		 checkBrownResnick, range_mpp, PREF_NOTHING,  
		 false, SUBMODEL_DEP, MAXMPPDIM, false, MISMATCH);
  addkappa(0, "xi", REALSXP, RandomType);
  addkappa(1, "mu", REALSXP, ShapeType);
  addkappa(2, "s",  REALSXP, ShapeType);
  subnames("phi", "tcf");
  //  addCov(BrownResnick, NULL, NULL);
  add_paramtype(ignoreall_paramtype);
  RandomShape(0, structBrownResnick, initBrownResnick, doBrownResnick); 

  GAUSSPROC = 
    IncludeModel("gauss.process", ProcessType, 1, 1, 2, NULL, 
		 XONLY, NO_ROTAT_INV,
		 checkgaussprocess, rangegaussprocess, PREF_NOTHING,
		 false, SUBMODEL_DEP, MAXSIMUDIM, false, MISMATCH);
  nickname("gauss");
  kappanames("loggauss", INTSXP, "stationary_only", INTSXP);  
  add_paramtype(ignoreall_paramtype);
  RandomShape(2, struct_gaussprocess, init_gaussprocess, do_gaussprocess);


  POISSONPROC =
    IncludeModel("poisson", ProcessType, 1, 1, 1, NULL,
		 XONLY, NO_ROTAT_INV,  
		 check_poisson, range_poisson, PREF_NOTHING,
		 false, SUBMODEL_DEP, MAXMPPDIM, false, MISMATCH);
  kappanames("intensity", REALSXP);
  add_paramtype(ignoreall_paramtype);
  RandomShape(0, struct_poisson, init_poisson, dompp);


  SCHLATHERPROC =
    IncludeModel("extremalgauss", ProcessType, 1, 2, 3, NULL, 
		 XONLY, NO_ROTAT_INV,  
		 check_schlather, range_mpp, PREF_NOTHING, 
		 false, SCALAR, MAXSIMUDIM, false, MISMATCH);
  nickname("schlather");
  addkappa(0, "xi", REALSXP, RandomType);
  addkappa(1, "mu", REALSXP, ShapeType);
  addkappa(2, "s",  REALSXP, ShapeType);
  subnames("phi", "tcf");
  add_paramtype(ignoreall_paramtype);
  addCov(extremalgaussian, NULL, NULL);
  RandomShape(0, struct_schlather, init_mpp, dompp);

  
  SMITHPROC =
    IncludeModel("smith", ProcessType, 1, 2, 3, NULL, XONLY, NO_ROTAT_INV,
		 check_smith, range_mpp, PREF_NOTHING,  
		 false, SCALAR, MAXMPPDIM, false, MISMATCH);
  addkappa(0, "xi", REALSXP, RandomType);
  addkappa(1, "mu", REALSXP, ShapeType);
  addkappa(2, "s",  REALSXP, ShapeType);
  subnames("shape", "tcf");
  add_paramtype(ignoreall_paramtype);
  RandomShape(0, struct_smith, init_mpp, dompp);
 
  CHI2PROC =
    IncludeModel("chi2", ProcessType, 1, 1, 1, NULL, XONLY, NO_ROTAT_INV,
		 checkchisqprocess, rangechisqprocess, PREF_NOTHING,
		 false, SCALAR, MAXSIMUDIM, false, MISMATCH);
  kappanames("f", INTSXP);  
  add_paramtype(ignoreall_paramtype);
  RandomShape(0, struct_chisqprocess, init_chisqprocess, do_chisqprocess);

 
  //printf("%d\n",   CHI2PROC);  assert(false);


  // addMarkov(MarkovhWittle);
  // Markovcircular, Markovexponential, Markovcubic, Markovdampedcosine, MarkovFD, 
  // Markovgauss, Markovgneiting, Markovgengneiting, Markovhyperbolic, Markovpenta, 
  // Markovpower, Markovqexponential, Markovspherical, Markovstable, 
 
  // addusersfunctions();

  int nr;
  for (nr=0; nr<currentNrCov; nr++) { 
    cov_fct *C = CovList + nr; // nicht gatternr
    // printf("name = %s %d\n", C->nick, nr);
    assert(C->Type != UndefinedType || C->TypeFct != NULL);

    assert({int k;
	for (k=0; k<C->kappas; k++) {
	  if (C->kappanames[k][0] == 'k' && C->kappanames[k][1] >= '0'
	      && C->kappanames[k][1] <= '9') {
	    PRINTF("%s %s\n", C->nick, C->kappanames[k]);
	    break;
	  }
	}; k >= C->kappas;
      });

    C->pref[Markov] = PREF_NONE;
    if (C->Type == RandomType) {    
       // printf("done = %s %d\n", CovList[nr].nick, nr);
      continue;
    }
    // int kappas = C->kappas;
    if (C->Type == NegDefType) 
      C->pref[TBM] = C->pref[CircEmbed] =C->pref[SpectralTBM] = PREF_NONE;
    if (nr != NUGGET){ // && (nr < DOLLAR || nr > LASTDOLLAR)) { 	
      for (i=0; i<Nothing; i++) {
	//	if (i == Specific) continue;


	//	if (nr == 36)
	//	printf("%s %s %d %d\n", C->nick, METHODNAMES[i], C->pref[i], C->implemented[i]);
	//
	//if (C->pref[i] > 0 &&  C->implemented[i] != IMPLEMENTED) {
	  //	printf("%d %d %s %s %d %d\n", nr, i, C->nick, METHODNAMES[i], 
	//	C->pref[i], C->implemented[i]);}
	C->pref[i] *= C->implemented[i] == IMPLEMENTED;
	//	printf(" %d\n", C->pref[i]);
	
      }
      //
      i = Specific;
      //printf("%d  %s %s %d %d\n", nr, C->nick, METHODNAMES[i], C->pref[i],
      //	     C->implemented[i]);
      i = Nothing;
      C->pref[i] *= (C->cov != ErrCov || C->nonstat_cov != ErrCovNonstat);
      // printf(" %d\n", C->pref[i]);
    }
  }
  
  //  PRINTF("simple checks auskommentiert\n");
  // assert(SimpleChecks());

}


bool isDollar(cov_model *cov) {
  int nr=cov->nr;
  return nr >= DOLLAR && nr <= LASTDOLLAR;
}


bool isDollarProc(cov_model *cov) {
  int nr=cov->nr;
  return nr == DOLLAR_PROC;
}

bool isAnyDollar(cov_model *cov) {
  int nr=cov->nr;
  return (nr >= DOLLAR && nr <= LASTDOLLAR) || nr == DOLLAR_PROC;
}

bool isTcf(Types type) {
  return type == TcfType;
}

bool isTcf(cov_model *cov) {
  Types type = CovList[cov->nr].Type;
  if (type != UndefinedType) return type == TcfType;

  assert(CovList[cov->nr].TypeFct != NULL);
  return CovList[cov->nr].TypeFct(TcfType, cov);
}

bool isPosDef(Types type) {
  return type == PosDefType || type == TcfType;
}

bool isPosDef(cov_model *cov) {
  Types type = CovList[cov->nr].Type;
  if (type != UndefinedType) return type == PosDefType|| type == TcfType;

  assert(CovList[cov->nr].TypeFct != NULL);
  return CovList[cov->nr].TypeFct(PosDefType, cov);
}

bool isNegDef(Types type) {
  return isPosDef(type) || type == NegDefType;
}

bool isNegDef(cov_model *cov) {
  Types type = CovList[cov->nr].Type;
  if (type != UndefinedType) return  isPosDef(type) || type == NegDefType;

  assert(CovList[cov->nr].TypeFct != NULL);
  return CovList[cov->nr].TypeFct(NegDefType, cov);
}

bool isProcess(Types type) {
  return 
    type == ProcessType || type == GaussMethodType || type == BrMethodType ;
}

bool isProcess(cov_model *cov) {
  Types type = CovList[cov->nr].Type;
  if (type != UndefinedType) 
    return type == ProcessType || type == GaussMethodType ||
      type == BrMethodType;

  assert(CovList[cov->nr].TypeFct != NULL);
  return CovList[cov->nr].TypeFct(ProcessType, cov);
}

bool isGaussMethod(Types type) {
  return type == GaussMethodType;
}

bool isGaussMethod(cov_model *cov) {
  Types type = CovList[cov->nr].Type;
  if (type != UndefinedType) return type == GaussMethodType;

  assert(CovList[cov->nr].TypeFct != NULL);
  return CovList[cov->nr].TypeFct(GaussMethodType, cov);
}

bool isPointShape(Types type) {
  return type == PointShapeType;
}

bool isPointShape(cov_model *cov) {
  Types type = CovList[cov->nr].Type;
  if (type != UndefinedType) return type == PointShapeType;

  assert(CovList[cov->nr].TypeFct != NULL);
  return CovList[cov->nr].TypeFct(PointShapeType, cov);
}

bool isRandom(Types type) {
  return type == RandomType;
}

bool isRandom(cov_model *cov) {
  Types type = CovList[cov->nr].Type;
  if (type != UndefinedType) return type == RandomType;

  assert(CovList[cov->nr].TypeFct != NULL);
  return CovList[cov->nr].TypeFct(RandomType, cov);
}

bool isShape(Types type) {
  return type == ShapeType || isNegDef(type);  
}

bool isShape(cov_model *cov) {
  Types type = CovList[cov->nr].Type;
 if (type != UndefinedType)  return type == ShapeType || isNegDef(type);  // || isTend(type) ?? to do ??

  assert(CovList[cov->nr].TypeFct != NULL);
  return CovList[cov->nr].TypeFct(ShapeType, cov);
}

bool isTrend(Types type) {
  return type == TrendType;
}

bool isTrend(cov_model *cov) {
  Types type = CovList[cov->nr].Type;
 if (type != UndefinedType)  return type == TrendType;


 assert(false);
  assert(CovList[cov->nr].TypeFct != NULL);
  return CovList[cov->nr].TypeFct(TrendType, cov);
}

bool isInterface(Types type) {
  return type == InterfaceType;
}

bool isInterface(cov_model *cov) {
  Types type = CovList[cov->nr].Type;
  if (type != UndefinedType) return type == InterfaceType;

  //PMI(cov, "isinterface");
  assert(CovList[cov->nr].TypeFct != NULL);
  return CovList[cov->nr].TypeFct(InterfaceType, cov);
}

bool isUndefinedType(Types type) {
  return type == UndefinedType;
}

bool isOtherType(Types type) {
  return type == OtherType;
}

bool isOtherType(cov_model *cov) {
  Types type = CovList[cov->nr].Type;
  if (type != UndefinedType) return type == OtherType;

  assert(CovList[cov->nr].TypeFct != NULL);
  return CovList[cov->nr].TypeFct(InterfaceType, cov);
}

////////////

bool isGaussProcess(cov_model *cov) {
  int nr = cov->nr;
  return nr == GAUSSPROC || isGaussMethod(cov);
}

bool isBernoulliProcess(cov_model *cov) {
  int nr = cov->nr;
  return nr == BINARYPROC;
}

bool isGaussBasedProcess(cov_model *cov) {
  int nr = cov->nr;
  return isGaussProcess(cov) || nr == CHI2PROC;
}

bool isBrownResnickProcess(cov_model *cov) {
  Types type = CovList[cov->nr].Type;
  int nr = cov->nr;
  return type == BrMethodType || nr == BROWNRESNICKPROC;
}

bool isBRuserProcess(cov_model *cov) {
  Types type = CovList[cov->nr].Type;
  return type == BrMethodType;
 }

bool isBRuserProcess(Types type) {
  return type == BrMethodType;
}

bool isMaxStable(cov_model *cov) {
  int nr = cov->nr;
  return isBrownResnickProcess(cov) || nr == SMITHPROC || nr == SCHLATHERPROC;
}


bool isCov(cov_model *cov) {
  int nr = cov->nr;
  return nr == COVFCTN || nr == COVMATRIX ;
}


bool hasNoRole(cov_model *cov) {
  int role = cov->role;
  return role == ROLE_BASE;
}

bool hasMaxStableRole(cov_model *cov) {
  int role = cov->role;
  return role == ROLE_MAXSTABLE || role == ROLE_BROWNRESNICK ||
    role == ROLE_SMITH || role == ROLE_SCHLATHER;
}

bool hasPoissonRole(cov_model *cov) {
  int role = cov->role;
  return role == ROLE_POISSON_GAUSS || role == ROLE_POISSON;
}


int role_of_process(int nr) {
  return 
    (nr == AVERAGE_USER || nr == AVERAGE_INTERN ||
     nr == RANDOMCOIN_USER) ? ROLE_POISSON
    : ((nr >= CIRCEMBED &&  nr <= CE_INTRINPROC_INTERN)
       || nr == DIRECT || nr == NUGGET || nr == NUGGET_INTERN
       || nr == SEQUENTIAL 
       || nr == SPECTRAL_PROC_USER || nr == SPECTRAL_PROC_INTERN
       || nr == TBM_PROC_USER || nr == TBM_PROC_INTERN 
       || nr == GAUSSPROC
       ) ? ROLE_GAUSS
    : nr == HYPERPLANE_USER || nr == HYPERPLANE_INTERN ? ROLE_GAUSS
    : nr == SPECIFIC ? ROLE_GAUSS
    : ( nr == BRSHIFTED_USER || nr == BRMIXED_USER || nr == BRORIGINAL_USER
	|| nr == BROWNRESNICKPROC) ? ROLE_BROWNRESNICK 
    : nr == BINARYPROC ? ROLE_BERNOULLI
    : nr == POISSONPROC ? ROLE_POISSON
    : nr == SCHLATHERPROC ? ROLE_SCHLATHER
    : nr == SMITHPROC ? ROLE_SMITH
    : ROLE_FAILED;
}


bool isMonotone(cov_model *cov) {
  int monotone = cov->monotone;
  return monotone >= MONOTONE && monotone <= NORMAL_MIXTURE ;
} 

bool isMonotone(int monotone) {
  return monotone >= MONOTONE && monotone <= NORMAL_MIXTURE;
}

bool isNormalMixture(cov_model *cov) {
  int monotone = cov->monotone;
  return monotone == NORMAL_MIXTURE || monotone == COMPLETELY_MON;
}


bool isNormalMixture(int monotone) {
  return monotone == NORMAL_MIXTURE || monotone == COMPLETELY_MON;
}


bool isBernstein(cov_model *cov) {
  int monotone = cov->monotone;
  return monotone == BERNSTEIN;
}

bool isBernstein(int monotone) {
  return monotone == BERNSTEIN;
}

bool isGneiting(cov_model *cov) {
  int monotone = cov->monotone;
  return monotone == GNEITING_MON || monotone == COMPLETELY_MON;
}

bool isGneiting(int monotone) {
  return monotone == GNEITING_MON || monotone == COMPLETELY_MON;
}



