#ifndef RandomShape_H
#define RandomShape_H 1

//-----------------------------------------------------------------
// MPP operator


void kappamppplus(int i, cov_model *cov, int *nr, int *nc);
int checkmppplus(cov_model *cov);
void rangempplus(cov_model *cov, range_type *range);
int struct_mppplus(cov_model *plus, cov_model **newmodel);
int init_mppplus(cov_model *cov, storage *s);
void do_mppplus(cov_model *cov, storage *s);
void mppplus(double *x, cov_model *cov, double *v); 
//void logmppplus(double *x, cov_model *cov, double *v, double *sign); 


void pts_given_shape(double *x, cov_model *cov, double *v); 
void logpts_given_shape(double *x, cov_model *cov, double *v, double *sign); 
int check_pts_given_shape(cov_model *cov);
int struct_pts_given_shape(cov_model *cov, cov_model **newmodel);
int init_pts_given_shape(cov_model *cov, storage *S);  
void do_pts_given_shape(cov_model *cov, storage *S);
void range_pts_given_shape(cov_model *cov, range_type *range);


void standard_shape(double *x, cov_model *cov, double *v); 
void logstandard_shape(double *x, cov_model *cov, double *v, double *sign); 
int check_standard_shape(cov_model *cov);
int struct_standard_shape(cov_model *cov, cov_model **newmodel);
int init_standard_shape(cov_model *cov, storage *S);  
void do_standard_shape(cov_model *cov, storage *S);

void stationary_shape(double *x, cov_model *cov, double *v); 
void logstationary_shape(double *x, cov_model *cov, double *v, double *sign); 
int check_stationary_shape(cov_model *cov);
int struct_stationary_shape(cov_model *cov, cov_model **newmodel);
int init_stationary_shape(cov_model *cov, storage *S);  
void do_stationary_shape(cov_model *cov, storage *S);


//-----------------------------------------------------------------
// Trends
void mixed(double *x,  cov_model *cov, double *v);
void MLEmixed(double *x,  cov_model *cov, double *v);
void mixed_nonstat(double *x,  double *y, cov_model *cov, double *v);
void MLEmixed_nonstat(double *x,  double *y, cov_model *cov, double *v);
void kappamixed(int i, cov_model *cov, int *nr, int *nc);
int  checkmixed(cov_model *cov);
void rangemixed(cov_model *cov, range_type* ra);
int  initmixed(cov_model *cov, storage *s);
void domixed(cov_model *cov, storage *s);
void covmatrix_mixed(cov_model *cov, double *v);
char iscovmatrix_mixed(cov_model *cov);

void trend(double *x, cov_model *cov, double *v);
void trend_nonstat(double *x, double *y, cov_model *cov, double *v);
int checktrend(cov_model *cov);
void rangetrend(cov_model *cov, range_type* ra);
int init_trend(cov_model *cov, storage *s);
void do_trend(cov_model *cov, storage *s);
double GetInternalMean(cov_model *cov);
void kappatrend(int i, cov_model *cov, int *nr, int *nc);


//-----------------------------------------------------------------
// shapes
void kappa_ave(int i, cov_model *cov, int *nr, int *nc);
void ave(double *x, cov_model *cov, double *v);
int checkave(cov_model *cov);
void rangeave(cov_model *cov, range_type* ra);
//void sd_standard(mpp_storage *s, cov_model *cov);
int structAve(cov_model *cov, cov_model **newmodel);

int check_shapeave(cov_model *cov);
int init_shapeave(cov_model *cov, storage *s);
void do_shapeave(cov_model *cov, storage *s);
void logshapeave(double *x, cov_model *cov, double *v, double *sign);

void ball(double *x, cov_model *cov, double *v);
int init_ball(cov_model *cov, storage *s);
void do_ball(cov_model *cov, storage *s); 
int struct_ball(cov_model *cov, cov_model **newmodel);
void Inverseball(double *x, cov_model *cov, double *v);

void kappa_EAxxA(int i, cov_model *cov, int *nr, int *nc);
void EAxxA(double *x, cov_model *cov, double *v);
int checkEAxxA(cov_model *cov);
void rangeEAxxA(cov_model *cov, range_type* ra);
void minmaxEigenEAxxA(cov_model *cov, double *mm);

void kappa_EtAxxA(int i, cov_model *cov, int *nr, int *nc);
void EtAxxA(double *x, cov_model *cov, double *v);
int checkEtAxxA(cov_model *cov);
void rangeEtAxxA(cov_model *cov, range_type* ra);
void minmaxEigenEtAxxA(cov_model *cov, double *mm);

void Polygon(double *x, cov_model *cov, double *v);
int check_polygon(cov_model *cov);
void range_polygon(cov_model *cov, range_type *range);
int init_polygon(cov_model *cov, storage *s);
void do_polygon(cov_model *cov, storage *s); 
int struct_polygon(cov_model *cov, cov_model **newmodel);
void Inversepolygon(double *x, cov_model *cov, double *v);
void InversepolygonNonstat(double *v, cov_model *cov, double *x, double *y);

void kappa_rotat(int i, cov_model *cov, int *nr, int *nc);
void rotat(double *x, cov_model *cov, double *v);
int checkrotat(cov_model *cov);
void rangerotat(cov_model *cov, range_type* ra);
void minmaxEigenrotat(cov_model *cov, double *mm);

void kappa_Rotat(int i, cov_model *cov, int *nr, int *nc);
void Rotat(double *x, cov_model *cov, double *v);
int checkRotat(cov_model *cov);
void rangeRotat(cov_model *cov, range_type* ra);


void kappa_rational(int i, cov_model *cov, int *nr, int *nc);
void minmaxEigenrational(cov_model *cov, double *mm);
void rational(double *x, cov_model *cov, double *v);
int checkrational(cov_model *cov);
void rangerational(cov_model *cov, range_type* ra);

void kappa_stp(int i, cov_model *cov, int *nr, int *nc);
void stp(double *x,  double *y, cov_model *cov, double *v);
int checkstp(cov_model *cov);
void rangestp(cov_model *cov, range_type* ra);
int structStp(cov_model *cov, cov_model **newmodel);

int check_shapestp(cov_model *cov);
int init_shapestp(cov_model *cov, storage *s);
void do_shapestp(cov_model *cov, storage *s);
void logshapestp(double *x, double *u, cov_model *cov, double *v, double *);

void truncsupport(double *x, cov_model *cov, double *v);
int checktruncsupport(cov_model *cov);
void truncsupportInverse(double *x, cov_model *cov, double *v);
void rangetruncsupport(cov_model *cov, range_type *range);
int struct_truncsupport(cov_model *cov, cov_model **newmodel);
int init_truncsupport(cov_model *cov, storage *s);
void do_truncsupport(cov_model *cov, storage *s);

void randomsign(double *x, cov_model *cov, double *v);
void lograndomsign(double *x, cov_model *cov, double *v, double *sign); 
void randomsignInverse(double *x, cov_model *cov, double *v);
void randomsignNonstatInverse(double *v, cov_model *cov, double *x, double *y);
int check_randomsign(cov_model *cov);
void range_randomsign(cov_model *cov, range_type *range);
int init_randomsign(cov_model *cov, storage *s);
int struct_randomsign(cov_model *cov, cov_model **newmodel);
void do_randomsign(cov_model *cov, storage *s);





//-----------------------------------------------------------------
// Gauss Methods
void kappa_ce(int i, cov_model *cov, int *nr, int *nc);
int check_ce(cov_model *cov); 
void range_ce(cov_model *cov, range_type *range);  
int init_circ_embed(cov_model *cov, storage *s);
void do_circ_embed(cov_model *cov, storage *s);

int struct_ce_approx(cov_model *cov, cov_model **newmodel);
int init_ce_approx(cov_model *cov, storage *S);
void do_ce_approx(cov_model *cov, storage *S);

void kappa_localproc(int i, cov_model *cov, int *nr, int *nc);

int check_co_proc(cov_model *cov);
int check_co_intern(cov_model *cov);
void range_co_proc(cov_model *cov, range_type* ra);

void distrD(double *x, cov_model *cov, double *v);
void distrP(double *x, cov_model *cov, double *v);
void distrQ(double *x, cov_model *cov, double *v);
void distrR(double *x, cov_model *cov, double *v);
void kappa_distr(int i, cov_model *cov, int *nr, int *nc);
int check_distr(cov_model *cov);
int init_distr(cov_model *cov, storage *s);
void do_distr(cov_model *cov, storage *s);
void range_distr(cov_model *cov, range_type *range);
void range_intrinCE(cov_model *cov, range_type* ra);
int check_local_proc(cov_model *cov);

int check_directGauss(cov_model *cov);
void range_direct(cov_model *cov, range_type *range);
int init_directGauss(cov_model *cov, storage *s);
void do_directGauss(cov_model *cov, storage *s);


int check_hyperplane(cov_model *cov);
int check_hyperplane_intern(cov_model *cov);
void range_hyperplane(cov_model *cov, range_type *range);
int struct_hyperplane(cov_model *cov, cov_model **newmodel);
int init_hyperplane(cov_model *cov, storage *s);
void do_hyperplane(cov_model *cov, storage *s);


int check_nugget_proc(cov_model *cov);
void range_nugget_proc(cov_model *cov, range_type *range);
int init_nugget(cov_model *cov, storage *s);
void do_nugget(cov_model *cov, storage *s );
int struct_nugget(cov_model *cov, cov_model **newmodel);


int check_randomcoin(cov_model *cov);
void range_randomcoin(cov_model *cov, range_type *range);
//void coin(double *x, cov_model *cov, double *v);
int init_randomcoin(cov_model *cov, storage *s);
int struct_randomcoin(cov_model *cov, cov_model **newmodel);
//void coinInverse(double *x, cov_model *cov, double *v);
int struct_randomcoin_extract(cov_model *cov, cov_model **newmodel);



int check_sequential(cov_model *cov) ;
void range_sequential(cov_model *cov, range_type *range) ;
int init_sequential(cov_model *cov, storage *s);
void do_sequential(cov_model *cov, storage *s);


int check_spectral(cov_model *cov) ;
void range_spectral(cov_model *cov, range_type *range) ;
int init_spectral(cov_model *cov, storage *s);
void do_spectral(cov_model *cov, storage *s );
int struct_spectral(cov_model *plus, cov_model **newmodel);


int check_specificGauss(cov_model *cov);
int struct_specificGauss(cov_model *cov, cov_model **newmodel);
int init_specificGauss(cov_model *cov, storage *S);
void do_specificGauss(cov_model *cov, storage *S);



void tbm_kappasproc(int i, cov_model *cov, int *nr, int *nc);
int checktbmproc(cov_model *cov);
void rangetbmproc(cov_model *cov, range_type* ra);
int struct_tbmproc(cov_model *cov, cov_model **newmodel);
int init_tbmproc(cov_model *cov, storage *s);
void do_tbmproc(cov_model *cov, storage *s);
//void tbm2num(double *x, cov_model *cov, double *v);


//-----------------------------------------------------------------
// MPP Methods

int init_BRorig(cov_model *cov, storage *s);
void do_BRorig(cov_model *cov, storage *s);

int init_BRshifted(cov_model *cov, storage *s);
void do_BRshifted(cov_model *cov, storage *s);

int check_BRmixed(cov_model *cov);
int init_BRmixed(cov_model *cov, storage *s);
void do_BRmixed(cov_model *cov, storage *s);
void range_BRmixed(cov_model *cov, range_type *range);
void kappaBRmixed(int i, cov_model *cov, int *nr, int *nc);
 
int initBRuser (cov_model *cov, storage *S);
int structBRuser(cov_model *cov, cov_model **newmodel);
int structBRintern(cov_model *cov, cov_model **newmodel);



//-----------------------------------------------------------------
// Processes

// checkSproc by checkS
int structSproc(cov_model *cov, cov_model **newmodel);
int initSproc(cov_model *cov, storage *s);
void doSproc(cov_model *cov, storage *s);

int checkplusproc(cov_model *cov);
 int structplusproc(cov_model *cov, cov_model **newmodel);
int initplusproc(cov_model *cov, storage *s);
void doplusproc(cov_model *cov, storage VARIABLE_IS_NOT_USED *s) ;
 
void binary(double *x, cov_model *cov, double *v);
int checkbinaryprocess(cov_model *cov);
void rangebinaryprocess(cov_model *cov, range_type *range);
int struct_binaryprocess(cov_model *cov, cov_model **newmodel);
int init_binaryprocess(cov_model *cov, storage *s);
void do_binaryprocess(cov_model *cov, storage *s);

void BrownResnick(double *x, cov_model *cov, double *v);
int checkBrownResnick(cov_model *cov);
int structBrownResnick(cov_model *cov, cov_model **newmodel);
int initBrownResnick(cov_model *cov, storage *s);
void doBrownResnick(cov_model *cov, storage *s);

int check_common_gauss(cov_model *cov);
void range_common_gauss(cov_model *cov, range_type *range);
int checkgaussprocess(cov_model *cov);
void rangegaussprocess(cov_model *cov, range_type *range);
int struct_gaussprocess(cov_model *cov, cov_model **newmodel);
int init_gaussprocess(cov_model *cov, storage *s);
void do_gaussprocess(cov_model *cov, storage *s);

int struct_extractdollar(cov_model *cov, cov_model **newmodel);

int check_poisson(cov_model *cov);
int struct_poisson(cov_model *cov, cov_model **newmodel);
void range_poisson(cov_model *cov, range_type *range);
int init_poisson(cov_model *cov, storage *S);

void extremalgaussian(double *x, cov_model *cov, double *v);
int check_schlather(cov_model *cov);
int struct_schlather(cov_model *cov, cov_model **newmodel);

int struct_smith(cov_model *cov, cov_model **newmodel);
int check_smith(cov_model *cov);


int checkchisqprocess(cov_model *cov) ;
void rangechisqprocess(cov_model *cov, range_type *range) ;
int struct_chisqprocess(cov_model *cov, cov_model **newmodel) ;
int init_chisqprocess(cov_model *cov, storage *s) ;
void do_chisqprocess(cov_model *cov, storage *s);

 
void dompp(cov_model *cov, storage *S); 
int init_mpp(cov_model *cov, storage *S);
void range_mpp(cov_model *cov, range_type *range);
int SetGEVetc(cov_model *cov, int type);


//-----------------------------------------------------------------
// Interfaces
void simulate(double *x, cov_model *cov, double *v);
int check_simulate(cov_model *cov); 
void range_simulate(cov_model *cov, range_type *range);  
int struct_simulate(cov_model *cov, cov_model **newmodel);
int init_simulate(cov_model *cov, storage *S);
// void do_simulate(cov_model *cov, storage *S);

void Cov(double *x, cov_model *cov, double *value) ;
int check_cov(cov_model *cov) ;
int struct_cov(cov_model *cov, cov_model **newmodel);

int check_fctn(cov_model *cov);

void CovMatrix(double *x, cov_model *cov, double *value) ;
int check_covmatrix(cov_model *cov) ;

void EvalDistr(double *x, cov_model *cov, double *v);
void kappa_EvalDistr(int i, cov_model *cov, int *nr, int *nc);
int check_EvalDistr(cov_model *cov); 
void range_EvalDistr(cov_model *cov, range_type *range);  
int struct_EvalDistr(cov_model *cov, cov_model **newmodel);
int init_EvalDistr(cov_model *cov, storage *S);
// void do_EvalDistr(cov_model *cov, storage *S);

void RFget(double *x, cov_model *cov, double *v);
int SearchParam(cov_model *cov, get_storage *s) ;
int check_RFget(cov_model *cov) ;
void range_RFget(cov_model *cov, range_type* range);
int struct_RFget(cov_model *cov, cov_model **newmodel);


void Variogram(double *x, cov_model *cov, double *value) ;
void Pseudovariogram(double *x, cov_model *cov, double *value) ;
int check_vario(cov_model *cov);
int struct_variogram(cov_model *cov, cov_model **newmodel);

void Dummy(double *x, cov_model *cov, double *value);
int check_dummy(cov_model *cov);
int struct_dummy(cov_model *cov, cov_model **newmodel);


//-----------------------------------------------------------------
// unsorted


#endif /* RandomShape_H */
 
