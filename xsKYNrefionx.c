/* KYNrefionx - relativistic disc reflection in lamp-post geometry (ionised disc)
 *              model subroutine for XSPEC
 * 
 * ref. Dovciak M., Karas V., Yaqoob T. (2004)
 * -----------------------------------------------------------------------------
 * OTHER REFERENCES:
 * 
 * Dovciak, M., Svoboda, J., Goosmann, R. W., et al.: 2014, in Proceedings
 * of RAGtime 14-16: Workshops on black holes and neutron stars (Silesian 
 * University in Opava). [arXiv:1412.8627]
 * 
 * Dovciak M., Karas V. & Yaqoob, T. (2004). An extended scheme for fitting 
 * X-ray data with accretion disk spectra in the strong gravity regime. 
 * ApJS, 153, 205.
 * 
 * Dovciak M., Karas V., Martocchia A., Matt G. & Yaqoob T. (2004). XSPEC model
 * to explore spectral features from black hole sources. In Proc. of the 
 * workshop on processes in the vicinity of black holes and neutron stars. 
 * S.Hledik & Z.Stuchlik, Opava. In press. [astro-ph/0407330]
 * 
 * Dovciak M. (2004). Radiation of accretion discs in strong gravity. Faculty of
 * Mathematics and Physics, Charles University, Prague. PhD thesis.
 * [astro-ph/0411605]
 * -----------------------------------------------------------------------------
 * 
 * This subroutine computes the emission from an acrretion disc that is
 * illuminated from the primary power-law source located on the axis above the
 * central black hole. All relativistic effects are taken into account (in all 
 * three parts of the light path - from the primary source to the observer, to 
 * disc and from the disc to the observer). This model calls subroutine ide() 
 * for integrating local emission over the disc and uses the FITS file 
 * 'KBHtablesNN.fits' defining the transfer functions needed for integration 
 * over disc as well as the FITS file 'KBHlamp_q.fits' defining the transfer 
 * functions between the source and the disc. For details on ide() and the FITS 
 * file 'KBHtablesNN.fits' see the subroutine ide() in xside.c, for details on 
 * the FITS file 'KBHlamp_q.fits' see the subroutine KYNrlpli() in xsKYNrlpli.c.
 * The reflection is taken from Ross & Fabian tables reflionx.mod.
 *
 * par1  ... a/M     - black hole angular momentum (-1 <= a/M <= 1)
 * par2  ... theta_o - observer inclination in degrees (0-pole, 90-disc)
 * par3  ... rin - inner edge of non-zero disc emissivity (in GM/c^2 or in 
 *                 r_mso)
 * par4  ... ms  - switch that defines the meaning/units of rin, rout
 *                 0: we integrate from inner edge = par3 
 *                 1: if the inner edge of the disc is below marginally stable
 *                    orbit then we integrate emission above MSO only
 *                 2: we integrate from inner edge given in units of MSO, i.e.
 *                    inner edge = par3 * r_mso (the same applies for outer 
 *                    edge)
 * par5  ... rout  - outer edge of non-zero disc emissivity (in GM/c^2 or in 
 *                   r_mso)
 * par6  ... phi   - lower azimuth of non-zero disc emissivity (deg)
 * par7  ... dphi  - (phi + dphi) is upper azimuth of non-zero disc emissivity
 *                   0 <= dphi <= 360  (deg)
 * par8  ... M/M8   - black hole mass in units of 10^8 solar masses
 * par9  ... height - height on the axis (measured from the center) at which
 *                    the primary source is located (GM/c^2)
 *                    DO NOT USE THE MODEL WITH NEGATIVE HEIGHT!!! 
 * par10 ... PhoIndex - power-law energy index of the primary flux
 * par11 ... L/Ledd   - dE/dt, the intrinsic local (if negative) or the 
 *                      observed (if positive) primary isotropic flux in the 
 *                      X-ray energy range 2-10keV in units of Ledd
 * par12 ... Np:Nr  - ratio of the primary to the reflected normalization
 *                    1 - self-consistent model for isotropic primary source
 *                    0 - only reflection, primary source is hidden
 *                  - if positive then L/Ledd (par11) means the luminosity 
 *                    towards the observer
 *                  - if negative then L/Ledd (par11) means the luminosity 
 *                    towards the disc
 * par13 ... density  - density profile normalization in 10^15 cm^(-3) 
 *                      if positive
 *                    - ionisation profile normalisation if it is negative
 * par14 ... den_prof - radial power-law density profile if par13 is positive
 *                    - radial ionisation profile if par13 is negative
 * par15 ... abun   - Fe abundance (in solar abundance)
 * par16 ... alpha  - position of the cloud centre in GM/c^2 in alpha coordinate
 *                    (alpha being the impact parameter in phi direction, 
 *                     positive for approaching side of the disc)
 * par17 ... beta   - position of the cloud centre in GM/c^2 in beta coordinate
 *                    (beta being the impact parameter in theta direction, 
 *                     positive in up direction, i.e. above the disc)
 * par18 ... rcloud - radius of the obscuring cloud (in GM/c^2)
 *                  - if negative, only the emission transmitted through
 *                    the cloud is taken into account
 * par19 ... zshift - overall Doppler shift
 * par20 ... limb   - limb darkening/brightening law (emission directionality)
 *                  - if = 0 the local emisivity is not multiplied by anything
 *                  - if = 1 the local emisivity is multiplied by 1+2.06*mu
 *                    (limb darkening)
 *                  - if = 2 the local emisivity is multiplied by ln(1+1/mu)
 *                    (limb brightening)
 * par21 ... tab - which reflion table to use 
 *                 1 -> reflion (the old one, lower cut-off energy at 1eV,
 *                      not good for PhoIndex > 2)
 *                 2 -> reflionx (the new one, lower cut-off energy at 100eV)
 * par22 ... sw  - swich for the way how to compute the refl. spectra
 *                 1 -> use the computed Xi for the interpolation in reflion,
 *                      i.e. use proper total incident intensity
 *                      with the shifted cut-offs
 *                 2 -> use the Xi correspondent to the computed normalization
 *                      of the incident flux, i.e. do not shift the cut-offs
 *                      when computing the total incident intensity
 * par23 ... ntable - table of relativistic transfer functions used in the model
 *                    (defines fits file with tables), 0<= ntable <= 99
 * par24 ... nrad   - number of grid points in radius
 * par25 ... division - type of division in r integration
 *                      0 -> equidistant radial grid (constant linear step)
 *                      1 -> exponential radial grid (constant logarithmic step)
 * par26 ... nphi   - number of grid points in azimuth
 * par27 ... smooth - whether to smooth the resulting spectrum (0-no, 1-yes)
 * par28 ... nthreads - number of threads to be used for computations
 * par29 ... norm     - has to be set to unity!
 *
 * NOTES:
 *  -> accuracy vs. speed trade off depends mainly on: nrad, nphi
 *
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "fitsio.h"

/*******************************************************************************
*******************************************************************************/
#ifdef OUTSIDE_XSPEC

#define IFL    1
#define NPARAM 29
#define NE     400
#define E_MIN  0.1
#define E_MAX  80.

int main() {
  
void KYNrefionx(const double *ear, int ne, const double *param, int ifl, 
                double *photar, double *photer, const char* init);

double ear[NE + 1], param[NPARAM], photar[NE], photer[NE];
char   initstr[0] = "";
int    ie;

param[ 0] = 1.;       // a/M
param[ 1] = 30.;      // theta_o
param[ 2] = 1.;       // rin
param[ 3] = 1.;       // ms
param[ 4] = 400.;     // rout
param[ 5] = 0.;       // phi0
param[ 6] = 360.;     // dphi
param[ 7] = 1.;       // M/M8
param[ 8] = 3.;       // height
param[ 9] = 2.;       // PhoIndex
param[10] = 0.001;    // L/Ledd
param[11] = 1.;       // Np:Nr
param[12] = 1.;       // density
param[13] = 0.;       // den_prof
param[14] = 1.;       // abun
param[15] = -6.;      // alpha
param[16] = 0.;       // beta
param[17] = 0.;       // rcloud
param[18] = 0.;       // zshift
param[19] = 0.;       // limb
param[20] = 2.;       // tab
param[21] = 1.;       // sw
param[22] = 80.;      // ntable
param[23] = 300.;     // nrad
param[24] = 1.;       // division
param[25] = 360.;     // nphi
param[26] = 1.;       // smooth
param[27] = 4.;       // nthreads
param[28] = 1.;       // norm

for (ie = 0; ie <= NE; ie++) {
//  ear[ie] = E_MIN + ie * (E_MAX-E_MIN) / NE;
  ear[ie] = E_MIN * pow(E_MAX / E_MIN, ((double) ie) / NE);
}

KYNrefionx(ear, NE, param, IFL, photar, photer, initstr);

return(0);
}

#endif
/*******************************************************************************
*******************************************************************************/

#define MPC_2 1.05e-49
#define ERG 6.241509e8
#define RG2 2.1819222e26
// Ledd is in erg (not W) and multiplied by 10^8 due to (M / (10^8*Msun)) scale
#define LEDD 1.26e46
#define HUBBLE 70.
#define CLIGHT 299792.458
#define LOGXI_NORM0 4.761609554
#define PI 3.14159265359
#define PI2 6.2831853071795865
#define LAMP "KBHlamp_q.fits"
#define REFLION1 "reflion.mod"
#define REFLION2 "reflionx.mod"
#define REFLIONX_NORM 1.e20
#define EC 300.

/* Let's declare variables that are common for the main and emissivity 
   subroutines */
static float  *radius, *xi, *logxi;
static double *gfac, *transf_d, *energy1, *flux1;
static double h, gamma0, nH0, qn, mass2, am2, r_plus, Np, logxi_norm;
static long   nrad, nxi;
static int    sw, limb;

extern char*  FGMODF(void);
extern char*  FGMSTR(char* dname);
extern void   FPMSTR(const char* value1, const char* value2);
extern int    xs_write(char* wrtstr, int idest);
extern double incgamma(double a, double x);
extern void   cutoffpl(double *ear, int ne, double *param, double *photar);

void KYNrefionx(const double *ear, int ne, const double *param, int ifl, 
                double *photar, double *photer, const char* init) {
  
extern int ide(const double *ear, const int ne, const int nt, double *far, 
               double *qar, double *uar, double *var, 
               const double *ide_param, void (*emissivity)(), 
               const int ne_loc);

void emis_KYNrefionx(double** ear_loc, const int ne_loc, const int nt, 
                     double *far_loc, double *qar_loc, double *uar_loc, 
                     double *var_loc, const double r, const double phi, 
                     const double cosmu, const double phiphoton, 
                     const double alpha_o, const double beta_o, 
                     const double delay, const double g);

/* Let's declare static variables whose values should be remembered for next
   run in XSPEC */
static char   kydir[255] = "";
static char   pname[128] = "KYDIR", pkyLxLamp[128] = "KYLxLamp";
static char   pkyxiin[128] = "KYXIin", pkyxiout[128] = "KYXIout";
static char   pkyRefl[128] = "KYRefl";
static long   nrh, nh, nincl, nener;
static float  *r_horizon, *height, *incl, *dWadWo, *q, *pr, *dWadSd, *abun,
              *gam, *emission, *energy0;
static long   ne_loc, nabun, ngam;
static double *energy2, *flux0, transf_o;
static double h_rh_old = -1., gam_old = -1., abun_old = -1., thetaO_old = -1.,
              am_old = -1.;
static int    rix_old = -1, first_h = 1, first_rix = 1;

FILE   *fw;
char   errortxt[80];
char   reflion[255], text[255];
char   kyxiin[32], kyxiout[32], kyLxLamp[32], kyRefl[32];
double ide_param[25];
double far[ne], qar[ne], uar[ne], var[ne], pd[ne], pa[ne], pa2[ne];
double ttmp, ttmp1, utmp, utmp1, vtmp, vtmp1, y1, y2, y3, y4, y5, y6, y7, y8;
double pr_final, pom, pom1, pom2, pom3;
double r, r2, delta, ULt, rms, tmp1, Ut, U_r, UrU_r, Lms, Ems;
//double q_final, U_phi, Ur;
double am, thetaO, rin, rout, h_rh, elow, ehigh;
double mass, Anorm, Dnorm, g_L, E0, Lx, flux_prim, flux_refl, refl_ratio;
double zzshift;
double pamin, pamax, pa2min, pa2max, NpNr;
double abundance, lensing, gfactor0, ionisation;
long   iabun, igam, ixi;
int    imin, imax, irh0, ih0, ith0, ir0, iabun0, igam0;
int    nspec, polar, stokes, rix;
int    i, ie, je, quit;
// the following are needed for cut-off power-law taken from XSPEC
double ear1[ne + 1], param1[2];
double photar1[ne]; 

// these are needed to work with a fits file...
fitsfile *fptr;
char     tables_file[255];
int      hdutype = 2;
int      colnum = 1;
long     frow = 1, felem = 1, nelems, nrow;
float    float_nulval = 0.;
int      nelements, nelements1, nelements2;
int      ihorizon, irow, anynul, status = 0;

// Let's initialize parameters for subroutine ide()
for (ie = 0; ie < ne; ie++) far[ie] = 0;
// am - black hole angular momentum
ide_param[0] = param[0];
am = ide_param[0];
am2 = am * am;
pom1 = pow(1. + am, 1. / 3.);
pom2 = pow(1. - am, 1. / 3.);
pom3 = pow(1. - am2, 1. / 3.);
pom = 1. + pom3 * (pom1 + pom2);
pom1 = sqrt(3. * am2 + pom * pom);
if (am >= 0) rms = 3. + pom1 - sqrt((3. - pom) * (3. + pom + 2. * pom1));
else rms = 3. + pom1 + sqrt((3. - pom) * (3. + pom + 2. * pom1));
r_plus = 1. + sqrt(1. - am2);
// thetaO - observer inclination
ide_param[1] = param[1];
thetaO = ide_param[1];
// rin - inner edge of non-zero disc emissivity
ide_param[2] = param[2];
if (param[3] == 1.)
  if (param[2] < rms) rin = rms;
  else rin = param[2];
else
  if (param[2] < r_plus) rin = r_plus;
  else rin = param[2];
// ms - whether to integrate from rin or rms
ide_param[3] = param[3];
// rout - outer edge of non-zero disc emissivity
ide_param[4] = param[4];
rout = param[4];
// phi  - lower azimuth of non-zero disc emissivity (deg)
ide_param[5] = param[5];
// dphi - (phi+dphi) is upper azimuth of non-zero disc emissivity (deg)
ide_param[6] = param[6];
// nrad - number of grid points in radius
ide_param[7] = param[23];
// division - type of division in r integration (0-equidistant, 1-exponential)
ide_param[8] = param[24];
// nphi - number of grid points in azimuth
ide_param[9] = param[25];
// smooth - whether to smooth the resulting spectrum (0-no, 1-yes)
ide_param[10] = param[26];
// normal - how to normalize the final spectrum
ide_param[11] = -1.;
// ntable - table model (defines fits file with tables)
ide_param[13] = param[22];
// M/M8 - black hole mass in units of 10^8 solar masses
mass = param[7];
mass2 = mass * mass;
logxi_norm = LOGXI_NORM0 + log10(mass);
// height - height on the axis (measured from the center) at which the primary
//          source is located (GM/c^2)
h = param[8];
if (h >= 0.)
  if (h < r_plus) h_rh = 0.;
  else h_rh = h - r_plus;
else {
  xs_write("kynrefionx: height has to be positive.", 5);
  for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
  return;
}
// PhoIndex - power-law energy index of the lamp emission
gamma0 = param[9];
// L/Ledd - dE/dt primary isotropic flux in Eddington luminosity
Np = param[10];
// Np:Nr - ratio of the primary to the reflected normalization
NpNr = param[11];
if( NpNr > 0. ) Np /= NpNr;
// nH0 - density/ionisation profile normalization in 10^15 cm^(-3)
nH0 = param[12];
if (nH0 == 0.) {
  xs_write("kynrefionx: density/ionisation must be non-zero!", 5);
  for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
  return;
}
// q_n - radial power-law density/ionisation profile
qn = param[13];
// Fe abundance (in solar abundance)
abundance = param[14];
// zshift - overall Doppler shift
if (param[18] > 0.) {
  ide_param[12] = param[18];
  Dnorm = pow(HUBBLE / CLIGHT / param[18], 2);
}else if (param[18] < 0.) {
  ide_param[12] = 0.;
  Dnorm = pow(-HUBBLE / CLIGHT / param[18], 2);
}else {
  ide_param[12] = 0.;
  Dnorm = 1.;
}
// zzshift - multiplication factor for gfac from zshift needed for primary
zzshift=1.0/(1.0+ide_param[12]);
// limb - table model (defines fits file with tables)
limb = (int) param[19];
if ((limb < 0) || (limb > 2)) {
  xs_write("kynrefionx: limb must be >= 0 and <= 2.", 5);
  for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
  return;
}
// tab - which reflion table to use 
rix = (int) param[20];
if (rix == 1) {
  E0 = 0.001;
  sprintf(reflion, REFLION1);
  nener = 500;
  elow = 0.01;
}
else{
  E0 = 0.1;
  sprintf(reflion, REFLION2);
  nener = 375;
  elow = 0.1;
}
ehigh = EC;
// sw  - swich for the way how to compute the refl. spectra
sw = (int) param[21];
// edivision - type of division in local energies (0-equidistant, 1-exponential)
ide_param[14] = 1.;
// periodic and dt are not needed for nt=1
// (ide_param[15], ide_param[16])
// polar - whether we need value of change in polarization angle (0-no,1-yes)
//stokes = (int) param[27];
stokes = 0;
if ((stokes < 0) || (stokes > 6)) {
  xs_write("kynrefionx: stokes has to be 0-6.", 5);
  for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
  return;
}
polar = 0;
if (stokes > 0) polar = 1;
ide_param[17] = polar;
//delay_r, delay_phi need not to be defined for non-timing computations
//ide_param[18]=-1.
//ide_param[19]=del_a
// number of threads for multithread computations
ide_param[20] = param[27];
// alpha - position of the cloud in alpha impact parameter (in GM/c^2)
ide_param[21] = param[15];
// beta - position of the cloud in beta impact parameter (in GM/c^2)
ide_param[22] = param[16];
// rcloud - radius of the cloud (in GM/c^2)
ide_param[23] = param[17];
//whether the flux defined in emissivity subroutine is local one (0) or the 
//observed one (1)
ide_param[24] = 0.;

// check if normalization parameter is equal to 1.
if ((param[18] != 0.) && (param[28] != 1.)) {
  xs_write("kynrefionx: the normalisation parameter par32 should be frozen to unity.", 5);
}
/******************************************************************************/
// Let's read the lamp post tables
if(first_h) {
// The status parameter must always be initialized.
  status = 0;
// Open the FITS file for readonly access
// - if set try KYDIR directory, otherwise look in the working directory
//   or in the xspec directory where tables are usually stored...
  sprintf(kydir, "%s", FGMSTR(pname));
  if (strlen(kydir) == 0) sprintf(tables_file, "./%s", LAMP);
  else if (kydir[strlen(kydir) - 1] == '/') sprintf(tables_file, "%s%s",
                                                    kydir, LAMP);
  else sprintf(tables_file, "%s/%s", kydir, LAMP);
// Let's read the 'KBHlamp_q' fits file
// The status parameter must always be initialized.
  status = 0;
  ffopen(&fptr, tables_file, READONLY, &status);
  if (status) {
    sprintf(tables_file, "%s%s", FGMODF(), LAMP);
    status = 0;
    ffopen(&fptr, tables_file, READONLY, &status);
  }
  if (status) {
    if (status) ffrprt(stderr, status);
    ffclos(fptr, &status);
    xs_write("\nkynrefionx: set the KYDIR to the directory with the KY tables.", 5);
    for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
    return;
  }
// Let's read tables (binary tables => hdutype=2)
// Move to the extension 'r_horizon' and read its values
  ffmrhd(fptr, 1, &hdutype, &status);
  ffgnrw(fptr, &nrh, &status);
//******************************************************************************
//  fprintf(stdout,"nrh = %ld\n",nrh);
//******************************************************************************   
// Allocate memory for r_horizon...
  if ((r_horizon = (float *) malloc(nrh * sizeof(float))) == NULL) {
    xs_write("kynrefionx: Failed to allocate memory for tmp arrays.", 5);
    for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
    return;
  }
// Read the data in the 'r_horizon' table
  nelems = nrh;
// FFGCV reads the VALUES from the first column.
  ffgcv(fptr, TFLOAT, colnum, frow, felem, nelems, &float_nulval, r_horizon,
        &anynul, &status);
//******************************************************************************
//  for ( i=0; i<nrh; i++)fprintf(stdout,"%f\n",r_horizon[i]);
//******************************************************************************   
// Move to the extension 'height' and read its values
  ffmrhd(fptr, 1, &hdutype, &status);
  ffgnrw(fptr, &nh, &status);
//******************************************************************************
//  fprintf(stdout,"nh = %ld\n",nh);
//******************************************************************************   
// Allocate memory for height...
  if ((height = (float *) malloc(nh * sizeof(float))) == NULL) {
    xs_write("kynrefionx: Failed to allocate memory for tmp arrays.", 5);
    for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
    return;
  }
// Read the data in the 'height' table
  nelems = nh;
// FFGCV reads the VALUES from the first column.
  ffgcv(fptr, TFLOAT, colnum, frow, felem, nelems, &float_nulval, height,
        &anynul, &status);
//******************************************************************************
//  for ( i=0; i<nh; i++)fprintf(stdout,"%f\n",height[i]);
//******************************************************************************   
// Move to the extension 'inclination' and read its values
  ffmrhd(fptr, 1, &hdutype, &status);
  ffgnrw(fptr, &nincl, &status);
//******************************************************************************
//  fprintf(stdout,"nincl = %ld\n",nincl);
//******************************************************************************   
// Allocate memory for inclination...
  if ((incl = (float *) malloc(nincl * sizeof(float))) == NULL) {
    xs_write("kynrefionx: Failed to allocate memory for tmp arrays.", 5);
    for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
    return;
  }
// Read the data in the 'inclination' table
  nelems = nincl;
// FFGCV reads the VALUES from the first column.
  ffgcv(fptr, TFLOAT, colnum, frow, felem, nelems, &float_nulval, incl,
        &anynul, &status);
//******************************************************************************
//  for ( i=0; i<nincl; i++)fprintf(stdout,"%f\n",incl[i]);
//******************************************************************************   
// Move to the extension 'r_rh' and read its values
  ffmrhd(fptr, 1, &hdutype, &status);
  ffgnrw(fptr, &nrad, &status);
//******************************************************************************
//  fprintf(stdout,"nrad = %ld\n",nrad);
//******************************************************************************   
// Allocate memory for r_rh...
  if ((radius = (float *) malloc(nrad * sizeof(float))) == NULL) {
    xs_write("kynrefionx: Failed to allocate memory for tmp arrays.", 5);
    for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
    return;
  }
// Read the data in the 'r_rh' table
  nelems = nrad;
// FFGCV reads the VALUES from the first column.
  ffgcv(fptr, TFLOAT, colnum, frow, felem, nelems, &float_nulval, radius,
        &anynul, &status);
//******************************************************************************
//  for ( i=0; i<nrad; i++)fprintf(stdout,"%f\n",radius[i]);
//******************************************************************************   
// Let's read the tables for dWadWo, q, p^r and dWadSd
// allocate memory for the arrays
  if ((dWadWo = (float *) malloc(nincl * nh * nrh * sizeof(float))) == NULL) {
    xs_write("kynrefionx: Failed to allocate memory for tmp arrays.", 5);
    for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
    return;
  }
  if ((q = (float *) malloc(nrad * nh * nrh * sizeof(float))) == NULL) {
    xs_write("kynrefionx: Failed to allocate memory for tmp arrays.", 5);
    for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
    return;
  }
  if ((pr = (float *) malloc(nrad * nh * nrh * sizeof(float))) == NULL) {
    xs_write("kynrefionx: Failed to allocate memory for tmp arrays.", 5);
    for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
    return;
  }
  if ((dWadSd = (float *) malloc(nrad * nh * nrh * sizeof(float))) == NULL) {
    xs_write("kynrefionx: Failed to allocate memory for tmp arrays.", 5);
    for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
    return;
  }
// read the tables
  for (ihorizon = 0; ihorizon < nrh; ihorizon++) {
    ffmrhd(fptr, 1, &hdutype, &status);
/*  to read the file only once we have to read in blocks (all columns
    from the extension are put to buffer together)
    let's find out how many rows are going to be read into the buffer */
    ffgrsz(fptr, &nrow, &status);
    nelements1 = nrow * nincl;
    nelements2 = nrow * nrad;
    for (irow = 0; irow < nh; irow += nrow) {
//    the last block to read may be smaller:
      if ((nrad - irow) < nrow) {
        nelements1 = (nrad - irow) * nincl;
        nelements2 = (nrad - irow) * nrad;
      }
      ffgcv(fptr, TFLOAT, 1, irow + 1, 1, nelements1, &float_nulval, 
            &dWadWo[irow * nincl + nh * nincl * ihorizon],
            &anynul, &status);
      ffgcv(fptr, TFLOAT, 2, irow + 1, 1, nelements2, &float_nulval, 
            &q[irow * nrad + nh * nrad * ihorizon],
            &anynul, &status);
      ffgcv(fptr, TFLOAT, 3, irow + 1, 1, nelements2, &float_nulval, 
            &pr[irow * nrad + nh * nrad * ihorizon],
            &anynul, &status);
      ffgcv(fptr, TFLOAT, 4, irow + 1, 1, nelements2, &float_nulval, 
            &dWadSd[irow * nrad + nh * nrad * ihorizon],
            &anynul, &status);
    }
  }
// The FITS file must always be closed before exiting the program.
  ffclos(fptr, &status);
/*******************************************************************************
  irh=20;
  ih=99;
  for ( i=0; i<nincl; i++ ) fprintf(stdout,"%d\t%f\t%f\n",i,incl[i],
    dWadWo[i+nincl*ih+nincl*nh*irh]);
  for ( i=0; i<nrad; i++) fprintf(stdout,"%d\t%f\t%f\t%f\t%f\n",i,radius[i],
    q[i+nrad*ih+nrad*nh*irh],pr[i+nrad*ih+nrad*nh*irh],
    dWadSd[i+nrad*ih+nrad*nh*irh]);
*******************************************************************************/
// Firstly we have to free allocated memory for the arrays gfac,
// cosin, phiph, transf_d
  if ((gfac = (double *) malloc(nrad * sizeof(double))) == NULL) {
    xs_write("kynrefionx: Failed to allocate memory for tmp arrays.", 5);
    for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
    return;
  }
/*if ((cosin = (double *) malloc(nrad * sizeof(double))) == NULL) {
    xs_write("kynrefionx: Failed to allocate memory for tmp arrays.", 5);
    for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
    return;
  }
  if ((phiph = (double *) malloc(nrad * sizeof(double))) == NULL) {
    xs_write("kynrefionx: Failed to allocate memory for tmp arrays.", 5);
    for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
    return;
  }*/
  if ((transf_d = (double *) malloc(nrad * sizeof(double))) == NULL) {
    xs_write("kynrefionx: Failed to allocate memory for tmp arrays.", 5);
    for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
    return;
  }
  first_h = 0;
}
/******************************************************************************/
if (h_rh > height[nh-1]) {
  sprintf(errortxt, "kynrefionx: the height must be lower than or equal to %f.",
          height[nh - 1] + r_plus);
  xs_write(errortxt, 5);
  for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
  return;
}
if (h_rh < height[0]) {
  sprintf(errortxt, "kynrefionx: the height is too low, we set it to %f.", 
          height[0] + r_plus);
  xs_write(errortxt, 5);
  h_rh = height[0];
  h = h_rh + r_plus;
}
// Let's interpolate the tables to desired spin and height
if ((am != am_old) || (h_rh != h_rh_old) || (thetaO != thetaO_old)) {
// given am->r_plus, find the corresponding index in r_horizon[]:
  imin = 0;
  imax = nrh;
  irh0 = nrh / 2;
  while ((imax - imin) > 1) {
    if (r_plus >= r_horizon[irh0 - 1]) imin = irh0;
    else imax = irh0;
    irh0 = (imin + imax) / 2;
  }
  if (irh0 == 0) irh0 = 1;
//if ((imax == nrh) && (r_plus > r_horizon[nrh - 1])) irh0 = nrh;
  ttmp = (r_plus - r_horizon[irh0 - 1]) / 
         (r_horizon[irh0] - r_horizon[irh0 - 1]);
  ttmp1 = 1. - ttmp;
// given h, find the corresponding index in height[]:
  imin = 0;
  imax = nh;
  ih0 = nh / 2;
  while ((imax - imin) > 1) {
    if (h_rh >= height[ih0 - 1]) imin = ih0;
    else imax = ih0;
    ih0 = (imin + imax) / 2;
  }
  if (ih0 == 0) ih0 = 1;
//if ((imax == nh) && (h_rh > height[nh - 1])) ih0 = nh;
  utmp = (h_rh - height[ih0 - 1]) / (height[ih0] - height[ih0 - 1]);
  utmp1 = 1. - utmp;
// given thetaO, find the corresponding index in incl[]:
  imin = 0;
  imax = nincl;
  ith0 = nincl / 2;
  while ((imax - imin) > 1) {
    if (thetaO >= incl[ith0 - 1]) imin = ith0;
    else imax = ith0;
    ith0 = (imin + imax) / 2;
  }
  if (ith0 == 0) ith0 = 1;
//if ((imax == nincl) && (thetaO > incl[nincl - 1])) ith0 = nincl;
  vtmp = (thetaO - incl[ith0 - 1]) / (incl[ith0] - incl[ith0 - 1]);
  vtmp1 = 1. - vtmp;
// transfer function from the axis to the observer
  y1 = dWadWo[ith0 - 1 + nincl * (ih0 - 1) + nincl * nh * (irh0 - 1)];
  y2 = dWadWo[ith0 - 1 + nincl * (ih0 - 1) + nincl * nh * irh0];
  y3 = dWadWo[ith0 - 1 + nincl * ih0 + nincl * nh * irh0];
  y4 = dWadWo[ith0 - 1 + nincl * ih0 + nincl * nh * (irh0 - 1)];
  y5 = dWadWo[ith0 + nincl * (ih0 - 1) + nincl * nh * (irh0 - 1)];
  y6 = dWadWo[ith0 + nincl * (ih0 - 1) + nincl * nh * irh0];
  y7 = dWadWo[ith0 + nincl * ih0 + nincl * nh * irh0];
  y8 = dWadWo[ith0 + nincl * ih0 + nincl * nh * (irh0 - 1)];
  transf_o = (vtmp1 * (utmp1 * (ttmp1 * y1 + ttmp * y2) + 
                       utmp *  (ttmp * y3 + ttmp1 * y4)) + 
              vtmp * (utmp1 * (ttmp1 * y5 + ttmp * y6) + 
                      utmp *  (ttmp * y7 + ttmp1 * y8)));
  if ((am != am_old) || (h_rh != h_rh_old)) {
    for (i = 0; i < nrad; i++) {
/*
// q from the axis to the disc
      y1 = q[i + nrad * (ih0 - 1) + nrad * nh * (irh0 - 1)];
      y2 = q[i + nrad * (ih0 - 1) + nrad * nh * irh0];
      y3 = q[i + nrad * ih0 + nrad * nh * irh0];
      y4 = q[i + nrad * ih0 + nrad * nh * (irh0 - 1)];
      q_final = utmp1 * (ttmp1 * y1 + ttmp * y2) + 
                utmp * (ttmp * y3 + ttmp1 * y4);
*/
// pr at the disc
      y1 = pr[i + nrad * (ih0 - 1) + nrad * nh * (irh0 - 1)];
      y2 = pr[i + nrad * (ih0 - 1) + nrad * nh * irh0];
      y3 = pr[i + nrad * ih0 + nrad * nh * irh0];
      y4 = pr[i + nrad * ih0 + nrad * nh * (irh0 - 1)];
      pr_final = utmp1 * (ttmp1 * y1 + ttmp * y2) + 
                 utmp * (ttmp * y3 + ttmp1 * y4);
// temporary variables
      r = r_plus + radius[i];
      r2 = r * r;
      delta = r2 - 2. * r + am2;
      ULt = sqrt((h * h + am2) / (h * h - 2. * h + am2));
      if (r >= rms) {
        tmp1 = sqrt(r2 - 3. * r + 2. * am * sqrt(r));
        Ut = (r2 + am * sqrt(r)) / r / tmp1;
//      U_phi = (r2 + am2 - 2. * am * sqrt(r)) / sqrt(r) / tmp1;
        U_r = 0.;
//      Ur = 0.;
        UrU_r = 0.;
      }
      else {
        tmp1 = sqrt(rms * (rms - 3.) + 2. * am * sqrt(rms));
        Lms = (rms * rms + am2 - 2. * am * sqrt(rms)) / sqrt(rms) / tmp1;
        Ems = (rms * (rms - 2.) + am * sqrt(rms)) / rms / tmp1;
        Ut = (Ems * (r * (r2 + am2) + 2. * am2) - 2. * am * Lms) / r / delta;
//      U_phi = Lms;
        UrU_r = -1. + ((r2 + am2 + 2. * am2 / r) * Ems * Ems - 4. * am / 
                r * Ems * Lms - (1. - 2. / r) * Lms * Lms) / delta;
        if (UrU_r < 0.) UrU_r = 0.;
        U_r = -sqrt(UrU_r / delta) * r;
//      Ur = -sqrt(delta * UrU_r) / r;
      }
      tmp1 = Ut - pr_final * U_r;
// gfactor  from the axis to the disc
      gfac[i] = tmp1 / ULt;
// cosin at the disc
//    cosin[i] = q_final / r / tmp1;
// phip_i at the disc
//    phiph[i] = atan2(-U_phi, r * (pr_final - Ur *tmp1));
// dWadSd from the axis to the disc
      y1 = dWadSd[i + nrad * (ih0 - 1) + nrad * nh * (irh0 - 1)];
      y2 = dWadSd[i + nrad * (ih0 - 1) + nrad * nh * irh0];
      y3 = dWadSd[i + nrad * ih0 + nrad * nh * irh0];
      y4 = dWadSd[i + nrad * ih0 + nrad * nh * (irh0 - 1)];
      transf_d[i] = utmp1 * (ttmp1 * y1 + ttmp * y2) + 
                    utmp * (ttmp * y3 + ttmp1 * y4);
    }
  }
//******************************************************************************
//    fprintf(stdout,"%f %f\n", thetaO, transf_o);
//    for(i = 0; i < nrad; i++) 
//      fprintf(stdout,"%d %f %f %f %f %f\n", i, radius[i], gfac[i], cosin[i], 
//              phiph[i], transf_d[i]);
//******************************************************************************
}
/******************************************************************************/
// Let's read the reflionx tables
if ((strcmp(kydir, FGMSTR(pname)) != 0) || (rix != rix_old) || first_rix) {
  sprintf(text, "kynrefionx: initializing %s tables...", reflion);
  xs_write(text, 5);
  xs_write("Ross & Fabian (2005), MNRAS, 358, 211",5);
// The status parameter must always be initialized.
  status = 0;
// Open the FITS file for readonly access
// - if set try KYDIR directory, otherwise look in the working directory
//   or in the xspec directory where tables are usually stored...
  sprintf(kydir, "%s", FGMSTR(pname));
  if (strlen(kydir) == 0) sprintf(tables_file, "./%s", reflion);
  else if (kydir[strlen(kydir) - 1] == '/') sprintf(tables_file, "%s%s",
                                                    kydir, reflion);
  else sprintf(tables_file, "%s/%s", kydir, reflion);
// Let's read the reflion(x) tables
// The status parameter must always be initialized.
  status = 0;
  ffopen(&fptr, tables_file, READONLY, &status);
  if (status) {
    sprintf(tables_file, "%s%s", FGMODF(), reflion);
    status = 0;
    ffopen(&fptr, tables_file, READONLY, &status);
  }
  if (status) {
    if (status) ffrprt(stderr, status);
    ffclos(fptr, &status);
    xs_write("\nkynrefionx: set the KYDIR to the directory with the KY tables.", 5);
    for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
    return;
  }
  if (((strcmp(kydir, FGMSTR(pname)) != 0) || (rix != rix_old)) && !first_rix) {
// Free memory from tmp arrays...
    free((void *) abun);
    abun = NULL;
    free((void *) gam);
    gam = NULL;
    free((void *) xi);
    xi = NULL;
    free((void *) logxi);
    logxi = NULL;
    free((void *) energy0);
    energy0 = NULL;
    free((void *) emission);
    emission = NULL;
    free((void *) flux0);
    flux0 = NULL;
    free((void *) flux1);
    flux1 = NULL;
    free((void *) energy1);
    energy1 = NULL;
    free((void *) energy2);
    energy2 = NULL;
  }
// Let's read tables (binary tables => hdutype=2)
// Move to the first extension ('PARAMETERS')
  ffmrhd(fptr, 1, &hdutype, &status);
// Read the values of parameters --> abundance, photon index and ionisation
  nelems = 1;
  colnum = 9;
  frow = 1;
  ffgcv(fptr, TLONG, colnum, frow, felem, nelems, &float_nulval, &nabun,
        &anynul, &status);
  frow = 2;
  ffgcv(fptr, TLONG, colnum, frow, felem, nelems, &float_nulval, &ngam,
        &anynul, &status);
  frow = 3;
  ffgcv(fptr, TLONG, colnum, frow, felem, nelems, &float_nulval, &nxi,
        &anynul, &status);
    nspec = nabun * ngam * nxi;
//******************************************************************************
//  fprintf(stdout,"%d %d %d\n", nabun, ngam, nxi);
//******************************************************************************
// Allocate memory for abun, gam and xi...
  if ((abun = (float *) malloc(nabun * sizeof(float))) == NULL) {
    xs_write("kynrefionx: Failed to allocate memory for tmp arrays.", 5);
    for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
    return;
  }
  if ((gam = (float *) malloc(ngam * sizeof(float))) == NULL) {
    xs_write("kynrefionx: Failed to allocate memory for tmp arrays.", 5);
    for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
    return;
  }
  if ((xi = (float *) malloc(nxi * sizeof(float))) == NULL) {
    xs_write("kynrefionx: Failed to allocate memory for tmp arrays.", 5);
    for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
    return;
  }
  if ((logxi = (float *) malloc(nxi * sizeof(float))) == NULL) {
    xs_write("kynrefionx: Failed to allocate memory for tmp arrays.", 5);
    for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
    return;
  }
  colnum = 10;
  frow = 1;
  nelems = nabun;
  ffgcv(fptr, TFLOAT, colnum, frow, felem, nelems, &float_nulval, abun,
        &anynul, &status);
  frow = 2;
  nelems = ngam;
  ffgcv(fptr, TFLOAT, colnum, frow, felem, nelems, &float_nulval, gam,
        &anynul, &status);
  frow = 3;
  nelems = nxi;
  ffgcv(fptr, TFLOAT, colnum, frow, felem, nelems, &float_nulval, xi,
        &anynul, &status);
  for (i = 0; i < nxi; i++) logxi[i] = log10(xi[i]);
//******************************************************************************
//  for(i = 0; i < 20; i++)
//    fprintf(stdout,"%4d %12.1f %12.1f %12.1f %12.1f\n", i, abun[i], gam[i], 
//            xi[i], logxi[i]);
//******************************************************************************
// Move to the second extension 'ENERGIES' and read energy values
  ffmrhd(fptr, 1, &hdutype, &status);
  ffgnrw(fptr, &ne_loc, &status);
  ne_loc++;
//******************************************************************************
//  fprintf(stdout,"%d\n", ne_loc);
//******************************************************************************
// Allocate memory for energy...
  if ((energy0 = (float *) malloc(ne_loc * sizeof(float))) == NULL) {
    xs_write("kynrefionx: Failed to allocate memory for tmp arrays.", 5);
    for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
    return;
  }
  nelems = ne_loc - 1;
  colnum = 1;
  frow = 1;
  ffgcv(fptr, TFLOAT, colnum, frow, felem, nelems, &float_nulval, energy0,
        &anynul, &status);
  nelems = 1;
  colnum = 2;
  frow = ne_loc - 1;
  ffgcv(fptr, TFLOAT, colnum, frow, felem, nelems, &float_nulval, &energy0[ne_loc - 1],
        &anynul, &status);
//******************************************************************************
//  for(i = 0; i < ne_loc; i++) fprintf(stdout,"%d %f\n", i, energy0[i]);
//******************************************************************************
// Allocate memory for emission...
  if ((emission = (float *) malloc(nspec * (ne_loc - 1) * sizeof(float))) == NULL) {
    xs_write("kynrefionx: Failed to allocate memory for tmp arrays.", 5);
    for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
    return;
  }
// Let's read the tables for emission
  ffmrhd(fptr, 1, &hdutype, &status);
// to read the file only once we have to read in blocks (all columns
// from the extension are put to buffer together)
// let's find out how many rows are going to be read into the buffer
  ffgrsz(fptr, &nrow, &status);
  nelements = nrow * (ne_loc - 1);
  for (irow = 0; irow < nspec; irow += nrow) {
// the last block to read may be smaller:
    if ((nspec - irow) < nrow) nelements = (nspec - irow) * (ne_loc - 1);
    iabun = irow / (ngam * nxi) + 1;
    igam = (irow - (iabun - 1) * ngam * nxi) / nxi + 1;
    ixi = irow + 1 - (iabun - 1) * ngam * nxi - (igam - 1) * nxi;
    ffgcv(fptr, TFLOAT, 2, irow + 1, 1, nelements, &float_nulval, 
          &emission[(ne_loc - 1) * (ixi - 1) + (ne_loc - 1) * nxi * (igam - 1) + 
                    (ne_loc - 1) * nxi * ngam * (iabun - 1)], &anynul, &status);
  }
// The FITS file must always be closed before exiting the program.
  ffclos(fptr, &status);
/*******************************************************************************
  iabun = 7;
  igam = 5;
  ixi = 3;
  fprintf(stdout, "%f %f %f %f\n", abun[iabun - 1], gam[igam - 1], 
                                   pow(10, xi[ixi - 1]), xi[ixi - 1]);
  for(ie = 0; ie < ne_loc - 1; ie++)
    fprintf(stdout, "%d %f %f\n", ie, energy0[ie], emission[ie + 
      (ne_loc - 1) * (ixi - 1) + (ne_loc - 1) * nxi * (igam - 1) + 
      (ne_loc - 1 ) * nxi * ngam * (iabun - 1)]);
*******************************************************************************/
// Allocate memory for local flux...
  if ((flux0 = (double *) malloc(ne_loc * nxi * sizeof(double))) == NULL) {
    xs_write("kynrefionx: Failed to allocate memory for tmp arrays.", 5);
    for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
    return;
  }
  if ((flux1 = (double *) malloc(nener * nxi * sizeof(double))) == NULL) {
    xs_write("kynrefionx: Failed to allocate memory for tmp arrays.", 5);
    for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
    return;
  }
  if ((energy1 = (double *) malloc((nener + 1) * sizeof(double))) == NULL) {
    xs_write("kynrefionx: Failed to allocate memory for tmp arrays.", 5);
    for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
    return;
  }
  if ((energy2 = (double *) malloc((nener + 1) * sizeof(double))) == NULL) {
    xs_write("kynrefionx: Failed to allocate memory for tmp arrays.", 5);
    for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
    return;
  }
  first_rix = 0;
  xs_write("------------------------------------------------", 5);
}
// end of reading reflionx fits file............................................
/******************************************************************************/
// Let's interpolate the emission to get rid of the non-radial dependences, i.e.
// interpolate in abundance and gamma0
if ((rix != rix_old) || (strcmp(kydir, FGMSTR(pname)) != 0) || ((abun_old == -1.) 
    || (abundance != abun_old)) || ((gam_old == -1) || (gamma0 != gam_old))) {
// given gamma0, find the corresponding index in gam():
  imin = 0;
  imax = ngam;
  igam0 = ngam / 2;
  while ((imax - imin) > 1) {
    if (gamma0 >= gam[igam0 - 1]) imin = igam0;
    else imax = igam0;
    igam0 = (imin + imax) / 2;
  }
  if (igam0 == 0) igam0 = 1;
  ttmp = (gamma0 - gam[igam0 - 1]) / (gam[igam0] - gam[igam0 - 1]);
  if (ttmp < 0.) ttmp = 0.;
  if (ttmp > 1.) ttmp = 1.;
  ttmp1 = 1. - ttmp;
// given abundance, find the corresponding index in abun[]:
  imin = 0;
  imax = nabun;
  iabun0 = nabun / 2;
  while ((imax - imin) > 1) {
    if (abundance >= abun[iabun0 - 1]) imin = iabun0;
    else imax = iabun0;
    iabun0 = (imin + imax) / 2;
  }
  if (iabun0 == 0) iabun0 = 1;
  utmp = (abundance - abun[iabun0 - 1]) / (abun[iabun0] - abun[iabun0 - 1]);
  if (utmp < 0.) utmp = 0.;
  if (utmp > 1.) utmp = 1.;
  utmp1 = 1. - utmp;
// In the following the last energy index is not used 
// (we have one less flux indices then in energy)
  for (ixi = 0; ixi < nxi; ixi++)
    for (ie = 0; ie < ne_loc - 1; ie++) {
      y1 = emission[ie + (ne_loc - 1) * ixi + (ne_loc - 1) * nxi * (igam0 - 1) + 
                    (ne_loc - 1) * nxi * ngam * (iabun0 - 1)];
      y2 = emission[ie + (ne_loc - 1) * ixi + (ne_loc - 1) * nxi * igam0 + 
                    (ne_loc - 1) * nxi * ngam * (iabun0 - 1)];
      y3 = emission[ie + (ne_loc - 1) * ixi + (ne_loc - 1) * nxi * igam0 + 
                    (ne_loc - 1) * nxi * ngam * iabun0];
      y4 = emission[ie + (ne_loc - 1) * ixi + (ne_loc - 1) * nxi * (igam0 - 1) + 
                    (ne_loc - 1) * nxi * ngam * iabun0];
      flux0[ie + ne_loc * ixi] = (utmp1 * (ttmp1 * y1 + ttmp * y2) + 
                                  utmp * (ttmp * y3 + ttmp1 * y4));
    }
/*******************************************************************************
  ixi = 8;
  fprintf(stdout, "%d %d %d\n", iabun0, igam0, ixi);
  fprintf(stdout, "%f %f %f %f\n", abun[iabun0 - 1], gam[igam0 - 1], 
          pow(10, logxi[ixi - 1]), xi[ixi - 1]);
  for(ie = 0; ie < ne_loc - 1; ie++) 
    fprintf(stdout, "%d %f %f\n", ie, energy0[ie], 
            flux0[ie + ne_loc * (ixi - 1)] / (energy0[i + 1] - energy0[i]));
*******************************************************************************/
// Let's rebin the spectra to the evenly log spaced energies
  if ((rix != rix_old) || (strcmp(kydir, FGMSTR(pname)) != 0)) {
    energy2[0] = elow / (1. + pow(ehigh / elow, 1. / (nener - 1.))) * 2.;
    for (ie = 1; ie <= nener; ie++) {
      energy2[ie] = energy2[0] * pow(ehigh / elow, ie / (nener - 1.));
      energy1[ie - 1] = elow * pow(ehigh / elow, (ie - 1.) / (nener - 1.));
    }
  }
  for (ixi = 0; ixi < nxi; ixi++)
    for (ie = 0; ie < nener; ie++) flux1[ie + nener * ixi] = 0.;
  ie = 1;
  while (energy0[ie] <= energy2[0]) ie++;
  je = 1;
  quit = 0;
  while ((ie <= (ne_loc - 1)) && (energy0[ie - 1] < energy2[nener])) {
    ttmp = energy0[ie] - energy0[ie - 1];
    while (energy2[je - 1] < energy0[ie]) {
      if (energy2[je - 1] < energy0[ie - 1]) utmp = energy0[ie - 1];
      else utmp = energy2[je - 1];
      if (energy2[je] < energy0[ie]) utmp1 = energy2[je];
      else utmp1 = energy0[ie];
      if (utmp1 > utmp)
        for (ixi = 0; ixi < nxi; ixi++) 
          flux1[je - 1 + nener * ixi] +=
            flux0[ie - 1 + ne_loc * ixi] * (utmp1 - utmp) / ttmp;
      if (je < nener) je++;
      else {
        quit = 1;
        break;
      }
    }
    if (quit) break;
    if (energy2[je - 1] > energy0[ie]) je--;
    ie++;
  }
  for (ixi = 0; ixi < nxi; ixi++)
    for (ie = 0; ie < nener; ie++) 
      flux1[ie + nener * ixi] /= (energy2[ie + 1] - energy2[ie]);
/*******************************************************************************
  ixi = 8;
  fprintf(stdout, "%d %d %d\n", iabun0, igam0, ixi);
  fprintf(stdout, "%f %f %f %f\n", abun[iabun0 - 1], gam[igam0 - 1], 
          pow(10, logxi[ixi - 1]), xi[ixi - 1]);
  for(ie = 0; ie < nener; ie++) 
    fprintf(stdout, "%d %f %f\n", ie, energy1[ie], flux1[ie + nener * (ixi - 1)]);
*******************************************************************************/
  abun_old = abundance;
  gam_old = gamma0;
}
sprintf(kydir, "%s", FGMSTR(pname));
rix_old = rix;
am_old = am;
thetaO_old = thetaO;
h_rh_old = h_rh;
/******************************************************************************/
//Let's compute the Np according to its definition, i.e. transform from 
//intrinsic or observed in 2-10keV to total intrinsic photon flux
//Lx is intrinsic photon flux in 2-10keV
//The flux is corrected for cosmological redshift, thus we do not use zzshift
//in luminosity computations!
g_L = sqrt(1. - 2. * h / (am2 + h * h));
if( Np < 0. ){
  Lx = -Np;
  Np *= - incgamma(2. - gamma0, E0 / EC) / 
        ( incgamma(2. - gamma0, 2. / EC) - incgamma(2. - gamma0, 10. / EC));
}else{
  Lx = Np / g_L / g_L / transf_o /
       ( incgamma(2. - gamma0, 2. / g_L / EC) - incgamma(2. - gamma0, 10. / g_L / EC));
  Np = Lx * incgamma(2. - gamma0, E0 / EC);        
  Lx *= ( incgamma(2. - gamma0, 2. / EC) - incgamma(2. - gamma0, 10. / EC));
}
//Let's write the intrinsic photon flux in 2-10keV into the xset XSPEC variable
//KYLxLamp
sprintf(kyLxLamp, "%e", Lx);
FPMSTR(pkyLxLamp, kyLxLamp);
// Let's compute the ionisation at rin and rout
for (i = 0; i < 2; i++) {
  if (i == 0)
    if (rin <= (radius[0] + r_plus)) r = radius[1] + r_plus;
    else r = rin;
  else r = rout;
  imin = 0;
  imax = nrad;
  ir0 = nrad / 2;
  while ((imax - imin) > 1) {
    if (r >= (radius[ir0 - 1] + r_plus)) imin = ir0;
    else imax = ir0;
    ir0 = (imin + imax) / 2;
  }
  if ((ir0 == nrad) && (r == (radius[nrad - 1] + r_plus))) ir0 = ir0 - 1;
  ttmp = (r - radius[ir0 - 1] - r_plus) / (radius[ir0] - radius[ir0 - 1]);
  ttmp1 = 1. - ttmp;
// Let's interpolate gfactor between two radii
  gfactor0 = ttmp * gfac[ir0] + ttmp1 * gfac[ir0 - 1];
// Let's interpolate lensing between two radii
  lensing = (ttmp * transf_d[ir0] + ttmp1 * transf_d[ir0 - 1]) * 
            h * sqrt(1. - 2. * h / (h * h + am2)) / (r * r + h * h) / r;
  if (lensing != 0.) {
    if(nH0 > 0.){
      if (qn != 0.) {
        if (sw == 1) ionisation = logxi_norm + log10(pow(r, -qn) * lensing * 
                                  gfactor0 * Np / mass2 / nH0);
        else ionisation = logxi_norm + log10(pow(r, -qn) * lensing * 
                          pow(gfactor0, gamma0 - 1.) * Np / mass2 / nH0);
      } else {
        if (sw == 1) ionisation = logxi_norm + 
                                  log10(lensing * gfactor0 * Np / mass2 / nH0);
        else ionisation = logxi_norm + 
               log10(lensing * pow(gfactor0, gamma0 - 1.) * Np / mass2 / nH0);
      }
    } else ionisation = log10(-nH0) + qn * log10(r);
    if (i == 0) {
      sprintf(kyxiin, "%e", pow(10, ionisation));
      FPMSTR(pkyxiin, kyxiin);
    }
    else {
      sprintf(kyxiout, "%e", pow(10, ionisation));
      FPMSTR(pkyxiout, kyxiout);
    }
  }
}
/******************************************************************************/
#ifdef OUTSIDE_XSPEC
// let's write the input parameters to a file
fw = fopen("kynrefionx.txt", "w");
fprintf(fw, "a/m         %12.6f\n", param[0]);
fprintf(fw, "thetaO      %12.6f\n", param[1]);
fprintf(fw, "rin         %12.6f\n", param[2]);
fprintf(fw, "ms          %12d\n", (int) param[3]);
fprintf(fw, "rout        %12.6f\n", param[4]);
fprintf(fw, "phi         %12.6f\n", param[5]);
fprintf(fw, "dphi        %12.6f\n", param[6]);
fprintf(fw, "M/M8        %12.6f\n", param[7]);
fprintf(fw, "height      %12.6f\n", param[8]);
fprintf(fw, "PhoIndex    %12.6f\n", param[9]);
fprintf(fw, "L/Ledd      %12.6f\n", param[10]);
fprintf(fw, "Np:Nr       %12.6f\n", param[11]);
fprintf(fw, "nH0         %12.6f\n", param[12]);
fprintf(fw, "qn          %12.6f\n", param[13]);
fprintf(fw, "abun        %12.6f\n", param[14]);
fprintf(fw, "limb        %12d\n", (int) param[19]);
fprintf(fw, "tab         %12d\n", (int) param[20]);
fprintf(fw, "sw          %12d\n", (int) param[21]);
fprintf(fw, "alpha       %12.6f\n", ide_param[21]);
fprintf(fw, "beta        %12.6f\n", ide_param[22]);
fprintf(fw, "rcloud      %12.6f\n", ide_param[23]);
fprintf(fw, "zshift      %12.6f\n", param[18]);
fprintf(fw, "ntable      %12d\n", (int) param[22]);
fprintf(fw, "nrad        %12d\n", (int) param[23]);
fprintf(fw, "division    %12d\n", (int) param[24]);
fprintf(fw, "nphi        %12d\n", (int) param[25]);
fprintf(fw, "smooth      %12d\n", (int) param[26]);
fprintf(fw, "Stokes      %12d\n", stokes);
fprintf(fw, "polar       %12d\n", polar);
fprintf(fw, "r_horizon   %12.6f\n", r_plus);
fprintf(fw, "r_ms        %12.6f\n", rms);
fprintf(fw, "edivision   %12d\n", (int) ide_param[14]);
fprintf(fw, "nener       %12ld\n", nener);
fprintf(fw, "norm        %12.6f\n", param[28]);
fprintf(fw, "nthreads    %12d\n", (int) param[27]);
fprintf(fw, "Xi_in       %s\n", kyxiin);
fprintf(fw, "Xi_out      %s\n", kyxiout);
fclose(fw);
#endif
/******************************************************************************/

// Let's integrate local emission over the accretion disc
if (ide(ear, ne, 1, far, qar, uar, var, ide_param, emis_KYNrefionx , nener)) {
  for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
  return;
}
// Let's normalize the reflected flux properly
for (ie = 0; ie < ne; ie++) {
  far[ie] *= Dnorm * mass2 * RG2 * MPC_2 * REFLIONX_NORM;
/*  if (polar == 1) {
    qar[ie] *= Dnorm * mass2 * RG2 * MPC_2 * REFLIONX_NORM;
    uar[ie] *= Dnorm * mass2 * RG2 * MPC_2 * REFLIONX_NORM;
    var[ie] *= Dnorm * mass2 * RG2 * MPC_2 * REFLIONX_NORM;
  }*/
}
// Let's add primary flux to the solution (note we multiply by dt later on)
refl_ratio=-1.;
if (NpNr != 0) {
// Let's compute the incomplete gamma function with the XSPEC incgamma 
// function
  Anorm = LEDD * mass * MPC_2 * ERG * Np / pow(EC, 2. - gamma0) / PI2 / 2. / 
          incgamma(2. - gamma0, E0 / EC);
// let's compute the cut-off powerlaw with the XSPEC routine cutoffPowerLaw
  for (ie = 0; ie <= ne; ie++) ear1[ie] = ear[ie];
  param1[0] = param[9];
  param1[1] = g_L * zzshift * EC;
  cutoffpl(ear1, ne, param1, photar1);
  Anorm *= Dnorm * fabs(NpNr) * transf_o * pow(g_L * zzshift, gamma0);
  flux_refl = flux_prim = 0.;
  for (ie = 0; ie < ne; ie++){
    flux_refl += far[ie];
    if (ear[ie] > g_L * zzshift * E0){
      flux_prim += Anorm * photar1[ie];
      far[ie] += Anorm * photar1[ie];
    }
  }
  refl_ratio = flux_refl / flux_prim;
}
sprintf(kyRefl, "%e", refl_ratio);
FPMSTR(pkyRefl, kyRefl);

// interface with XSPEC
if (!stokes) for (ie = 0; ie < ne; ie++) photar[ie] = far[ie];
else {
  pamin = 1e30;
  pamax = -1e30;
  pa2min = 1e30;
  pa2max = -1e30;
  for (ie = ne - 1; ie >= 0; ie--) {
    pd[ie] = sqrt(qar[ie] * qar[ie] + uar[ie] * uar[ie] + var[ie] * var[ie]) /
             (far[ie] + 1e-30);
    pa[ie] = 0.5 * atan2(uar[ie], qar[ie]) / PI * 180.;
    if (ie < (ne - 1)) {
      while ((pa[ie] - pa[ie + 1]) > 90.) pa[ie] -= 180.;
      while ((pa[ie + 1] - pa[ie]) > 90.) pa[ie] += 180.;
    }
    if (pa[ie] < pamin) pamin = pa[ie];
    if (pa[ie] > pamax) pamax = pa[ie];
    pa2[ie] = 0.5 * asin(var[ie] / sqrt(qar[ie] * qar[ie] + uar[ie] * uar[ie] +
                                  var[ie] * var[ie] + 1e-30)) / PI * 180.;
    if (ie < (ne - 1)) {
      while ((pa2[ie] - pa2[ie + 1]) > 90.) pa2[ie] -= 180.;
      while ((pa2[ie + 1] - pa2[ie]) > 90.) pa2[ie] += 180.;
    }
    if (pa2[ie] < pa2min) pa2min = pa2[ie];
    if (pa2[ie] > pa2max) pa2max = pa2[ie];
  }
  fw = fopen("stokes.dat", "w");
  for (ie = 0; ie < ne; ie++) {
    if ((pamax + pamin) > 180.) pa[ie] -= 180.;
    if ((pamax + pamin) < -180.) pa[ie] += 180.;
    if ((pa2max + pa2min) > 180.) pa2[ie] -= 180.;
    if ((pa2max + pa2min) < -180.) pa2[ie] += 180.;
    fprintf(fw,
      "%14.6f\t%14.6f\t%14.6f\t%14.6f\t%14.6f\t%14.6f\t%14.6f\t%14.6f\n", 
      0.5 * (ear[ie] + ear[ie+1]), far[ie] / (ear[ie+1] - ear[ie]), 
      qar[ie] / (ear[ie+1] - ear[ie]), uar[ie] / (ear[ie+1] - ear[ie]), 
      var[ie] / (ear[ie+1] - ear[ie]), pd[ie], pa[ie], pa2[ie]);
//interface with XSPEC..........................................................
    if (stokes == 1) photar[ie] = qar[ie];
    if (stokes == 2) photar[ie] = uar[ie];
    if (stokes == 3) photar[ie] = var[ie];
    if (stokes == 4) photar[ie] = pd[ie] * (ear[ie + 1] - ear[ie]);
    if (stokes == 5) photar[ie] = pa[ie] * (ear[ie + 1] - ear[ie]);
    if (stokes == 6) photar[ie] = pa2[ie] * (ear[ie + 1] - ear[ie]);
  }
  fclose(fw);
}

/******************************************************************************/
#ifdef OUTSIDE_XSPEC
// final spectrum output -- write ear[] and photar[] into file:
fw = fopen("kynrefionx_photar.dat", "w");

if( NpNr != 0. )
  for (ie = 0; ie < ne; ie++) fprintf(fw, "%14.6f\t%E\t%E\n", 
    0.5*(ear[ie]+ear[ie+1]), 
    (photar[ie]-Anorm*photar1[ie]) / (ear[ie+1] - ear[ie]),
    (Anorm*photar1[ie]) / (ear[ie+1] - ear[ie]));  
else
  for (ie = 0; ie < ne; ie++) fprintf(fw, "%14.6f\t%E\n", 
    0.5*(ear[ie]+ear[ie+1]), photar[ie] / (ear[ie+1] - ear[ie]));
fclose(fw);
#endif
/******************************************************************************/

return;
}

/*******************************************************************************
*******************************************************************************/
void emis_KYNrefionx(double** ear_loc, const int ne_loc, const int nt, 
                     double *far_loc, double *qar_loc, double *uar_loc, 
                     double *var_loc, const double r, const double phi, 
                     const double cosmu, const double phiphoton, 
                     const double alpha_o, const double beta_o, 
                     const double delay, const double g) {
// local emissivity --> far_loc(:) array
// local energy array --> ear_loc()
// ne_loc --> number of points in local energy where the local photon flux
// density in keV is defined;
// this is a steady model (nt = 1);
// disc surface in polar coords r, phi;
// cosine of local emission angle --> cosmu

double factor, factor1, factor2, gfactor, lensing, ionisation, fluxe[ne_loc];
double ttmp, ttmp1, y1, y2;
long   ixi0;
int    ie, imin, imax, ir0;

*ear_loc = energy1;
// Normalization due to an imposed emissivity law
//if (limb == 0) factor = 1. / PI2;
//if (limb == 1) factor = 1. / PI2 / (1. + 2.06 / 2.) * (1. + 2.06 * cosmu);
//if (limb == 2) factor = 1. / PI2 / log(4.) * log(1. + 1. / cosmu);
if (limb == 0) factor = 1. / PI;
if (limb == 1) factor = 1. / PI / (1. + 2.06 * 2. / 3.) * (1. + 2.06 * cosmu);
if (limb == 2) factor = 1. / PI * log(1. + 1. / cosmu);
if(nH0 > 0.) factor *= nH0;
// given r, find corresponding indices in radius:
imin = 0;
imax = nrad;
ir0 = nrad / 2;
while ((imax - imin) > 1) {
  if (r >= (radius[ir0 - 1] + r_plus)) imin = ir0;
  else imax = ir0;
  ir0 = (imin + imax) / 2;
}
//if (ir0 == 0) ir0 = 1;
//if ((imax == nrad) && (r > (radius[nrad - 1] + r_plus)) ir0 = nrad;
if ((ir0 == nrad) && (r == (radius[nrad - 1] + r_plus))) ir0 -= 1;
if ((ir0 == 0) || (ir0 >= nrad)) {
  for (ie = 0; ie < ne_loc; ie++) far_loc[ie] = 0.;
/*if (polar) {
    for (ie = 0; ie < ne_loc; ie++) {
      qar_loc[ie] = 0.;
      uar_loc[ie] = 0.;
      var_loc[ie] = 0.;
    }
  }*/ 
}else  {
  ttmp = (r - radius[ir0 - 1] - r_plus) / (radius[ir0] - radius[ir0 - 1]);
  ttmp1 = 1. - ttmp;
// Let's interpolate gfactor between two radii
  gfactor = ttmp * gfac[ir0] + ttmp1 * gfac[ir0 - 1];
// Let's interpolate cosmu0 between two radii
//  cosmu0 = ttmp * cosin[ir0] + ttmp1 * cosin[ir0 - 1];
// Let's interpolate lensing between two radii
  lensing = (ttmp * transf_d[ir0] + ttmp1 * transf_d[ir0 - 1]) * 
            h * sqrt(1. - 2. * h / (h * h + am2)) / (r * r + h * h) / r;
// Let's compute the emitted flux at the particular radius
  if (lensing != 0.) {
    if(nH0 > 0.){
      if (qn != 0.) {
        if (sw == 1) ionisation = logxi_norm + log10(pow(r, -qn) * lensing * 
                                  gfactor * Np / mass2 / nH0);
        else ionisation = logxi_norm + log10(pow(r, -qn) * lensing * 
                          pow(gfactor, gamma0 - 1.) * Np / mass2 / nH0);
      } else {
        if (sw == 1) ionisation = logxi_norm + log10(lensing * gfactor * 
                                  Np / mass2 / nH0);
        else ionisation = logxi_norm + log10(lensing * 
                          pow(gfactor, gamma0 - 1.) * Np / mass2 / nH0);
      }
      factor1 = 1.;
    } else{
      ionisation = log10(-nH0) + qn * log10(r);
      if(sw==1) factor1 = pow(10, logxi_norm - qn * log10(r)) * 
                          lensing * gfactor * Np / mass2 / (-nH0);
      else factor1 = pow10(logxi_norm - qn * log10(r)) *
                     lensing * pow(gfactor, gamma0 - 1.) * Np / mass2 / (-nH0);
    }
// give ionisation, find the corresponding index in logXi():
    imin = 0;
    imax = nxi;
    ixi0 = nxi / 2;
    while ((imax - imin) > 1) {
      if (ionisation >= logxi[ixi0 - 1]) imin = ixi0;
      else imax = ixi0;
      ixi0 = (imin + imax) / 2;
    }
    if (ixi0 == 0) ixi0 = 1;
    ttmp = (ionisation - logxi[ixi0 - 1]) / (logxi[ixi0] - logxi[ixi0 - 1]);
    if (ttmp < 0.) {
      ttmp = 0.;
      factor1 *= pow(10, ionisation - logxi[0]);
    }
    if (ttmp > 1.) {
      ttmp = 1.;
      factor1 *= pow(10, ionisation - logxi[nxi - 1]);
    }
    ttmp1 = 1. - ttmp;
    for (ie = 0; ie < ne_loc; ie++) {
      y1 = flux1[ie + ne_loc * (ixi0 - 1)];
      y2 = flux1[ie + ne_loc * ixi0];
      factor2 = exp( energy1[ie] / EC * (1.-1./gfactor) );
      fluxe[ie] = (ttmp1 * y1 + ttmp * y2) * factor * factor1 * factor2;
      if (sw == 1) fluxe[ie] *= pow(gfactor, gamma0 - 2.);
      if (nH0 > 0. && qn != 0.) fluxe[ie] *= pow(r, qn);
    }
  }else for (ie = 0; ie < ne_loc; ie++) fluxe[ie] = 0;
  for (ie = 0; ie < ne_loc; ie++) far_loc[ie] = fluxe[ie];
}
/*******************************************************************************
// local spectrum output -- write energy1[] and far_local[] into file
  fw = fopen("kynrefionx_photar_loc.dat", "w");
  for (ie = 0; ie < ne_loc; ie++)
    fprintf(fw, "%14.6f\t%E\n", energy1[ie], fluxe[ie]);
  fclose(fw)
*******************************************************************************/
return;
}
