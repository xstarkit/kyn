/* KYNphebb - phenomenological black-body emission - non-axisymmetric version
 * 
 * ref. Dovciak M., Karas V., Yaqoob T. (2004)
 * -----------------------------------------------------------------------------
 * OTHER REFERENCES:
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
 * This model computes phenomenological black-body emission from an accretion
 * disc around a black hole. It is phenomenological in the sense that the radial
 * dependence of the disc temperature is a simple powerlaw. The local flux is
 * defined as flux ~ E^2/(exp(E/kT)-1), where T=Tin*(r/rin)^(-BBindex).
 * All relativistic effects that change properties of the light on its way from
 * the disc to the observer are taken into account. This model calls subroutine
 * ide() for integrating local emission over the disc and uses the FITS file
 * 'KBHtablesNN.fits' defining the transfer functions needed for integration.
 * For details on ide() and the FITS file see the subroutine ide() in xside.c.
 *
 * The polarisation is computed from Chandrasekhar's formular for infinite 
 * optical depth (parameter tau=11) or by the STOKES Monte Carlo code for 
 * optical depths tau = 0.2, 0.5, 1.0, 2.0, 5.0, and 10.0. stored in the tables
 * 'goosmann.fits' (computed by the author of STOKES code Rene Goosmann). The
 * tau=10 results are the same as tau=infinity. The tables correspond 
 * to a model setup as follows: a plane-parallel electron scattering disc is 
 * irradiated from its midplane and evaluated for the Stokes parameters I and Q 
 * at 40 different viewing angles, which are given by their cosine values. 
 * The values of I and Q are normalized by the total number of photons sampled. 
 * A positive value of Q denotes a polarization vector that is parallel with 
 * the disk's symmetry axis, a negative value stands for a vector being 
 * perpendicular to this axis. The U values are basically zero in all cases. 
 * The intrinsic irradiation at the midplane is assumed to be isotropic at 
 * every point. The electron scattering is realized by Thomson scattering and 
 * therefore wavelength-independent.
 *
 * par1  ... a/M     - black hole angular momentum (-1 <= a/M <= 1)
 * par2  ... theta_o - observer inclination in degrees (0-pole, 90-disc)
 *                     if negative, the results are computed for the opposite 
 *                     (and up-side down) observer located on the other side of 
 *                     the disc (with the same inclination angle), i.e. 
 *                     the whole system (both the disc and black hole) is
 *                     rotating counter-clockwise direction for the positive 
 *                     inclination, and clockwise direction for the negative 
 *                     inclination; this is important only for polarisation
 *                     properties (Stokes parameter, par19, larger than 1)
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
 *                 - if outer edge is equal or larger than 100000 GM/c^2 then 
 *                   the emission from above this radius is added
 * par6  ... phi   - lower azimuth of non-zero disc emissivity (deg)
 * par7  ... dphi  - (phi + dphi) is upper azimuth of non-zero disc emissivity
 *                   0 <= dphi <= 360  (deg)
 * par8  ... Tin     - temperature in keV at the inner edge of the disc
 * par9  ... BBindex - radial power-law index for radial dependence of the
 *                     black-body temperature
 * par10 ... alpha  - position of the cloud centre in GM/c^2 in alpha coordinate
 *                    (alpha being the impact parameter in phi direction, 
 *                     positive for approaching side of the disc)
 * par11 ... beta   - position of the cloud centre in GM/c^2 in beta coordinate
 *                    (beta being the impact parameter in theta direction, 
 *                     positive in up direction, i.e. above the disc)
 * par12 ... rcloud - radius of the obscuring cloud (in GM/c^2)
 *                  - if negative, only the emission transmitted through
 *                    the cloud is taken into account
 * par13 ... zshift - overall Doppler shift
 * par14 ... ntable - table of relativistic transfer functions used in the model
 *                    (defines fits file with tables), 0<= ntable <= 99
 * par15 ... nrad   - number of grid points in radius
 * par16 ... division - type of division in r integration
 *                      0 -> equidistant radial grid (constant linear step)
 *                      1 -> exponential radial grid (constant logarithmic step)
 * par17 ... nphi   - number of grid points in azimuth
 * par18 ... smooth - whether to smooth the resulting spectrum (0-no, 1-yes)
 * par19 ... Stokes - what should be stored in photar() array, i.e. as output
 *                    = -1 - the output is defined according to the XFLT0001 
 *                           keyword of the SPECTRUM extension of the data file,
 *                           where "Stokes:0" means photon number density flux,
 *                           "Stokes:1" means Stokes parameter Q devided by 
 *                           energy and "Stokes:2" means Stokes parameter U 
 *                           devided by energy
 *                    =  0 - array of photon number density flux per bin
 *                          (array of Stokes parameter I devided by energy)
 *                           with the polarisation computations switched off
 *                    =  1 - array of photon number density flux per bin
 *                          (array of Stokes parameter I devided by energy),
 *                           here, the polarisation computations are switched on
 *                           and different approximation for computed flux is 
 *                           used with non-isotropic emission directionality
 *                    =  2 - array of Stokes parameter Q devided by energy
 *                    =  3 - array of Stokes parameter U devided by energy
 *                    =  4 - array of Stokes parameter V devided by energy
 *                    =  5 - array of degree of polarization
 *                    =  6 - array of polarization angle psi=0.5*atan(U/Q)
 *                    =  7 - array of "Stokes" angle
 *                           beta=0.5*asin(V/sqrt(Q*Q+U*U+V*V))
 *                    =  8 - array of Stokes parameter Q devided by I
 *                    =  9 - array of Stokes parameter U devided by I
 *                    = 10 - array of Stokes parameter V devided by I
 * par20 ... chi0     - orientation of the system (-90 < chi0 < 90), 
 *                      the orientation angle (in degrees) of the system 
 *                      rotation axis with direction up, this angle is added to 
 *                      the computed polarisation angle at infinity, 
 *                      the orientation is degenarate by 180 degrees
 * par21 ... tau      - tau of the disc atmosphere, 
 *                    - tables created by Monte Carlo code Stokes for tau = 0.2, 
 *                      0.5, 1., 2., 5., 10.
 *                    - Chandrasekhar's relations for infinite optical depth for
 *                      tau > 10 (I in tables is actually the same already for 
 *                      tau=5. and Q for tau=10.)
 * par22 ... nthreads - number of threads to be used for computations
 *
 * NOTES:
 *  -> accuracy vs. speed trade off depends mainly on: nrad, nphi
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "fitsio.h"

/*******************************************************************************
*******************************************************************************/
#define NPARAM 22

#ifdef OUTSIDE_XSPEC

#define NE     200
#define E_MIN  0.01
#define E_MAX  20.
#define IFL    1

int main() {

void KYNBBphen(const double *ear, int ne, const double *param, int ifl, 
               double *photar, double *photer, const char* init);

double ear[NE + 1], param[NPARAM], photar[NE], photer[NE];
char   initstr[0] = "";
int    ie;

param[ 0] = 1.;         // a/M
param[ 1] = 30.;        // theta_o
param[ 2] = 1.;         // rin
param[ 3] = 1.;         // ms
param[ 4] = 1000.;       // rout
param[ 5] = 0.;         // phi
param[ 6] = 360.;       // dphi
param[ 7] = 1.;         // Tin
param[ 8] = 0.75;       // BBindex
param[ 9] = -3.;        // alpha
param[10] = 0.;         // beta
param[11] = 0.;         // rcloud
param[12] = 0.;         // zshift
param[13] = 80.;        // ntable
param[14] = 150.;       // nrad
param[15] = 1.;         // division
param[16] = 180.;       // nphi
param[17] = 0.;         // smooth
param[18] = 1.;         // Stokes
param[19] = 0.;         // chi0
param[20] = 11.;        // tau
param[21] = 4.;         // nthreads


for (ie = 0; ie <= NE; ie++) {
//  ear[ie] = E_MIN + ie * (E_MAX-E_MIN) / NE;
  ear[ie] = E_MIN * pow(E_MAX/E_MIN, ((double) ie) / NE);
}

KYNBBphen(ear, NE, param, IFL, photar, photer, initstr);

return(0);
}
#endif
/*******************************************************************************
*******************************************************************************/

#define GOOSMANN "goosmann.fits"
#define PI     3.14159265358979
#define H_KEVS 4.13566743e-18
#define C_MS   2.99792458e8
#define MSOLAR 1.989e+30
#define G      6.6743e-11
// the "kpc" below is 10kpc in cm
#define KPC    3.0857e+22
#define NCOSE0 21
#define ROUTMAX 1000.

/* Let's declare variables that are common for the main and emissivity 
   subroutines */
static double   *ener_loc, *flx, *I_loc, *Q_loc;
static float    *tau, *cose;
static double   Tin, BBindex, theta_o, rin, rout, tau0;
static int      polar;
static long int ntau, ncose;
// the following values are for mu_e = 0, 0.05, 0.1, ..., 0.9, 0.95, 1
static double I_l0[NCOSE0] = {0.18294,0.21613,0.24247,0.26702,0.29057,0.31350,
                              0.33599,0.35817,0.38010,0.40184,0.42343,0.44489,
                              0.46624,0.48750,0.50869,0.52981,0.55087,0.57189,
                              0.59286,0.61379,0.63469};
static double I_r0[NCOSE0] = {0.23147,0.25877,0.28150,0.30299,0.32381,0.34420,
                              0.36429,0.38417,0.40388,0.42346,0.44294,0.46233,
                              0.48165,0.50092,0.52013,0.53930,0.55844,0.57754,
                              0.59661,0.61566,0.63469};

extern char* FGMODF(void);
extern char* FGMSTR(char* dname);
extern int   xs_write(char* wrtstr, int idest);
extern float DGFILT(int ifl, const char* key);

void KYNBBphen(const double *ear, int ne, const double *param, int ifl, 
               double *photar, double *photer, const char* init) {

extern int ide(const double *ear, const int ne, const int nt, double *far, 
               double *qar, double *uar, double *var, 
               const double *ide_param, void (*emissivity)(), 
               const int ne_loc);

void emis_BBphen(double** ear_loc, const int ne_loc, const int nt, 
                 double *far_loc, double *qar_loc, double *uar_loc, 
                 double *var_loc, const double r, const double phi, 
                 const double cosmu, const double phiphoton, 
                 const double alpha_o, const double beta_o, 
                 const double delay, const double g);

void outer_disc_phen(const double *ear, const int ne, double *flux);

/* Let's declare static variables whose values should be remembered for next
   run in XSPEC */
static char   kydir[255] = "";
static char   pname[128] = "KYDIR";
static int    polar_old = -1, param_unchanged = -1, ne_old = -1; 
static float  *IQ;
static double param_old[NPARAM];
static double *far, *qar, *uar, *var;

FILE *fw;
double ide_param[25], flux_out[ne + 1], qar_final[ne], uar_final[ne],
       pd[ne], pa[ne], pa2[ne];
double I_l, I_r, Q_l, cose0, ttmp, ttmp1, y1, y2;
double pamin, pamax, pa2min, pa2max;
double am, am2, pom, pom1, pom2, pom3, rms, r_plus, chi0;
int    ne_loc, stokes, ie, irow, imin, imax, i0, itau0, orientation, iparam;
float  data_type;
char   data_type_c[8] = "Stokes";
int    outerdisc = 0;

// these are needed to work with a fits file...
fitsfile *fptr;
char     tables_file[255];
int      hdutype = 2;
int      colnum = 1;
long     frow = 1, felem = 1, nelems, nrow;
float    float_nulval = 0.;
int      nelements;
int      itau, icose, anynul, status = 0;//, maxdim=1000, naxis;

// Free memory for far, qar, uar and var if they has already been created and
// if their dimensions has changed...
if( ne_old != -1 && ne != ne_old ){
  free((void *) far); far = NULL;
  free((void *) qar); qar = NULL;
  free((void *) uar); uar = NULL;
  free((void *) var); var = NULL;
}
// Allocate memory for far, qar, uar and var...
if (ne_old == -1 || ne != ne_old )
if ((far = (double *) malloc( ne * sizeof(double))) == NULL ||
    (qar = (double *) malloc( ne * sizeof(double))) == NULL ||
    (uar = (double *) malloc( ne * sizeof(double))) == NULL ||
    (var = (double *) malloc( ne * sizeof(double))) == NULL ) {
  xs_write("kynphebb: Failed to allocate memory for Stokes arrays.", 5);
  for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
  return;
}
ne_old = ne;

// polar - whether we need value of change in polarization angle (0-no,1-yes)
stokes = (int) param[18];
if ((stokes < -1) || (stokes > 10)) {
  xs_write("kynphebb: Stokes has to be -1, 0-10", 5);
  for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
  return;
}
if(stokes == -1){
  data_type = DGFILT(ifl, data_type_c);
  if (data_type == 0. || data_type == 1. || data_type == 2.){
    stokes = 1 + (int) data_type;
  }
  else {
    xs_write("kynphebb: no or wrong information on data type (counts, q, u)", 5);
    xs_write("kynphebb: stokes = par19 = 1 (i.e. counts) will be used", 5);
    stokes=1;
  }
}

polar = 0;
if (stokes != 0) polar = 1;
ide_param[17] = polar;
chi0 = param[19]/180.*PI;
if (((chi0 < -90.) || (chi0 > 90.)) && polar) {
  xs_write("kynphebb: chi0 has to be between -90 and 90 degrees", 5);
  for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
  return;
}
//we can skip main part of computations if we change only Stokes parameter or 
//orientation of the system or number of threads for computations
if(param_unchanged != -1){
  iparam=0;
  param_unchanged=1;
  while( param_unchanged && iparam < NPARAM){
    if( ( iparam !=1 && iparam != 18 && iparam != 19 && iparam != 21
          && param[iparam] != param_old[iparam] ) ||
        ( iparam == 1 && fabs(param[1]) != fabs(param_old[1]) ) ||
        ( iparam == 18 && polar != polar_old && polar == 0 || polar_old == 0 ) )
      param_unchanged = 0;
    iparam++;
  }
  if( param_unchanged ) goto skip_computation;
}
param_unchanged=0;

// Let's initialize parameters for subroutine ide()
// a/M - black hole angular momentum
ide_param[0] = param[0];
am = param[0];
am2 = am * am;
r_plus= 1. + sqrt(1. - am2);
pom1 = pow(1. + am, 1. / 3.);
pom2 = pow(1. - am, 1. / 3.);
pom3 = pow(1. - am2, 1. / 3.);
pom = 1. + pom3 * (pom1 + pom2);
pom1 = sqrt(3. * am2 + pom * pom);
if (am >= 0) rms= 3. + pom1 - sqrt((3. - pom) * (3. + pom + 2. * pom1));
else rms = 3. + pom1 + sqrt((3. - pom) * (3. + pom + 2. * pom1));
// theta_o - observer inclination
ide_param[1] = fabs(param[1]);
theta_o = fabs(param[1]);
// rin - inner edge of non-zero disc emissivity
ide_param[2] = param[2];
// ms - whether to integrate from rin or rms
ide_param[3] = param[3];
// rout - outer edge of non-zero disc emissivity
ide_param[4] = param[4];
// rin, rout - inner, outer edge of non-zero disc emissivity
if( param[3] == 1. ){
  if( param[2] < rms ) rin = rms;
  else rin = param[2];
  rout = param[4];
}else if( param[3] == 2. ){
  rin  = param[2] * rms;
  rout = param[4] * rms;
}else if( param[3] == 0. ){
  rin  = param[2];
  rout = param[4];
}
if(rin  < r_plus) rin  = r_plus;
if(rout < r_plus) rout = r_plus;
//if rout >= 100*ROUTMAX than we assume rout to be infinity
//we then re-define it to ROUTMAX (or rin if rin>ROUTMAX)
//and compute the flux above ROUTMAX with outer_disc subroutine
if( rout >= 100*ROUTMAX ){
  rout = ROUTMAX;
  outerdisc = 1;
  if( param[3] == 1 || param[3] == 0 ){
    if(rout < param[2]) rout = param[2];
    ide_param[4] = rout;
  }
  else if( param[3] == 2 ){
    if(rout < param[2] * rms) rout = param[2] * rms;
    ide_param[4] = rout / rms;    
  }
}
// phi - lower azimuth of non-zero disc emissivity (deg)
ide_param[5] = param[5];
// dphi - (phi+dphi) is upper azimuth of non-zero disc emissivity (deg)
ide_param[6] = param[6];
// nrad - number of grid points in radius
ide_param[7] = param[14];
// division - type of division in r integration (0-equidistant, 1-exponential)
ide_param[8] = param[15];
// nphi - number of grid points in azimuth
ide_param[9] = param[16];
// smooth - whether to smooth the resulting spectrum (0-no, 1-yes)
ide_param[10] = param[17];
// normal - how to normalize the final spectrum
ide_param[11] = -1.;
// Tin - temperature in keV at the inner edge of the disc
Tin = param[7];
// BBindex - radial power-law index for radial dependence
// of black-body temperature
BBindex = param[8];
// zshift - overall Doppler shift
ide_param[12] = param[12];
// ntable - table model (defines fits file with tables)
ide_param[13] = param[13];
ne_loc = ne;
// edivision - type of division in local energies (0-equidistant, 1-exponential)
ide_param[14] = 1.;
// periodic and dt are not needed for nt = 1
// (ide_param[15], ide_param[16])
tau0 = param[20];
if ((tau0 < 0.2) && polar) {
  xs_write("kynphebb: tau has to be larger or equal to 0.2", 5);
  for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
  return;
}
// delay_r and delay_phi are not used
// (ide_param[18], ide_param[19])
// number of threads for multithread computations
ide_param[20] = param[21];
// alpha - position of the cloud in alpha impact parameter (in GM/c^2)
ide_param[21] = param[9];
// beta - position of the cloud in beta impact parameter (in GM/c^2)
ide_param[22] = param[10];
// rcloud - radius of the cloud (in GM/c^2)
ide_param[23] = param[11];
//whether the flux defined in emissivity subroutine is local one (0) or the 
//observed one (1)
ide_param[24] = 1.;

/******************************************************************************/
#ifdef OUTSIDE_XSPEC
// Let's write input parameters to a text file
fw = fopen("kynphebb.txt", "w");
fprintf(fw, "a/M          %12.6f\n", param[0]);
fprintf(fw, "theta_o      %12.6f\n", param[1]);
fprintf(fw, "rin          %12.6f\n", param[2]);
fprintf(fw, "ms           %12d\n", (int) param[3]);
fprintf(fw, "rout         %12.6f\n", param[4]);
fprintf(fw, "phi          %12.6f\n", param[5]);
fprintf(fw, "dphi         %12.6f\n", param[6]);
fprintf(fw, "Tin          %12.6f\n", param[7]);
fprintf(fw, "BBindex      %12.6f\n", param[8]);
fprintf(fw, "alpha        %12.6f\n", ide_param[21]);
fprintf(fw, "beta         %12.6f\n", ide_param[22]);
fprintf(fw, "rcloud       %12.6f\n", ide_param[23]);
fprintf(fw, "zshift       %12.6f\n", param[12]);
fprintf(fw, "ntable       %12d\n", (int) param[13]);
fprintf(fw, "nrad         %12d\n", (int) param[14]);
fprintf(fw, "division     %12d\n", (int) param[15]);
fprintf(fw, "nphi         %12d\n", (int) param[16]);
fprintf(fw, "smooth       %12d\n", (int) param[17]);
fprintf(fw, "Stokes       %12d\n", (int) param[18]);
fprintf(fw, "chi0         %12.6f\n", param[19]);
fprintf(fw, "tau          %12.6f\n", tau0);
fprintf(fw, "r_horizon    %12.6f\n", r_plus);
fprintf(fw, "r_ms         %12.6f\n", rms);
fprintf(fw, "edivision    %12d\n", (int) ide_param[14]);
fprintf(fw, "normal       %12.6f\n", ide_param[11]);
fprintf(fw, "nthreads     %12d\n", (int) ide_param[20]);
fclose(fw);
#endif
/******************************************************************************/

// initialize some variables needed for local flux defined in local energies
// Allocate memory for ener_loc and flx...
if ((ener_loc = (double *) malloc((ne_loc + 1) * sizeof(double))) == NULL) {
  xs_write("kynphebb: Failed to allocate memory for tmp arrays.", 5);
  for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
  return;
}
if ((flx = (double *) malloc((ne_loc + 1) * sizeof(double))) == NULL) {
  xs_write("kynphebb: Failed to allocate memory for tmp arrays.", 5);
  for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
  return;
}
for (ie = 0; ie <= ne_loc; ie++) {
  ener_loc[ie] = ear[ie];
  flx[ie] = 2 * pow(ener_loc[ie] / H_KEVS / C_MS / KPC, 2.) /
                H_KEVS * pow(G * MSOLAR / (C_MS * C_MS), 2.);
}

/******************************************************************************
// local spectrum output -- write ener_loc[] and flx[] into file:
fw = fopen("kynphebb_photar_loc.dat", "w");
for (ie = 0; ie < ne_loc; ie++)
  fprintf(fw, "%14.6f\t%E\n", ener_loc[ie], flx[ie]);
fclose(fw);
******************************************************************************/

/******************************************************************************/
if ( polar && polar_old == -1 ){
// Let's read the local polarisation tables
// The status parameter must always be initialized.
  status = 0;
// Open the FITS file for readonly access
// - if set try KYDIR directory, otherwise look in the working directory
//   or in the xspec directory where tables are usually stored...
  sprintf(kydir, "%s", FGMSTR(pname));
  if (strlen(kydir) == 0) sprintf(tables_file, "./%s", GOOSMANN);
  else if (kydir[strlen(kydir) - 1] == '/') sprintf(tables_file, "%s%s",
                                                    kydir, GOOSMANN);
  else sprintf(tables_file, "%s/%s", kydir, GOOSMANN);
// Let's read the 'goosmann.fits' file
// The status parameter must always be initialized.
  status = 0;
  ffopen(&fptr, tables_file, READONLY, &status);
  if (status) {
    sprintf(tables_file, "%s%s", FGMODF(), GOOSMANN);
    status = 0;
    ffopen(&fptr, tables_file, READONLY, &status);
  }
  if (status) {
    if (status) ffrprt(stderr, status);
    ffclos(fptr, &status);
    xs_write("\nkynphebb: set the KYDIR to the directory with the KY tables",5);
    for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
    return;
  }
// Let's read tables (binary tables => hdutype=2)
// Move to the extension 'r_horizon' and read its values
  ffmrhd(fptr, 1, &hdutype, &status);
  ffgnrw(fptr, &ntau, &status);
/******************************************************************************/
//  fprintf(stdout,"ntau = %ld\n",ntau);
/******************************************************************************/
// Allocate memory for tau...
  if ((tau = (float *) malloc(ntau * sizeof(float))) == NULL) {
    xs_write("kynphebb: Failed to allocate memory for tmp arrays.", 5);
    for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
    return;
  }
// Read the data in the 'tau' table
  nelems = ntau;
// FTGCV reads the VALUES from the first column.
  ffgcv(fptr, TFLOAT, colnum, frow, felem, nelems, &float_nulval, tau,
        &anynul, &status);
/******************************************************************************/
//  for ( itau=0; itau<ntau; itau++)fprintf(stdout,"%f\n",tau[itau]);
/******************************************************************************/
// Move to the extension 'cose' and read its values
  ffmrhd(fptr, 1, &hdutype, &status);
  ffgnrw(fptr, &ncose, &status);
/******************************************************************************/
//  fprintf(stdout,"ncose = %ld\n",ncose);
/******************************************************************************/
// Allocate memory for height...
  if ((cose = (float *) malloc(ncose * sizeof(float))) == NULL) {
    xs_write("kynphebb: Failed to allocate memory for tmp arrays.", 5);
    for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
    return;
  }
// Read the data in the 'cose' table
  nelems = ncose;
// FTGCV reads the VALUES from the first column.
  ffgcv(fptr, TFLOAT, colnum, frow, felem, nelems, &float_nulval, cose,
        &anynul, &status);
/******************************************************************************/
//  for ( icose=0; icose<ncose; icose ++)fprintf(stdout,"%f\n",cose[icose]);
/******************************************************************************/
// Let's read the tables for I and Q
// allocate memory for the arrays
  if ((IQ = (float *) malloc(ntau * ncose * 2 * sizeof(float))) == NULL) {
    xs_write("kynphebb: Failed to allocate memory for tmp arrays.", 5);
    for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
    return;
  }
// read the tables
  ffmrhd(fptr, 1, &hdutype, &status);
/* to read the file only once we have to read in blocks (all columns
   from the extension are put to buffer together)
   let's find out how many rows are going to be read into the buffer */
  ffgrsz(fptr, &nrow, &status);
//  if( nrow > ncose ) nrow = ncose;
  nelements = nrow * 2;
  for (irow = 0; irow < ncose; irow += nrow) {
//  the last block to read may be smaller:
    if ((ncose - irow) < nrow) nelements = (ncose - irow) * 2;
    for( itau=0; itau < ntau; itau++ )
      ffgcv(fptr, TFLOAT, itau+1, irow + 1, 1, nelements, &float_nulval, 
            &(IQ[itau*ncose*2 + irow*2]), &anynul, &status);
  }
// The FITS file must always be closed before exiting the program.
  ffclos(fptr, &status);
/*******************************************************************************
  itau=3;
  for ( icose=0; icose<ncose; icose++ ) 
    fprintf(stdout,"%d\t%f\t%e\t%e\n",icose+1,cose[icose],
    IQ[itau*ncose*2+icose*2], IQ[itau*ncose*2+icose*2+1]);
*******************************************************************************/
// We have to allocate memory for the arrays I_loc[] and Q_loc[]
  if ((I_loc = (double *) malloc(ncose * sizeof(double))) == NULL) {
    xs_write("kynphebb: Failed to allocate memory for tmp arrays.", 5);
    for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
    return;
  }
  if ((Q_loc = (double *) malloc(ncose * sizeof(double))) == NULL) {
    xs_write("kynphebb: Failed to allocate memory for tmp arrays.", 5);
    for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
    return;
  }
  polar_old = polar;
}
/******************************************************************************/

// given tau0, find the corresponding index in tau[] and compute I and Q:
if(polar && tau0 <= tau[ntau - 1]){
  imin = 0;
  imax = ntau;
  itau0 = ntau / 2;
  while ((imax - imin) > 1) {
    if (tau0 >= tau[itau0 - 1]) imin = itau0;
    else imax = itau0;
    itau0 = (imin + imax) / 2;
  }
  if (itau0 == 0) itau0 = 1;
//if ((imax == ntau) && (tau0 > tau[ntau - 1])) itau0 = ntau;
  ttmp = (tau0 - tau[itau0 - 1]) / (tau[itau0] - tau[itau0 - 1]);
  ttmp1 = 1. - ttmp;
  for( icose = 0; icose < ncose; icose++ ){
    y1 = IQ[(itau0-1)*ncose*2+icose*2];
    y2 = IQ[itau0*ncose*2+icose*2];
    I_loc[icose] = ttmp1 * y1 + ttmp * y2;
    y1 = IQ[(itau0-1)*ncose*2+icose*2+1];
    y2 = IQ[itau0*ncose*2+icose*2+1];
    Q_loc[icose] = ttmp1 * y1 + ttmp * y2;
  }
}

if (ide(ear, ne, 1, far, qar, uar, var, ide_param, emis_BBphen, ne_loc)) {
  for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
  return;
}
if( outerdisc ) outer_disc_phen(ear, ne, flux_out);

// interface with XSPEC
// final spectrum output -- write ear[] and photar[] into file:
if (!stokes){
  if( outerdisc )
    for (ie = 0; ie < ne; ie++) far[ie] += flux_out[ie];
} else {
  if( outerdisc ){
    cose0 = cos( theta_o / 180. * PI );
    if( tau0 <= tau[ntau - 1] ){
// Rene's tables
      imin = 1;
      imax = ncose;
      i0 = ( 1 + ncose) / 2;
      while ( ( imax - imin ) > 1 ){
        if( cose0 >= cose[i0] ) imin = i0;
        else imax = i0;
        i0 = ( imin + imax ) / 2;
      }
      I_l = ( I_loc[imin+1] - I_loc[imin] ) / 
            ( cose[imin+1] - cose[imin] ) * ( cose0 - cose[imin] ) + 
              I_loc[imin];
      Q_l = ( Q_loc[imin+1] - Q_loc[imin] ) / 
            ( cose[imin+1] - cose[imin] ) * ( cose0 - cose[imin] ) + 
              Q_loc[imin];
      for( ie=0; ie<ne; ie++ ){
        far[ie] += I_l * flux_out[ie];
        qar[ie] += Q_l * flux_out[ie];
      }
    } else{
// Chandrasekhar's formulae
      imin =  0;
      imax = 20;
      i0   = 10;
      while ( ( imax - imin ) > 1 ){
       if( cose0 >= i0 * 0.05 ) imin = i0;
       else imax = i0;
       i0 = ( imin + imax ) / 2;
      }
// let's interpolate the I_l and I_r between cos(theta_o)
      I_l = ( I_l0[imin+1] - I_l0[imin] ) / 0.05 * ( cose0 - imin * 0.05 ) +
              I_l0[imin];
      I_r = ( I_r0[imin+1] - I_r0[imin] ) / 0.05 * ( cose0 - imin * 0.05 ) +
              I_r0[imin];
      for( ie=0; ie<ne; ie++ ){
/* model with Rayleigh scattering according to the Chandrasekhar table XXIV
   overall flux integrated in emission angles (eq.96 in 68.6 chapter X in 
   Chandrasekhar) */
        far[ie] += ( I_l + I_r ) * flux_out[ie];
        qar[ie] += ( I_l - I_r ) * flux_out[ie];
      }
    }
  }
}

skip_computation:
if (!stokes){
  for (ie = 0; ie < ne; ie++) photar[ie] = far[ie];
} else{
// let's change the angle to opposite due to opposite system rotation
  if(param[1] >= 0.)orientation=1.;
  else orientation=-1.;
// let's change the orientation of the system 
  if(chi0 != 0.)
    for( ie=0; ie<ne; ie++ ){
      qar_final[ie] = qar[ie]*cos(2*chi0)-orientation*uar[ie]*sin(2*chi0);
      uar_final[ie] = orientation*uar[ie]*cos(2*chi0)+qar[ie]*sin(2*chi0);
    }
  else
    for( ie=0; ie<ne; ie++ ){
      qar_final[ie] = qar[ie];
      uar_final[ie] = orientation*uar[ie];
    }

  pamin = 1e30;
  pamax = -1e30;
  pa2min = 1e30;
  pa2max = -1e30;
  for (ie = ne - 1; ie >= 0; ie--) {
    pd[ie] = sqrt(qar_final[ie] * qar_final[ie] + uar_final[ie] * uar_final[ie] 
                  + var[ie] * var[ie]) / (far[ie] + 1e-99);
    pa[ie] = 0.5 * atan2(uar_final[ie], qar_final[ie]) / PI * 180.;
    if (ie < (ne - 1)) {
      while ((pa[ie] - pa[ie + 1]) > 90.) pa[ie] -= 180.;
      while ((pa[ie + 1] - pa[ie]) > 90.) pa[ie] += 180.;
    }
    if (pa[ie] < pamin) pamin = pa[ie];
    if (pa[ie] > pamax) pamax = pa[ie];
    pa2[ie] = 0.5 * asin(var[ie] / sqrt(qar_final[ie] * qar_final[ie] 
                         + uar_final[ie] * uar_final[ie] + var[ie] * var[ie] 
                         + 1e-99)) / PI * 180.;
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
      "%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\n", 
      0.5 * (ear[ie] + ear[ie+1]), far[ie] / (ear[ie+1] - ear[ie]), 
      qar_final[ie] / (ear[ie+1] - ear[ie]), 
      uar_final[ie] / (ear[ie+1] - ear[ie]), 
      var[ie] / (ear[ie+1] - ear[ie]), pd[ie], pa[ie], pa2[ie]);
//interface with XSPEC..........................................................
    if (stokes ==  1) photar[ie] = far[ie];
    if (stokes ==  2) photar[ie] = qar_final[ie];
    if (stokes ==  3) photar[ie] = uar_final[ie];
    if (stokes ==  4) photar[ie] = var[ie];
    if (stokes ==  5) photar[ie] = pd[ie] * (ear[ie + 1] - ear[ie]);
    if (stokes ==  6) photar[ie] = pa[ie] * (ear[ie + 1] - ear[ie]);
    if (stokes ==  7) photar[ie] = pa2[ie] * (ear[ie + 1] - ear[ie]);
    if (stokes ==  8) photar[ie] = qar_final[ie] / (far[ie]+1e-99) * (ear[ie + 1] - ear[ie]);
    if (stokes ==  9) photar[ie] = uar_final[ie] / (far[ie]+1e-99) * (ear[ie + 1] - ear[ie]);
    if (stokes == 10) photar[ie] = var[ie] / (far[ie]+1e-99) * (ear[ie + 1] - ear[ie]);
  }
  fclose(fw);
}

/******************************************************************************/
#ifdef OUTSIDE_XSPEC
// final spectrum output -- write ear[] and photar[] into file:
fw = fopen("kynphebb_photar.dat", "w");
if( outerdisc )
  for (ie = 0; ie < ne; ie++) {
    fprintf(fw, "%14.6f\t%E\t%E\t%E\n", 0.5 * (ear[ie] + ear[ie+1]),
      photar[ie] / (ear[ie + 1] - ear[ie]),
      far[ie] / (ear[ie+1] - ear[ie]), flux_out[ie] / (ear[ie + 1] - ear[ie]));
  }
else
  for (ie = 0; ie < ne; ie++) {
    fprintf(fw, "%14.6f\t%E\t%E\t%E\n", 0.5 * (ear[ie] + ear[ie+1]),
      photar[ie] / (ear[ie + 1] - ear[ie]));
  }
fclose(fw);
#endif
/******************************************************************************/

// Free memory from tmp arrays...
free((void *) ener_loc);
ener_loc = NULL;
free((void *) flx);
flx = NULL;

// Store old parameter values
for (iparam = 0; iparam < NPARAM; iparam++) param_old[iparam] = param[iparam];

return;
}

/*******************************************************************************
*******************************************************************************/

void emis_BBphen(double** ear_loc, const int ne_loc, const int nt, 
                 double *far_loc, double *qar_loc, double *uar_loc, 
                 double *var_loc, const double r, const double phi, 
                 const double cosmu, const double phiphoton, 
                 const double alpha_o, const double beta_o, 
                 const double delay, const double g) {

double flx0, flx1, g2, I_l, Q_l, I_r;
int    ie, imin, imax, i0;

*ear_loc = ener_loc;
g2 = g * g;
flx0 = flx[0] / (exp(*(*ear_loc) /
                (g * Tin * pow(r / rin, -BBindex))) - 1.);
for (ie = 0; ie < ne_loc; ie++) {
  flx1 = flx[ie+1] / (exp(*(*ear_loc + ie + 1) /
                (g * Tin * pow(r / rin, -BBindex))) - 1.);
  far_loc[ie] = (flx0 + flx1) / 2. *
                  ( *(*ear_loc + ie +1) - *(*ear_loc + ie) ) / g2;
  flx0 = flx1;
}
if (polar) {
  if( tau0 <= tau[ntau - 1] ){
// Rene's tables
    imin = 1;
    imax = ncose;
    i0 = ( 1 + ncose) / 2;
    while ( ( imax - imin ) > 1 ){
      if( cosmu >= cose[i0] ) imin = i0;
      else imax = i0;
      i0 = ( imin + imax ) / 2;
    }
    I_l = ( I_loc[imin+1] - I_loc[imin] ) / 
          ( cose[imin+1] - cose[imin] ) * ( cosmu - cose[imin] ) + 
            I_loc[imin];
    Q_l = ( Q_loc[imin+1] - Q_loc[imin] ) / 
          ( cose[imin+1] - cose[imin] ) * ( cosmu - cose[imin] ) + 
            Q_loc[imin];
    for( ie = 0; ie < ne_loc; ie++ ){
      qar_loc[ie] = Q_l * far_loc[ie];
      uar_loc[ie] = 0.;
      var_loc[ie] = 0.;
      far_loc[ie] *= I_l;
    }
  } else{
// Chandrasekhar's formulae
    imin =  0;
    imax = 20;
    i0   = 10;
    while ( ( imax - imin ) > 1 ){
     if( cosmu >= i0 * 0.05 ) imin = i0;
     else imax = i0;
     i0 = ( imin + imax ) / 2;
    }
// let's interpolate the I_l and I_r between cos(theta_o)
    I_l = ( I_l0[imin+1] - I_l0[imin] ) / 0.05 * ( cosmu - imin * 0.05 ) +
            I_l0[imin];
    I_r = ( I_r0[imin+1] - I_r0[imin] ) / 0.05 * ( cosmu - imin * 0.05 ) +
            I_r0[imin];
    for( ie = 0; ie < ne_loc; ie++ ){
/* model with Rayleigh scattering according to the Chandrasekhar table XXIV
   overall flux integrated in emission angles (eq.96 in 68.6 chapter X in 
   Chandrasekhar) */
      qar_loc[ie] = ( I_l - I_r ) * far_loc[ie];
      uar_loc[ie] = 0.;
      var_loc[ie] = 0.;
      far_loc[ie] *= ( I_l + I_r );
    }
  }
}
return;
}

/*******************************************************************************
*******************************************************************************/

void outer_disc_phen(const double *ear, const int ne, double *flux) {

double y0[23] = {1.932161609, 1.770917995, 1.497175205, 1.207475655,
                 0.9408507085, 0.7133969015, 0.5289516748, 0.3848869662,
                 0.2756142188, 0.1946729206, 0.1358796849, 0.0938688542,
                 0.0642642844, 0.0436488382, 0.02943951894, 0.01973255119,
                 0.01315284209, 0.008723435269, 0.00575970436, 0.003787405908,
                 0.002481262304, 0.00162006534, 0.00105449408};
double Tout, norm, x, y;
int    ie, imin, imax, i0;

// the factor rg^2 is due to rin which is in geometrical units!!!
norm = 4. * PI * cos(theta_o / 180. * PI) * rin * rin * pow(Tin, 2. / BBindex) /
       H_KEVS / pow(H_KEVS * C_MS * KPC, 2.) / BBindex * pow(G * MSOLAR /
       (C_MS * C_MS), 2.);
Tout = Tin * pow(rout / rin, -BBindex);
for (ie = 0; ie <= ne; ie++) {
  x = ear[ie] / Tout;
  if (x > 11.) flux[ie] = 0.;
  else {
    imin = 0;
    imax = 22;
    i0 = 11;
    while ((imax - imin) > 1) {
      if (x >= i0 * 0.5) imin = i0;
      else imax = i0;
      i0 = (imin + imax) / 2;
    }
    y = (y0[imin + 1] - y0[imin]) / 0.5 * (x - imin * 0.5) + y0[imin];
    flux[ie] = pow(ear[ie], 2. * (BBindex - 1) / BBindex) * y;
  }
}
for (ie = 0; ie < ne; ie++) flux[ie] = norm * (flux[ie] + flux[ie+1]) /
                                       2. * (ear[ie+1] - ear[ie]);
return;
}
