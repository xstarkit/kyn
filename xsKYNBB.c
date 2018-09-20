/* KYNBB - black-body emission - non-axisymmetric version
 *         model subroutine for XSPEC
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
 * This model computes the black-body emission from an accretion disc around
 * a black hole. The local flux is defined as flux ~ E^2/(exp(E/kT)-1),
 * where the temperature T is defined according to Novikov-Thorne.
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
 *                     properties (Stokes parameter, par20, larger than 1)
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
 *                 - if outer edge is equal or larger than 1000 GM/c^2 then 
 *                   the emission from above this radius is added
 * par6  ... phi   - lower azimuth of non-zero disc emissivity (deg)
 * par7  ... dphi  - (phi + dphi) is upper azimuth of non-zero disc emissivity
 *                   0 <= dphi <= 360  (deg)
 * par8  ... BHmass  - the black hole mass in units of Solar mass
 * par9  ... arate  - accretion rate in units of Solar mass per Julian year
 *                    (365.25days)
 * par10 ... f_col  - spectral hardening factor
 * par11 ... alpha  - position of the cloud centre in GM/c^2 in alpha coordinate
 *                    (alpha being the impact parameter in phi direction, 
 *                     positive for approaching side of the disc)
 * par12 ... beta   - position of the cloud centre in GM/c^2 in beta coordinate
 *                    (beta being the impact parameter in theta direction, 
 *                     positive in up direction, i.e. above the disc)
 * par13 ... rcloud - radius of the obscuring cloud (in GM/c^2)
 *                  - if negative, only the emission transmitted through
 *                    the cloud is taken into account
 * par14 ... zshift - overall Doppler shift
 * par15 ... ntable - table of relativistic transfer functions used in the model
 *                    (defines fits file with tables), 0<= ntable <= 99
 * par16 ... nrad   - number of grid points in radius
 * par17 ... division - type of division in r integration
 *                      0 -> equidistant radial grid (constant linear step)
 *                      1 -> exponential radial grid (constant logarithmic step)
 * par18 ... nphi   - number of grid points in azimuth
 * par19 ... smooth - whether to smooth the resulting spectrum (0-no, 1-yes)
 * par20 ... Stokes - what should be stored in photar() array, i.e. as output
 *                    = 0 - array of photon number density flux per bin
 *                         (array of Stokes parameter I devided by energy)
 *                          with the polarisation computations switched off
 *                    = 1 - array of photon number density flux per bin
 *                         (array of Stokes parameter I devided by energy),
 *                          here, the polarisation computations are switched on
 *                          and different approximation for computed flux is 
 *                          used with non-isotropic emission directionality
 *                    = 2 - array of Stokes parameter Q devided by energy
 *                    = 3 - array of Stokes parameter U devided by energy
 *                    = 4 - array of Stokes parameter V devided by energy
 *                    = 5 - array of degree of polarization
 *                    = 6 - array of polarization angle psi=0.5*atan(U/Q)
 *                    = 7 - array of "Stokes" angle
 *                          beta=0.5*asin(V/sqrt(Q*Q+U*U+V*V))
 * par21 ... chi0     - orientation of the system (-90 < chi0 < 90), 
 *                      the orientation angle (in degrees) of the system 
 *                      rotation axis with direction up, this angle is added to 
 *                      the computed polarisation angle at infinity, 
 *                      the orientation is degenarate by 180 degrees
 * par22 ... tau      - tau of the disc atmosphere, 
 *                    - tables created by Monte Carlo code Stokes for tau = 0.2, 
 *                      0.5, 1., 2., 5., 10.
 *                    - Chandrasekhar's relations for infinite optical depth for
 *                      tau > 10 (I in tables is actually the same already for 
 *                      tau=5. and Q for tau=10.)
 * par23 ... nthreads - number of threads to be used for computations
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

#define NE     200
#define E_MIN  0.01
#define E_MAX  20.
#define NPARAM 22
#define IFL    1

int main() {

void KYNBB(const double *ear, int ne, const double *param, int ifl, 
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
param[ 7] = 3.;         // BHmass
param[ 8] = 1.e-9;       // arate
param[ 9] = 1.7;        // f_col
param[10] = -3.;        // alpha
param[11] = 0.;         // beta
param[12] = 0.;         // rcloud
param[13] = 0.;         // zshift
param[14] = 80.;        // ntable
param[15] = 150.;       // nrad
param[16] = 1.;         // division
param[17] = 180.;       // nphi
param[18] = 0.;         // smooth
param[19] = 1.;         // Stokes
param[20] = 0.;         // chi0
param[21] = 11.;        // tau
param[22] = 4.;         // nthreads

for (ie = 0; ie <= NE; ie++) {
//  ear[ie] = E_MIN + ie * (E_MAX-E_MIN) / NE;
  ear[ie] = E_MIN * pow(E_MAX/E_MIN, ((double) ie) / NE);
}

KYNBB(ear, NE, param, IFL, photar, photer, initstr);

return(0);
}

#endif
/*******************************************************************************
*******************************************************************************/

#define GOOSMANN "goosmann.fits"
#define PI     3.14159265358979
#define H_KEVS 4.13566743e-18
#define K_KEVK 8.6174e-8
#define C_MS   2.99792458e8
#define MSOLAR 1.989e+30
#define G      6.6743e-11
#define SIGMA  5.6704e-8
#define YEAR   31557600.0
// the "kpc" below is 10kpc in cm
#define KPC    3.0857e+22
#define NCOSE0 21
#define ROUTMAX 1000.

/* Let's declare variables that are common for the main and emissivity 
   subroutines */
static double   *ener_loc, *flx;
static float    *tau, *cose, *I_loc, *Q_loc;
static double   am, x0 , x1, x2, x3;
static double   arate, BHmass, f_col, rout, Tout, Tnorm, theta_o, tau0;
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

extern int xs_write(char* wrtstr, int idest);

void KYNBB(const double *ear, int ne, const double *param, int ifl, 
           double *photar, double *photer, const char* init) {

extern int ide(const double *ear, const int ne, const int nt, double *far, 
               double *qar, double *uar, double *var, 
               const double *ide_param, void (*emissivity)(), 
               const int ne_loc);

void emis_BB(double** ear_loc, const int ne_loc, const int nt, 
             double *far_loc, double *qar_loc, double *uar_loc, 
             double *var_loc, const double r, const double phi, 
             const double cosmu, const double phiphoton, 
             const double alpha_o, const double beta_o, 
             const double delay, const double g);

void outer_disc(const double *ear, const int ne, double *flux);

/* Let's declare static variables whose values should be remembered for next
   run in XSPEC */
static char   kydir[255] = "";
static char   pname[128] = "KYDIR";
static int    polar_old = -1; 
static float  *IQ;

FILE *fw;
double ide_param[25], flux_out[ne + 1];
double far[ne], qar[ne], uar[ne], var[ne], pd[ne], pa[ne], pa2[ne];
double I_l, I_r, Q_l, cose0, ttmp, ttmp1, y1, y2, tmp;
double pamin, pamax, pa2min, pa2max;
double am2, pom, pom1, pom2, pom3, rms, r_plus, x, Ccal, Lcal, arcosa3,
       orientation, chi0;
int    ne_loc, stokes, ie, irow, imin, imax, i0, itau0;

// these are needed to work with a fits file...
fitsfile *fptr;
char     tables_file[255];
int      hdutype = 2;
int      colnum = 1;
long     frow = 1, felem = 1, nelems, nrow;
float    float_nulval = 0.;
int      nelements;
int      itau, icose, anynul, status = 0;//, maxdim=1000, naxis;

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
x0 = sqrt(rms);
arcosa3 = acos(am) / 3.;
x1 = 2 * cos(arcosa3 - PI / 3.);
x2 = 2 * cos(arcosa3 + PI / 3.);
x3 = -2 * cos(arcosa3);
// theta_o - observer inclination
orientation = 1.;
if(param[1] < 0.) orientation = -1.;
ide_param[1] = fabs(param[1]);
theta_o = fabs(param[1]);
// rin - inner edge of non-zero disc emissivity
ide_param[2] = param[2];
// ms - whether to integrate from rin or rms
ide_param[3] = param[3];
// rout - outer edge of non-zero disc emissivity
ide_param[4] = param[4];
if( param[3] == 1 || param[3] == 0 ) rout = param[4];
else if( param[3] == 2 ) rout = param[4] * rms;
else rout=0.;
if( rout > ROUTMAX ) rout = ROUTMAX;
// phi - lower azimuth of non-zero disc emissivity (deg)
ide_param[5] = param[5];
// dphi - (phi+dphi) is upper azimuth of non-zero disc emissivity (deg)
ide_param[6] = param[6];
// nrad - number of grid points in radius
ide_param[7] = param[15];
// division - type of division in r integration (0-equidistant, 1-exponential)
ide_param[8] = param[16];
// nphi - number of grid points in azimuth
ide_param[9] = param[17];
// smooth - whether to smooth the resulting spectrum (0-no, 1-yes)
ide_param[10] = param[18];
// normal - how to normalize the final spectrum
ide_param[11] = -1.;
// BHmass - the black hole mass in Msolar
BHmass = param[7];
// arate - accretion rate in units of 1 solar mass per Julian year (365.25days)
arate = param[8];
// f_col - spectral hardening factor
f_col = param[9];
// zshift - overall Doppler shift
ide_param[12] = param[13];
// ntable - table model (defines fits file with tables)
ide_param[13] = param[14];
// set the outer temperature Tout (its normalization will be set later)
x = sqrt(rout);
// Bcal = 1. + am / (x*x*x);
Ccal = 1. - 3. / (x*x) + 2. * am / (x*x*x);
if (am < 1.)
  Lcal = 1. / x * (x - x0 - 1.5 * am * log(x / x0) -
         3. * pow(x1 - am,2.)/x1/(x1 - x2)/(x1 - x3) * log((x - x1)/(x0 - x1))-
         3. * pow(x2 - am,2.)/x2/(x2 - x1)/(x2 - x3) * log((x - x2)/(x0 - x2))-
         3. * pow(x3 - am,2.)/x3/(x3 - x1)/(x3 - x2) * log((x - x3)/(x0 - x3)));
else Lcal = 1. / x * (x - 1 + 1.5 * log((x + 2.) / 3. / x));
Tout = pow(x, -1.5) * pow(arate, 0.25) *
       pow(BHmass, -0.5) * pow(Lcal / Ccal, 0.25);
// Tout = pow(x, -1.5) * pow(1 - x0 / x, 0.25) *
//          pow(arate, 0.25) * pow(BHmass, -0.5);
// Let's define local energy and local flux
ne_loc = ne;
// edivision - type of division in local energies (0-equidistant, 1-exponential)
ide_param[14] = 1.;
// periodic and dt are not needed for nt = 1
// (ide_param[15], ide_param[16])
// polar - whether we need value of change in polarization angle (0-no,1-yes)
stokes = (int) param[19];
if ((stokes < 0) || (stokes > 7)) {
  xs_write("kynbb: Stokes has to be 0-7", 5);
  for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
  return;
}
polar = 0;
if (stokes > 0) polar = 1;
ide_param[17] = polar;
chi0 = param[20]/180.*PI;
if (((chi0 < -90.) || (chi0 > 90.)) && polar) {
  xs_write("kynbb: chi0 has to be between -90 and 90 degrees", 5);
  for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
  return;
}
tau0 = param[21];
if ((tau0 < 0.2) && polar) {
  xs_write("kynbb: tau has to be larger or equal to 0.2", 5);
  for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
  return;
}
// delay_r and delay_phi are not used
// (ide_param[18], ide_param[19])
// number of threads for multithread computations
ide_param[20] = param[22];
// alpha - position of the cloud in alpha impact parameter (in GM/c^2)
ide_param[21] = param[10];
// beta - position of the cloud in beta impact parameter (in GM/c^2)
ide_param[22] = param[11];
// rcloud - radius of the cloud (in GM/c^2)
ide_param[23] = param[12];
//whether the flux defined in emissivity subroutine is local one (0) or the 
//observed one (1)
ide_param[24] = 1.;

/******************************************************************************/
#ifdef OUTSIDE_XSPEC
// Let's write input parameters to a text file
fw = fopen("kynbb.txt", "w");
fprintf(fw, "a/M          %12.6f\n", param[0]);
fprintf(fw, "theta_o      %12.6f\n", param[1]);
fprintf(fw, "rin          %12.6f\n", param[2]);
fprintf(fw, "ms           %12d\n", (int) param[3]);
fprintf(fw, "rout         %12.6f\n", param[4]);
fprintf(fw, "phi          %12.6f\n", param[5]);
fprintf(fw, "dphi         %12.6f\n", param[6]);
fprintf(fw, "BBmass       %12.6f\n", param[7]);
fprintf(fw, "arate        %12.6f\n", param[8]);
fprintf(fw, "f_col        %12.6f\n", param[9]);
fprintf(fw, "alpha      %12.6f\n", ide_param[21]);
fprintf(fw, "beta       %12.6f\n", ide_param[22]);
fprintf(fw, "rcloud       %12.6f\n", ide_param[23]);
fprintf(fw, "zshift       %12.6f\n", param[13]);
fprintf(fw, "ntable       %12d\n", (int) param[14]);
fprintf(fw, "nrad         %12d\n", (int) param[15]);
fprintf(fw, "division     %12d\n", (int) param[16]);
fprintf(fw, "nphi         %12d\n", (int) param[17]);
fprintf(fw, "smooth       %12d\n", (int) param[18]);
fprintf(fw, "Stokes       %12d\n", (int) param[19]);
fprintf(fw, "chi0         %12.6f\n", param[20]);
fprintf(fw, "tau          %12.6f\n", tau0);
fprintf(fw, "r_horizon    %12.6f\n", r_plus);
fprintf(fw, "r_ms         %12.6f\n", rms);
fprintf(fw, "edivision    %12d\n", (int) ide_param[14]);
fprintf(fw, "ne_loc       %12d\n", ne_loc);
fprintf(fw, "normal       %12.6f\n", ide_param[11]);
fprintf(fw, "nthreads     %12d\n", (int) ide_param[20]);
fclose(fw);
#endif
/******************************************************************************/

// initialize some variables needed for local flux defined in local energies
// Allocate memory for ener_loc and flx...
if ((ener_loc = (double *) malloc((ne_loc + 1) * sizeof(double))) == NULL) {
  xs_write("kynbb: Failed to allocate memory for tmp arrays.", 5);
  for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
  return;
}
if ((flx = (double *) malloc((ne_loc + 1) * sizeof(double))) == NULL) {
  xs_write("kynbb: Failed to allocate memory for tmp arrays.", 5);
  for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
  return;
}
for (ie = 0; ie <= ne_loc; ie++) {
  ener_loc[ie] = ear[ie];
// the factor rg^2 is due to the fact that we integrate in dS = r * dr * dphi
// which we do in radius in geometrical units!!!
  flx[ie] = 2 * pow(ener_loc[ie] / H_KEVS / C_MS / KPC, 2.) / H_KEVS /
            pow(f_col, 4.) * pow(G * MSOLAR * BHmass / (C_MS * C_MS), 2.);
}
Tnorm = f_col * K_KEVK * pow(C_MS, 1.5) *
        pow(3. / (8. * PI * G * G * MSOLAR * YEAR * SIGMA), 0.25);

/******************************************************************************
// local spectrum output -- write ener_loc[] and flx[] into file:
fw = fopen("kynbb_photar_loc.dat", "w");
for (ie = 0; ie < ne_loc; ie++)
  fprintf(fw, "%14.6f\t%E\n", ener_loc[ie], flx[ie]);
fclose(fw);
******************************************************************************/

/******************************************************************************/
if ( polar & polar_old == -1 ){
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
    xs_write("\nkynbb: set the KYDIR to the directory with the KY tables",5);
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
    xs_write("kynbb: Failed to allocate memory for tmp arrays.", 5);
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
    xs_write("kynbb: Failed to allocate memory for tmp arrays.", 5);
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
    xs_write("kynbb: Failed to allocate memory for tmp arrays.", 5);
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
    xs_write("kynbb: Failed to allocate memory for tmp arrays.", 5);
    for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
    return;
  }
  if ((Q_loc = (double *) malloc(ncose * sizeof(double))) == NULL) {
    xs_write("kynbb: Failed to allocate memory for tmp arrays.", 5);
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

if (ide(ear, ne, 1, far, qar, uar, var, ide_param, emis_BB, ne_loc)) {
  for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
  return;
}
if( rout == ROUTMAX ) outer_disc(ear, ne, flux_out);

// interface with XSPEC
// final spectrum output -- write ear[] and photar[] into file:
if (!stokes){
  if( rout == ROUTMAX ) 
    for (ie = 0; ie < ne; ie++) photar[ie] = far[ie] + flux_out[ie];
  else for (ie = 0; ie < ne; ie++) photar[ie] = far[ie];
} else {
  if( rout == ROUTMAX ){
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

// let's change the angle to opposite due to opposite system rotation
  if(orientation == -1.) 
    for( ie=0; ie<ne; ie++ )
      uar[ie] = -uar[ie];
// let's add the angle due to the orientation of the system 
  if(chi0 != 0.)
    for( ie=0; ie<ne; ie++ ){
      tmp = qar[ie];
      qar[ie] = qar[ie]*cos(2*chi0)-uar[ie]*sin(2*chi0);
      uar[ie] = uar[ie]*cos(2*chi0)+tmp*sin(2*chi0);
    }

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
      "%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\n", 
      0.5 * (ear[ie] + ear[ie+1]), far[ie] / (ear[ie+1] - ear[ie]), 
      qar[ie] / (ear[ie+1] - ear[ie]), uar[ie] / (ear[ie+1] - ear[ie]), 
      var[ie] / (ear[ie+1] - ear[ie]), pd[ie], pa[ie], pa2[ie]);
//interface with XSPEC..........................................................
    if (stokes == 1) photar[ie] = far[ie];
    if (stokes == 2) photar[ie] = qar[ie];
    if (stokes == 3) photar[ie] = uar[ie];
    if (stokes == 4) photar[ie] = var[ie];
    if (stokes == 5) photar[ie] = pd[ie] * (ear[ie + 1] - ear[ie]);
    if (stokes == 6) photar[ie] = pa[ie] * (ear[ie + 1] - ear[ie]);
    if (stokes == 7) photar[ie] = pa2[ie] * (ear[ie + 1] - ear[ie]);
  }
  fclose(fw);
}

/******************************************************************************/
#ifdef OUTSIDE_XSPEC
// final spectrum output -- write ear[] and photar[] into file:
fw = fopen("kynbb_photar.dat", "w");
if( rout == ROUTMAX )
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

return;
}

/*******************************************************************************
*******************************************************************************/

void emis_BB(double** ear_loc, const int ne_loc, const int nt, 
             double *far_loc, double *qar_loc, double *uar_loc, 
             double *var_loc, const double r, const double phi, 
             const double cosmu, const double phiphoton, 
             const double alpha_o, const double beta_o, 
             const double delay, const double g) {

double temp, x, Ccal, Lcal, flx0, flx1, g2, I_l, Q_l, I_r;
int    ie, imin, imax, i0;

*ear_loc = ener_loc;
x = sqrt(r);
// Bcal = 1. + am / pow(x, 3.);
Ccal = 1. - 3. / (x*x) + 2. * am / (x*x*x);
if (am < 1.) 
  Lcal = 1. / x * (x - x0 - 1.5 * am * log(x / x0) - 
         3. * pow(x1 - am,2.)/x1/(x1 - x2)/(x1 - x3) * log((x - x1)/(x0 - x1))-
         3. * pow(x2 - am,2.)/x2/(x2 - x1)/(x2 - x3) * log((x - x2)/(x0 - x2))-
         3. * pow(x3 - am,2.)/x3/(x3 - x1)/(x3 - x2) * log((x - x3)/(x0 - x3)));
else Lcal = 1. / x * (x - 1 + 1.5 * log((x + 2.) / 3. / x));
temp = Tnorm * pow(x, -1.5) * pow(arate, 0.25) *
       pow(BHmass, -0.5) * sqrt(sqrt(Lcal / Ccal));
temp *= g;
g2 = g * g;
flx0 = flx[0] / ( exp( *(*ear_loc) / temp ) - 1. );
for (ie = 0; ie < ne_loc; ie++) {
  flx1 = flx[ie+1] / ( exp( *(*ear_loc + ie + 1) / temp ) - 1. );
  far_loc[ie] = ( flx0 + flx1 ) / 2. * 
                  ( *(*ear_loc + ie + 1) - *(*ear_loc + ie) ) / g2 ;
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

void outer_disc(const double *ear, const int ne, double *flux) {

double y0[23] = {1.932161609, 1.770917995, 1.497175205, 1.207475655,
                 0.9408507085, 0.7133969015, 0.5289516748, 0.3848869662,
                 0.2756142188, 0.1946729206, 0.1358796849, 0.0938688542,
                 0.0642642844, 0.0436488382, 0.02943951894, 0.01973255119,
                 0.01315284209, 0.008723435269, 0.00575970436, 0.003787405908,
                 0.002481262304, 0.00162006534, 0.00105449408};
double norm, x, y;
int    ie, imin, imax, i0;

Tnorm = f_col * K_KEVK * pow(C_MS, 1.5) *
        pow(3. / (8. * PI * G * G * MSOLAR * YEAR * SIGMA), 0.25);
Tout = Tnorm * Tout;
// the factor rg^2 is due to rout which is in geometrical units!!!
norm = 16. * PI * cos(theta_o / 180. * PI) * rout * rout * pow(Tout, 8. / 3.) /
       H_KEVS / pow(H_KEVS * C_MS* KPC, 2.) / 3. / pow(f_col, 4.) *
       pow(G * MSOLAR * BHmass / (C_MS * C_MS), 2.);
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
    flux[ie] = pow(ear[ie], -2./3.) * y;
  }
}
for (ie = 0; ie < ne; ie++) flux[ie] = norm * (flux[ie] + flux[ie+1]) /
                                       2. * (ear[ie+1] - ear[ie]);
return;
}
