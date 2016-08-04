/* KYNlpcr - reflection
 * C version of Fortran77 model subroutine for xspec (supplemented by a main
 * driver).
 *
 * ref. Dovciak M., Karas V., Martocchia A., Matt G., Yaqoob T. (2004)
 * -------------------------------------------------------------------
 * REFERENCES:
 *
 * Dovciak M., Karas V., Martocchia A., Matt G. & Yaqoob T. (2004). XSPEC model
 * to explore spectral features from black hole sources.
 * In Proc. of the workshop on processes in the vicinity of black holes
 * and neutron stars. S.Hledik & Z.Stuchlik, Opava. In press. [astro-ph/0407330]
 *
 * Dovciak M., Karas V. & Yaqoob, T. (2004). An extended scheme
 * for fitting X-ray data with accretion disk spectra in the strong gravity
 * regime. ApJS, 153, 205.
 *
 * Dovciak M. (2004). Radiation of accretion discs in strong gravity. Faculty of
 * Mathematics and Physics, Charles University, Prague. PhD thesis.
 * [astro-ph/0411605]
 * -------------------------------------------------------------------
 *
 * This subroutine computes the emission from an acrretion disc that is
 * illuminated from the primary power-law source located on the axis above the
 * central black hole. All relativistic effects are taken into account (in both
 * parts of the light path - from the primary source to the disc and from the
 * disc to the observer). This subroutine calls subroutine ide() for integrating
 * local emission over the disc and uses the fits file 'KBHtablesNN.fits'
 * defining the transfer functions needed for the integration. For details on
 * ide() and the fits file see the subroutine ide(). This subroutine uses
 * KBHlamp.fits file with transfer functions for the light coming from primary
 * source and hitting the accretion disc (see the description of this file
 * below). The reflection is modelled by Monte Carlo simulations of Compton
 * scattering which is precalculated in the reflspectra.fits file (see the
 * description of this file below).
 *
 * par1  ... a/M     - black hole angular momentum (-1 <= a/M <= 1)
 * par2  ... theta_o - observer inclination in degrees (0-pole, 90-disc)
 * par3  ... rin - inner edge of non-zero disc emissivity (in GM/c^2 or in 
 *                 r_mso)
 * par4  ... ms  - 0 - we integrate from inner edge = par3 
 *                 1 - if the inner edge of the disc is below marginally stable
 *                     orbit then we integrate emission above MSO only
 *                 2 - we integrate from inner edge given in units of MSO, i.e.
 *                     inner edge = par3 * r_mso (the same applies for outer 
 *                     edge)
 * par5  ... rout  - outer edge of non-zero disc emissivity (in GM/c^2 or in 
 *                   r_mso)
 * par6  ... phi   - lower azimuth of non-zero disc emissivity (deg)
 * par7  ... dphi  - (phi + dphi) is upper azimuth of non-zero disc emissivity
 *                   0 <= dphi <= 360  (deg)
 * par8 ... M/M8   - black hole mass in units of 10^8 solar masses
 * par9  ... height   - height on the axis (measured from the center) at which
 *                      the primary source is located (GM/c^2)
 *                      DO NOT USE THE MODEL WITH NEGATIVE HEIGHT!!! 
 * par10 ... PhoIndex - power-law energy index of the primary flux
 * par11 ... L/Ledd   - dE/dt primary isotropic flux in Ledd
 * par12 ... NpNr   - ratio of the primary to the reflected normalization
 *                    1 - self-consistent model for isotropic primary source
 *                    0 - only reflection, primary source is hidden
 * par13 ... line  - whether to include lines and/or reflection continuum in
 *                   the spectra
 *                   0 - only continuum
 *                   1 - Kalpha Fe line with continuum
 *                   2 - Kalpha and Kbeta lines with continuum
 *                   3 - all lines computed by NOAR for neutral disc
 *                  -1 - only Kalpha Fe line without the reflection continuum
 *                  -2 - Kalpha and Kbeta lines without the reflection continuum
 *                  -3 - all lines computed by NOAR for neutral disc without
 *                       the reflection continuum
 * par14 ... Ec     - energy cut-off
 * par15 ... alpha  - position of the cloud centre in GM/c^2 in alpha coordinate
 *                    (alpha being the impact parameter in phi direction, 
 *                     positive for approaching side of the disc)
 * par16 ... beta   - position of the cloud centre in GM/c^2 in beta coordinate
 *                    (beta being the impact parameter in theta direction, 
 *                     positive in up direction, i.e. away from the disc)
 * par17 ... rcloud - radius of the obscuring cloud
 * par18 ... zshift - overall Doppler shift
 * par19 ... ntable - defines fits file with tables (0 <= ntable <= 99)
 *                    by default only the tables with ntable=80 are correct
 *                    for this model
 * par20 ... nrad   - number of grid points in radius
 * par21 ... division - type of division in r integration
 *                      (0-equidistant, 1-exponential)
 * par22 ... nphi   - number of grid points in azimuth
 * par23 ... smooth - whether to smooth the resulting spectrum (0-no, 1-yes)
 * par24 ... Stokes - what should be stored in photar() array
 *                    = 0 - photon number density flux per bin
 *                         (Stokes parameter I devided by energy)
 *                    = 1 - Stokes parameter Q devided by energy
 *                    = 2 - Stokes parameter U devided by energy
 *                    = 3 - Stokes parameter V devided by energy
 *                    = 4 - degree of polarization
 *                    = 5 - polarization angle psi=0.5*atan(U/Q)
 *                    = 6 - the other pol. angle
 *                          beta=0.5*asin(V/sqrt(Q*Q+U*U+V*V))
 * par25 ... nthreads - how many threads should be used for computations
 *
 * If the height of the primary source above the black hole is negative
 * then this model assumes radial power-law dependence of local emissivity
 * (tables in the KBHlamp.fits file are not used) and it is assumed 
 * that the disc is illuminated from all directions isotropically
 * (we integrate the tables in reflspectra.fits file over azimuthal angles
 * and cosines of incident angles).
 *
 * NOTES:
 *  -> accuracy vs. speed trade off depends mainly on: nrad, nphi
 *
 *  -> the normalization is in the same physical units as the local flux
 *     i.e. we do not renormalize the spectra
 *
 *  -> for positive values of height this model includes a physical model of
 *     polarization based on Rayleigh scattering in single scattering
 *     approximation, for negative values of height the local emission is
 *     completely linearly polarized in the direction perpendicular to the disc
 *
 * -------------------------------------------
 * KBHlamp.fits
 *
 * This file contains pre-calculated values of the functions needed for the
 * lamp-post models. It is supposed that a primary source of emission is placed
 * on the axis at a height h from the centre above the Kerr black hole.
 * The matter in the disc rotates on stable circular (free) orbits above the
 * marginally stable orbit and it is freely falling below this orbit
 * where it has the same energy and angular momentum as the matter 
 * which is on the marginally stable orbit. It is assumed that the corona
 * between the source and the disc is optically thin, therefore ray-tracing
 * in the vacuum Kerr space-time could be used for computing the functions.
 *
 * There are five functions stored in the KBHlamp.fits file as columns. They are
 * parametrized by the value of the horizon of the black hole (FITS table
 * extensions), height of the primary source (rows) and either the
 * inclination angles of the observer (elements) or the radius 
 * in GM/c^2 at which a photon strikes the disc (elements). All these are
 * defined as vectors at the beginning of the FITS file. The tables are 
 * defined for
 * r_horizon = 1.00, 1.05, ..., 1.95, 2.00,
 * h-r_horizon = 0.1 - 100 (100 values with exponentially growing step),
 * inclination = 0.1, 1, 5, 10, 15, ..., 80, 85, 89 and
 * r-r_horizon = 0.01 - 1000 (100 values with exponentially growing step).
 *
 * The functions included are:
 * dWadWo - amplification of the primary source:
 *        = sin(theta_axis_local) / sin(theta_observer) *
 *          dtheta_axis_local/dtheta_observer
 * g-factor - the ratio of the energy of a photon hitting the disc to the energy
 *            of the same photon when emitted from a primary source,
 * cosine of the incident angle - an absolute value of the cosine of the local
 *                                incident angle between the incident light ray
 *                                and local disc normal
 * azimuthal incident angle - the angle (in radians) between the projection of
 *                            the three-momentum of the incident photon into the
 *                            disc (in the local rest frame co-moving with the
 *                            disc) and the radial tetrad vector.
 * dWadSd - part of the amplification of the incident flux on the disc from the
 *          primary source:
 *        = sin(theta_axis_local) * dtheta_axis_local / dtheta_fake, with the
 *          theta_fake defined as tan(theta_fake) = radius_incident / height
 *          one has to multiply by h/(r^2+h^2) to get the full amplification
 *
 * The definition of the file KBHlamp.fits:
 * 0. All of the extensions defined below are binary.
 * 1. The first extension contains a vector of the horizon values in GM/c^2
 *    (1.00 <= r_horizon <= 2.00).
 * 2. The second extension contains a vector of the values of heights h of
 *    a primary source in GM/c^2 (0.1 <= h-r_horizon <= 100).
 * 3. The third extension contains a vector of the values of the observer
 *    inclinations (0.1 <= inclination <= 89).
 * 4. The fourth extension contains a vector of the values of the incident
 *    radius (0.01 <= radius-r_horizon <= 1000).
 * 5. In the following extensions the functions are defined, each extension is
 *    for a particular value of r_horizon, each row is for a particular value of
 *    height and each element is either for a particular value of inclination
 *    (dWadWo) or incident radiation (the rest of the functions).
 * 6. Each of these extensions has five columns. In each column, a particular
 *    function is stored - dWadWo, g-factor, cosine of the local incident angle,
 *    azimuthal incident angle and dWadSd, respectively.
 * -------------------------------------------
 *
 * reflspectra.fits
 *
 * Dependence of the reflection coefficients on the angle of incidence, angle
 * of emission and emission azimuthal angle is stored in this FITS file.
 * The emission is induced by a power-law incident radiation. Values were
 * computed by the Monte Carlo simulations of Compton scattering,
 * for details see Matt, Perola & Piro (1991). The new tables were computed 
 * by NOAR code (see ???). The reflected radiation depends on the photon index
 * of the incident radiation.
 * There are several binary extensions in this fits file:
 * - the first extension contains energy values in keV where the coefficients
 *   are computed, currently the interval from 0.8 to 100 keV is covered
 * - the second extension contains the values of the relative azimuthal angle of
 *   the incident and emitted photon
 * - the third extension contains the absolute values of the cosine of the
 *   incident angles
 * - the fourth extension contains the values of the cosine of the emission
 *   angles
 * - the fifth extension contains the values of the photon indices of the
 *   incident powerlaw, currently there are tables only for Gamma  = 2.0
 * - in the following extensions the reflection coefficients are stored, each
 *   extension is for a particular value of PhoIndex, here values of the
 *   coefficients are stored as a 4D array for different incident angles,
 *   and for different angles of emission, azimuthal emission angle and energy
 *   in fortran 77:
 *     memr(f+ie-1+ne*(iazim-1)+ne*nazim*(icosi-1)+ne*nazim*ncosi*(icose-1))
 *   in c: f[icose][icosi][iazim][ie]
 *
 *==============================================================================
 *
 * NB: accuracy vs. speed trade off depends mainly on: nrad, nphi
 *
 *
 * 13.12.2007. changed to work with rene's tables of cold disc
 *             the result is not much smooth though
 * 27. 8.2009  changed to work with new tables with azimuthal dependence
 * 12.11.2009  changed to work with new lamp-post tables, we can fit for height
 *             as well now
 *
 *==============================================================================
 */

/*******************************************************************************
*******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "fitsio.h"

/*******************************************************************************
*******************************************************************************/
#ifdef OUTSIDE_XSPEC

#define IFL    1
#define NPARAM 25
#define NE     200
#define E_MIN  1.2
#define E_MAX  100.

int main() {

void KYNlpcr(const double *ear, int ne, const double *param, int ifl, 
            double *photar, double *photer, const char* init);

double ear[NE+1], photar[NE], photer[NE], param[NPARAM];
char   initstr[0] = "";
int    ie;

param[ 0] = 1.;         // a/M
param[ 1] = 80.;        // thetaO
param[ 2] = 1.;         // rin
param[ 3] = 1.;         // ms
param[ 4] = 400.;       // rout
param[ 5] = 0.;         // phi0
param[ 6] = 360.;       // dphi
param[ 7] = 0.05;         // M/M8
param[ 8] = 2.5;         // height
param[ 9] = 2.;         // PhoIndex
param[10] = 0.001;      // L/Ledd
param[11] = 1.;         // Np:Nr
param[12] = 3.;        // line
param[13] = 300.;       // Ec
param[14] = -6.;        // alpha
param[15] = 0.;         // beta
param[16] = 0.;         // rcloud
param[17] = 0.;         // zshift
param[18] = 80.;        // ntable
param[19] = 300.;       // nrad
param[20] = 1.;         // division
param[21] = 360.;       // nphi
param[22] = 1.;         // smooth
param[23] = 0.;         // Stokes
param[24] = 4.;         // nthreads

for(ie = 0; ie <= NE; ie++) {
//  ear[ie] = E_MIN + ie * (E_MAX-E_MIN) / NE;
  ear[ie] = E_MIN * pow(E_MAX / E_MIN, ((double) ie) / NE);
}

KYNlpcr(ear, NE, param, IFL, photar, photer, initstr);
return(0);
}

#endif
/*******************************************************************************
*******************************************************************************/

#define LAMP "KBHlamp_q.fits"
#define REFSPECTRA "reflspectra.fits"
#define NLINES 6
#define PI 3.14159265358979
#define PI2 6.2831853071795865
#define MPC_2 1.05e-49
#define ERG 6.241509e8
#define E0 0.1
// Ledd is in erg (not W) and multiplied by 10^8 due to (M / (10^8*Msun)) scale
#define LEDD 1.26e46
#define HUBBLE 70.
#define CLIGHT 299792.458

/* Let's declare variables that are common for the main and emissivity 
   subroutines */
static float  *cosi, *cose, *azimuth, *radius;
static double *energy, *flux_c, *flux_l;
static double *gfac, *cosin, *phiph, *transf_d;
static double h, h_rh, am2, r_plus, gam;
static int    polar, line;
static long   ncosi, ncose, nazim, nrad;

extern char*  FGMODF(void);
extern char*  FGMSTR(char* dname);
extern int    xs_write(char* wrtstr, int idest);
extern double incgamma(double a, double x);
extern void   cutoffpl(double *ear, const int ne, double *param, double *photar);

void KYNlpcr(const double *ear, int ne, const double *param, int ifl, 
             double *photar, double *photer, const char* init) {

extern int ide(const double *ear, const int ne, const int nt, double *far, 
               double *qar, double *uar, double *var, 
               const double *ide_param, void (*emissivity)(), 
               const int ne_loc);

void emis_KYNlpcr(double** ear_loc, const int ne_loc, const int nt, 
                  double *far_loc, double *qar_loc, double *uar_loc, 
                  double *var_loc, const double r, const double phi, 
                  const double cosmu, double phiphoton, 
                  const double alpha_o, const double beta_o, 
                  const double delay, const double g);

/* Let's declare static variables whose values should be remembered for next
   run in XSPEC */
static char   kydir[255]="";
static char   pname[128]="KYDIR", pkyLxLamp[128] = "KYLxLamp";
static char   pkyRefl[128] = "KYRefl";
static long   nrh, ngamma, nh, nincl, nener;
static float  *r_horizon, *gamma, *height, *incl, *dWadWo, *dWadSd, *q, *pr, 
              *flux0;
static double transf_o;
static double h_rh_old = -1., gam_old = -1., am_old = -1., thetaO_old = -1.;
static int    first = 1, first_h = 1, line_old = -99;
static const int npoints[6] = {3, 3, 3, 3, 10, 2};
static const int fpoint[6] = {13, 69, 116, 157, 302, 318};

FILE   *fw;
char   errortxt[80];
char   kyLxLamp[32], kyRefl[32];
double ide_param[25];
double far[ne], qar[ne], uar[ne], var[ne], pd[ne], pa[ne], pa2[ne];
float  *energy0;
double ttmp, ttmp1, utmp, utmp1, vtmp, vtmp1, y1, y2, y3, y4, y5, y6, y7, y8;
double pr_final, q_final, pom, pom1, pom2, pom3;
double r, r2, delta, ULt, rms, tmp1, Ut, U_phi, U_r, Ur, UrU_r, Lms, Ems;
double am, thetaO, Ec, f0, f1, ratio;
double mass, Np, Anorm, Dnorm, g_L, Lx, flux_prim, flux_refl, refl_ratio;
double zzshift;
double pamin, pamax, pa2min, pa2max, NpNr;
int    imin, imax, irh0, ih0, ith0, igamma, igam0, icosi, icose, iazim, iline;
int    i, ie, je, stokes;
int    ipoints, iener;
// the following are needed for cut-off power-law taken from XSPEC
double ear1[ne + 1], param1[2];
double photar1[ne];
//int    irh, ih;
//char   gg[3], tables1[34]="cold-disk-azimuth-reflection-gamma", 
//       tables2[5]="-Inc.", tables3[4]="-Azi", cazim[3], ccosi[3],
//       file_dat[52];

// these are needed to work with a fits file...
fitsfile *fptr;
char     tables_file[255];
int      hdutype = 2;
int      colnum = 1;
long     frow = 1, felem = 1, nelems, nrow;
float    float_nulval = 0.;
int      nelements1, nelements2;
int      ihorizon, irow, anynul, status = 0;//, maxdim=1000, naxis;

// Let's initialize parameters for subroutine ide()
for (ie = 0; ie < ne; ie ++) far[ie] = 0.;
// am - black hole angular momentum
ide_param[0] = param[0];
am = ide_param[0];
am2 = am * am;
pom1 = pow(1. + am, 1./3.);
pom2 = pow(1. - am, 1./3.);
pom3 = pow(1. - am2, 1./3.);
pom = 1. + pom3 * (pom1 + pom2);
pom1 = sqrt(3. * am2 + pom * pom);
if (am >= 0.) rms = 3. + pom1 -sqrt((3. - pom) * (3. + pom + 2. *pom1));
else rms = 3. + pom1 + sqrt((3. - pom) * (3. + pom + 2. *pom1));
r_plus = 1. + sqrt(1. - am2);
// thetaO - observer inclination
ide_param[1] = param[1];
thetaO = ide_param[1];
// rin - inner edge of non-zero disc emissivity
ide_param[2] = param[2];
// ms - whether to integrate from rin or rms
ide_param[3] = param[3];
// rout - outer edge of non-zero disc emissivity
ide_param[4] = param[4];
// phi  - lower azimuth of non-zero disc emissivity (deg)
ide_param[5] = param[5];
// dphi - (phi+dphi) is upper azimuth of non-zero disc emissivity (deg)
ide_param[6] = param[6];
// nrad - number of grid points in radius
ide_param[7] = param[19];
// division - type of division in r integration (0-equidistant, 1-exponential)
ide_param[8] = param[20];
// nphi - number of grid points in azimuth
ide_param[9] = param[21];
// smooth - whether to smooth the resulting spectrum (0-no, 1-yes)
ide_param[10] = param[22];
// normal - how to normalize the final spectrum
ide_param[11] = -1.;
// zshift - overall Doppler shift
ide_param[12] = param[17];
if (param[17] > 0.) {
  ide_param[12] = param[17];
  Dnorm = pow(HUBBLE / CLIGHT / param[17], 2.);
}else if(param[17] < 0.) {
  ide_param[12] = 0.;
  Dnorm = pow(-HUBBLE / CLIGHT /param[17], 2.);
}else{
  ide_param[12] = 0;
  Dnorm = 1.;
}
// zzshift - multiplication factor for gfac from zshift needed for primary
zzshift=1.0/(1.0+ide_param[12]);
// ntable - table model (defines fits file with tables)
ide_param[13] = param[18];
// M/M8 - black hole mass in units of 10^8 solar masses
mass = param[7];
// height - height of the lamp above the black hole
h = param[8];
if (h >= 0.) if (h < r_plus) h_rh = 0.;
             else h_rh = h - r_plus;
else h_rh = h;
// PhoIndex - power-law energy index of the lamp emission
gam = param[9];
// L/Ledd - dE/dt primary isotropic flux in Ledd
Np = param[10];
// NpNr - ratio of the primary normalization to the reflected normalization
NpNr = param[11];
// line - whether to include line in the spectra
line = (int) param[12];
// Ec - energy cut-off
Ec = param[13];
// edivision - type of division in local energies (0-equidistant, 1-exponential)
ide_param[14] = 1.;
// periodic and dt need not to be set for nt=1
// (ide_param[15], ide_param[16])
// polar - whether we need value of change in polarization angle (0-no,1-yes)
stokes = (int) param[23];
if ((stokes < 0.) || (stokes > 6)) {
  xs_write("kynlpcr: Stokes has to be 0-6", 5);
  for(ie = 0; ie < ne; ie++) photar[ie] = 0.;
  return;  
}
polar = 0;
if (stokes > 0) polar = 1;
ide_param[17] = polar;
// delay_r and delay_phi are not used
// (ide_param[18], ide_param[19])
// number of threads for multithread computations
ide_param[20] = param[24];
// alpha - position of the cloud in alpha impact parameter (in GM/c^2)
ide_param[21] = param[14];
// beta - position of the cloud in beta impact parameter (in GM/c^2)
ide_param[22] = param[15];
// rcloud - radius of the cloud (in GM/c^2)
ide_param[23] = param[16];
//whether the flux defined in emissivity subroutine is local one (0) or the 
//observed one (1)
ide_param[24] = 0.;

/******************************************************************************/
// Let's read the lamp post tables
if (first_h && (h_rh >= 0.)) {
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
    xs_write("\nkynlpcr: set the KYDIR to the directory with the KY tables",5);
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
    xs_write("kynlpcr: Failed to allocate memory for tmp arrays.", 5);
    for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
    return;
  }
// Read the data in the 'r_horizon' table
  nelems = nrh;
// FTGCV reads the VALUES from the first column.
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
    xs_write("kynlpcr: Failed to allocate memory for tmp arrays.", 5);
    for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
    return;
  }
// Read the data in the 'height' table
  nelems = nh;
// FTGCV reads the VALUES from the first column.
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
    xs_write("kynlpcr: Failed to allocate memory for tmp arrays.", 5);
    for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
    return;
  }
// Read the data in the 'inclination' table
  nelems = nincl;
// FTGCV reads the VALUES from the first column.
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
// Allocate memory for radius...
  if ((radius = (float *) malloc(nrad * sizeof(float))) == NULL) {
    xs_write("kynlpcr: Failed to allocate memory for tmp arrays.", 5);
    for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
    return;
  }
// Read the data in the 'r_rh' table
  nelems = nrad;
// FTGCV reads the VALUES from the first column.
  ffgcv(fptr, TFLOAT, colnum, frow, felem, nelems, &float_nulval, radius,
        &anynul, &status);
//******************************************************************************
//  for ( i=0; i<nrad; i++)fprintf(stdout,"%f\n",radius[i]);
//******************************************************************************   
// Let's read the tables for dWadWo, q, p^r and dWadSd
// allocate memory for the arrays
  if ((dWadWo = (float *) malloc(nincl * nh * nrh * sizeof(float))) == NULL) {
    xs_write("kynlpcr: Failed to allocate memory for tmp arrays.", 5);
    for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
    return;
  }
  if ((q = (float *) malloc(nrad * nh * nrh * sizeof(float))) == NULL) {
    xs_write("kynlpcr: Failed to allocate memory for tmp arrays.", 5);
    for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
    return;
  }
  if ((pr = (float *) malloc(nrad * nh * nrh * sizeof(float))) == NULL) {
    xs_write("kynlpcr: Failed to allocate memory for tmp arrays.", 5);
    for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
    return;
  }
  if ((dWadSd = (float *) malloc(nrad * nh * nrh * sizeof(float))) == NULL) {
    xs_write("kynlpcr: Failed to allocate memory for tmp arrays.", 5);
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
  irh=0;
  ih=30;
  for ( i=0; i<nincl; i++ ) fprintf(stdout,"%d\t%f\t%f\n",i,incl[i],
    dWadWo[i+nincl*ih+nincl*nh*irh]);
  for ( i=0; i<nrad; i++) fprintf(stdout,"%d\t%f\t%f\t%f\t%f\n",i,radius[i],
    q[i+nrad*ih+nrad*nh*irh],pr[i+nrad*ih+nrad*nh*irh],
    dWadSd[i+nrad*ih+nrad*nh*irh]);
*******************************************************************************/
// Firstly we have to free allocated memory for the arrays gfac,
// cosin, phiph, transf_d
  if ((gfac = (double *) malloc(nrad * sizeof(double))) == NULL) {
    xs_write("kynlpcr: Failed to allocate memory for tmp arrays.", 5);
    for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
    return;
  }
  if ((cosin = (double *) malloc(nrad * sizeof(double))) == NULL) {
    xs_write("kynlpcr: Failed to allocate memory for tmp arrays.", 5);
    for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
    return;    
  }
  if ((phiph = (double *) malloc(nrad * sizeof(double))) == NULL) {
    xs_write("kynlpcr: Failed to allocate memory for tmp arrays.", 5);
    for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
    return;
  }
  if ((transf_d = (double *) malloc(nrad * sizeof(double))) == NULL) {
    xs_write("kynlpcr: Failed to allocate memory for tmp arrays.", 5);
    for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
    return;    
  }
  first_h = 0;
}
/******************************************************************************/
if (h >= 0.) {
  if (h_rh > height[nh - 1]) {
    sprintf(errortxt, "kynlpcr: the height must be lower than or equal to %f.",
            height[nh - 1] + r_plus);
    xs_write(errortxt, 5);
    for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
    return;
  }
  if (h_rh < height[0]) {
    sprintf(errortxt, "kynlpcr: the height is too low, we set it to %f.",
	    height[0] + r_plus);
    xs_write(errortxt, 5);
    h_rh = height[0];
    h = h_rh + r_plus;
  }
}
// Let's interpolate the tables to desired spin and height
if (((am != am_old) || (h_rh != h_rh_old) || (thetaO != thetaO_old))
   && (h_rh >= 0.)) {
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
  ttmp = (r_plus - r_horizon[irh0 - 1]) / (r_horizon[irh0] - r_horizon[irh0 - 1]);
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
  transf_o = (vtmp1 * (utmp1 * (ttmp1 * y1 + ttmp * y2) + utmp *
             (ttmp * y3 + ttmp1 * y4)) + vtmp * (utmp1 *
             (ttmp1 * y5 + ttmp * y6) + utmp * (ttmp * y7 + ttmp1 * y8)));
  if ((am != am_old) || (h_rh != h_rh_old)) {
    for (i = 0; i < nrad; i++) {
// q from the axis to the disc
      y1 = q[i + nrad * (ih0 - 1) + nrad * nh * (irh0 - 1)];
      y2 = q[i + nrad * (ih0 - 1) + nrad * nh * irh0];
      y3 = q[i + nrad * ih0 + nrad * nh * irh0];
      y4 = q[i + nrad * ih0 + nrad * nh * (irh0 - 1)];
      q_final = utmp1 * (ttmp1 * y1 + ttmp * y2) + 
                utmp * (ttmp * y3 + ttmp1 * y4);
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
        U_phi = (r2 + am2 - 2. * am * sqrt(r)) / sqrt(r) / tmp1;
        U_r = 0.;
        Ur = 0.;
        UrU_r = 0.;
      }
      else {
        tmp1 = sqrt(rms * (rms - 3.) + 2. * am * sqrt(rms));
        Lms = (rms * rms + am2 - 2. * am * sqrt(rms)) / sqrt(rms) / tmp1;
        Ems = (rms * (rms - 2.) + am * sqrt(rms)) / rms / tmp1;
        Ut = (Ems * (r * (r2 + am2) + 2. * am2) - 2. * am * Lms) / r / delta;
        U_phi = Lms;
        UrU_r = -1. + ((r2 + am2 + 2. * am2 / r) * Ems * Ems - 4. * am /
                r * Ems * Lms - (1. - 2. / r) * Lms * Lms) / delta;
        if (UrU_r < 0.) UrU_r = 0.;
        U_r = -sqrt(UrU_r / delta) * r;
        Ur = -sqrt(delta * UrU_r) / r;
      }
      tmp1 = Ut - pr_final * U_r;
// gfactor  from the axis to the disc
      gfac[i] = tmp1 / ULt;
// cosin at the disc
      cosin[i] = q_final / r / tmp1;
// phip_i at the disc
      phiph[i] = atan2(-U_phi, r * (pr_final - Ur * tmp1));
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
if (first) {
// Let's read the 'refspectra' fits file
// The status parameter must always be initialized.
  status = 0;
// Set up the directory and the fits file name of the tables
// and open the fits file with the tables
  sprintf(kydir, "%s", FGMSTR(pname));
  if (strlen(kydir) == 0) sprintf(tables_file, "./%s", REFSPECTRA);
  else if (kydir[strlen(kydir) - 1] == '/') sprintf(tables_file, "%s%s",
                                                    kydir, REFSPECTRA);
  else sprintf(tables_file, "%s/%s", kydir, REFSPECTRA);
// Let's read the 'fluorescent_line' fits file
// The status parameter must always be initialized.
  status = 0;
  ffopen(&fptr, tables_file, READONLY, &status);
  if(status) {
    sprintf(tables_file, "%s%s", FGMODF(), REFSPECTRA);
    status = 0;
    ffopen(&fptr, tables_file, READONLY, &status);
  }
  if (status){
    if (status) ffrprt(stderr, status);
    ffclos(fptr, &status);
    xs_write("\nkynlpcr: set the KYDIR to the directory with the KY tables.",5);
    for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
    return;
  }
// Let's read tables (binary tables => hdutype=2)
// Move to the extension 'energy' and read its values
  ffmrhd(fptr, 1, &hdutype, &status);
  ffgnrw(fptr, &nener, &status);
//******************************************************************************
//  fprintf(stdout,"nener = %ld\n",nener);
//******************************************************************************   
// Allocate memory for energy...
  if ((energy0 = (float *) malloc((nener + 1) * sizeof(float))) == NULL) {
    xs_write("kynlpcr: Failed to allocate memory for tmp arrays.", 5);
    for(ie = 0; ie < ne; ie++) photar[ie] = 0.;
    return;
  }
  if ((energy = (double *) malloc((nener + 1) * sizeof(double))) == NULL) {
    xs_write("kynlpcr: Failed to allocate memory for tmp arrays.", 5);
    for(ie = 0; ie < ne; ie++) photar[ie] = 0.;
    return;
  }
// Read the data in the 'energy' table
  nelems = nener;
// FTGCVE reads the VALUES from the first column.    
  ffgcv(fptr, TFLOAT, colnum, frow, felem, nelems, &float_nulval, energy0,
        &anynul, &status);
  for(ie = 0; ie < nener; ie++) energy[ie] = (double) energy0[ie];
//******************************************************************************
//  for ( ie=0; ie<nener; ie++)fprintf(stdout,"%f\n",energy[ie]);
//******************************************************************************   
// Move to the extension 'azimuthal_angle' and read its values
  ffmrhd(fptr, 1, &hdutype, &status);
  ffgnrw(fptr, &nazim, &status);
//******************************************************************************
//  fprintf(stdout,"nazim = %d\n",nazim);
//******************************************************************************   
// Allocate memory for azimuth...
  if ((azimuth = (float *) malloc(nazim * sizeof(float))) == NULL) {
    xs_write("kynlpcr: Failed to allocate memory for tmp arrays.", 5);
    for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
    return;  
  }
// Read the data in the 'azimuthal_angle' table
  nelems = nazim;
// FTGCVE reads the VALUES from the first column.    
  ffgcv(fptr, TFLOAT, colnum, frow, felem, nelems, &float_nulval, azimuth,
        &anynul, &status);
//******************************************************************************
//  for ( i=0; i<nazim; i++)fprintf(stdout,"%f\n",azimuth[i]);
//******************************************************************************   
// Move to the extension 'cos_theta_i' and read its values
  ffmrhd(fptr, 1, &hdutype, &status);
  ffgnrw(fptr, &ncosi, &status);
//******************************************************************************
//  fprintf(stdout,"ncosi = %ld\n",ncosi);
//******************************************************************************   
// Allocate memory for cosi...
  if ((cosi = (float *) malloc(ncosi * sizeof(float))) == NULL) {
    xs_write("kynlpcr: Failed to allocate memory for tmp arrays.", 5);
    for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
    return;
  }
// Read the data in the 'cos_theta_i' table
  nelems = ncosi;
// FTGCVE reads the VALUES from the first column.    
  ffgcv(fptr, TFLOAT, colnum, frow, felem, nelems, &float_nulval, cosi,
        &anynul, &status);
//******************************************************************************
//  for ( i=0; i<ncosi; i++)fprintf(stdout,"%f\n",cosi[i]);
//******************************************************************************   
// Move to the extension 'cos_theta_e' and read its values
  ffmrhd(fptr, 1, &hdutype, &status);
  ffgnrw(fptr, &ncose, &status);
//******************************************************************************
//  fprintf(stdout,"ncose = %d\n",ncose);
//******************************************************************************   
// Allocate memory for cose...
  if ((cose = (float *) malloc(ncose * sizeof(float))) == NULL) {
    xs_write("kynlpcr: Failed to allocate memory for tmp arrays.", 5);
    for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
    return;      
  }
// Read the data in the 'cos_theta_e' table
  nelems = ncose;
// FTGCVE reads the VALUES from the first column.    
  ffgcv(fptr, TFLOAT, colnum, frow, felem, nelems, &float_nulval, cose,
        &anynul, &status);
//******************************************************************************
//  for ( i=0; i<ncose; i++)fprintf(stdout,"%f\n",cose[i]);
//******************************************************************************   
// Move to the extension 'Gamma' and read its values
  ffmrhd(fptr, 1, &hdutype, &status);
  ffgnrw(fptr, &ngamma, &status);
//******************************************************************************
//  fprintf(stdout,"ngamma = %d\n",ngamma);
//******************************************************************************   
// Allocate memory for gamma...
  if ((gamma = (float *) malloc(ngamma * sizeof(float))) == NULL) {
    xs_write("kynlpcr: Failed to allocate memory for tmp arrays.", 5);
    for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
    return;    
  }
// Read the data in the 'Gamma' table
  nelems = ngamma;
// FTGCVE reads the VALUES from the first column.    
  ffgcv(fptr, TFLOAT, colnum, frow, felem, nelems, &float_nulval, gamma,
        &anynul, &status);
//******************************************************************************
//  for ( i=0; i<ngamma; i++)fprintf(stdout,"%f\n",gamma[i]);
//******************************************************************************   
// Allocate memory for flux0...
  if ((flux0 = (float *) malloc(ncosi * ncose * nazim * nener * ngamma * 
                                sizeof(float))) == NULL) {
    xs_write("kynlpcr: Failed to allocate memory for tmp arrays.", 5);
    for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
    return;
  }
  for (igamma = 0; igamma < ngamma; igamma++) {
// Move to the extension for the particular table
    ffmrhd(fptr, 1, &hdutype, &status);
// Read the data in the table
    nelems = nener * nazim * ncosi * ncose;
// FTGCVE reads the VALUES from the first column.
    ffgcv(fptr, TFLOAT, colnum, frow, felem, nelems, &float_nulval, 
          &flux0[ncosi * ncose * nazim * nener * igamma], &anynul, &status);
  }
// The FITS file must always be closed before exiting the program.
  ffclos(fptr, &status);
/*******************************************************************************
//  for (igamma = 0; igamma < ngamma; igamma++) {
  for (igamma = 5; igamma <= 5; igamma++) {
    if (gamma[igamma] != 1.)
      Anorm = (pow(200., 1.-gamma[igamma])-pow(0.8, 1.-gamma[igamma])) /
            (1.-gamma[igamma]);
    else Anorm = log(250.);
    Anorm = 0.05 * 2. * 5. * 0.003 * log(10.) / 180. * PI / Anorm;
//    for (iazim = 0; iazim < nazim; iazim++) {
    for (iazim = 0; iazim <= 0; iazim++) {
//      for (icosi = 0; icosi < ncosi; icosi++) {
      for (icosi = 0; icosi <= 0; icosi++) {
//      sprintf(gg, "%3.1f", gamma[igamma]);
//      if (azimuth[iazim] < 180.)
//        sprintf(cazim,"%03d", (int) (azimuth[iazim]+2.5));
//      else
//        sprintf(cazim,"%03d",(int) (azimuth[iazim]-2.5));
//      sprintf(ccosi, "%03d", (int) (1000.*cosi[icosi]));
//      sprintf(file_dat, "%s%s%s%s%s%s", tables1, gg, tables2, ccosi, tables3, cazim);
//      fprintf(stdout,"|%s|\n",file_dat);
//      fw = fopen(file_dat, "w");
        fprintf(stdout, "%3.1f %03d %03d\n", gamma[igamma], 
          (int) (1000.*cosi[icosi]), (int) (azimuth[iazim]+2.5));
        fw = fopen("test.dat", "w");
        for (ie = 0; ie < nener; ie++) {
          fprintf(fw, "%14.6f", energy[ie]);
          for (icose = 0; icose < ncose; icose++) {
            fprintf(fw, "\t%E", flux0[ie+nener*iazim+nener*nazim*icosi+
             nener*nazim*ncosi*icose+nener*nazim*ncosi*ncose*igamma]*Anorm*
             energy[ie]);
          }
          fprintf(fw,"\n");
        }
        fclose(fw);
      }
    }
  }
*******************************************************************************/
}

/******************************************************************************/
// Let us check if our model is valid for the value of cut-off energy Ecut
if (Ec < energy[nener-1]) {
  sprintf(errortxt, "kynlpcr: WARNING - Ecut = %f\n", Ec);
  xs_write(errortxt, 5);
  sprintf(errortxt, "The reflection model kynlpcr cannot account for Ecut < %f keV\n", 
          energy[nener-1]);
  xs_write(errortxt, 5);
}
if (first || (h_rh * h_rh_old <= 0.) || (gam != gam_old) || (line != line_old)){
// given gam, find corresponding index in gamma():
  imin = 0;
  imax = ngamma;
  igam0 = ngamma / 2;
  while ((imax - imin) > 1) {
    if (gam >= gamma[igam0-1]) imin = igam0;
    else imax = igam0;
    igam0 = (imin + imax) / 2;
  }
  if (igam0 == 0) igam0 = 1;
  if (igam0 == ngamma) igam0 = ngamma - 1;
  ttmp = (gam - gamma[igam0 - 1]) / (gamma[igam0] - gamma[igam0 - 1]);
  ttmp1 = 1. - ttmp;
// Firstly we have to free allocated memory for these arrays if we have
// already allocated it and the dimensions of the arrays have changed
  if ((!first) && (h_rh * h_rh_old <= 0.)) {
// Free memory from tmp arrays...
    free((void *) flux_c);
    free((void *) flux_l);
  }
  if (h_rh >= 0.) {
// Allocate memory for flux...
    if (first || (h_rh * h_rh_old <= 0.)) {
      if ((flux_c = (double *) malloc(ncosi * ncose * nazim * nener * sizeof(double))) == NULL) {
        xs_write("kynlpcr: Failed to allocate memory for tmp arrays.", 5);
        for(ie = 0; ie < ne; ie++) photar[ie] = 0.;
        return;
      }
      if ((flux_l = (double *) malloc(ncosi * ncose * nazim * nener * sizeof(double))) == NULL) {
        xs_write("kynlpcr: Failed to allocate memory for tmp arrays.", 5);
        for(ie = 0; ie < ne; ie++) photar[ie] = 0.;
        return;
      }
    }
    for(icosi = 0; icosi < ncosi; icosi++) {
      for(icose = 0; icose < ncose; icose++) {
        for(iazim = 0; iazim < nazim; iazim++) {
          for(je = 0; je < nener; je++) {
            ie = je + nener * iazim + nener * nazim * icosi + nener * nazim * ncosi * icose;
            y1 = flux0[ie + nener * nazim * ncosi * ncose * (igam0-1)] * 
                 pow(energy[je], gamma[igam0 - 1]);
            y2 = flux0[ie + nener * nazim * ncosi * ncose * igam0] * 
                 pow(energy[je], gamma[igam0]);
            flux_c[ie] = (ttmp1 * y1 + ttmp * y2) / pow(energy[je], gam);
            flux_l[ie] = 0.;
          }
// Let us divide the spectra to continuum and line(s)
          ie = nener * iazim + nener * nazim * icosi + nener * nazim * ncosi * icose;
// NOTE: the lines at 0.87, 1.28, 1.77, 2.35, 6.4 and 7.132 keV
//       are substracted for flux_c and we separate the desired lines in flux_l
          for(iline = 0; iline < NLINES; iline++) {
            je = fpoint[iline];
            ipoints = npoints[iline];
            f0 = flux_c[je - ipoints + ie];
            ratio = (flux_c[je + 1 + ie] - f0) /
                    (energy[je + 1] - energy[je - ipoints]);
            for (iener = 1; iener <= ipoints; iener++) {
              f1 = f0 + ratio * (energy[je - ipoints + iener] - 
                                 energy[je - ipoints]);
              if (f1 < flux_c[je - ipoints + iener + ie]) {
// for lines at 0.87, 1.28, 1.77, 2.35 keV
                if ((abs(line) == 3) && (iline <= 3)) 
                  flux_l[je - ipoints + iener + ie] = 
                    flux_c[je - ipoints + iener + ie] - f1;
// for the line at 6.4 keV
                if ((abs(line) >= 1) && (iline == 4)) 
                  flux_l[je - ipoints + iener + ie] = 
                    flux_c[je - ipoints + iener + ie] - f1;
// for the line at 7.132 keV
                if ((abs(line) >= 2) && (iline == 5)) 
                  flux_l[je - ipoints + iener + ie] = 
                    flux_c[je - ipoints + iener + ie] - f1;
                flux_c[je - ipoints + iener + ie] = f1;
              }
            }
          }
        }
      }
    }
  } else {
// Allocate memory for flux...
    if (first || (h_rh * h_rh_old <= 0.)) {
      if ((flux_c = (double *) malloc(ncose * nener * sizeof(double))) == NULL) {
        xs_write("kynlpcr: Failed to allocate memory for tmp arrays.", 5);
        for(ie = 0; ie < ne; ie++) photar[ie] = 0.;
        return;
      }
      if ((flux_l = (double *) malloc(ncose * nener * sizeof(double))) == NULL) {
        xs_write("kynlpcr: Failed to allocate memory for tmp arrays.", 5);
        for(ie = 0; ie < ne; ie++) photar[ie] = 0.;
        return;
      }
    }
    for(icose = 0; icose < ncose; icose++) {
      for(je = 0; je < nener; je++) {
// let's interpolate the values in flux0() table to desired gamma and integrate
// in incident angles and azimuthal angles (isotropic illumination):
        flux_c[je + nener * icose] = 0.;
        flux_l[je + nener * icose] = 0.;
        utmp = pow(energy[je], gamma[igam0]);
        utmp1 = pow(energy[je], gamma[igam0-1]);
        for(iazim = 0; iazim < nazim; iazim++) {
          for(icosi = 0; icosi < ncosi; icosi++) {
            ie = je + nener * iazim + nener * nazim * icosi + nener * nazim * ncosi * icose;
            y1 = flux0[ie + nener * nazim * ncosi * ncose * (igam0 - 1)] * utmp1;
            y2 = flux0[ie + nener * nazim * ncosi * ncose * igam0] * utmp;
            flux_c[je + nener * icose] += (ttmp1 * y1 + ttmp * y2);
          }
        }
        flux_c[je + nener * icose] = flux_c[je + nener * icose]/ncosi/nazim/ 
                                     pow(energy[je], gam);
      }
// Let us divide the spectra to continuum and line(s)
      ie = nener * icose;
// NOTE: the lines at 0.87, 1.28, 1.77, 2.35, 6.4 and 7.132 keV
//       are substracted for flux_c and we separate the desired lines in flux_l
      for(iline = 0; iline < NLINES; iline++) {
        je = fpoint[iline];
        ipoints = npoints[iline];
        f0 = flux_c[je - ipoints + ie];
        ratio = (flux_c[je + 1 + ie] - f0) / 
                (energy[je + 1] - energy[je - ipoints]);
        for(iener = 1; iener <= ipoints; iener++) {
          f1 = f0 + ratio * (energy[je - ipoints + iener] - 
                             energy[je - ipoints]);
          if (f1 < flux_c[je - ipoints + iener + ie]) {
// for lines at 0.87, 1.28, 1.77, 2.35 keV
            if ((abs(line) == 3) && (iline <= 3)) 
              flux_l[je - ipoints + iener + ie] = 
                flux_c[je - ipoints + iener + ie] - f1;
// for the line at 6.4 keV
            if ((abs(line) >= 1) && (iline == 4))
              flux_l[je - ipoints + iener + ie] = 
                flux_c[je - ipoints + iener + ie] - f1;
// for the line at 7.132 keV
            if ((abs(line) >= 2) && (iline == 5)) 
              flux_l[je - ipoints + iener + ie] = 
                flux_c[je - ipoints + iener + ie] - f1;
            flux_c[je - ipoints + iener + ie] = f1;
          }
        }
      }
    }

/*******************************************************************************
// local spectrum output -- write energy0() and flux() into file:
    fw = fopen("kynlpcr_photar_loc.dat", "w");
    icose=0;
    for (ie = 0; ie < nener; ie++)
      fprintf(fw, "%14.6f\t%E\t%E\n", energy[ie],flux_c[ie+nener*icose],
                                      flux_l[ie+nener*icose]);
    fclose(fw);
*******************************************************************************/
  }
}
first = 0;
am_old = am;
thetaO_old = thetaO;
h_rh_old = h_rh;
gam_old = gam;
line_old = line;

/******************************************************************************/
#ifdef OUTSIDE_XSPEC
// let's write the input parameters to a file
fw = fopen("kynlpcr.txt", "w");
fprintf(fw, "a/M         %12.6f\n", param[0]);
fprintf(fw, "theta_o     %12.6f\n", param[1]);
fprintf(fw, "rin         %12.6f\n", param[2]);
fprintf(fw, "ms          %12d\n", (int) param[3]);
fprintf(fw, "rout        %12.6f\n", param[4]);
fprintf(fw, "phi         %12.6f\n", param[5]);
fprintf(fw, "dphi        %12.6f\n", param[6]);
fprintf(fw, "mass        %12.6f\n", param[7]);
fprintf(fw, "height      %12.6f\n", param[8]);
fprintf(fw, "PhoIndex    %12.6f\n", param[9]);
fprintf(fw, "L/Ledd      %12.6f\n", param[10]);
fprintf(fw, "Np:Nr       %12.6f\n", param[11]);
fprintf(fw, "line        %12d\n", (int) param[12]);
fprintf(fw, "Ecut        %12.6f\n", param[13]);
fprintf(fw, "alpha       %12.6f\n", ide_param[21]);
fprintf(fw, "beta        %12.6f\n", ide_param[22]);
fprintf(fw, "rcloud      %12.6f\n", ide_param[23]);
fprintf(fw, "zshift      %12.6f\n", param[17]);
fprintf(fw, "ntable      %12d\n", (int) param[18]);
fprintf(fw, "nrad        %12d\n", (int) param[19]);
fprintf(fw, "division    %12d\n", (int) param[20]);
fprintf(fw, "nphi        %12d\n", (int) param[21]);
fprintf(fw, "smooth      %12d\n", (int) param[22]);
fprintf(fw, "Stokes      %12d\n", (int) param[23]);
fprintf(fw, "polar       %12d\n", polar);
fprintf(fw, "r_horizon   %12.6f\n", r_plus);
fprintf(fw, "r_ms        %12.6f\n", rms);
fprintf(fw, "edivision   %12d\n", (int) ide_param[14]);
fprintf(fw, "nener       %12ld\n", nener);
fprintf(fw, "normal      %12.6f\n", ide_param[11]);
fprintf(fw, "nthreads    %12d\n", (int) ide_param[20]);
fclose(fw);
#endif
/******************************************************************************/
// Let's integrate local emission over the accretion disk
if (ide(ear, ne, 1, far, qar, uar, var, ide_param, emis_KYNlpcr, nener)) {
  for(ie = 0; ie < ne; ie++) photar[ie] = 0.;
  return;
}

//Let's compute the Np according to its definition, i.e. transform from 
//intrinsic or observed in 2-10keV to total intrinsic photon flux
//Lx is intrinsic photon flux in 2-10keV
g_L = sqrt(1. - 2. * h / (am2 + h * h));
if( Np < 0. ){
  Lx = -Np;
  Np *= - incgamma(2. - gam, E0 / Ec) / 
        ( incgamma(2. - gam, 2. / Ec) - incgamma(2. - gam, 10. / Ec));
}else{  
  Lx = Np / g_L / g_L / transf_o /
       ( incgamma(2. - gam, 2. / g_L / Ec) - incgamma(2. - gam, 10. / g_L / Ec));
  Np = Lx * incgamma(2. - gam, E0 / Ec);        
  Lx *= ( incgamma(2. - gam, 2. / Ec) - incgamma(2. - gam, 10. / Ec));
}
//Let's write the intrinsic photon flux in 2-10keV into the xset XSPEC variable
//KYLxLamp
sprintf(kyLxLamp, "%e", Lx);
FPMSTR(pkyLxLamp, kyLxLamp);

// Let's compute the incomplete gamma function with the XSPEC incGamma 
// function
Anorm = LEDD * mass * MPC_2 * ERG * Np / pow(Ec, 2. - gam) / PI2 / 2. / 
        incgamma(2. - gam, E0 / Ec) * Dnorm;
//if the primary is power-law without the cut-offs
/*if(gam!=2.) Anorm = LEDD * mass * MPC_2 * ERG * Np / PI2 / 2. /
                      ((pow(Ec,2.-gam) - pow(E0,2.-gam)) / (2.-gam)) * Dnorm
else Anorm = LEDD * mass * MPC_2 * ERG * Np / PI2 / 2. / log(Ec/E0) * Dnorm
*/
for(ie = 0; ie < ne; ie++){
  far[ie] *= Anorm;
  if (polar == 1) {
    qar[ie] *= Anorm;
    uar[ie] *= Anorm;
    var[ie] *= Anorm;
  }
}
// Let's add primary flux to the solution
refl_ratio=-1.;
if (NpNr > 0.) {  
// let's compute the cut-off powerlaw with the XSPEC routine cutoffPowerLaw
  for(ie = 0; ie <= ne; ie++) ear1[ie] = (double) ear[ie];
  param1[0] = (double) gam;
  param1[1] = (double) g_L * zzshift * Ec;
  cutoffpl(ear1, ne, param1, photar1);
  Anorm *= NpNr * transf_o * pow(g_L * zzshift, gam);
  flux_refl = flux_prim = 0.;
  for(ie = 0; ie < ne; ie++) 
    if (ear[ie] > g_L * zzshift * E0){
      flux_refl += far[ie];
      flux_prim += Anorm * photar1[ie];
      far[ie] += Anorm * photar1[ie];
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
// final spectrum output -- write ear() and photar() into file:
fw = fopen("kynlpcr_photar.dat", "w");
if( NpNr > 0. )
  for (ie = 0; ie < ne; ie++) fprintf(fw, "%14.6f\t%E\t%E\n", 
    0.5*(ear[ie]+ear[ie+1]), 
    (photar[ie]-Anorm * photar1[ie]) / (ear[ie+1] - ear[ie]),
    Anorm * photar1[ie] / (ear[ie+1] - ear[ie]));  
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
void emis_KYNlpcr(double** ear_loc, const int ne_loc, const int nt, 
                  double *far_loc, double *qar_loc, double *uar_loc, 
                  double *var_loc, const double r, const double phi, 
                  const double cosmu, double phiphoton, 
                  const double alpha_o, const double beta_o, 
                  const double delay, const double g) {
// local emissivity --> far_loc(:) array
// local energy array --> ear_loc()
// ne_loc --> number of points in local energy where the local photon flux
// density in keV is defined;
// this is a steady model (nt = 1);
// disc surface in polar coords r, phi;
// cosine of local emission angle --> cosmu

double azim, rq, factor, gfactor, cosmu0, phiphoton0, lensing, f_c, f_l;
double utmp, utmp1, ttmp, ttmp1, wtmp, wtmp1, vtmp, vtmp1;
double y1, y2, y3, y4, y5, y6, y7, y8;
double m0, m02, m, m2, fil, fir, fiu, poldeg, chi2;
int    iazim0, icosi0, icose0, imin, imax, ir0;
int    i, ie;

*ear_loc = energy;
//given cosmu, find corresponding indices in cose:
icose0 = 0;
for(i = 1; i <= ncose; i++) {
  if (cosmu >= cose[i-1]) icose0 = i;
  else break;
}
if (h_rh >= 0) {
  factor = h * sqrt(1. - 2. * h / (h * h + am2)) / (r * r + h * h) / r / cosmu;
// given r, find corresponding indices in radius:
  imin = 0;
  imax = nrad;
  ir0 = nrad / 2;
  while ((imax - imin) > 1) {
    if (r >= (radius[ir0 - 1] + r_plus)) imin = ir0;
    else imax = ir0;
    ir0 = (imin + imax) / 2;
  }
  if((ir0 == nrad) && (r == radius[nrad-1] + r_plus)) ir0--;
  if ((ir0 == 0) || (ir0 >= nrad)) {
    for(ie = 0; ie < ne_loc; ie++) far_loc[ie] = 0.;
    if (polar)
      for(ie = 0; ie < ne_loc; ie++) {
        qar_loc[ie] = 0.;
        uar_loc[ie] = 0.;
        var_loc[ie] = 0.;
      }
  }
  else {
    ttmp = (r - radius[ir0 - 1] - r_plus) / (radius[ir0] - radius[ir0 - 1]);
    ttmp1 = 1. - ttmp;
// Let's interpolate gfactor between two radii
    gfactor = ttmp * gfac[ir0] + ttmp1 * gfac[ir0-1];
// Let's interpolate cosmu0 between two radii
    cosmu0 = ttmp * cosin[ir0] + ttmp1 * cosin[ir0 - 1];
// Let's interpolate photon0 between two radii
    y1 = phiph[ir0 - 1];
    y2 = phiph[ir0];
    if (fabs(y2 - y1) > PI) {
      if (y1 < 0.) y1+=PI2;
      if (y2 < 0.) y2+=PI2;
    }
    phiphoton0 = (ttmp1 * y1 + ttmp * y2);
    if (phiphoton0 > PI) phiphoton0-=PI2;
    if (phiphoton0 < -PI) phiphoton0+=PI2;
// Let's interpolate lensing between two radii
    lensing = ttmp * transf_d[ir0] + ttmp1 * transf_d[ir0 - 1];
// Let's find the rel. azimuthal angle between inc. and refl. rays
    if (fabs(phiphoton - phiphoton0) > PI) {
      if (phiphoton < 0.) phiphoton+=PI2;
      if (phiphoton0 < 0.) phiphoton0+=PI2;
    }
    azim = phiphoton - phiphoton0;
    if (azim > PI) azim-=PI2;
    if (azim < -PI) azim+=PI2;
// given azimuth, find corresponding indices in azim:
    iazim0 = 0;
    for (i = 1; i <= nazim; i++) {
      if ((fabs(azim) / PI * 180.) >= azimuth[i-1]) iazim0 = i;
      else break;
    }
// given cosmu0, find corresponding indices in cosi:
    icosi0 = 0;
    for (i = 1; i <= ncosi; i++) {
      if (cosmu0 >= cosi[i-1]) icosi0 = i;
      else break;
    }
    factor = factor * pow(gfactor, gam - 1.) * lensing;
// cose interpolation
    if (icose0 == 0) {
      icose0 = 1;
      utmp = 0.;
      utmp1 = 1.;
    }
    else {
      if (icose0 == ncose) {
        icose0 = ncose - 1;
        utmp = 1.;
        utmp1 = 0.;
      }
      else {
        utmp= (cosmu - cose[icose0-1]) / (cose[icose0] - cose[icose0 - 1]);
        utmp1 = 1. - utmp;
      }
    }
/*
// cosi interpolation
    if (icosi0 == 0) {
      icosi0 = 1;
      vtmp = 0.;
      vtmp1 = 1.;
    }
    else {
      if (icosi0 == ncosi) {
        icosi0 = ncosi - 1;
        vtmp = 1.;
        vtmp1 = 0.;
      }
      else {
        vtmp = (cosmu0 - cosi[icosi0 - 1]) / (cosi[icosi0] - cosi[icosi0 - 1]);
        vtmp1 = 1. - vtmp;
      }
    }
*/
    vtmp = (cosmu0 - cosi[icosi0 - 1]) / (cosi[icosi0] - cosi[icosi0 - 1]);
    vtmp1 = 1. - vtmp;
    if (iazim0 == 0) {
      iazim0 = 1;
      wtmp = 0.;
      wtmp1 = 1.;
    }
    else {
      if (iazim0 == nazim) {
        iazim0 = nazim - 1;
        wtmp = 1.;
        wtmp1 = 0.;
      }
      else {
        wtmp = (fabs(azim) / PI * 180. - azimuth[iazim0 - 1]) / 
               (azimuth[iazim0] - azimuth[iazim0 - 1]);
        wtmp1 = 1. - wtmp;
      }
    }
    for (ie = 0; ie < ne_loc; ie++) {
      if (line >= 0) {
        y1 = flux_c[ie + ne_loc * (iazim0 - 1) + ne_loc * nazim * (icosi0 - 1) + 
                    ne_loc * nazim * ncosi * (icose0 - 1)];
        y2 = flux_c[ie + ne_loc * (iazim0 - 1) + ne_loc * nazim * icosi0 + 
                    ne_loc * nazim * ncosi * (icose0 - 1)];
        y3 = flux_c[ie + ne_loc * (iazim0 - 1) + ne_loc * nazim * icosi0 + 
                    ne_loc * nazim * ncosi * icose0];
        y4 = flux_c[ie + ne_loc * (iazim0 - 1) + ne_loc * nazim * (icosi0 - 1) + 
                    ne_loc * nazim * ncosi * icose0];
        y5 = flux_c[ie + ne_loc * iazim0 + ne_loc * nazim * (icosi0 - 1) + 
                    ne_loc * nazim * ncosi * (icose0 - 1)];
        y6 = flux_c[ie + ne_loc * iazim0 + ne_loc * nazim * icosi0 + 
                    ne_loc * nazim * ncosi * (icose0 - 1)];
        y7 = flux_c[ie + ne_loc * iazim0 + ne_loc * nazim * icosi0 + 
                    ne_loc * nazim * ncosi * icose0];
        y8 = flux_c[ie + ne_loc * iazim0 + ne_loc * nazim * (icosi0 - 1) + 
                    ne_loc * nazim * ncosi * icose0];
        f_c = factor * (wtmp1 * (utmp1 * (vtmp1 * y1 + vtmp * y2) + 
                                 utmp * (vtmp * y3 + vtmp1 * y4)) + 
                        wtmp * (utmp1 * (vtmp1 * y5 + vtmp * y6) + 
                                utmp * (vtmp * y7 + vtmp1 * y8)));
      }
      else f_c = 0.;
      if (line != 0) {
        y1 = flux_l[ie + ne_loc * (iazim0 - 1) + ne_loc * nazim * (icosi0 - 1) + 
                    ne_loc * nazim * ncosi * (icose0 - 1)];
        y2 = flux_l[ie + ne_loc * (iazim0 - 1) + ne_loc * nazim * icosi0 + 
                    ne_loc * nazim * ncosi * (icose0 - 1)];
        y3 = flux_l[ie + ne_loc * (iazim0 - 1) + ne_loc * nazim * icosi0 + 
                    ne_loc * nazim * ncosi * icose0];
        y4 = flux_l[ie + ne_loc * (iazim0 - 1) + ne_loc * nazim * (icosi0 - 1) + 
                    ne_loc * nazim * ncosi * icose0];
        y5 = flux_l[ie + ne_loc * iazim0 + ne_loc * nazim * (icosi0 - 1) + 
                    ne_loc * nazim * ncosi * (icose0 - 1)];
        y6 = flux_l[ie + ne_loc * iazim0 + ne_loc * nazim * icosi0 + 
                    ne_loc * nazim *ncosi * (icose0 - 1)];
        y7 = flux_l[ie + ne_loc * iazim0 + ne_loc * nazim * icosi0 + 
                    ne_loc * nazim * ncosi * icose0];
        y8 = flux_l[ie + ne_loc * iazim0 + ne_loc * nazim * (icosi0 - 1) + 
                    ne_loc * nazim * ncosi * icose0];
        f_l = factor * (wtmp1 * (utmp1 * (vtmp1 * y1 + vtmp * y2) + 
                                 utmp * (vtmp * y3 + vtmp1 * y4)) + 
                        wtmp * (utmp1 * (vtmp1 * y5 + vtmp * y6) + 
                                         utmp * (vtmp * y7 + vtmp1 * y8)));
      }
      else f_l = 0.;
      far_loc[ie] = f_c + f_l;
      if (polar) {
        m0 = cosmu0;
        m02 = cosmu0 * cosmu0;
        m = cosmu;
        m2 = cosmu * cosmu;
// full Chandrasekhar's formulae
        fil = (m2 * (1. + m02) + 2. * (1. - m2) * (1. - m02) - 
              4. * m * m0 * sqrt((1. - m2) * (1. - m02)) * cos(azim) - 
              m2 * (1. - m02) * cos(2. * azim));
        fir = (1. + m02 + (1. - m02) * cos(2. * azim));
        fiu = (-4. * m0 * sqrt((1. - m2) * (1. - m02)) * sin(azim) - 
              2. * m * (1. - m02) * sin(2. * azim));
        poldeg = sqrt((fil - fir) * (fil - fir) + fiu * fiu) / (fil + fir);
        chi2 = atan2(fiu, fil - fir);
        qar_loc[ie] = poldeg * f_c * cos(chi2);
        uar_loc[ie] = poldeg * f_c * sin(chi2);
        var_loc[ie] = 0.;
      }
    }
  }
}
else {
  rq = pow(r, h) / cosmu;
/*
// cose interpolation
    if (icose0 == 0) {
      icose0 = 1;
      utmp = 0.;
      utmp1 = 1.;
    }
    else {
      if (icose0 == ncose) {
        icose0 = ncose - 1;
        utmp = 1.;
        utmp1 = 0.;
      }
      else {
        utmp= (cosmu - cose[icose0-1]) / (cose[icose0] - cose[icose0 - 1]);
        utmp1 = 1. - utmp;
      }
    }
*/
  utmp = (cosmu - cose[icose0 - 1]) / (cose[icose0] - cose[icose0 - 1]);
  utmp1 = 1. - utmp;
  for(ie = 0; ie < ne_loc; ie++) {
    y1 = flux_c[ie + ne_loc * (icose0 - 1)];
    y2 = flux_c[ie + ne_loc * icose0];
    f_c = (utmp1 * y1 + utmp * y2) * rq;
    if (line > 0) {
      y1 = flux_l[ie + ne_loc * (icose0 - 1)];
      y2 = flux_l[ie + ne_loc * icose0];
      f_l = (utmp1 * y1 + utmp * y2) * rq;
    }
    else f_l = 0.;
    far_loc[ie] = f_c + f_l;
  }
  if (polar)
    for(ie = 0; ie < ne_loc; ie++) {
      qar_loc[ie] = far_loc[ie];
      uar_loc[ie] = 0.;
      var_loc[ie] = 0.;
    }
}
return;
}
