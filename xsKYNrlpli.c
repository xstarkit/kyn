/* KYNrlpli - relativistic line in lamp-post geometry - non-axisymmetric version
 *            model subroutine for XSPEC
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
 * S.Hledik & Z.Stuchlik, Opava. [astro-ph/0407330]
 * 
 * Dovciak M. (2004). Radiation of accretion discs in strong gravity. Faculty of
 * Mathematics and Physics, Charles University, Prague. PhD thesis.
 * [astro-ph/0411605]
 * -----------------------------------------------------------------------------
 * 
 * This subroutine computes the emission from an acrretion disc that is
 * illuminated from the primary power-law source located on the axis above the
 * central black hole. All relativistic effects are taken into account (in both
 * parts of the light path - from the primary source to the disc and from the
 * disc to the observer). This subroutine calls subroutine ide() for integrating
 * local emission over the disc and uses the fits file 'KBHtablesNN.fits'
 * defining the transfer functions needed for the integration. For details on
 * ide() and the fits file see the subroutine ide(). This subroutine uses the
 * FITS 'KBHlamp80.fits' file with transfer functions for the light coming from
 * primary source and hitting the accretion disc (see the description of this 
 * file below). The fluorescent iron line is modelled by Monte Carlo code NOAR 
 * and is stored in the FITS file 'fluorescent_line.fits' (see the description 
 * of this file below).
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
 * par8  ... height - height on the axis (measured from the center) at which
 *                    the primary source is located (GM/c^2)
 * par9  ... PhoIndex - power-law energy index of the primary flux
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
 *                    = 0 - array of photon number density flux per bin
 *                         (array of Stokes parameter I devided by energy)
 *                          with the polarisation computations switched off
 *                    = 1 - array of photon number density flux per bin
 *                         (array of Stokes parameter I devided by energy),
 *                          with the polarisation computations switched on
 *                    = 2 - array of Stokes parameter Q devided by energy
 *                    = 3 - array of Stokes parameter U devided by energy
 *                    = 4 - array of Stokes parameter V devided by energy
 *                    = 5 - array of degree of polarization
 *                    = 6 - array of polarization angle psi=0.5*atan(U/Q)
 *                    = 7 - array of "Stokes" angle
 *                          beta=0.5*asin(V/sqrt(Q*Q+U*U+V*V))
 * par20 ... nthreads - number of threads to be used for computations
 * par21 ... normtype - how to normalize the spectra
 *                      = 0: normalization to the total photon flux
 *                      > 0: normalization to the photon flux at 'par21' keV
 *                      = -1: the photon flux is not re-normalized,
 *                      = -2: normalization to the maximum of the photon flux
 *
 * NOTES:
 *  -> accuracy vs. speed trade off depends mainly on: nrad, nphi
 *
 *  -> in this model it is assumed that local emission is completely
 *     linearly polarized in the direction perpendicular to the disc
 *
 * -----------------------------------------------------------------------------
 * 
 * KBHlamp80.fits
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
 * There are seven functions stored in the KBHlamp80.fits file as columns. 
 * They are parametrized by the value of the horizon of the black hole (FITS 
 * table extensions), height of the primary source (rows) and either the
 * inclination angles of the observer (elements) or the radius in GM/c^2
 * at which a photon strikes the disc (elements). All these are defined
 * as vectors at the beginning of the FITS file. The tables are defined for
 * horizon radius r_horizon = 1.00, 1.05, ..., 1.95, 2.00,
 * height h-r_horizon = 0.1 - 100 (100 values with exponentially growing step),
 * inclination = 0.1, 1, 5, 10, 15, ..., 80, 85, 89 and
 * radius r-r_horizon = 0.01 - 1000 (100 values with exponentially growing step).
 * 
 * The functions included are:
 * q2_a   - constant of motion, q^2, defining the photon angular momentum 
 *          between the axis and the observer
 * dWadWo - amplification of the primary source flux:
 *        = sin(theta_axis_local) / sin(theta_observer) *
 *          dtheta_axis_local/dtheta_observer
 * delay_a - delay between the emission of the photon by the lamp and detecting
 *           it by the observer
 * q  - constant of motion defining the photon angular momentum between the axis
 *      and the disc
 * pr - the radial component of the photon momentum at the disc
 * dWadSd - part of the amplification of the incident flux on the disc from the
 *          primary source:
 *        = sin(theta_axis_local) * dtheta_axis_local / dtheta_fake, with the
 *          theta_fake defined as tan(theta_fake) = radius_incident / height
 *          one has to multiply by h/(r^2+h^2) to get the full amplification
 * delay - delay between the emission of the photon by the lamp and falling
 *         down at the disc
 * 
 * The rest of the functions below are computed from the above ones:
 * g-factor - the ratio of the energy of a photon hitting the disc to the energy
 *            of the same photon when emitted from a primary source,
 * cosine of the incident angle - an absolute value of the cosine of the local
 *                                incident angle between the incident light ray
 *                                and local disc normal
 * azimuthal incident angle - the angle (in radians) between the projection of
 *                            the three-momentum of the incident photon into the
 *                            disc (in the local rest frame co-moving with the
 *                            disc) and the radial tetrad vector.
 *
 * The definition of the file KBHlamp80.fits:
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
 *    (q2_a, dWadWo, delay_a) or incident radius (the rest of the functions).
 * 6. Each of these extensions has seven columns. In each column, a particular
 *    function is stored - q2_a, dWadWo, delay_a, q, pr, dWadSd and delay, 
 *    respectively.
 * 
 * -----------------------------------------------------------------------------
 *
 * fluorescent_line.fits
 *
 * Dependence of the reflection coefficients on the angle of incidence, angle
 * of emission is stored in this FITS file. The emission is induced by a 
 * power-law incident radiation. Emission line was computed with NOAR code.
 * 
 * There are several binary extensions in this fits file:
 * - the first extension contains primary source Photon index values 
 *   (1.1, 1.2, ..., 3.0)
 * - the second extension contains the absolute values of the cosine of the
 *   incident angles
 * - the third extension contains the values of the cosine of the emission
 *   angles
 * - in the forth and second extensions the two multiplicative parts of the 
 *   emission line are stored (the first one depends only on incident and 
 *   emission angles, the second one depends on photon index and emission angle)
 *
 *******************************************************************************
 *
 * 12.11.2009  changed to work with new lamp-post tables, we can fit for height
 *             as well now
 * 28. 3.2011  changed to work with lamp-post tables in q and p^r
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
#define NPARAM 21
#define NE     200
#define E_MIN  0.
#define E_MAX  10.

int main() {

void KYNrlpli(const double *ear, int ne, const double *param, int ifl, 
              double *photar, double *photer, const char* init);

double ear[NE+1], photar[NE], photer[NE], param[NPARAM];
char   initstr[0] = "";
int    ie;

param[ 0] = 1.;       // a/M
param[ 1] = 30.;      // theta_o
param[ 2] = 1.;       // rin
param[ 3] = 1.;       // ms
param[ 4] = 1000.;    // rout
param[ 5] = 0.;       // phi0
param[ 6] = 360.;     // dphi
param[ 7] = 3.;       // height
param[ 8] = 2.;       // PhoIndex
param[ 9] = -3.;      // alpha
param[10] = 0.;       // beta
param[11] = 0.;       // rcloud
param[12] = 0.;       // zshift
param[13] = 80.;      // ntable
param[14] = 1000.;    // nrad
param[15] = 1.;       // division
param[16] = 720.;     // nphi
param[17] = 1.;       // smooth
param[18] = 0.;       // Stokes
param[19] = 2.;       // nthreads
param[20] = 0.;       // normtype

for (ie = 0; ie <= NE; ie++) {
  ear[ie] = E_MIN + ie * (E_MAX-E_MIN) / NE;
//  ear[ie] = E_MIN * pow(E_MAX / E_MIN, ((double) ie) / NE);
}

KYNrlpli(ear, NE, param, IFL, photar, photer, initstr);

return(0);
}
#endif
/*******************************************************************************
*******************************************************************************/

#define NE_LOC 3
#define LAMP "KBHlamp80.fits"
#define FLUORESCENT_LINE "fluorescent_line.fits"
#define PI 3.14159265359

/* Let's declare variables that are common for the main and emissivity 
   subroutines */
static double *energy, *cosin, *phiph, *transf_d, *gfac;
static float  *radius, *cosi, *cose, *flux1, *flux2;
static double h, h_rh, transf_o, am2, r_plus, gam, gtmp, gtmp1;
static long   ncose, ncosi, nrad;
static int    polar, igam0;

extern char* FGMODF(void);
extern char* FGMSTR(char* dname);
extern int xs_write(char* wrtstr, int idest);

void KYNrlpli(const double *ear, int ne, const double *param, int ifl, 
              double *photar, double *photer, const char* init) {

extern int ide(const double *ear, const int ne, const int nt, double *far, 
               double *qar, double *uar, double *var, 
               const double *ide_param, void (*emissivity)(), 
               const int ne_loc);

void emis_KYNrlpli(double** ear_loc, const int ne_loc, const int nt, 
                   double *far_loc, double *qar_loc, double *uar_loc, 
                   double *var_loc, const double r, const double phi, 
                   const double cosmu, const double phiphoton, 
                   const double alpha_o, const double beta_o, 
                   const double delay, const double g);

/* Let's declare static variables whose values should be remembered for next
   run in XSPEC */
static char   kydir[255] = "";
static char   pname[128] = "KYDIR";
static long   nrh, ngamma, nh, nincl;
static float  *r_horizon, *gamma, *height, *incl, *dWadWo, *dWadSd, *q, *pr;
static double h_rh_old = -1., am_old = -1., thetaO_old = -1.;
static int    first = 1, first_h = 1;

FILE   *fw;
char   errortxt[80];
double ide_param[25];
double far[ne], qar[ne], uar[ne], var[ne], pd[ne], pa[ne], pa2[ne];
double ttmp, ttmp1, utmp, utmp1, vtmp, vtmp1;
double pr_final, q_final, pom, pom1, pom2, pom3, y1, y2, y3, y4, y5, y6, y7, y8;
double r, r2, delta, ULt, rms, tmp1, Ut, U_phi, U_r, Ur, UrU_r, Lms, Ems;
double am, thetaO;
double pamin, pamax, pa2min, pa2max;
int    imin, imax, irh0, ih0, ith0;
int    ie, i, stokes;
//int    j, irh, ih, ir;

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
// a/M - black hole angular momentum
ide_param[0] = param[0];
am = ide_param[0];
am2 = am * am;
pom1 = pow(1. + am, 1. / 3.);
pom2 = pow(1. - am, 1. / 3.);
pom3 = pow(1. - am2, 1. / 3.);
pom = 1. + pom3 * (pom1 + pom2);
pom1 = sqrt(3. * am2 + pom * pom);
if (am >= 0) rms= 3. + pom1 - sqrt((3. - pom) * (3. + pom + 2. * pom1));
else rms = 3. + pom1 + sqrt((3. - pom) * (3. + pom + 2. * pom1));
r_plus= 1. + sqrt(1. - am2);
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
ide_param[7] = param[14];
// division - type of division in r integration (0-equidistant, 1-exponential)
ide_param[8] = param[15];
// nphi - number of grid points in azimuth
ide_param[9] = param[16];
// smooth - whether to smooth the resulting spectrum (0-no, 1-yes)
ide_param[10] = param[17];
// normal - how to normalize the final spectrum
ide_param[11] = param[20];
//ide_param[11] = -1.;
//ide_param[11] = -2.;
// zshift - overall Doppler shift
ide_param[12] = param[12];
// ntable - table model (defines fits file with tables)
ide_param[13] = param[13];
// PhoIndex - power-law energy index of the lamp emission
gam = param[8];
// height - height of the lamp above the black hole
h = param[7];
if (h >= 0.) {
  if (h < r_plus) h_rh = 0.;
  else h_rh = h - r_plus;
} else h_rh = h;
// edivision - type of division in local energies (0-equidistant, 1-exponential)
ide_param[14] = 1.;
// periodic and dt need not to be set for nt=1
// (ide_param[15], ide_param[16])
// polar - whether we need value of change in polarization angle (0-no,1-yes)
stokes = (int) param[18];
if ((stokes < 0) || (stokes > 7)) {
  xs_write("kynrlpli: Stokes has to be 0-7", 5);
  for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
  return;
}
polar = 0;
if (stokes > 0) polar = 1;
ide_param[17] = polar;
// delay_r and delay_phi are not used
// (ide_param[18], ide_param[19])
// number of threads for multithread computations
ide_param[20] = param[19];
// alpha_cloud - position of the cloud in alpha impact parameter (in GM/c^2)
ide_param[21] = param[9];
// beta_cloud - position of the cloud in beta impact parameter (in GM/c^2)
ide_param[22] = param[10];
// rcloud - radius of the cloud (in GM/c^2)
ide_param[23] = param[11];
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
// Let's read the 'KBHlamp80' fits file
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
    xs_write("\nkynrlpli: set the KYDIR to the directory with the KY tables",5);
    for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
    return;
  }
// Let's read tables (binary tables => hdutype=2)
// Move to the extension 'r_horizon' and read its values
  ffmrhd(fptr, 1, &hdutype, &status);
  ffgnrw(fptr, &nrh, &status);
//******************************************************************************
//  fprintf(stdout,"nrh = %d\n",nrh);
//******************************************************************************   
// Allocate memory for r_horizon...
  if ((r_horizon = (float *) malloc(nrh * sizeof(float))) == NULL) {
    xs_write("kynrlpli: Failed to allocate memory for tmp arrays.", 5);
    for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
    return;
  }
// Read the data in the 'r_horizon' table
  nelems = nrh;
// FTGCVE reads the VALUES from the first column.
  ffgcv(fptr, TFLOAT, colnum, frow, felem, nelems, &float_nulval, r_horizon,
        &anynul, &status);
//******************************************************************************
//  for ( i=0; i<nrh; i++)fprintf(stdout,"%f\n",r_horizon[i]);
//******************************************************************************   
// Move to the extension 'height' and read its values
  ffmrhd(fptr, 1, &hdutype, &status);
  ffgnrw(fptr, &nh, &status);
//******************************************************************************
//  fprintf(stdout,"nh = %d\n",nh);
//******************************************************************************   
// Allocate memory for height...
  if ((height = (float *) malloc(nh * sizeof(float))) == NULL) {
    xs_write("kynrlpli: Failed to allocate memory for tmp arrays.", 5);
    for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
    return;
  }
// Read the data in the 'height' table
  nelems = nh;
// FTGCVE reads the VALUES from the first column.
  ffgcv(fptr, TFLOAT, colnum, frow, felem, nelems, &float_nulval, height,
        &anynul, &status);
//******************************************************************************
//  for ( i=0; i<nh; i++)fprintf(stdout,"%f\n",height[i]);
//******************************************************************************   
// Move to the extension 'inclination' and read its values
  ffmrhd(fptr, 1, &hdutype, &status);
  ffgnrw(fptr, &nincl, &status);
//******************************************************************************
//  fprintf(stdout,"nincl = %d\n",nincl);
//******************************************************************************   
// Allocate memory for inclination...
  if ((incl = (float *) malloc(nincl * sizeof(float))) == NULL) {
    xs_write("kynrlpli: Failed to allocate memory for tmp arrays.", 5);
    for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
    return;
  }
// Read the data in the 'inclination' table
  nelems = nincl;
// FTGCVE reads the VALUES from the first column.
  ffgcv(fptr, TFLOAT, colnum, frow, felem, nelems, &float_nulval, incl,
        &anynul, &status);
//******************************************************************************
//  for ( i=0; i<nincl; i++)fprintf(stdout,"%f\n",incl[i]);
//******************************************************************************   
// Move to the extension 'r_rh' and read its values
  ffmrhd(fptr, 1, &hdutype, &status);
  ffgnrw(fptr, &nrad, &status);
//******************************************************************************
//  fprintf(stdout,"nrad = %d\n",nrad);
//******************************************************************************   
// Allocate memory for radius...
  if ((radius = (float *) malloc(nrad * sizeof(float))) == NULL) {
    xs_write("kynrlpli: Failed to allocate memory for tmp arrays.", 5);
    for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
    return;
  }
// Read the data in the 'r_rh' table
  nelems = nrad;
// FTGCVE reads the VALUES from the first column.
  ffgcv(fptr, TFLOAT, colnum, frow, felem, nelems, &float_nulval, radius,
        &anynul, &status);
//******************************************************************************
//  for ( i=0; i<nrad; i++)fprintf(stdout,"%f\n",radius[i]);
//******************************************************************************   
// Let's read the tables for dWadWo, q, p^r and dWadSd
// allocate memory for the arrays
  if ((dWadWo = (float *) malloc(nincl * nh * nrh * sizeof(float))) == NULL) {
    xs_write("kynrlpli: Failed to allocate memory for tmp arrays.", 5);
    for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
    return;
  }
  if ((q = (float *) malloc(nrad * nh * nrh * sizeof(float))) == NULL) {
    xs_write("kynrlpli: Failed to allocate memory for tmp arrays.", 5);
    for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
    return;
  }
  if ((pr = (float *) malloc(nrad * nh * nrh * sizeof(float))) == NULL) {
    xs_write("kynrlpli: Failed to allocate memory for tmp arrays.", 5);
    for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
    return;
  }
  if ((dWadSd = (float *) malloc(nrad * nh * nrh * sizeof(float))) == NULL) {
    xs_write("kynrlpli: Failed to allocate memory for tmp arrays.", 5);
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
      if ((nh - irow) < nrow) {
        nelements1 = (nh - irow) * nincl;
        nelements2 = (nh - irow) * nrad;
      }
      ffgcv(fptr, TFLOAT, 2, irow+1, 1, nelements1, &float_nulval, 
            &dWadWo[irow * nincl + nh * nincl * ihorizon],
            &anynul, &status);
      ffgcv(fptr, TFLOAT, 4, irow+1, 1, nelements2, &float_nulval, 
            &q[irow * nrad + nh * nrad * ihorizon],
            &anynul, &status);
      ffgcv(fptr, TFLOAT, 5, irow+1, 1, nelements2, &float_nulval, 
            &pr[irow * nrad + nh * nrad * ihorizon],
            &anynul, &status);
      ffgcv(fptr, TFLOAT, 6, irow+1, 1, nelements2, &float_nulval, 
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
    q[i+nrad*ih+nrad*nh*irh],pr[i+nrad*ih+nrad*nh*ir],
    dWadSd[i+nrad*ih+nrad*nh*irh]);
*******************************************************************************/
// Firstly we have to free allocated memory for the arrays transf_o, gfac,
// cosin, phiph, transf_d
  if ((gfac = (double *) malloc(nrad * sizeof(double))) == NULL) {
    xs_write("kynrlpli: Failed to allocate memory for tmp arrays.", 5);
    for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
    return;
  }
  if ((cosin = (double *) malloc(nrad * sizeof(double))) == NULL) {
    xs_write("kynrlpli: Failed to allocate memory for tmp arrays.", 5);
    for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
    return;
  }
  if ((phiph = (double *) malloc(nrad * sizeof(double))) == NULL) {
    xs_write("kynrlpli: Failed to allocate memory for tmp arrays.", 5);
    for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
    return;
  }
  if ((transf_d = (double *) malloc(nrad * sizeof(double))) == NULL) {
    xs_write("kynrlpli: Failed to allocate memory for tmp arrays.", 5);
    for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
    return;
  }
  first_h = 0;
}
/******************************************************************************/
if (h >= 0.) {
  if (h_rh > height[nh - 1]) {
    sprintf(errortxt, "The height must be lower than or equal to %f\n",
      height[nh - 1] + r_plus);
    xs_write(errortxt, 5);
    for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
    return;
  }
  if (h_rh < height[0]) {
    sprintf(errortxt, "The height is too low, we set it to %f\n",
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
    if (r_plus >= r_horizon[irh0-1]) imin = irh0;
    else imax = irh0;
    irh0 = (imin + imax) / 2;
  }
  if (irh0 == 0) irh0 = 1;
//if ((imax == nrh) && (r_plus > r_horizon[nrh-1])) irh0 = nrh;
  ttmp = (r_plus - r_horizon[irh0 - 1])
         / (r_horizon[irh0] - r_horizon[irh0 - 1]);
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
//if ((imax == nh) && (h_rh > height[nh-1])) ih0 = nh;
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
      q_final = utmp1 * (ttmp1 * y1 + ttmp *y2) +
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
      transf_d[i] = utmp1 * (ttmp1 * y1 + ttmp * y2) + utmp
                    * (ttmp * y3 + ttmp1 * y4);
    }
  }
}
/******************************************************************************/
if (first) {
// Let's create energy array
  if ((energy = (double *) malloc((NE_LOC + 1) * sizeof(double))) == NULL) {
    xs_write("kynrlpli: Failed to allocate memory for tmp arrays.", 5);
    for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
    return;
  }
  energy[0] = 6.4 - 0.5e-3;
  energy[1] = 6.4;
  energy[2] = 6.4 + 0.5e-3;
// Let's read the 'fluorescent_line' fits file
// The status parameter must always be initialized.
  status = 0;
// Set up the directory and the fits file name of the tables
// and open the fits file with the tables
  sprintf(kydir, "%s", FGMSTR(pname));
  if (strlen(kydir) == 0) sprintf(tables_file, "./%s", FLUORESCENT_LINE);
  else if (kydir[strlen(kydir) - 1] == '/') sprintf(tables_file, "%s%s",
                                                    kydir, FLUORESCENT_LINE);
  else sprintf(tables_file, "%s/%s", kydir, FLUORESCENT_LINE);
// Let's read the 'fluorescent_line' fits file
// The status parameter must always be initialized.
  status = 0;
  ffopen(&fptr, tables_file, READONLY, &status);
  if(status) {
    sprintf(tables_file, "%s%s", FGMODF(), FLUORESCENT_LINE);
    status = 0;
    ffopen(&fptr, tables_file, READONLY, &status);
  }
  if (status){
    if (status) ffrprt(stderr, status);
    ffclos(fptr, &status);
    xs_write("\nkynrlpli: set the KYDIR to the directory with the KY tables",5);
    for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
    return;
  }
// Let's read tables (binary tables => hdutype=2)
// Move to the extension 'Gamma' and read its values
  ffmrhd(fptr, 1, &hdutype, &status);
  ffgnrw(fptr, &ngamma, &status);
//******************************************************************************
//  fprintf(stdout,"ngamma = %d\n",ngamma);
//******************************************************************************   
// Allocate memory for gamma...
  if ((gamma = (float *) malloc(ngamma * sizeof(float))) == NULL) {
    xs_write("kynrlpli: Failed to allocate memory for tmp arrays.", 5);
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
// Move to the extension 'cos_theta_i' and read its values
  ffmrhd(fptr, 1, &hdutype, &status);
  ffgnrw(fptr, &ncosi, &status);
//******************************************************************************
//  fprintf(stdout,"ncosi = %d\n",ncosi);
//******************************************************************************   
// Allocate memory for cosi...
  if ((cosi = (float *) malloc(ncosi * sizeof(float))) == NULL) {
    xs_write("kynrlpli: Failed to allocate memory for tmp arrays.", 5);
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
    xs_write("kynrlpli: Failed to allocate memory for tmp arrays.", 5);
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
// Move to the extension for flux part 1
  ffmrhd(fptr, 1, &hdutype, &status);
// Allocate memory for flux1
  if ((flux1 = (float *) malloc(ncosi * ncose * sizeof(float))) == NULL) {
    xs_write("kynrlpli: Failed to allocate memory for tmp arrays.", 5);
    for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
    return;
  }
// Read the data in the table
  nelems = ncosi * ncose;
// FTGCVE reads the VALUES from the first column.
  ffgcv(fptr, TFLOAT, 1, frow, felem, nelems, &float_nulval, flux1,
        &anynul, &status);
//******************************************************************************
//  for ( j=0; j<ncosi; j++){
//    for ( i=0; i<ncose; i++) fprintf(stdout,"%f\t",flux1[i+ncose*j]);
//    fprintf(stdout,"\n");
//  }
//******************************************************************************   
// Move to the extension for flux part 2
  ffmrhd(fptr, 1, &hdutype, &status);
// Allocate memory for and flux2...
  if ((flux2 = (float *) malloc(ngamma * ncose * sizeof(float))) == NULL) {
    xs_write("kynrlpli: Failed to allocate memory for tmp arrays.", 5);
    for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
    return;
  }
// Read the data in the table
  nelems = ngamma * ncose;
// FTGCVE reads the VALUES from the first column.
  ffgcv(fptr, TFLOAT, 1, frow, felem, nelems, &float_nulval, flux2,
        &anynul, &status);
//******************************************************************************
//  for ( j=0; j<ngamma; j++){
//    for ( i=0; i<ncose; i++) fprintf(stdout,"%f\t",flux2[i+ncose*j]);
//    fprintf(stdout,"\n");
//  }
//******************************************************************************   
// The FITS file must always be closed before exiting the program.
  ffclos(fptr, &status);
}
/******************************************************************************/
// given gam, find corresponding index in gamma[]:
imin = 0;
imax = ngamma;
igam0 = ngamma / 2;
while ((imax - imin) > 1) {
  if (gam >= gamma[igam0 - 1]) imin = igam0;
  else imax = igam0;
  igam0 = (imin + imax) / 2;
}
if ((igam0 == ngamma) && (gam == gamma[ngamma - 1])) igam0 = igam0 - 1;
if ((igam0 == 0) || (igam0 >= ngamma)) {
  xs_write("Wrong Gamma given, following values are allowed:", 5);
  sprintf(errortxt, "%f <= Gamma <= %f", gamma[0], gamma[ngamma - 1]);
  xs_write(errortxt, 5);
  for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
  return;
}
gtmp = (gam - gamma[igam0 - 1]) / (gamma[igam0] - gamma[igam0 - 1]);
gtmp1 = 1. - gtmp;
first = 0;
am_old = am;
thetaO_old = thetaO;
h_rh_old = h_rh;

/******************************************************************************/
#ifdef OUTSIDE_XSPEC
// let's write the input parameters to a file
fw = fopen("kynrlpli.txt", "w");
fprintf(fw, "a/M         %12.6f\n", param[0]);
fprintf(fw, "theta_o     %12.6f\n", param[1]);
fprintf(fw, "rin         %12.6f\n", param[2]);
fprintf(fw, "ms          %12d\n", (int) param[3]);
fprintf(fw, "rout        %12.6f\n", param[4]);
fprintf(fw, "phi         %12.6f\n", param[5]);
fprintf(fw, "dphi        %12.6f\n", param[6]);
fprintf(fw, "height      %12.6f\n", param[7]);
fprintf(fw, "PhoIndex    %12.6f\n", param[8]);
fprintf(fw, "alpha_c     %12.6f\n", ide_param[21]);
fprintf(fw, "beta_c      %12.6f\n", ide_param[22]);
fprintf(fw, "rcloud      %12.6f\n", ide_param[23]);
fprintf(fw, "zshift      %12.6f\n", param[12]);
fprintf(fw, "ntable      %12d\n", (int) param[13]);
fprintf(fw, "nrad        %12d\n", (int) param[14]);
fprintf(fw, "division    %12d\n", (int) param[15]);
fprintf(fw, "nphi        %12d\n", (int) param[16]);
fprintf(fw, "smooth      %12d\n", (int) param[17]);
fprintf(fw, "Stokes      %12d\n", (int) param[18]);
fprintf(fw, "polar       %12d\n", polar);
fprintf(fw, "r_horizon   %12.6f\n", r_plus);
fprintf(fw, "r_ms        %12.6f\n", rms);
fprintf(fw, "edivision   %12d\n", (int) ide_param[14]);
fprintf(fw, "ne_loc      %12d\n", NE_LOC);
fprintf(fw, "normal      %12.6f\n", ide_param[11]);
fprintf(fw, "nthreads    %12d\n", (int) ide_param[20]);
fclose(fw);
#endif
/******************************************************************************/

// Let's integrate local emission over the accretion disk
if (ide(ear, ne, 1, far, qar, uar, var, ide_param, emis_KYNrlpli, NE_LOC)) {
  for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
  return;
}
// interface with XSPEC
if (!stokes) for (ie = 0; ie < ne; ie++) photar[ie] = far[ie];
else {
// final spectrum output -- write ear[] and photar[] into file:
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
// final spectrum output -- write ear() and photar() into file:
fw = fopen("kynrlpli_photar.dat", "w");
for (ie = 0; ie < ne; ie++) fprintf(fw, "%14.6f\t%E\n", 0.5*(ear[ie]+ear[ie+1]), 
  photar[ie] / (ear[ie+1] - ear[ie]));
fclose(fw);
#endif
/******************************************************************************/

return;
}

/*******************************************************************************
*******************************************************************************/
void emis_KYNrlpli(double** ear_loc, const int ne_loc, const int nt, 
                   double *far_loc, double *qar_loc, double *uar_loc, 
                   double *var_loc, const double r, const double phi, 
                   const double cosmu, const double phiphoton, 
                   const double alpha_o, const double beta_o, 
                   const double delay, const double g) {

// the height must be positive in this model (broken-powerlaw is not implemented 
// here, one has to use KYNrline XSPEC model)
  
// local emissivity --> far_loc(:) array
// local energy array --> ear_loc()
// this is a steady model (nt = 1);
// disc surface in polar coords r, phi;
// cosine of local emission angle --> cosmu
  
double y1, y2, y3, y4, factor, gfactor, cosmu0, lensing;
double ttmp, ttmp1, utmp, utmp1, vtmp, vtmp1;
int imin, imax, ir0, icosi0, icose0, ie, i;
//double ralpha

*ear_loc = energy;
// given cosmu, find corresponding indices in cose:
icose0 = 0;
for (i = 1; i <= ncose; i++) {
  if (cosmu >= cose[i-1]) icose0 = i;
  else break;
}
//if (h_rh >= 0) {
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
  if ((ir0 == nrad) && (r == (radius[nrad - 1] + r_plus))) ir0--;
  if ((ir0 == 0) || (ir0 >= nrad)) {
    for (ie = 0; ie < ne_loc; ie++) far_loc[ie] = 0.;
    if (polar)
      for (ie = 0; ie < ne_loc; ie++) {
        qar_loc[ie] = 0.;
        uar_loc[ie] = 0.;
        var_loc[ie] = 0.;
      }
  }
  else {
    ttmp = (r - radius[ir0 - 1] - r_plus) / (radius[ir0] - radius[ir0 - 1]);
    ttmp1 = 1. - ttmp;
// Let's interpolate gfactor between two radii
    gfactor = ttmp * gfac[ir0] + ttmp1 * gfac[ir0 - 1];
// Let's interpolate cosmu0 between two radii
    cosmu0 = ttmp * cosin[ir0] + ttmp1 * cosin[ir0 - 1];
// Let's interpolate lensing between two radii
    lensing = ttmp * transf_d[ir0] + ttmp1 * transf_d[ir0 - 1];
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
        utmp = (cosmu - cose[icose0 - 1]) / (cose[icose0] - cose[icose0-1]);
        utmp1 = 1. - utmp;
      }
    }
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
    y1 = flux1[(icosi0 - 1) * ncose + (icose0 - 1)];
    y2 = flux1[icosi0 * ncose + (icose0 - 1)];
    y3 = flux1[icosi0 * ncose + icose0];
    y4 = flux1[(icosi0 - 1) * ncose + icose0];
    far_loc[1] = factor * ((utmp1 * (vtmp1 * y1 + vtmp * y2) +
                 utmp * (vtmp * y3 + vtmp1 * y4)));
    y1 = flux2[(igam0 - 1) * ncose + (icose0 - 1)];
    y2 = flux2[igam0 * ncose + (icose0 - 1)];
    y3 = flux2[igam0 * ncose + icose0];
    y4 = flux2[(igam0 - 1) * ncose + icose0];
    far_loc[1] = far_loc[1] * ((utmp1 * (gtmp1 * y1 + gtmp * y2) +
                 utmp * (gtmp * y3 + gtmp1 * y4)));
    far_loc[1] = far_loc[1] / ((*ear_loc)[2] - (*ear_loc)[0]) * 2.;    
    far_loc[2]=far_loc[0]=0.;
    if (polar) {
      for (ie = 0; ie < ne_loc; ie++) {
        qar_loc[ie] = far_loc[ie];
        uar_loc[ie] = 0.;
        var_loc[ie] = 0.;
      }
    }
  }
/*}
else {
  ralpha = pow(r, h) / cosmu;
// cose interpolation
    if (icose0 == 0) {
      icose0 = 1;
      utmp = 0;
      utmp1 = 1;
    }
    else {
      if (icose0 == ncose) {
        icose0 = ncose - 1;
        utmp = 1;
        utmp1 = 0;
      }
      else {
        utmp = (cosmu - cose[icose0 - 1]) / (cose[icose0] - cose[icose0 - 1]);
        utmp1 = 1. - utmp;
      }
    }
// The following is totally wrong, the code does not work for power-law radial
// emisson !!!
    for (ie = 0; ie < ne_loc; ie++) {
      y1 = flux1[ie + ne_loc * (icose0 - 1)];
      y2 = flux1[ie + ne_loc * icose0];
      far_loc[ie] = (utmp1 * y1 + utmp * y2) * ralpha;
    }
    if (polar) for (ie = 0; ie < ne_loc; ie++) {
                 qar_loc[ie] = far_loc[ie];
                 uar_loc[ie] = 0.;
                 var_loc[ie] = 0.;
               }
}*/
return;
}
/*******************************************************************************
*******************************************************************************/
