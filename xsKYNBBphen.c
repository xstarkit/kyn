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
 * par19 ... nthreads - number of threads to be used for computations
 *
 * NOTES:
 *  -> accuracy vs. speed trade off depends mainly on: nrad, nphi
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/*******************************************************************************
*******************************************************************************/
#ifdef OUTSIDE_XSPEC

#define NE     200
#define E_MIN  0.1
#define E_MAX  20.
#define NPARAM 19
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
param[ 4] = 400.;       // rout
param[ 5] = 0.;         // phi
param[ 6] = 360.;       // dphi
param[ 7] = 1.;         // Tin
param[ 8] = 0.75;       // BBindex
param[ 9] = -3.;        // alpha
param[10] = 0.;         // beta
param[11] = 0.;         // rcloud
param[12] = 0.;         // zshift
param[13] = 80.;        // ntable
param[14] = 200.;       // nrad
param[15] = 1.;         // division
param[16] = 180.;       // nphi
param[17] = 0.;         // smooth
param[18] = 2.;         // nthreads


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

#define PI     3.14159265358979
#define H_KEVS 4.13566743e-18
#define C_MS   2.99792458e8
#define MSOLAR 1.989e+30
#define G      6.6743e-11
// the "kpc" below is 10kpc in cm
#define KPC    3.0857e+22

static double *ener_loc, *flx;
static double Tin, BBindex, rin, theta_o, rout;
static int    polar;

extern int xs_write(char* wrtstr, int idest);

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

FILE *fw;
double ide_param[25], flux_out[ne + 1];
double far[ne], qar[ne], uar[ne], var[ne];
double am, am2, pom, pom1, pom2, pom3, rms, r_plus;
int    ne_loc, ie;


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
ide_param[1] = param[1];
theta_o = param[1];
// rin - inner edge of non-zero disc emissivity
ide_param[2] = param[2];
if (param[3] == 1.) rin = rms;
else if (param[3] < r_plus) rin = r_plus;
else rin = param[2];
// ms - whether to integrate from rin or rms
ide_param[3] = param[3];
// rout - outer edge of non-zero disc emissivity
ide_param[4] = param[4];
rout = param[4];
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
// polar - whether we need value of change in polarization angle (0-no,1-yes)
ide_param[17] = 0;
polar = 0;
// delay_r and delay_phi are not used
// (ide_param[18], ide_param[19])
// number of threads for multithread computations
ide_param[20] = param[18];
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
fw = fopen("kynbbphen.txt", "w");
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
  xs_write("kynbbphen: Failed to allocate memory for tmp arrays.", 5);
  for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
  return;
}
if ((flx = (double *) malloc((ne_loc + 1) * sizeof(double))) == NULL) {
  xs_write("kynbbphen: Failed to allocate memory for tmp arrays.", 5);
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
fw = fopen("kynbbphen_photar_loc.dat", "w");
for (ie = 0; ie < ne_loc; ie++)
  fprintf(fw, "%14.6f\t%E\n", ener_loc[ie], flx[ie]);
fclose(fw);
******************************************************************************/

if (ide(ear, ne, 1, far, qar, uar, var, ide_param, emis_BBphen, ne_loc)) {
  for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
  return;
}
outer_disc_phen(ear, ne, flux_out);

//interface with XSPEC..........................................................
for (ie = 0; ie < ne; ie++) photar[ie] = far[ie] + flux_out[ie];

/******************************************************************************/
#ifdef OUTSIDE_XSPEC
// final spectrum output -- write ear[] and photar[] into file:
fw = fopen("kynbbphen_photar.dat", "w");
for (ie = 0; ie < ne; ie++) {
  fprintf(fw, "%14.6f\t%E\t%E\t%E\n", 0.5 * (ear[ie] + ear[ie+1]), 
    photar[ie] / (ear[ie + 1] - ear[ie]),
    far[ie] / (ear[ie+1] - ear[ie]), flux_out[ie] / (ear[ie + 1] - ear[ie]));
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

void emis_BBphen(double** ear_loc, const int ne_loc, const int nt, 
                 double *far_loc, double *qar_loc, double *uar_loc, 
                 double *var_loc, const double r, const double phi, 
                 const double cosmu, const double phiphoton, 
                 const double alpha_o, const double beta_o, 
                 const double delay, const double g) {

int ie;
double flx0, flx1, g2;

*ear_loc = ener_loc;
g2 = g * g;
flx0 = flx[0] / (exp(*(*ear_loc) /
                (g * Tin * pow(r / rin, -BBindex))) - 1.);
for (ie = 1; ie <= ne_loc; ie++) {
  flx1 = flx[ie] / (exp(*(*ear_loc + ie) /
                (g * Tin * pow(r / rin, -BBindex))) - 1.);
  far_loc[ie-1] = (flx0 + flx1) / 2. *
                  ( *(*ear_loc + ie) - *(*ear_loc + ie -1) ) / g2;
  flx0 = flx1;
  if (polar) {
    qar_loc[ie] = 0.;
    uar_loc[ie] = 0.;
    var_loc[ie] = 0.;
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
norm = 4. * PI * cos(theta_o / 180. * PI) * rin* rin * pow(Tin, 2. / BBindex) /
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
