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
 *                    = 1 - array of Stokes parameter Q devided by energy
 *                    = 2 - array of Stokes parameter U devided by energy
 *                    = 3 - array of Stokes parameter V devided by energy
 *                    = 4 - array of degree of polarization
 *                    = 5 - array of polarization angle psi=0.5*atan(U/Q)
 *                    = 6 - array of "Stokes" angle
 *                          beta=0.5*asin(V/sqrt(Q*Q+U*U+V*V))
 * par21 ... nthreads - number of threads to be used for computations
 *
 * NOTES:
 *  -> accuracy vs. speed trade off depends mainly on: nrad, nphi
 * 
 *  -> in this model it is assumed that local emission is completely
 *     linearly polarized in the direction perpendicular to the disc
 * 
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
#define NPARAM 21
#define IFL    1

int main() {

void KYNBB(const double *ear, int ne, const double *param, int ifl, 
           double *photar, double *photer, const char* init);

double ear[NE + 1], param[NPARAM], photar[NE], photer[NE];
char   initstr[0] = "";
int    ie;

param[ 0] = 0.9982;     // a/M
param[ 1] = 30.;        // theta_o
param[ 2] = 1.;         // rin
param[ 3] = 1.;         // ms
param[ 4] = 400.;       // rout
param[ 5] = 0.;         // phi
param[ 6] = 360.;       // dphi
param[ 7] = 3.;         // BHmass
param[ 8] = 1e-9;       // arate
param[ 9] = 1.;         // f_col
param[10] = -3.;        // alpha
param[11] = 0.;         // beta
param[12] = 0.;         // rcloud
param[13] = 0.;         // zshift
param[14] = 80.;        // ntable
param[15] = 200.;       // nrad
param[16] = 1.;         // division
param[17] = 180.;       // nphi
param[18] = 0.;         // smooth
param[19] = 0.;         // Stokes
param[20] = 2.;         // nthreads

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

static double *ener_loc, *flx;
static double am, x0 , x1, x2, x3;
static double arate, BHmass, f_col, rout, Tout, Tnorm, theta_o;
static int    polar;

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

FILE *fw;
double ide_param[25], flux_out[ne + 1];
double far[ne], qar[ne], uar[ne], var[ne], pd[ne], pa[ne], pa2[ne];
double pamin, pamax, pa2min, pa2max;
double am2, pom, pom1, pom2, pom3, rms, r_plus, x, Ccal, Lcal, arcosa3;
int    ne_loc, stokes, ie;

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
ide_param[1] = param[1];
theta_o = param[1];
// rin - inner edge of non-zero disc emissivity
ide_param[2] = param[2];
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
if ((stokes < 0) || (stokes > 6)) {
  xs_write("kynbb: Stokes has to be 0-6", 5);
  for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
  return;
}  
polar = 0;
if (stokes > 0) polar = 1;
ide_param[17] = polar;
// delay_r and delay_phi are not used
// (ide_param[18], ide_param[19])
// number of threads for multithread computations
ide_param[20] = param[20];
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

if (ide(ear, ne, 1, far, qar, uar, var, ide_param, emis_BB, ne_loc)) {
  for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
  return;
}
outer_disc(ear, ne, flux_out);

// interface with XSPEC
if (!stokes) for (ie = 0; ie < ne; ie++) photar[ie] = far[ie] + flux_out[ie];
else {
// final spectrum output -- write ear[] and photar[] into file:
  for (ie = 0; ie < ne; ie++){
    far[ie] = far[ie] + flux_out[ie];
    qar[ie] = qar[ie] + flux_out[ie];
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
fw = fopen("kynbb_photar.dat", "w");
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

void emis_BB(double** ear_loc, const int ne_loc, const int nt, 
             double *far_loc, double *qar_loc, double *uar_loc, 
             double *var_loc, const double r, const double phi, 
             const double cosmu, const double phiphoton, 
             const double alpha_o, const double beta_o, 
             const double delay, const double g) {

double temp, x, Ccal, Lcal, flx0, flx1, g2;
int    ie;

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
for (ie = 1; ie <= ne_loc; ie++) {
  flx1 = flx[ie] / ( exp( *(*ear_loc + ie) / temp ) - 1. );
  far_loc[ie-1] = ( flx0 + flx1 ) / 2. * 
                  ( *(*ear_loc + ie) - *(*ear_loc + ie - 1) ) / g2 ;
  flx0 = flx1;
  if (polar) {
    qar_loc[ie] = far_loc[ie];
    uar_loc[ie] = 0.;
    var_loc[ie] = 0.;
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
