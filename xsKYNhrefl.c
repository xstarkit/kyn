/* KYNhrefl - general relativistic extension of the XSPEC model hrefl(powerlaw)
 *            model subroutine for XSPEC
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
 * Compton reflection model for an accretion disc around a Kerr black hole.
 * This model is based on an existing multiplicative HREFL model in combination
 * with the POWERLAW model. Local emission is the same as in HREFL*POWERLAW with
 * the parameters thetamin = 0 and thetamax = 90 (i.e. it is assumed that the
 * disc is illuminated from all directions isotropically) and with a broken
 * power-law radial dependence added. This model can be interpreted as a
 * Compton-reflection model for which the source of primary irradiation is near
 * above the disc, in contrast to the lamp-post scheme with the source on the
 * axis. The approximations for Compton reflection used in HREFL (and therefore
 * also in this model) are valid below ~15keV (in the disc rest frame).
 *
 * All relativistic effects are taken into account. This subroutine calls
 * subroutine ide() for integrating local emission over the disc and uses the
 * FITS file 'KBHtablesNN.fits' defining the transfer functions needed for the
 * integration. For details on ide() and the FITS file see the subroutine ide()
 * in xside.c.
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
 * par8  ... PhoIndex  - photon index of primary power-law illumination
 * par9  ... q_out - power-law index for radial dependence of emissivity for
 *                   outer region, scales as r^(-q_out)
 * par10 ... q_in  - power-law index for radial dependence of emissivity for
 *                   inner region, scales as rb^(q_in-q_out)*r^(-q_in)
 * par11 ... rb    - boundary between the region with power-law index q_out and
 *                   q_in
 *                 - if > 0 then the boundary is in units of MSO, i.e.
 *                   boundary = rb * r_mso
 *                 - if <= 0 then the boundary is equal to -(rb-r_horizon) where 
 *                   rb is in GM/c^2
 * par12 ... jump  - ratio of local flux in inner region to local flux in outer
 *                   region at boundary radius defined by rb
 * par13 ... Feabun  - iron abundance relative to Solar
 * par14 ... FeKedge - iron K-edge energy
 * par15 ... Escfrac - normalization of the original powerlaw emission
 * par16 ... covfac - normalization of the reflected emission
 * par17 ... alpha  - position of the cloud centre in GM/c^2 in alpha coordinate
 *                    (alpha being the impact parameter in phi direction, 
 *                     positive for approaching side of the disc)
 * par18 ... beta   - position of the cloud centre in GM/c^2 in beta coordinate
 *                    (beta being the impact parameter in theta direction, 
 *                     positive in up direction, i.e. above the disc)
 * par19 ... rcloud - radius of the obscuring cloud (in GM/c^2)
 *                  - if negative, only the emission transmitted through
 *                    the cloud is taken into account
 * par20 ... zshift - overall Doppler shift
 * par21 ... ntable - table of relativistic transfer functions used in the model
 *                    (defines fits file with tables), 0<= ntable <= 99
 * par22 ... nrad   - number of grid points in radius
 * par23 ... division - type of division in r integration
 *                      0 -> equidistant radial grid (constant linear step)
 *                      1 -> exponential radial grid (constant logarithmic step)
 * par24 ... nphi   - number of grid points in azimuth
 * par25 ... smooth - whether to smooth the resulting spectrum (0-no, 1-yes)
 * par26 ... Stokes - what should be stored in photar() array, i.e. as output
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
 * par27 ... nthreads - number of threads to be used for computations
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
#define E_MAX  15.
#define NPARAM 27
#define IFL    1

int main() {

void KYNhrefl(const double *ear, int ne, const double *param, int ifl, 
              double *photar, double *photer, const char* init);

double ear[NE + 1], param[NPARAM], photar[NE], photer[NE];
char   initstr[0] = "";
int    ie;

param[ 0] = 0.9982;     // a/M
param[ 1] = 30.;        // theta_o
param[ 2] = 1.;         // rin
param[ 3] = 0.;         // ms
param[ 4] = 10.;        // rout
param[ 5] = 0.;         // phi
param[ 6] = 360.;       // dphi
param[ 7] = 1.7;        // PhoIndex
param[ 8] = 3;          // q_out
param[ 9] = 2.;         // q_in
param[10] = 0.;         // rb
param[11] = 0.;         // jump
param[12] = 1.;         // Feabun
param[13] = 7.11;       // FeKedge
param[14] = 0.;         // Escfrac
param[15] = 1.;         // covfac
param[16] = -3.;        // alpha
param[17] = 0.;         // beta
param[18] = 0.;         // rcloud
param[19] = 0.;         // zshift
param[20] = 80.;        // ntable
param[21] = 500.;       // nrad
param[22] = 1.;         // division
param[23] = 300.;       // nphi
param[24] = 1.;         // smooth
param[25] = 0.;         // Stokes
param[26] = 2.;         // nthreads

for (ie = 0; ie <= NE; ie++) {
//  ear[ie] = E_MIN + ie * (E_MAX-E_MIN) / NE;
  ear[ie] = E_MIN * pow(E_MAX/E_MIN, ((double) ie) / NE);
}

KYNhrefl(ear, NE, param, IFL, photar, photer, initstr);
return(0);
}

#endif
/*******************************************************************************
*******************************************************************************/

#define PI 3.14159265358979

static double *ear_local, *far_local;
static double qout, qin, rb, fjump, afe, ekedge, escfrac, covfac;
static int    polar;

extern int xs_write(char* wrtstr, int idest);

void KYNhrefl(const double *ear, int ne, const double *param, int ifl, 
              double *photar, double *photer, const char* init) {

extern int ide(const double *ear, const int ne, const int nt, double *far, 
               double *qar, double *uar, double *var, 
               const double *ide_param, void (*emissivity)(), 
               const int ne_loc);

void emis_hrfl(double** ear_loc, const int ne_loc, const int nt, 
               double *far_loc, double *qar_loc, double *uar_loc, 
               double *var_loc, const double r, const double phi, 
               const double cosmu, const double phiphoton, 
               const double alpha_o, const double beta_o, 
               const double delay, const double g);

static int ne_loc_old = 1;

FILE *fw;
double ide_param[25];
double far[ne], qar[ne], uar[ne], var[ne], pd[ne], pa[ne], pa2[ne];
double emin, emax, dener, phoindex, pamin, pamax, pa2min, pa2max;
double am, am2, r_plus, rms, pom, pom1, pom2, pom3;
int    ie, ne_loc, stokes;

// Let's initialize parameters for subroutine ide()
// a/M - black hole angular momentum
ide_param[0] = param[0];
am = param[0];
am2 = am * am;
pom1 = pow(1. + am, 1. / 3.);
pom2 = pow(1. - am, 1. / 3.);
pom3 = pow(1. - am2, 1. / 3.);
pom = 1. + pom3 * (pom1 + pom2);
pom1 = sqrt(3. * am2 + pom * pom);
if (am >= 0) rms= 3. + pom1 - sqrt((3. - pom) * (3. + pom + 2. * pom1));
else rms = 3. + pom1 + sqrt((3. - pom) * (3. + pom + 2. * pom1));
r_plus= 1. + sqrt(1. - am2);
// theta_o - observer inclination
ide_param[1] = param[1];
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
ide_param[7] = param[21];
// division - type of division in r integration (0-equidistant, 1-exponential)
ide_param[8] = param[22];
// nphi - number of grid points in azimuth
ide_param[9] = param[23];
// smooth - whether to smooth the resulting spectrum (0-no, 1-yes)
ide_param[10] = param[24];
// normal - how to normalize the final spectrum
ide_param[11] = 1.;
// zshift - overall Doppler shift
ide_param[12] = param[19];
// ntable - table model (defines fits file with tables)
ide_param[13] = param[20];
// dener - width of the interval of local energies
emin = 0.5 * ear[0];
emax = 3 * ear[ne];
dener = 0.05;
ne_loc = (int) (round(log10(emax / emin) / dener)) + 1;
dener = log10(emax / emin) / (ne_loc - 1);
// edivision - type of division in local energies (0-equidistant, 1-exponential)
ide_param[14] = 1.;
// PhoIndex - power-law energy index of the incident radiation
phoindex = param[7];
// q_out, q_in, rb and jump
qout = param[8];
qin = param[9];
rb = param[10];
if (rb > 0.) rb *= rms;
else rb = -rb + r_plus;
fjump = param[11];
if (fjump < 0.) fjump = 1.;
// Feabun - iron abundance relative to Solar
afe = param[12];
// FeKedge - iron K-edge energy
ekedge = param[13];
// Escfrac - fraction of the direct flux seen by the observer
escfrac = param[14];
// covfac -  normalization of the reflected continuum
covfac = param[15];
if ((ne_loc_old == -1) || (ne_loc != ne_loc_old)) {
  if ((ne_loc_old != -1) && (ne_loc != ne_loc_old)) {
    free((void *) ear_local);
    ear_local = NULL;
    free((void *) far_local);
    far_local = NULL;
  }
  if ((ear_local = (double *) malloc((ne_loc + 1) * sizeof(double))) == NULL) {
    xs_write("kynhrefl: Failed to allocate memory for tmp arrays.", 5);
    for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
    return;
  }
  if ((far_local = (double *) malloc((ne_loc) * sizeof(double))) == NULL) {
    xs_write("kynhrefl: Failed to allocate memory for tmp arrays.", 5);
    for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
    return;
  }
}
for (ie = 0; ie < ne_loc; ie++) {
  ear_local[ie] = emin * pow(emax / emin, ((double) ie) / (ne_loc - 1));
  far_local[ie] = pow(ear_local[ie], -phoindex);
}
// periodic and dt are not needed for nt = 1
// (ide_param[15], ide_param[16])
// polar - whether we need value of change in polarization angle (0-no,1-yes)
stokes = (int) param[25];
if ((stokes < 0) || (stokes > 7)) {
  xs_write("kynhrefl: Stokes has to be 0-7", 5);
  for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
  return;
}
polar = 0;
if (stokes > 0) polar = 1;
ide_param[17] = polar;
// delay_r and delay_phi are not used
// (ide_param[18], ide_param[19])
// number of threads for multithread computations
ide_param[20] = param[26];
// alpha - position of the cloud in alpha impact parameter (in GM/c^2)
ide_param[21] = param[16];
// beta - position of the cloud in beta impact parameter (in GM/c^2)
ide_param[22] = param[17];
// rcloud - radius of the cloud (in GM/c^2)
ide_param[23] = param[18];
//whether the flux defined in emissivity subroutine is local one (0) or the 
//observed one (1)
ide_param[24] = 0.;

/******************************************************************************/
#ifdef OUTSIDE_XSPEC
// Let's write input parameters to a text file
fw = fopen("kynhrefl.txt", "w");
fprintf(fw, "a/M          %12.6f\n", param[0]);
fprintf(fw, "theta_o      %12.6f\n", param[1]);
fprintf(fw, "rin          %12.6f\n", param[2]);
fprintf(fw, "ms           %12d\n", (int) param[3]);
fprintf(fw, "rout         %12.6f\n", param[4]);
fprintf(fw, "phi          %12.6f\n", param[5]);
fprintf(fw, "dphi         %12.6f\n", param[6]);
fprintf(fw, "PhoIndex     %12.6f\n", param[7]);
fprintf(fw, "q_out        %12.6f\n", param[8]);
fprintf(fw, "q_in         %12.6f\n", param[9]);
fprintf(fw, "rb           %12.6f\n", param[10]);
fprintf(fw, "jump         %12.6f\n", param[11]);
fprintf(fw, "Feabun       %12.6f\n", param[12]);
fprintf(fw, "FeKedge      %12.6f\n", param[13]);
fprintf(fw, "Escfrac      %12.6f\n", param[14]);
fprintf(fw, "covfac       %12.6f\n", param[15]);
fprintf(fw, "alpha        %12.6f\n", ide_param[21]);
fprintf(fw, "beta         %12.6f\n", ide_param[22]);
fprintf(fw, "rcloud       %12.6f\n", ide_param[23]);
fprintf(fw, "zshift       %12.6f\n", param[19]);
fprintf(fw, "ntable       %12d\n", (int) param[20]);
fprintf(fw, "nrad         %12d\n", (int) param[21]);
fprintf(fw, "division     %12d\n", (int) param[22]);
fprintf(fw, "nphi         %12d\n", (int) param[23]);
fprintf(fw, "smooth       %12d\n", (int) param[24]);
fprintf(fw, "Stokes       %12d\n", (int) param[25]);
fprintf(fw, "polar        %12d\n", polar);
fprintf(fw, "r_horizon    %12.6f\n", r_plus);
fprintf(fw, "r_ms         %12.6f\n", rms);
fprintf(fw, "e_min        %12.6f\n", emin);
fprintf(fw, "e_max        %12.6f\n", emax);
fprintf(fw, "edivision    %12d\n", (int) ide_param[14]);
fprintf(fw, "ne_loc       %12d\n", ne_loc);
fprintf(fw, "dener        %12.6f\n", dener);
fprintf(fw, "normal       %12.6f\n", ide_param[11]);
fprintf(fw, "nthreads     %12d\n", (int) ide_param[20]);
fclose(fw);
#endif
/******************************************************************************/

/******************************************************************************
// local spectrum output -- write ear_local[] and far_local[] into file:
fw = fopen("kynhrefl_photar_loc.dat", "w");
for (ie = 0; ie < ne_loc; ie++)
  fprintf(fw, "%14.6f\t%E\n", ear_local[ie], far_local[ie]);
fclose(fw);
******************************************************************************/

if (ide(ear, ne, 1, far, qar, uar, var, ide_param, emis_hrfl, ne_loc)) {
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
// final spectrum output -- write ear[] and photar[] into file:
fw = fopen("kynhrefl_photar.dat", "w");
for (ie = 0; ie < ne; ie++)
  fprintf(fw, "%14.6f\t%E\n", 0.5 * (ear[ie] + ear[ie+1]), 
          photar[ie] / (ear[ie+1] - ear[ie]));
fclose(fw);
#endif
/******************************************************************************/

return;
}

/*******************************************************************************
*******************************************************************************/

void emis_hrfl(double** ear_loc, const int ne_loc, const int nt, 
               double *far_loc, double *qar_loc, double *uar_loc, 
               double *var_loc, const double r, const double phi, 
               const double cosmu, const double phiphoton, 
               const double alpha_o, const double beta_o, 
               const double delay, const double g) {

double KYhrefl(double xmu, double e, double afe, double ekedge);

double rq;
int    ie;

*ear_loc = ear_local;
if (r >= rb) rq = 1. / pow(r, qout);
else rq = fjump * pow(rb, qin - qout) / pow(r, qin);
for (ie = 0; ie < ne_loc; ie++) {
  far_loc[ie] = far_local[ie] * rq *
             (escfrac + covfac * KYhrefl(cosmu, *(*ear_loc + ie), afe, ekedge));
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

// Here functions from the HREFL model are defined

double KYhrefl(double xmu, double e, double afe, double ekedge) {

/* simple reflection model (see below) - valid for E<15 keV in the rest frame
 * T. Yaqoob sometime between 1991-1993; cleaned up a lot by C. Day.
 *
 *    Function computes the reflected spectrum from a cold semi-infinite
 *    slab in the approximation of elastic electron scattering, i.e.
 *    incident photon energy << mc^2 in the electron rest-frame. Should
 *    be good for photon energies up to ~15 keV.
 *
 *    USEAGE: Suppose direct spectrum is: spec (photons/cm/cm/s/keV)
 *    then the TOTAL observed spectrum is spec*(1+KYHREFL) in the same
 *    units.
 *
 *    PARAMETERS:
 *       XMU ..... cos(Angle) Angle is between the observer's line of sight
 *                     and the slab normal (in local fluid's frame).
 *       E ........ Energy of incident photon.
 *       AFE ...... Iron abundance relative to Solar.
 *       EKEDGE ... Iron K-edge energy.
 *    Recommended that D1,D2,FEBUN and EKEDGE be fixed.
 *
 *    The routine uses the solutions to the transfer equations utilizing
 *    the H-functions [see Basko, 1978]. Analytic approxiamtions to the
 *    H-functions are used, as well as an average of the H-functions for
 *    the incident rays. Errors from these approximations are typically
 *    less than a few percent, snf the results agree well with Monte-
 *    Carlo simulations.
 */

double albd(double afe, double ekedge, double e);

double hmean(double albedo);

double hfuna(double albedo, double xmu);

double albedo, f2;

if (xmu < -1.) xmu = -1.;
if (xmu > 1.) xmu = 1.;
albedo = albd(afe, ekedge, e);
if (xmu == 0.) f2 = 1.;
else f2 = log((1. + xmu) / xmu);
return(albedo * hmean(albedo) * hfuna(albedo, xmu) * f2 * 0.5);
// return(albedo * hmean(albedo) * hfuna(albedo, xmu) * xmu * 0.5);
}

//******************************************************************************

double albd(double afe, double ekedge, double e) {
// Computes the ratio of the Thompson Cross-section to the total
// (scattering + absorption) cross-section.

double abscrs(double e, double cd, double afe, double ekedge);

double thom = 0.665e-3;

return(thom / (thom + (abscrs(e, 1., afe, ekedge) / 1.2)));
}

//******************************************************************************

double hfuna(double albedo, double xmu) {

return((1. + sqrt(3.) * xmu) / (1. + sqrt(3. * (1. - albedo)) * xmu));
}

//******************************************************************************

double hmean(double albedo) {

double a, b, r;

a = sqrt(3.);
if (albedo == 1.) return(1. + a / 2.);
else {
  b = sqrt(3. * (1. - albedo));
  r = a / b;
  return(1. + (r - 1.) * (1. - (log(1. + b) / b)));
}
}

//******************************************************************************

double abscrs(double e, double cd, double afe, double ekedge) {

/*     Computes photoelectric absorption.
 *
 *     Uses piecewise polynomial fit of Morrison and McCammon Ap.J. 270,
 *     119 for range 0.03 to 10 keV. Below 0.03 keV uses power law fit to
 *     hydrogen and helium edge profiles interpolated/extrapolated from
 *     Henke data (1982) also used by Morrisom and McMammon. Above 10 keV
 *     crude "eyeball" fit provided by Gordon Stewart!? Corrected error
 *     in low energy stuff (CD missing!) RW 1988-May-27. Author Dick
 *     Willingale 1986-Sep-4
 *
 *     PARAMETERS
 *        E ........ Energy (keV)
 *        CD ....... Column density 10^21 cm-2
 *        EKEDGE ... Energy (keV)
 */

double feabs(double e, double cd, double afe, double ekedge);

double c0[14] = {17.3, 34.6, 78.1, 71.4, 95.5, 308.9, 120.6, 141.3, 202.7,
                 342.7, 352.2, 433.9, 629.0, 701.2};
double c1[14] = {608.1, 267.9, 18.8, 66.8, 145.8, -380.6, 169.3, 146.8, 104.7,
                 18.7, 18.7, -2.4, 30.9, 25.2};
double c2[14] = {-2150., -476.1, 4.3, -51.4, -61.1, 294., -47.7, -31.5, -17.,
                 0., 0., 0.75, 0., 0.};
double et[15] = {.03, .1, .284, .4, .532, .707, .867, 1.303, 1.84, 2.471, 3.21,
                 4.038, 7.111, 8.331, 10.};
double tau, aeff;
int    ke;

if (e < 0.0136) tau = 0.;
else if ((e >= 0.0136) && (e < 0.0246)) tau = cd * 9.49e-3 / pow(e, 3.26);
else if ((e >= 0.0246) && (e < 0.03)) tau = cd * 8.63e-3 / pow(e, 3.26) + 
                                            289.e-3 / pow(e, 2.11);
else if (e >= 10.) tau = (0.3 * cd) / (e * e * sqrt(e));
else for (ke = 0; ke < 14; ke++) {
       if ((e >= et[ke]) && (e < et[ke + 1])) {
         tau = (c0[ke] + c1[ke] * e + c2[ke] * e * e) / e / e / e;
         tau = 0.001 * tau * cd;
       }
     }
aeff = 1. - afe;
return(tau - feabs(e, cd, aeff, ekedge));
}

//******************************************************************************

double feabs(double e, double cd, double afe, double ekedge) {

  /*     Computes iron absorption
 *     PARAMETERS
 *       EKEDGE ... Iron K-edge energy (keV)
 *       E ........ Energy (keV)
 *       AFE ...... Iron abundance relative to Solar
 *       CD ....... Column density in units of 10^16.52 cm-2 so that
 *                     values are equivalent to hydrogen column in units
 *                     10^21 cm-2, assuming abundance ratio Fe/H of
 *                     10^(7.52-12).
 */

if ((e >= 0.73) && (e < ekedge)) return(0.043 * afe * cd * pow(0.73 / e, 2.2));
else if (e >= ekedge) return(0.0012 * afe * cd * pow(ekedge / e, 3.11));
else return(0.);
}
