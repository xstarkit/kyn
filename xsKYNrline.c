/* KYNrline - relativistic line with radial broken power-law emissivity
 *            - non-axisymmetric version
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
 * This subroutine takes local Gaussian line emission and gives total spectrum
 * of an accretion disc around a black hole. All relativistic effects are taken
 * into account. It calls subroutine ide() for integrating local emission over
 * the disc and uses the fits file 'KBHtablesNN.fits' defining the transfer
 * functions needed for integration. For details on ide() and the fits file see
 * the subroutine ide() in xside.c.
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
 * par8  ... Erest - rest energy of the line (keV)
 * par9  ... sigma - width of the line - Gaussian sigma (eV)
 * par10 ... q_out - power-law index for radial dependence of emissivity for
 *                   outer region, scales as r^(-q_out)
 * par11 ... q_in  - power-law index for radial dependence of emissivity for
 *                   inner region, scales as rb^(q_in-q_out)*r^(-q_in)
 * par12 ... rb    - boundary between the region with power-law index q_out and
 *                   q_in
 *                 - if > 0 then the boundary is in units of MSO, i.e.
 *                   boundary = rb * r_mso
 *                 - if <= 0 then the boundary is equal to -rb+r_horizon where 
 *                   rb is in GM/c^2
 * par13 ... jump  - ratio of local flux in inner region to local flux in outer
 *                   region at boundary radius defined by rb
 * par14 ... limb  - limb darkening/brightening law (emission directionality)
 *                 - if =  0 the local emisivity is not multiplied by anything
 *                 - if = -1 the local emisivity is multiplied by 1+2.06*mu
 *                   (limb darkening)
 *                 - if = -2 the local emisivity is multiplied by ln(1+1/mu)
 *                   (limb brightening)
 *                 - if different from 0, -1 and -2 then the local emisivity
 *                   is multiplied by mu^(limb)
 * par15 ... alpha  - position of the cloud centre in GM/c^2 in alpha coordinate
 *                    (alpha being the impact parameter in phi direction, 
 *                     positive for approaching side of the disc)
 * par16 ... beta   - position of the cloud centre in GM/c^2 in beta coordinate
 *                    (beta being the impact parameter in theta direction, 
 *                     positive in up direction, i.e. above the disc)
 * par17 ... rcloud - radius of the obscuring cloud (in GM/c^2)
 *                  - if negative, only the emission transmitted through
 *                    the cloud is taken into account
 * par18 ... zshift - overall Doppler shift
 * par19 ... ntable - table of relativistic transfer functions used in the model
 *                    (defines fits file with tables), 0<= ntable <= 99
 * par20 ... nrad   - number of grid points in radius
 * par21 ... division - type of division in r integration
 *                      0 -> equidistant radial grid (constant linear step)
 *                      1 -> exponential radial grid (constant logarithmic step)
 * par22 ... nphi   - number of grid points in azimuth
 * par23 ... smooth - whether to smooth the resulting spectrum (0-no, 1-yes)
 * par24 ... Stokes - what should be stored in photar() array, i.e. as output
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
 * par25 ... nthreads - number of threads to be used for computations
 * par26 ... normtype - how to normalize the spectra
 *                      = 0: normalization to the total photon flux
 *                      > 0: normalization to the photon flux at 'par26' keV
 *                      = -1: the photon flux is not re-normalized,
 *                      = -2: normalization to the maximum of the photon flux
 * 
 *  NOTES:
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

#define IFL    1
#define NPARAM 26
#define NE     400
#define E_MIN  0.1
#define E_MAX  10.

int main() {
  
void KYNrline(const double *ear, int ne, const double *param, int ifl, 
              double *photar, double *photer, const char* init);

double ear[NE + 1], param[NPARAM], photar[NE], photer[NE];
char   initstr[0] = "";
int    ie;

param[ 0] = 1.;   // a/M
param[ 1] = 30.;  // theta_o
param[ 2] = 1.;   // rin
param[ 3] = 1.;   // ms
param[ 4] = 400.; // rout
param[ 5] = 0.;   // phi
param[ 6] = 360.; // dphi
param[ 7] = 6.4;  // Erest
param[ 8] = 2.;   // sigma
param[ 9] = 3.;   // q_out
param[10] = 4.;   // q_in
param[11] = 0.;   // rb
param[12] = 1.;   // jump
param[13] = 0.;   // limb
param[14] = 100.; // alpha
param[15] = 0.;   // beta
param[16] = 0.;   // rcloud
param[17] = 0.;   // zshift
param[18] = 80.;  // ntable
param[19] = 500.; // nrad
param[20] = 1.;   // division
param[21] = 720.; // nphi
param[22] = 1.;   // smooth
param[23] = 0.;   // Stokes
param[24] = 2.;   // nthreads
param[25] = 0.;   // normtype

for (ie = 0; ie <= NE; ie++) {
//  ear[ie] = E_MIN + ie * (E_MAX-E_MIN) / NE;
  ear[ie] = E_MIN * pow(E_MAX/E_MIN, ((double) ie) / NE);
}

KYNrline(ear, NE, param, IFL, photar, photer, initstr);
return(0);
}
#endif
/*******************************************************************************
*******************************************************************************/

#define SQRPI  1.77245385091
#define PI     3.14159265358979
#define NE_LOC 9

/* Let's declare variables that are common for the main and emissivity 
   subroutines */
static double ener_loc[NE_LOC + 1], flx[NE_LOC];
static double qout, qin, rb, fjump, cosin;
static int    polar;

extern int xs_write(char* wrtstr, int idest);

void KYNrline(const double *ear, int ne, const double *param, int ifl,
              double *photar, double *photer, const char* init) {

extern int ide(const double *ear, const int ne, const int nt, double *far, 
               double *qar, double *uar, double *var, 
               const double *ide_param, void (*emissivity)(), 
               const int ne_loc);

void emis_plaw_line(double **ear_loc, const int ne_loc, const int nt, 
                    double *far_loc, double *qar_loc, double *uar_loc, 
                    double *var_loc, const double r, const double phi, 
                    const double cosmu, const double phiphoton, 
                    const double alpha_o, const double beta_o, 
                    const double delay, const double g);
 
FILE *fw;
double ide_param[25];
double far[ne], qar[ne], uar[ne], var[ne], pd[ne], pa[ne], pa2[ne];
double am, am2, r_plus, rms, pom, pom1, pom2, pom3;
double dener, erest, ssigma;
double pamin, pamax, pa2min, pa2max;
int    stokes, ie;

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
// phi - lower azimuth of non-zero disc emissivity (deg)
ide_param[5] = param[5];
// dphi - (phi+dphi) is upper azimuth of non-zero disc emissivity (deg)
ide_param[6] = param[6];
// parameters of the local Gaussian line:
// Erest and Gaussian ssigma=\sqrt{2}\sigma [keV]:
erest = param[7];
if (erest <= 0.) {
  xs_write("kynrline: Erest has to be larger than 0.", 5);
  for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
  return;
}
ssigma = sqrt(2.) * param[8] / 1000.;
if (ssigma <= 0.) {
  xs_write("kynrline: sigma has to be larger than 0.", 5);
  for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
  return;
}
// dener - width of the interval of local energies
dener = 5.2 * ssigma / ((double) NE_LOC - 1);
// q_out, q_in, rb and jump
qout = param[9];
qin = param[10];
rb = param[11];
if (rb > 0.) rb *= rms;
else rb = -rb + r_plus;
fjump = param[12];
if (fjump < 0.) fjump = 1.;
// limb darkening/brightening law
cosin = param[13];
// nrad - number of grid points in radius
ide_param[7] = param[19];
// division - type of division in r integration (0-equidistant, 1-exponential)
ide_param[8] = param[20];
// nphi - number of grid points in azimuth
ide_param[9] = param[21];
// smooth - whether to smooth the resulting spectrum (0-no, 1-yes)
ide_param[10] = param[22];
// normal - how to normalize the final spectrum
ide_param[11] = param[25];
// zshift - overall Doppler shift
ide_param[12] = param[17];
// ntable - table model (defines fits file with tables)
ide_param[13] = param[18];
// edivision - type of division in local energies (0-equidistant, 1-exponential)
ide_param[14] = 0.;
// periodic and dt are not needed for nt = 1
// (ide_param[15], ide_param[16])
// polar - whether we need value of change in polarization angle (0-no, 1-yes)
stokes = (int) param[23];
if ((stokes < 0) || (stokes > 7)) {
  xs_write("kynrline: Stokes has to be 0-7", 5);
  for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
  return;
}
polar = 0;
if (stokes > 0) polar = 1;
ide_param[17] = polar;
// delay_r and delay_phi are not used
// (ide_param[18], ide_param[10])
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
#ifdef OUTSIDE_XSPEC
// Let's write input parameters to a text file
fw = fopen("kynrline.txt", "w");
fprintf(fw, "a/M       %12.6f\n", param[0]);
fprintf(fw, "theta_o   %12.6f\n", param[1]);
fprintf(fw, "rin       %12.6f\n", param[2]);
fprintf(fw, "ms        %12d\n", (int) param[3]);
fprintf(fw, "rout      %12.6f\n", param[4]);
fprintf(fw, "phi       %12.6f\n", param[5]);
fprintf(fw, "dphi      %12.6f\n", param[6]);
fprintf(fw, "Erest     %12.6f\n", param[7]);
fprintf(fw, "sigma     %12.6f\n", param[8]);
fprintf(fw, "q_out     %12.6f\n", param[9]);
fprintf(fw, "q_in      %12.6f\n", param[10]);
fprintf(fw, "rb        %12.6f\n", param[11]);
fprintf(fw, "jump      %12.6f\n", param[12]);
fprintf(fw, "limb      %12.6f\n", param[13]);
fprintf(fw, "alpha     %12.6f\n", ide_param[21]);
fprintf(fw, "beta      %12.6f\n", ide_param[22]);
fprintf(fw, "rcloud    %12.6f\n", ide_param[23]);
fprintf(fw, "zshift    %12.6f\n", param[17]);
fprintf(fw, "ntable    %12d\n", (int) param[18]);
fprintf(fw, "nrad      %12d\n", (int) param[19]);
fprintf(fw, "division  %12d\n", (int) param[20]);
fprintf(fw, "nphi      %12d\n", (int) param[21]);
fprintf(fw, "smooth    %12d\n", (int) param[22]);
fprintf(fw, "Stokes    %12d\n", (int) param[23]);
fprintf(fw, "polar     %12d\n", polar);
fprintf(fw, "r_horizon %12.6f\n", r_plus);
fprintf(fw, "r_ms      %12.6f\n", rms);
fprintf(fw, "edivision %12d\n", (int) ide_param[14]);
fprintf(fw, "ne_loc    %12d\n", NE_LOC);
fprintf(fw, "dener     %12.6f\n", dener);
fprintf(fw, "ssigma    %12.6f\n", ssigma);
fprintf(fw, "normal    %12.6f\n", ide_param[11]);
fprintf(fw, "nthreads  %12d\n", (int) ide_param[20]);
fclose(fw);
#endif
/******************************************************************************/
      
// local flux defined in local energies
for (ie = 0; ie < NE_LOC; ie++) {
  ener_loc[ie] = erest - dener * (NE_LOC - 1) / 2.0 + dener * ie;
  flx[ie] = exp(-pow((ener_loc[ie] - erest) / ssigma, 2.)) / ssigma / SQRPI;
}

/*******************************************************************************
// local spectrum output -- write ener_loc[] and flx[] into file:
fw = fopen("kynrline_photar_loc.dat", "w");
for (ie = 0; ie < NE_LOC; ie++) {
  fprintf(fw, "%14.6f\t%E\n", ener_loc[ie], flx[ie]);
}
fclose(fw);
*******************************************************************************/
      
// Let's integrate local emission over the accretion disc
if (ide(ear, ne, 1, far, qar, uar, var, ide_param, emis_plaw_line, NE_LOC)) {
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
fw = fopen("kynrline_photar.dat", "w");
for (ie = 0; ie < ne; ie++) fprintf(fw, "%14.6f\t%E\n", 
                                    0.5 * (ear[ie] + ear[ie + 1]),
                                    photar[ie] / (ear[ie + 1] - ear[ie]));
fclose(fw);
#endif
/******************************************************************************/

return;
}
/*******************************************************************************
*******************************************************************************/

void emis_plaw_line(double** ear_loc, const int ne_loc, const int nt, 
                    double *far_loc, double *qar_loc, double *uar_loc, 
                    double *var_loc, const double r, const double phi, 
                    const double cosmu, const double phiphoton, 
                    const double alpha_o, const double beta_o, 
                    const double delay, const double g) {
  
double rq;
int ie;

*ear_loc = ener_loc;
if (r >= rb) rq = 1. / pow(r, qout);
else rq = fjump * pow(rb, qin - qout) / pow(r, qin);
if ((cosin != 0.) && (cosin != -1.) && (cosin != -2.)) 
  rq = rq * pow(cosmu, cosin);
if (cosin == -1.) rq = rq * (1. + 2.06 * cosmu);
if (cosin == -2.) rq = rq * log(1. + 1. / cosmu);
for (ie = 0; ie < ne_loc; ie++) {
  far_loc[ie] = flx[ie] * rq;
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
