/* KYNconv - convolution general relativistic model - non-axisymmetric version
 *          model subroutine for XSPEC
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
 * This convolution subroutine takes input photar[] array as a definition of the
 * local flux across the accretion disc around a black hole and adds a broken
 * power-law radial dependence and limb darkening/brightening law to it. The
 * output is total spectrum of an accretion disc. All relativistic effects are
 * taken into account. This model calls subroutine ide() for integrating local
 * emission over the disc and uses the fits file 'KBHtablesNN.fits' defining the
 * transfer functions needed for integration. For details on ide() and the fits
 * file see the subroutine ide() in xside.c.
 * 
 * There are several restrictions that arise from the fact that we use existing
 * XSPEC models for definition of the local flux:
 * - only the energy dependence of the photon flux can be defined by local XSPEC
 *   models,
 * - only a certain type of radial dependence of the local photon flux can be
 *   imposed - we have chosen to use a broken power-law radial dependence,
 * - there is no intrinsic azimuthal dependence of the local photon flux, the
 *   only azimuthal dependence comes through limb darkening/brightening law
 *   (emission angle depends on azimuth)
 * - local flux can highly depend on the energy resolution, i.e. on the energy
 *   binning used, if the energy resolution is not high enough. This is because
 *   the flux is defined in the centre of each bin. A large number of bins is
 *   needed for highly varying local flux with energy.
 * 
 * For emissivities that cannot be defined by existing XSPEC models, or where 
 * the limitations mentioned above are too restrictive, one has to add a new
 * user-defined model to XSPEC (by adding a new subroutine to XSPEC). This 
 * method is more flexible and faster than when using this convolution model, 
 * and hence it is recommended even for cases when this model could be used. 
 * In any new model for XSPEC one can use the common ray-tracing driver for 
 * relativistic smearing of the local emission: ide() for non-axisymmetric 
 * models and idre() for axisymmetric ones.
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
 * par8  ... q_out - power-law index for radial dependence of emissivity for
 *                   outer region, scales as r^(-q_out)
 * par9  ... q_in  - power-law index for radial dependence of emissivity for
 *                   inner region, scales as rb^(q_in-q_out)*r^(-q_in)
 * par10 ... rb    - boundary between the region with power-law index q_out and
 *                   q_in
 *                 - if > 0 then the boundary is in units of MSO, i.e.
 *                   boundary = rb * r_mso
 *                 - if <= 0 then the boundary is equal to -(rb-r_horizon) where 
 *                   rb is in GM/c^2
 * par11 ... jump  - ratio of local flux in inner region to local flux in outer
 *                   region at boundary radius defined by rb
 * par12 ... limb  - limb darkening/brightening law (emission directionality)
 *                 - if =  0 the local emisivity is not multiplied by anything
 *                 - if = -1 the local emisivity is multiplied by 1+2.06*mu
 *                   (limb darkening)
 *                 - if = -2 the local emisivity is multiplied by ln(1+1/mu)
 *                   (limb brightening)
 *                 - if different from 0, -1 and -2 then the local emisivity
 *                   is multiplied by mu^(limb)
 * par13 ... alpha  - position of the cloud centre in GM/c^2 in alpha coordinate
 *                    (alpha being the impact parameter in phi direction, 
 *                     positive for approaching side of the disc)
 * par14 ... beta   - position of the cloud centre in GM/c^2 in beta coordinate
 *                    (beta being the impact parameter in theta direction, 
 *                     positive in up direction, i.e. above the disc)
 * par15 ... rcloud - radius of the obscuring cloud (in GM/c^2)
 *                  - if negative, only the emission transmitted through
 *                    the cloud is taken into account
 * par16 ... zshift - overall Doppler shift
 * par17 ... ntable - table of relativistic transfer functions used in the model
 *                    (defines fits file with tables), 0<= ntable <= 99
 * par18 ... nrad   - number of grid points in radius
 * par19 ... division - type of division in r integration
 *                      0 -> equidistant radial grid (constant linear step)
 *                      1 -> exponential radial grid (constant logarithmic step)
 * par20 ... nphi   - number of grid points in azimuth
 * par21 ... ne_loc - number of grid points in local energy
 *                    (energy resolution of local flux, the grid is equidistant
 *                    in logarithmic scale)
 * par22 ... smooth - whether to smooth the resulting spectrum (0-no, 1-yes)
 * par23 ... Stokes - what should be stored in photar() array, i.e. as output
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
 * par24 ... nthreads - number of threads to be used for computations
 * par25 ... normtype - how to normalize the spectra
 *                      =  0: the total photon flux is the same as the total 
 *                            flux before convolution,
 *                      >  0: the photon flux at 'par25' keV is the same as the 
 *                            photon flux at that energy before convolution,
 *                      = -1: the photon flux is not re-normalized,
 *                      = -2: the maximum of the photon flux is the same as the 
 *                            maximum of the photon flux before convolution,
 * 
 * NOTES:
 *  -> accuracy vs. speed trade off depends mainly on: nrad, nphi, ne_loc
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
#define NPARAM 25
#define NE     200
#define E_MIN  0.01
#define E_MAX  15.

int main() {

void KYNconv(const double *ear, int ne, const double *param, int ifl, 
             double *photar, double *photer, const char* init);
  
double ear[NE + 1], param[NPARAM], photar[NE], photer[NE];
char   initstr[0] = "";
int    ie;
  
param[ 0] = 1.;     // a/M
param[ 1] = 30.;    // theta_o
param[ 2] = 1.;     // rin
param[ 3] = 1.;     // ms
param[ 4] = 400.;   // rout
param[ 5] = 0.;     // phi
param[ 6] = 360.;   // dphi
param[ 7] = 3.;     // q_out
param[ 8] = 2.;     // q_in
param[ 9] = 0.;     // rb
param[10] = 1.;     // jump
param[11] = 0.;     // limb
param[12] = -3.;    // alpha
param[13] = 0.;     // beta
param[14] = 0.;     // rcloud
param[15] = 0.;     // zshift
param[16] = 80.;    // ntable
param[17] = 100.;   // nrad
param[18] = 1.;     // division
param[19] = 180.;   // nphi
param[20] = 100.;   // ne_loc
param[21] = 1.;     // smooth
param[22] = 0.;     // Stokes
param[23] = 2.;     // nthreads
param[24] = 0.;     // normtype

for (ie = 0; ie <= NE; ie++) {
//  ear[ie] = E_MIN + ie * (E_MAX - E_MIN) / NE;
  ear[ie] = E_MIN * pow(E_MAX / E_MIN, ((double) ie) / NE);
}

for (ie = 0; ie < NE; ie++)
  photar[ie] = 2. / (ear[ie + 1] + ear[ie]);

KYNconv(ear, NE, param, IFL, photar, photer, initstr);
return(0);
}
#endif
/*******************************************************************************
*******************************************************************************/

#define PI     3.14159265358979

static double *ear_local, *far_local;
static double qout, qin, rb, fjump, cosin;
static int    polar;

extern int xs_write(char* wrtstr, int idest);

void KYNconv(const double *ear, int ne, const double *param, int ifl, 
             double *photar, double *photer, const char* init) {
  
extern int ide(const double *ear, const int ne, const int nt, double *far, 
               double *qar, double *uar, double *var, 
               const double *ide_param, void (*emissivity)(), 
               const int ne_loc);

void emis_KYNconv(double **ear_loc, const int ne_loc, const int nt, 
                  double *far_loc, double *qar_loc, double *uar_loc, 
                  double *var_loc, const double r, const double phi, 
                  const double cosmu, const double phiphoton, 
                  const double alpha_o, const double beta_o, 
                  const double delay, const double g);

static int ne_loc_old = -1;

FILE *fw;
double ide_param[25];
double far[ne], qar[ne], uar[ne], var[ne], pd[ne], pa[ne], pa2[ne];
double am, am2, r_plus, rms, pom, pom1, pom2, pom3;
double pamin, pamax, pa2min, pa2max, normal, norm0, norm1, ttmp, ttmp1;
int    ie, ie0, je, ne_loc, stokes, imax, imin;

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
// nrad - number of grid points in radius
ide_param[7] = param[17];
// division - type of division in r integration (0-equidistant, 1-exponential)
ide_param[8] = param[18];
// nphi - number of grid points in azimuth
ide_param[9] = param[19];
// smooth - whether to smooth the resulting spectrum (0-no, 1-yes)
ide_param[10] = param[21];
// normal - how to normalize the final spectrum
normal = param[24];
// zshift - overall Doppler shift
ide_param[12] = param[15];
// ntable - table model (defines fits file with tables)
ide_param[13] = param[16];
// number of points in local energy (resolution of the local flux)
ne_loc = (int) param[20];
if ((ne_loc_old == -1) || (ne_loc != ne_loc_old)) {
  if ((ne_loc_old != -1) && (ne_loc != ne_loc_old)) {
    free((void *) ear_local);
    ear_local = NULL;
    free((void *) far_local);
    far_local = NULL;
  }
  if ((ear_local = (double *) malloc((ne_loc + 1) * sizeof(double))) == NULL) {
    xs_write("kynconv: Failed to allocate memory for tmp arrays.", 5);
    for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
    return;
  }
  if ((far_local = (double *) malloc(ne_loc * sizeof(double))) == NULL) {
    xs_write("kynconv: Failed to allocate memory for tmp arrays.", 5);
    for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
    return;
  }
  ne_loc_old = ne_loc;
}
for (ie = 0; ie < ne_loc; ie++) {
  ear_local[ie] = (ear[0] + ear[1]) / 2. * pow((ear[ne - 1] + ear[ne]) /
                  (ear[0] + ear[1]), ((double) ie) / (ne_loc - 1.));
}
far_local[0] = photar[0] / (ear[1] - ear[0]);
far_local[ne_loc - 1] = photar[ne - 1] / (ear[ne] - ear[ne - 1]);
for (ie = 1; ie < ne_loc - 1; ie++){
  for (je = 1; je < ne; je++) {
    if ((ear_local[ie] > (ear[je] + ear[je - 1]) / 2.) && 
        (ear_local[ie] <= (ear[je + 1] + ear[je]) / 2.)) {
      far_local[ie] = photar[je - 1] / (ear[je] - ear[je - 1]) + 
      (2. * ear_local[ie] - (ear[je] + ear[je - 1])) /
      (ear[je + 1] - ear[je - 1]) * (photar[je] /
      (ear[je + 1] - ear[je]) - photar[je - 1] / (ear[je] - ear[je - 1]));
      break;
    }
  }
}
norm0=0.;
for (ie = 0; ie < ne; ie++) norm0+=photar[ie];
if(normal>=0. || normal==-2){
  if((normal >0. && (ear_local[0] > normal || ear_local[ne_loc-1] < normal)) ||
      normal==0.){
    ide_param[11]=0.;
  }else if(normal == -2.){
    for (ie = 0; ie < ne_loc; ie++) if(norm0 < far_local[ie])norm0=far_local[ie];
    ide_param[11]=-2.;    
  }else{
// given normal, find the corresponding index in ear_local[]:
    imin = 0;
    imax = ne_loc-1;
    ie0 = ne_loc / 2;
    while ((imax - imin) > 1) {
      if (normal >= ear_local[ie0 - 1]) imin = ie0;
      else imax = ie0;
      ie0 = (imin + imax) / 2;
    }
    if (ie0 == 0) ie0 = 1;
    ttmp = (normal - ear_local[ie0 - 1]) /
           (ear_local[ie0] - ear_local[ie0 - 1]);
    ttmp1 = 1. - ttmp;
    norm1=ttmp1*far_local[ie0-1]+ttmp*far_local[ie0];
    if(norm1 != 0.){
      ide_param[11]=normal;
      norm0=norm1;
    }else ide_param[11]=0.;
  }
}else{
  norm0=1.;
  ide_param[11]=-1.;
}
// edivision - type of division in local energies (0-equidistant, 1-exponential)
ide_param[14] = 1.;
// q_out, q_in, rb and jump
qout = param[7];
qin = param[8];
rb = param[9];
if (rb > 0.) rb *= rms;
else rb = -rb + r_plus;
fjump = param[10];
if (fjump < 0.) fjump = 1.;
// limb darkening/brightening law
cosin = param[11];
// periodic and dt are not needed for nt = 1
// (ide_param[15], ide_param[16])
// polar - whether we need value of change in polarization angle (0-no, 1-yes)
stokes = (int) param[22];
if ((stokes < 0.) || (stokes > 7)) {
  xs_write("kynconv: Stokes has to be 0-7", 5);
  for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
  return;
}
polar = 0;
if (stokes > 0) polar = 1;
ide_param[17] = polar;
// delay_r and delay_phi are not used
// (ide_param[18], ide_param[19])
// number of threads for multithread computations
ide_param[20] = param[23];
// alpha - position of the cloud in alpha impact parameter (in GM/c^2)
ide_param[21] = param[12];
// beta - position of the cloud in beta impact parameter (in GM/c^2)
ide_param[22] = param[13];
// rcloud - radius of the cloud (in GM/c^2)
ide_param[23] = param[14];
//whether the flux defined in emissivity subroutine is local one (0) or the 
//observed one (1)
ide_param[24] = 0.;

/******************************************************************************/
#ifdef OUTSIDE_XSPEC
// Let's write input parameters to a text file
fw = fopen("kynconv.txt", "w");
fprintf(fw, "a/M          %12.6f\n", param[0]);
fprintf(fw, "theta_o      %12.6f\n", param[1]);
fprintf(fw, "rin          %12.6f\n", param[2]);
fprintf(fw, "ms           %12d\n", (int) param[3]);
fprintf(fw, "rout         %12.6f\n", param[4]);
fprintf(fw, "phi          %12.6f\n", param[5]);
fprintf(fw, "dphi         %12.6f\n", param[6]);
fprintf(fw, "q_out        %12.6f\n", param[7]);
fprintf(fw, "q_in         %12.6f\n", param[8]);
fprintf(fw, "rb           %12.6f\n", param[9]);
fprintf(fw, "jump         %12.6f\n", param[10]);
fprintf(fw, "limb         %12.6f\n", param[11]);
fprintf(fw, "alpha        %12.6f\n", ide_param[21]);
fprintf(fw, "beta         %12.6f\n", ide_param[22]);
fprintf(fw, "rcloud       %12.6f\n", ide_param[23]);
fprintf(fw, "zshift       %12.6f\n", param[15]);
fprintf(fw, "ntable       %12d\n", (int) param[16]);
fprintf(fw, "nrad         %12d\n", (int) param[17]);
fprintf(fw, "division     %12d\n", (int) param[18]);
fprintf(fw, "nphi         %12d\n", (int) param[19]);
fprintf(fw, "ne_loc       %12d\n", (int) param[20]);
fprintf(fw, "smooth       %12d\n", (int) param[21]);
fprintf(fw, "Stokes       %12d\n", (int) param[22]);
fprintf(fw, "norm type    %12.6f\n", param[24]);
fprintf(fw, "polar        %12d\n", polar);
fprintf(fw, "r_horizon    %12.6f\n", r_plus);
fprintf(fw, "r_ms         %12.6f\n", rms);
fprintf(fw, "edivision    %12d\n", (int) ide_param[14]);
fprintf(fw, "nthreads     %12d\n", (int) ide_param[20]);
fclose(fw);
#endif
/******************************************************************************/

/*******************************************************************************
// local spectrum output -- write ear_local[] and far_local[] into file:
fw = fopen("kynconv_photar_loc.dat", "w");
for (ie = 0; ie < ne_loc; ie++)
  fprintf(fw, "%11.6f\t%11.6f\n", ear_local[ie], far_local[ie]);
fclose(fw);
*******************************************************************************/

// Let's integrate local emission over the accretion disc
if (ide(ear, ne, 1, far, qar, uar, var, ide_param, emis_KYNconv, ne_loc)) {
  for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
  return;
}

// interface with XSPEC
// let's renormalize the output to match the "input" flux
for (ie = 0; ie < ne; ie++) far[ie]*=norm0;

if (!stokes) for (ie = 0; ie < ne; ie++) photar[ie] = far[ie];
else {
// final spectrum output -- write ear[] and photar[] into file:
  fw = fopen("stokes.dat", "w");
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
fw = fopen("kynconv_photar.dat", "w");
for (ie = 0; ie < ne; ie++) fprintf(fw, "%14.6f\t%E\n",
                                    0.5 * (ear[ie] + ear[ie+1]),
                                    photar[ie] / (ear[ie+1] - ear[ie]));
fclose(fw);
#endif
/******************************************************************************/

return;
}

/*******************************************************************************
*******************************************************************************/

void emis_KYNconv(double** ear_loc, const int ne_loc, const int nt, 
                  double *far_loc, double *qar_loc, double *uar_loc, 
                  double *var_loc, const double r, const double phi, 
                  const double cosmu, const double phiphoton, 
                  const double alpha_o, const double beta_o, 
                  const double delay, const double g) {
  
double rq;
int ie;

*ear_loc = ear_local;
if (r >= rb) rq = 1. / pow(r,qout);
else rq = fjump * pow(rb, qin - qout) / pow(r, qin);
if ((cosin != 0.) && (cosin != -1.) && (cosin != -2.)) 
  rq = rq * pow(cosmu, cosin);
if (cosin == -1.) rq *= (1. + 2.06 * cosmu);
if (cosin == -2.) rq *= log(1. + 1. / cosmu);
for (ie = 0; ie < ne_loc; ie++) {
  far_loc[ie] = far_local[ie] * rq;
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
