/*
 * ide - Integrating Disc Emission for KYN package of XSPEC models
 *     - C integration subroutine for Keplerian dics emission
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
 * This subroutine integrates local emission and local Stokes parameters
 * for (partially) polarized emission of the accretion disc near rotating
 * (Kerr) black hole (characterized by angular momentum a/M) for an observer
 * with inclination angle theta_o.
 * 
 * Six transfer functions are needed for the integration:
 * a) delay (only for non-stationary disc emission, e.g. orbiting spot),
 * b) g-factor (energy shift),
 * c) cosine of the local emission angle (this is needed only when disc
 *    emission depends on it),
 * d) lensing,
 * e) change of polarization angle (only needed for Stokes parameters)
 * f) azimuthal emission angle (only needed when emission depends on this angle)
 * 
 * These functions differ for different spin, a/M, and inclination, theta_o. 
 * Some are read and interpolated for particular a/M and theta_o from the fits 
 * file called 'KBHtables80.fits', others are computed from other functions
 * stored in this FITS file. The format of this fits file is described in detail
 * below.
 * 
 * By 'local' it is meant 'with respect to the local inertial frame
 * connected with the fluid in the accretion disc' everywhere in this code
 * (local emission, local angle of emitted ray, local disc normal, ...).
 * The local frame is defined in the code not in the tables in KBHtables80.fits
 * file.
 * 
 * The subroutine 'ide' has 10 parameters:
 * ear - array of energy bins (same as 'ear' for local XSPEC models)
 * ne  - number of energy bins (same as 'ne' for local XSPEC models)
 * nt  - number of grid points in time (nt=1 means stationary model)
 * far[ne,nt] - array of photon number density flux per bin
 *              (same as 'photar' for local XSPEC models but with time
 *               resolution)
 * qar[ne,nt] - array of Stokes parameter Q divided by energy
 *              ("photon number density flux per bin for Stokes parameter Q")
 *            - it characterizes a linear polarization perpendicular or
 *              parallel to the disc
 * uar[ne,nt] - array of Stokes parameter U divided by energy
 *              ("photon number density flux per bin for Stokes parameter U")
 *            - it characterizes a linear polarization in the direction +45
 *              (if positive) or -45 degrees (if negative) with the direction
 *              perpendicular to the disc (the angle +45 means that when we are
 *              looking towards the coming emitted light beam, the direction of
 *              polarization is 45 degrees counter-clockwise from the "up"
 *              direction...)
 * var[ne,nt] - array of Stokes parameter V divided by energy
 *              ("photon number density flux per bin for Stokes parameter V")
 *            - it characterizes circular polarization - counter-clockwise
 *              (right-handed) if positive and clockwise (left-handed) if
 *              negative
 * ide_param  - 25 more parameters needed for integration (explained below)
 * emissivity - name of an external emissivity subroutine, where local emission
 *              of the disc is defined (explained in detail below)
 * ne_loc - number of points (in energies) where local photon number density 
 *          flux (per keV) is defined in emissivity subroutine
 * 
 * For an example on how to use 'ide' and 'emissivity' subroutines in local 
 * models see e.g. relativistic line model "xsKYNrline.c".
 *
 * -----------------------------------------------------------------------------
 * 
 * ide_param:
 * 
 * ide_param[0]: am     - black hole angular momentum (-1 <= am <= 1)
 * ide_param[1]: thetaO - observer inclination in degrees (0-pole, 90-equator)
 * ide_param[2]: rin    - inner edge of non-zero disc emissivity (in GM/c^2 or 
 *                        in r_mso)
 * (int) ide_param[3]:  ms - switch that defines the meaning/units of rin, rout
 *                           0 - we integrate from inner edge = rin
 *                           1 - if the inner edge of the disc is below 
 *                               marginally stable orbit then we integrate
 *                               emission above MSO only
 *                           2 - we integrate from inner edge given in units of 
 *                               MSO, i.e. inner edge = rin * r_mso (the same 
 *                               applies for outer edge)
 * ide_param[4]: rout - outer edge of non-zero disc emissivity (in GM/c^2 or in 
 *                      r_mso)
 *                    - if rin > rout => zero total flux
 * ide_param[5]: phi0 - lower azimuth of non-zero disc emissivity (deg)
 *                      -180 <= phi0 <= 180 (but may be any other angle too)
 * ide_param[6]: dphi - (phi0+dphi) = upper azimuth of non-zero disc emissivity
 *                      0 <= dphi <= 360 (deg)
 * (int) ide_param[7]:  nrad      - number of grid points in radius
 * (int) ide_param[8]:  rdivision - type of division in radial integration
 *                      0 -> equidistant radial grid (constant linear step)
 *                      1 -> exponential radial grid (constant logarithmic step)
 * (int) ide_param[9]:  nphi    - number of grid points in azimuth
 * (int) ide_param[10]: smooth  - whether to smooth the resulting spectrum
 *                                (0-no, 1-yes)
 * ide_param[11]: normal - how to normalize the final spectrum
 *                         = 0.  - normalization to unity (total flux=1.)
 *                                 (e.g. for line)
 *                         > 0.  - normalization to the flux at 'normal' keV
 *                                 (usually used for continuum)
 *                         = -1. - final spectrum is not normalized in ide
 *                         = -2. - final spectrum is normalized to have maximum
 *                                 photon number density flux equal unity
 * ide_param[12]:       zshift - overall Doppler shift
 * (int) ide_param[13]: ntable - table of transfer functions used in the model
 *                               (defines fits file with tables),
 *                               0<= ntable <= 99
 * (int) ide_param[14]: edivision - type of division in local energies
 *                                  (0-equidistant, 1-exponential, i.e. more
 *                                   points with lower energy)
 * (int) ide_param[15]: periodic  - has to be set to 1 for periodic emission
 *                                  (e.g. orbiting spot) and 0 for all other 
 *                                  cases
 *                                - need not to be set if nt=1
 * ide_param[16]:       dt - time step (in GM/c^2)
 *                         - need not to be set if nt=1
 * (int) ide_param[17]: polar - whether polarization from fits tables is needed
 *                              1 - yes, 0 - no
 *                              (we usually do not need polarization in XSPEC
 *                              and so we need not to read this information from
 *                              FITS file - the tables will occupy less space in
 *                              memory)
 * ide_param[18]: delay_r   - (r-r_plus) - delay at this radial coord. will be 
 *                                         subtracted in non-stationary 
 *                                         calculations
 * ide_param[19]: delay_phi - phi - delay at this azimuthal BL. coord. 
 *                                  (given in degrees) will be subtracted in 
 *                                  non-stationary calculations
 * ---> the last two parameters are used if we want to set the beginning of
 *      time axis, e.g. for the moment when photons emitted from the spot
 *      that was at the closest approach to the observer come to the observer
 * ide_param[20]: nthreads  - how many threads should be used for computations
 * ide_param[21]: alpha  - position of the cloud in alpha impact parameter 
 *                         (in GM/c^2)
 * ide_param[22]: beta   - position of the cloud in beta impact parameter 
 *                         (in GM/c^2)
 * ide_param[23]: rcloud - radius of the cloud (in GM/c^2)
 *                       - if negative, only the emission transmitted through
 *                         the cloud is taken into account
 * (int) ide_param[24]: observed_flux - whether the flux defined in emissivity 
 *                                      subroutine is local one (0) or the
 *                                      observed one (1)
 * 
 * NOTE: accuracy vs. speed trade off depends mainly on: nrad, nphi
 * 
 * -----------------------------------------------------------------------------
 * 
 * emissivity subroutine:
 * 
 * emissivity(ear_loc,ne_loc,nt,far_loc,qar_loc,uar_loc,var_loc,
 *            rnow,pnow,cosine,phiphoton,alpha_o,beta_o,delay,g)
 * 
 * External emissivity subroutine has 14 parameters:
 * ear_loc[ne_loc] - array of local energies where local photon number density
 *                     flux far_loc is defined, ear_loc[ne_loc] is not used
 * ne_loc - number of points where local photon number density flux is defined
 *          (in energies)
 * nt     - number of grid points in time (nt=1 means stationary model)
 * far_loc[ne_loc+1,nt] - array of local photon number density flux (per keV)
 *                        for 'each time'
 *                      - if the local emissivity is zero for a certain time
 *                        (jt) for all energies we can put far_loc[ne_loc,jt]=0
 *                        and then we exclude adding zeros to the
 *                        total flux (speeding up the computations), if local 
 *                        flux for jt is NOT zero for all energies then 
 *                        far_loc[ne_loc,jt] MUST BE NONZERO!!!
 * qar_loc[ne_loc,nt] - array of local Stokes parameter Q divided by local
 *                      energy for 'each time' ("local photon number density
 *                      flux per keV for Stokes parameter Q")
 *                    - it characterizes a linear polarization perpendicular
 *                      or parallel to the disc in local frame
 * uar_loc[ne_loc,nt] - array of local Stokes parameter U divided by local
 *                      energy for 'each time' ("local photon number density
 *                      flux per keV for Stokes parameter U")
 *                    - it characterizes a linear polarization in the direction
 *                      +45 (if positive) or -45 (if negative) degrees with the
 *                      direction perpendicular to the disc in local frame (the
 *                      angle +45 means that when we are looking toward coming
 *                      emitted light beam, the direction of polarization is 45
 *                      degrees counter-clockwise from the "up" direction...)
 * var_loc[ne_loc,nt] - array of local Stokes parameter V divided by local
 *                      energy for 'each time' ("local photon number density
 *                      flux per keV for Stokes parameter V")
 *                    - it characterizes circular polarization -
 *                      counter-clockwise (right-handed) if positive and
 *                      clockwise (left-handed) if negative in local frame
 * rnow   - radius where the local photon flux far_loc (or Stokes parameters
 *          qar_loc, uar_loc, var_loc) at local energy ear_loc is wanted
 * pnow   - azimuth where the local photon flux far_loc (or Stokes parameters
 *          qar_loc, uar_loc, var_loc)at local energy ear_loc is wanted
 * cosine - cosine of the local angle between emitted ray and local disc normal
 * phiphoton  - angle between emitted ray projected onto the plane of the disc
 *              (in the local frame of the moving disc) and radial component of
 *              the local tetrade (in rad)
 * alpha_o - impact parameter alpha (x-coordinate on the sky of the observer), 
 *           correspondent to light ray emitted from the rnow, pnow coordinates
 *           on the disc, positive for approaching side of the disc (in GM/c^2)
 * beta_o  - impact parameter beta  (y-coordinate on the sky of the observer), 
 *           correspondent to light ray emitted from the rnow, pnow coordinates
 *           on the disc, positive in up direction, i.e. above the disc 
 *           (in GM/c^2)
 * delay   - delay with which photon arrives at the observer from the rnow, pnow
 *           coordinates on the disc
 * g - g-factor, overall redshift between the observed energy and local energy
 *     of the photon when emitted from the rnow, pnow coordinates on the disc, 
 *     g=E_observed/E_local
 *
 * IMPORTANT NOTE: while integrated arrays far, qar and var are evaluated per
 *                 energy bin (i.e. these quantities are integrated over energy
 *                 bin), local arrays far_loc, qar_loc and var_loc are evaluated
 *                 per keV (i.e. they are NOT integrated over energy bin)!!!
 * 
 * -----------------------------------------------------------------------------
 * 
 * Functions in the fits file KBHtablesNN.fits
 * 
 * Functions stored in this file are defined here as tables for different values
 * of horizon (r_horizon, not a/M!) and observer inclination angle theta_o. 
 * Each table consists of values of particular function for different r (rows) 
 * and phi (columns) coordinates. Particular r_horizon, theta_o, r and phi where
 * the functions are given are defined at the beginning of the fits file as 
 * vectors.
 * 
 * Definition of KBHtablesNN.fits:
 * 0. All of the extensions defined below are binary.
 * 1. The first extension contains 5 integers defining which of the functions
 *    are present in the tables. The integers correspond to alpha,
 *    beta, pr, delay and lensing, respectively. Value 0 means that the 
 *    function is not present in the tables, value 1 means it is.
 * 2. The second extension contains vector of horizon values.
 *    (1.00 <= r_horizon <= 2.00)
 * 3. The third extension contains vector of values of the observer's
 *    inclination angle in degrees. (0 <= theta_o <= 90, 0-pole, 90-disc)
 * 4. The fourth extension contains vector of (r-r_horizon) values.
 *    Note it is not r-coordinate itself but relative 'coordinate distance' from
 *    the horizon!!! (0 <= r < infinity, 0 - horizon)
 * 5. The fifth extension contains vector of phi values. Subroutine ide expects
 *    phi to be Kerr ingoing coordinate, not Boyer-Lindquist one!
 *    (0 <= phi <= 2*PI)
 * 6. All previous vectors have to have values sorted in increasing order!
 * 7. In following extensions the functions are defined, each extension
 *    is for particular value of r_horizon and theta_o. The values of r_horizon
 *    and theta_o are changing with each extension in the following order:
 *    r_horizon[1] x theta_o[1],
 *    r_horizon[1] x theta_o[2],
 *    r_horizon[1] x theta_o[3],
 *    ...
 *    ...
 *    r_horizon[2] x theta_o[1],
 *    r_horizon[2] x theta_o[2],
 *    r_horizon[2] x theta_o[3],
 *    ...
 *    ...
 * 8. Each of these extensions has the same number of columns (up to 5). In each
 *    column a particular function is defined - alpha, beta, pr, delay and 
 *    lensing, repsectively.
 *    The order of the functions is important but some of the functions may be
 *    missing as defined in the first extension (see 1. above).
 * 
 *    alpha - impact parameter alpha (x-coordinate on the sky of the observer), 
 *            correspondent to light ray emitted from the r, phi coordinates
 *            on the disc, positive for approaching side of the disc (in GM/c^2)
 *    beta  - impact parameter beta  (y-coordinate on the sky of the observer), 
 *            correspondent to light ray emitted from the r, phi coordinates
 *            on the disc, positive in up direction, i.e. above the disc 
 *            (in GM/c^2)
 *    pr    - the radial component of the photon momentum at the disc
 *    delay - Boyer-Lindquist time that elapses between emission of a photon 
 *            from the disc and absorption of the photon by the observers eye at
 *            infinity minus a particular constant (so that delay = zero close 
 *            to the x-axis in the equatorial plane, observer is at infinity in 
 *            the direction of y-axis with x=0, the black hole is at x=y=z=0)
 *    lensing - ratio of the flux tube cross-section at infinity to the flux 
 *              tube cross-section at the disc
 * 
 * The following functions are computed from these in the code:
 * 
 * g-factor - ratio of the energy of a photon received by the observer at
 *            infinity to the local energy of the same photon when emitted
 *            from an accretion disc
 * cosine - cosine of the local angle between the emitted light ray and local
 *          disc normal
 * phi_photon - angle between emitted ray projected onto the plane of the disc
 *              (in the local frame of the moving disc) and radial component of
 *               the local tetrade (in rad)
 * polarization angle
 *  - if the light emitted from the disc is linearly polarized then the
 *    direction of polarization will be changed by this angle in infinity -
 *    counter-clockwise if positive, clockwise if negative (we are looking
 *    toward the coming emitted beam), on the disc we measure the angle of
 *    polarization with respect to "up" direction perpendicular to the disc
 *    with respect to the local rest frame, in infinity we also measure the
 *    angle of polarization with respect to "up" direction perpendicular to
 *    the disc - "polarization angle" is the difference between these two
 *    angles
 * 
 * 9. Each row corresponds to a particular value of r (see 4. above).
 * 10. Each element of an extension is a vector. Each element of this vector
 *     corresponds to a particular value of phi (see 5. above).
 * 
 * For an example of a FITS file defined in this way, see KBHtables80.fits.
 * 
 * -----------------------------------------------------------------------------
 * 
 * KBHtables80.fits
 * 
 * These tables are computed for an optically thick and geometrically thin
 * accretion disc near Kerr black hole. The medium between the disc and the
 * observer is supposed to be optically thin for the wavelengths one is
 * interested in. Therefore ray tracing in vacuum Kerr space-time could be used
 * for calculating the transfer functions.
 * 
 * When calculating the transfer functions, it was supposed that the fluid in
 * the disc rotates on stable circular (free) orbits above marginally stable
 * orbit (MSO). Below MSO the fluid is freely falling and has the same energy
 * and angular momentum as the matter which is on the MSO.
 * 
 * The observer is placed in the direction phi = PI/2. The black hole rotates
 * counter-clockwise.
 * 
 * - all of the functions are present in these tables
 * - values of r_horizon are: 1.00, 1.05, 1.10, 1.15, ..., 1.90, 1.95, 2.00
 *   (21 elements)
 * - values of theta_o are: 0.1, 1, 5, 10, 15, 20, 25, ..., 75, 80, 85, 89
 *   (20 elements)
 * - values of r are exponentially increasing from 0 to 999 (150 elements)
 * - values of phi are equidistantly spread from 0 to 2*PI rad with much denser
 *   cover 'behind' the black hole, i.e. near phi = 3/2*PI (because some of the
 *   functions, e.g. alpha, beta, lensing, are changing heavily in this area
 *   for higher inclination angles (theta_o > 70 deg)). (200 elements)
 * 
 *******************************************************************************
 *
 * 18.5.2005 interpolation of g-factor changed from bilinear to bicubic
 * 26.9.2006 better time integration
 * 12.3.2010 possibility to have a<0 added
 * 22.3.2010 computing in threads added
 * 25.1.-4.2.2014 conversion from f77 to c
 *
 ******************************************************************************/

//#include <stdio.h>

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "fitsio.h"

typedef double Real;

#define KBHTABLES "KBHtables"
#define PI  3.1415926535897932
#define PI2 6.2831853071795865
#define TRUE 1
#define FALSE 0

extern char* FGMODF(void);
extern char* FGMSTR(char* dname);
extern void  FPMSTR(const char* value1, const char* value2);
extern int   xs_write(char* wrtstr, int idest);

static double am,rms,am2,rms2, r_plus, r_minus;
static double thetaO;
static double zzshift;
static double ddel, dt;
static int rdivision, edivision, polar, periodic, observed_flux;
static int nrad, nphi;
static double rin, rout, phi0, phi1;
/* number of (r,phi) integration loops that one thread in multithread version 
   will compute during one call */
static int nloops;
static double alphac, betac, rcloud, rcloud2;
//static FILE  *outfile;

//Should we initialize these???
static char kydir_old[255]="";
static char kydir[255]="";
static int column[5]={0, 0, 0, 0, 0};
static double *eard;
static double *dc;
static long nr, nph;
static double *r_vec;
static float *phi_vec;

/*******************************************************************************
*******************************************************************************/
int ide(const double *ear, const int ne, const int nt, double *far, 
        double *qar, double *uar, double *var, const double *ide_param, 
        void (*emissivity)(), const int ne_loc){
  
extern void create_ide_threads(const int nthreads, const int nthreads_old, 
                               const int ne, const int ne_old,
                               const int nt, const int nt_old, 
                               const int polar, const int polar_old);
void element(int iloop, const int ne_loc, const int ne, const int nt, 
             double *fce, double *qce, double *uce, double *vce, 
             void (*emissivity)());
extern void ide_threads(const int nthreads, int nloops, const int polar, 
                        const int ne_loc, const int ne, const int nt, 
                        double *fc, double *qc, double *uc, double *vc, 
                        void (*emissivity)(), void (*element)());
int ide_ini(int ntable, const int nt);

static double am_old=-1.;
static double thetaO_old=-1.;
static int    ntable_old=-1;
static int    nt_old=1;
static int    nthreads_old=1;
static int    ne_old=0;
static int    ntt_old=1;
static int    polarr_old=0;
static char pname[128]="KYDIR";
static const char pkyrh[128]="KYRH";
static const char pkyrin[128]="KYRIN";
static const char pkyrms[128]="KYRMS";

double far_int, fc[ne*nt], qc[ne*nt], uc[ne*nt], vc[ne*nt];
double fce[ne*nt], qce[ne*nt], uce[ne*nt], vce[ne*nt];
int    smooth, ntable, ms;
double pom, pom1, pom2, pom3, normal, zshift, dphi, dlograd;
int    nthreads, nnloops, ie, it, imin, imax, ir0, iphi0, iloop;
double delay_r, delay_phi;
double ttmp, ttmp1, utmp, utmp1, y1, y2, y3, y4;
char   text[80], kyrh[32], kyrin[32], kyrms[32];

//let's check if input parameters have reasonable values
am=ide_param[0];
if(fabs(am)>1.){
 xs_write("ide: a/M must be >= -1 and <= 1",5);
 return(1);
}
am2=am*am;
//outer and inner horizons
r_plus=1.+sqrt(1.-am2);
r_minus=1.-sqrt(1.-am2);
// marginally stable orbit in the Kerr geometry:
pom1=pow(1.+am,0.33333333);
pom2=pow(1.-am,0.33333333);
pom3=pow(1.-am2,0.33333333);
pom=1.+pom3*(pom1+pom2);
pom1=sqrt(3.*am2+pom*pom);
if(am>=0)rms=3.+pom1-sqrt((3.-pom)*(3.+pom+2.*pom1));
else rms=3.+pom1+sqrt((3.-pom)*(3.+pom+2.*pom1));
rms2=rms*rms;
thetaO=ide_param[1];
if(thetaO<0. || thetaO>90.){
 xs_write("ide: theta_o must be >= 0 and <= 90",5);
 return(1);
}
ms=(int) ide_param[3];
if(ms!=0 && ms!=1 && ms!=2){
  xs_write("ide: ms must be 0, 1 or 2",5);
  return(1);
}
if(ms==1){
 if(ide_param[2]<rms)rin=rms-r_plus;
 else rin=ide_param[2]-r_plus;
 rout=ide_param[4]-r_plus;
}else if(ms==2){
 rin=ide_param[2]*rms-r_plus;
 rout=ide_param[4]*rms-r_plus;
}else if(ms==0){
 rin=ide_param[2]-r_plus;
 rout=ide_param[4]-r_plus;
}
if(rin<0.)rin=0.;
if(rout<0.)rout=0.;
phi0=(ide_param[5]/180.*PI);
phi1=phi0+ide_param[6]/180.*PI;
nrad=(int) ide_param[7];
if(nrad<1){
  xs_write("ide: nrad must be larger than 0",5);
  return(1);
}
rdivision=(int) ide_param[8];
if(rdivision!=0 && rdivision!=1){
  xs_write("ide: division must be 0 or 1",5);
  return(1);
}
nphi=(int) ide_param[9];
if(nphi<1){
  xs_write("ide: nphi must be larger than 0",5);
  return(1);
}
smooth=(int) ide_param[10];
if(smooth!=0 && smooth!=1){
  xs_write("ide: smooth must be 0 or 1",5);
  return(1);
}
normal=ide_param[11];
zshift=ide_param[12];
if(zshift<=-1.){
  xs_write("ide: zshift must be larger than -1",5);
  return(1);
}
ntable=(int) ide_param[13];
if(ntable<0 || ntable>99){
  xs_write("ide: ntable must be >= 0 and <= 99",5);
  return(1);
}
edivision=(int) ide_param[14];
if(edivision!=0 && edivision!=1){
  xs_write("ide: edivision must be 0 or 1",5);
  return(1);
}
periodic=(int) ide_param[15];
if(periodic!=0 && periodic!=1 && nt>1){
  xs_write("ide: periodic must be 0 or 1",5);
  return(1);
}
dt=ide_param[16];
if(dt<=0 && nt>1){
  xs_write("ide: dt must be larger than 0",5);
  return(1);
}
polar=(int) ide_param[17];
if(polar!=0 && polar!=1){
  xs_write("ide: polar must be 0 or 1",5);
  return(1);
}
delay_r=ide_param[18];
if(delay_r<0 && nt>1){
  xs_write("ide: delay_r must be positive",5);
  return(1);
}
delay_phi=ide_param[19]/180.*PI;
if(ne<1){
  xs_write("ide: ne must be larger than 0",5);
  return(1);
}
if(ne_loc<3){
  xs_write("ide: ne_loc must be larger than 2",5);
  return(1);
}
if(nt<=0){
  xs_write("ide: nt must be larger than 0",5);
  return(1);
}
//number of threads for multithread computations
nthreads=(int) ide_param[20];
if(nthreads<1){
  xs_write("ide: nthreads must be larger than 0",5);
  return(1);
}
//position of obscuring cloud centre and its radius
alphac=ide_param[21];
betac=ide_param[22];
rcloud=ide_param[23];
rcloud2=rcloud*rcloud;
//whether the flux defined in emissivity subroutine is local one (0) or the 
//observed one (1)
observed_flux = (int) ide_param[24];

//initialize some variables:
/* number of (r,phi) integration loops that one thread in multithread
   version will compute during one call */
nloops=100;
// zzshift - multiplication factor for gfac from zshift
zzshift=1.0/(1.0+zshift);
// let's check whether rin, rout, phi0 and phi1 are reasonable
if(rin>=rout)return(1);
if(phi0>=phi1)return(1);
//somehow we have to call FPMSTR to "initialize" it, otherwise we get weird
//results
FPMSTR(pkyrh,kyrh);
sprintf(kyrh,"%g",r_plus);
FPMSTR(pkyrh,kyrh);
sprintf(kyrin,"%g",rin+r_plus);
FPMSTR(pkyrin,kyrin);
sprintf(kyrms,"%g",rms);
FPMSTR(pkyrms,kyrms);

sprintf(kydir,"%s",FGMSTR(pname));
/* initialize data tables dc() (interpolated for given am, thetaO):
   we read the tables if we haven't yet or if we have changed them or am,thetaO
   or nt or polar or if we have changed the KYDIR directory...
   so that we can get new dc() */
if(ntable_old==-1 || strcmp(kydir_old,kydir) || ntable_old!=ntable || 
   am_old!=am || thetaO_old!=thetaO || (nt>1 && nt_old==1 && column[3]==1)){
  if(ide_ini(ntable, nt))return(1);
  am_old=am;
  thetaO_old=thetaO;
  if(nt>1)nt_old=nt;
  if(nt==1 && ntable!=ntable_old)nt_old=1;
  ntable_old=ntable;
  sprintf(kydir_old,"%s",kydir);
}
// now let's check whether rin and rout fall in the range covered by data tables
if(rin<r_vec[0]){
  xs_write("ide: r out of range covered by data tables",5);
  rin=r_vec[0];
}
if(rout>(r_vec[nr-1]+1e-6)){
  xs_write("ide: r out of range covered by data tables",5);
//  rout=r_vec[nr-1];
}
// let's check whether the step in phi is not too small
dphi=(phi1-phi0)/((double) nphi);
if((dphi<1e-6*fabs(phi1) && fabs(phi1)>fabs(phi0)) || 
   (dphi<1e-6*fabs(phi0) && fabs(phi0)>fabs(phi1))){
  if(fabs(phi1)>fabs(phi0))dphi=1e-6*fabs(phi1);
  else dphi=1e-6*fabs(phi0);
  nphi=(int) ((phi1-phi0)/dphi);
  xs_write("ide: too small step in phi, because nphi is too large",5);
  sprintf(text,"    ...changing the step in phi (nphi = %d)",nphi);
  xs_write(text,5);
  if(nphi<1)return(1);
}
/* let's check whether the step in r is not too small
   and initialize some variables needed in integration */
//...for equidistant grid in r
if(rdivision==0){
  if((rout-rin)/((double) nrad)<1e-6*(rout+r_plus)){    
    nrad=(int) ((rout-rin)/(1e-6*(rout+r_plus)));
    xs_write("ide: too small step in r, because nrad is too large",5);
    sprintf(text,"    ...changing the step in r (nrad = %d)",nrad);
    xs_write(text,5);
    if(nrad<1)return(1);
  }
}
//...for equidistant grid in log(r)
if(rdivision==1){
  dlograd=(log10(rout+r_plus)-log10(rin+r_plus))/((double) nrad);
  if(dlograd<2e-6*log10(rout+r_plus)){
    dlograd=2e-6*log10(rout+r_plus);
    nrad=(int) ((log10(rout+r_plus)-log10(rin+r_plus))/dlograd);
    xs_write("ide: too small step in log(r), because nrad is too large",5);
    sprintf(text,"    ...changing the step in log(r) (nrad = %d)",nrad);
    xs_write(text,5);
    if(nrad<1)return(1);
  }
}
//Let's calculate the delay at delay_r, delay_phi in non-stationary calculations
if(nt>1){
// given (r,phi), find corresponding indices in r_vec,phi_vec():
  imin=0;
  imax=nr;
  ir0=nr/2;
  while((imax-imin)>1){
    if(delay_r>=r_vec[ir0-1])imin=ir0;
    else imax=ir0;
    ir0=(imin+imax)/2;
  }
  if(imax==nr && delay_r>r_vec[nr-1])ir0=nr;
//... we must convert phi from Boyer-Lindquist coordinate to Kerr coordinate
  if(fabs(am)==1.)delay_phi+=am/(delay_r+r_plus-1.);
  else delay_phi-=am/(r_plus-r_minus)*log(delay_r/(delay_r+r_plus-r_minus));
//phi in the data tables is >= 0 and <= 2*PI
//...let's modify delay_phi so that it lies in this range
  if(am<0.)delay_phi=PI-delay_phi;
  delay_phi=fmod(delay_phi,PI2);
  if(delay_phi<0.)delay_phi+=PI2;
  imin=0;
  imax=nph;
  iphi0=nph/2;
  while((imax-imin)>1){
    if(delay_phi>=phi_vec[iphi0-1])imin=iphi0;
    else imax=iphi0;
    iphi0=(imin+imax)/2;
  }
  if(imax==nph && delay_phi>phi_vec[nph-1])iphi0=nph;
//let's interpolate the values in dc() table to desired (delay_r,delay_phi):
  ttmp=(delay_r-r_vec[ir0-1])/(r_vec[ir0]-r_vec[ir0-1]);
  utmp=(delay_phi-phi_vec[iphi0-1])/(phi_vec[iphi0]-phi_vec[iphi0-1]);
  ttmp1=1.-ttmp;
  utmp1=1.-utmp;
//time delay at delay_r+r_plus and delay_phi
  y1=dc[ir0-1+nr*(iphi0-1)];
  y2=dc[ir0+nr*(iphi0-1)];
  y3=dc[ir0+nr*iphi0];
  y4=dc[ir0-1+nr*iphi0];
  ddel=utmp1*(ttmp1*y1+ttmp*y2)+utmp*(ttmp*y3+ttmp1*y4);
}
/*******************************************************************************
  Let's write the functions to the file for TESTing purposes...
   outfile=fopen("test2.dat","w");
*******************************************************************************/
//let's integrate...loop over the disc;
//(rnow,pnow) are polar coords (radius in units of GM/c^2):
for(it=0;it<nt;it++){
  for(ie=0;ie<ne;ie++){
    fc[ie+ne*it]=0.;
    if(polar){
      qc[ie+ne*it]=0.;
      uc[ie+ne*it]=0.;
      vc[ie+ne*it]=0.;
    }
  }
}
eard=ear;
//let's create, destroy or modify the threads
if(nthreads!=nthreads_old || ne!=ne_old || nt!=ntt_old || polar!=polarr_old){
  create_ide_threads(nthreads, nthreads_old, ne, ne_old, nt, ntt_old,
                     polar, polarr_old);
  nthreads_old=nthreads;
  ne_old=ne;
  ntt_old=nt;
  polarr_old=polar;
}
if(nthreads==1){
//let's clear the element flux contribution arrays
  for(ie=0;ie<ne;ie++){
    for(it=0;it<nt;it++){
      fce[ie+ne*it]=0.;
      if(polar){
        qce[ie+ne*it]=0.;
        uce[ie+ne*it]=0.;
        vce[ie+ne*it]=0.;
      }
    }
  }
  iloop=1;
  nloops=nrad*nphi;
  element(iloop,ne_loc,ne,nt,fce,qce,uce,vce,emissivity);
  for(it=0;it<nt;it++){
    for(ie=0;ie<ne;ie++){
      fc[ie+ne*it]+=fce[ie+ne*it];
      if(polar){
        qc[ie+ne*it]+=qce[ie+ne*it];
        uc[ie+ne*it]+=uce[ie+ne*it];
        vc[ie+ne*it]+=vce[ie+ne*it];
      }
    }
  }
}else{
//if we have more CPU or CPU cores we will use more threads for the computation
  nnloops=(nrad*nphi)/nloops+1;
  ide_threads(nthreads, nnloops, polar, ne_loc, ne, nt, fc, qc, uc, vc, 
              emissivity, element);
}
//let's smooth the spectra a little bit if it is required
if(smooth==0){
  for(it=0;it<nt;it++){
    for(ie=0;ie<ne;ie++){
      far[ie+ne*it]=fc[ie+ne*it];
      if(polar){
        qar[ie+ne*it]=qc[ie+ne*it];
        uar[ie+ne*it]=uc[ie+ne*it];
        var[ie+ne*it]=vc[ie+ne*it];
      }
    }
  }
}else{
  for(it=0;it<nt;it++){
    far[ne*it]=fc[ne*it];
    far[ne-1+ne*it]=fc[ne-1+ne*it];
    for(ie=1;ie<ne-1;ie++)
      far[ie+ne*it]=(fc[ie-1+ne*it]/(ear[ie]-ear[ie-1])+
      2.0*fc[ie+ne*it]/(ear[ie+1]-ear[ie])+
      fc[ie+1+ne*it]/(ear[ie+2]-ear[ie+1]))/4.0*(ear[ie+1]-ear[ie]);
//let's smooth the other Stokes parameters if they are being calculated
    if(polar){
      qar[ne*it]=qc[ne*it];
      uar[ne*it]=uc[ne*it];
      var[ne*it]=vc[ne*it];
      qar[ne-1+ne*it]=qc[ne-1+ne*it];
      uar[ne-1+ne*it]=uc[ne-1+ne*it];
      var[ne-1+ne*it]=vc[ne-1+ne*it];
      for(ie=1;ie<ne-1;ie++){
        qar[ie+ne*it]=(qc[ie-1+ne*it]/(ear[ie]-ear[ie-1])+
          2.0*qc[ie+ne*it]/(ear[ie+1]-ear[ie])+
          qc[ie+1+ne*it]/(ear[ie+2]-ear[ie+1]))/4.0*(ear[ie+1]-ear[ie]);
        uar[ie+ne*it]=(uc[ie-1+ne*it]/(ear[ie]-ear[ie-1])+
          2.0*uc[ie+ne*it]/(ear[ie+1]-ear[ie])+
          uc[ie+1+ne*it]/(ear[ie+2]-ear[ie+1]))/4.0*(ear[ie+1]-ear[ie]);
        var[ie+ne*it]=(vc[ie-1+ne*it]/(ear[ie]-ear[ie-1])+
          2.0*vc[ie+ne*it]/(ear[ie+1]-ear[ie])+
          vc[ie+1+ne*it]/(ear[ie+2]-ear[ie+1]))/4.0*(ear[ie+1]-ear[ie]);
      }
    }
  }
}
//let's normalize the spectrum if it is required (i.e. normal >= 0.)
if(normal>=0. || normal==-2.){
  for(it=0;it<nt;it++){
//  if normal is 0 or if normal is out of energy range ear[0]-ear[ne]
//  we normalize the spectrum to total flux equal to unity
    if((normal>0. && (ear[0]>normal || ear[ne]<normal)) || normal==0.){
      far_int=0.;
      for(ie=0;ie<ne;ie++)far_int+=far[ie+ne*it];
    }else if(normal==-2.){
      far_int=0.;
      for(ie=0;ie<ne;ie++){
        if(far[ie+ne*it]/(ear[ie+1]-ear[ie])>far_int)
          far_int=far[ie+ne*it]/(ear[ie+1]-ear[ie]);}
    }else{
/*  otherwise we normalize to the maximum at normal keV (i.e. the photon flux
    density per keV is 1 at normal keV) */
//  let's find the energy bin where normal falls
      for(ie=0;ie<ne;ie++)
        if(ear[ie]<=normal && ear[ie+1]>=normal)break;
/*    let's calculate flux at normal keV - we interpolate two neighbouring
      values (that are at the middle of energy bins) */
      if(((ear[ie+1]+ear[ie])/2.>=normal && ie>0) || 
         ((ear[ie+1]+ear[ie])/2.<normal && ie==ne-1))ie=ie-1;
      far_int=far[ie+ne*it]/(ear[ie+1]-ear[ie])+
        (2.*normal-(ear[ie+1]+ear[ie]))/(ear[ie+2]-ear[ie])*
        (far[ie+1+ne*it]/(ear[ie+2]-ear[ie+1])-
         far[ie+ne*it]/(ear[ie+1]-ear[ie]));
//    if the flux at normal keV is too small (zero) we normalize the spectrum
//    to total flux equal to unity
      if(far_int<1e-30){
        far_int=0.;
        for(ie=0;ie<ne;ie++)far_int+=far[ie+ne*it];
      }
    }
    if(far_int>1e-30){
      for(ie=0;ie<ne;ie++){
        far[ie+ne*it]=far[ie+ne*it]/far_int;
        if(polar){
          qar[ie+ne*it]=qar[ie+ne*it]/far_int;
          uar[ie+ne*it]=uar[ie+ne*it]/far_int;
          var[ie+ne*it]=var[ie+ne*it]/far_int;
        }
      }
    }
  }
}
/*******************************************************************************
  Let's write the functions to the file for TESTing purposes...
   fclose(outfile);
*******************************************************************************/
return(0);
}
/*******************************************************************************
*******************************************************************************/
/*
   read the data tables from the fits file and
   initialize data table dc(:,:,:) for given ntable and parameters am and thetaO
*/
int ide_ini(int ntable, const int nt){

static int column_old[5]={0,0,0,0,0};
static int ntable_old=-1;
static int nt_old=1;
static int nhorizon_old=0;
static int nincl_old=0;
static int nr_old=0;
static int nph_old=0;
static long   nhorizon, nincl;
static int    ncolumns;
static float  *rmvec, *thvec, *alpha0, *beta0, *pr0, *delay, *lensing;
static double *amvec;
static int  read_alpha=FALSE, read_beta=FALSE, read_pr=FALSE, read_del=FALSE, 
            read_lens=FALSE;
double alpha, beta, pr;
double gfac, cosine, polarization, phiphoton, thetaOO;
int  ihorizon, iincl, irow, icol;
int  iph, ir, imin, imax, iam0, ith0, i;
double ttmp, ttmp1, utmp, utmp1, y1, y2, y3, y4;
double kr, kphi, deltar, deltaphi, dg1, dg2, dg3, dg4, dg5, dg6, dg7, dg8;
double l, q2, p_th, sqr_ms, Ems, Lms;
double rnow, rnow2, delta, sqr, p_r, pt, pth, pph, AA;
double Ut, Ur, Uph, kappa1, kappa2, ft, fr, fth, fph;
//FILE   *outfile1;
// these are needed to work with a fits file...
char     tables_file[255];
fitsfile *fptr;
int      hdutype=2;
int      colnum=1;
long     frow=1, felem=1, nelems, nrow;
float    float_nulval=0.;
int      anynul, status=0;//, maxdim=1000, naxis;

/* reading table files including vectors r_vec and phi_vec -- execute on the
   first pass only or when we choose different tables or if we change nt=1 to
   nt>1 ...
   if we change KYDIR directory and we have already read tables (i.e.
   ntables_old!=-1) we will change ntable_old so that we read new ones
   and re-allocations of the arrays are performed properly... */
if(strcmp(kydir_old,kydir) && ntable_old!=-1)ntable_old=-2;
if(ntable_old==-1 || ntable!=ntable_old || (nt>1 && nt_old==1 && 
   column[3]==1)){
  xs_write("ide: initializing data tables, please wait...",5);
  xs_write("     Ref.: Dovciak M., Karas V. & Yaqoob T.",5);
  xs_write("           ApJS July 2004, Volume 153, Issue 1, pp. 205-221",5);
  xs_write("     ------------------------------------------------------",5);
// Set up the directory and the fits file name of the tables
// and open the fits file with the tables
  if(strlen(kydir)==0)sprintf(tables_file,"./%s%2.2d.fits",KBHTABLES,ntable);
  else if(kydir[strlen(kydir)-1]=='/')sprintf(tables_file,"%s%s%2.2d.fits",
                                              kydir,KBHTABLES,ntable);
  else sprintf(tables_file,"%s/%s%2.2d.fits",kydir,KBHTABLES,ntable);
//The status parameter must always be initialized.
  status=0;
  ffopen(&fptr, tables_file, READONLY, &status);
  if(status!=0){
    sprintf(tables_file,"%s%s",FGMODF(),KBHTABLES);
    status=0;
    ffopen(&fptr, tables_file, READONLY, &status);
  }
  if(status!=0){
    if(status)ffrprt(stderr, status);
    ffclos(fptr, &status);
    xs_write("\nide: set the KYDIR to the directory with the KY tables",5);
    ntable=-1;
    return(1);
  }
//Let's read tables (binary tables => hdutype=2)
// Move to the extension 'tables' and read its values
  ffmrhd(fptr, 1, &hdutype, &status);
  if(ntable_old==-1 || ntable!=ntable_old){
    nelems=5;
//  ffgcv reads the VALUES from the first column.
    ffgcv(fptr, TINT, colnum, frow, felem, nelems, &float_nulval, column,
          &anynul, &status);
    ncolumns=0;
    for(i=0;i<5;i++)ncolumns+=column[i];
/*******************************************************************************
   for(i=0;i<5;i++)fprintf(stdout, "%d ", column[i]);
   fprintf(stdout, "\n%d\n", ncolumns);
*******************************************************************************/
  }
// Move to the extension 'r_horizon' and read r_horizon values
  ffmrhd(fptr, 1, &hdutype, &status);
  if(ntable_old==-1 || ntable!=ntable_old){
//  read number of r_horizon values
    ffgnrw(fptr, &nhorizon, &status);
/*******************************************************************************
    fprintf(stdout, "\nNumber of r_horizon values: %d\n", nhorizon);
*******************************************************************************/
/*  We allocate memory for arrays only if we read the tables for the first
    time or if we have changed the dimensions of the arrays (in case the
    tables have changed) */
    if(ntable_old==-1 || nhorizon!=nhorizon_old){
/*    Firstly we have to free allocated memory for these arrays if we have
      already allocated it and the dimensions of the arrays have changed */
      if((ntable_old!=-1) && (nhorizon!=nhorizon_old)){
//      Free memory from tmp arrays...
        free(amvec);
        free(rmvec);
      }
//    Allocate memory for amvec and rmvec...
      amvec=(double *) malloc((size_t) nhorizon*sizeof(double*));
      if(!amvec){xs_write("failed to allocate memory for tmp arrays...",5);
        ntable=-1;return(1);}
      rmvec=(float *) malloc((size_t) nhorizon*sizeof(float*));
      if(!rmvec){xs_write("failed to allocate memory for tmp arrays...",5);
        ntable=-1;return(1);}
    }
// Read the data in the 'r_horizon' table
    nelems=nhorizon;        
//  ffgcv reads the VALUES from the first column.
    ffgcv(fptr, TFLOAT, colnum, frow, felem, nelems, &float_nulval, rmvec,
          &anynul, &status);
//  let's calculate am values for r_horizon values
    for(i=0;i<nhorizon;i++){
      amvec[i]=sqrt(1.-(rmvec[i]-1.)*(rmvec[i]-1));
/*******************************************************************************
      fprintf(stdout, "%d %g %g\n", i, rmvec[i], amvec[i]);
*******************************************************************************/
    }
//  check if current value of am is in the range covered by data tables
    if(nhorizon>1 && ((fabs(am)<amvec[nhorizon-1]) || (fabs(am)>amvec[0]))){
      xs_write("ide: a out of table range",5);
      ntable=-1;return(1);}
    if(nhorizon==1 && ((fabs(fabs(am)-amvec[0]))>1e-4))
      xs_write("ide: a/M out of table range",5);
  }
//Move to the extension 'inclination' and read inclination values
  ffmrhd(fptr, 1, &hdutype, &status);
  if(ntable_old==-1 || ntable!=ntable_old){
//  read number of inclination values
    ffgnrw(fptr, &nincl, &status);
/*******************************************************************************
    fprintf(stdout, "\nNumber of inclination values: %d\n", nincl);
*******************************************************************************/
/*  We allocate memory for arrays only if we read the tables for the first
    time or if we have changed the dimensions of the arrays (in case the
    tables have changed) */
    if(ntable_old==-1 || nincl!=nincl_old){
/*    Firstly we have to free allocated memory for these arrays if we have
      already allocated it and the dimensions of the arrays have changed */
//    Free memory from tmp arrays...
      if((ntable_old!=-1) && (nincl!=nincl_old))free(thvec);
//    Allocate memory for thvec...
      thvec=(float *) malloc((size_t) nincl*sizeof(float*));
      if(!thvec){xs_write("failed to allocate memory for tmp arrays...",5);
        ntable=-1;return(1);}
    }
//  Read the data in the 'inclination' table
    nelems=nincl;
//  ffgcv reads the VALUES from the first column.
    ffgcv(fptr, TFLOAT, colnum, frow, felem, nelems, &float_nulval, thvec,
          &anynul, &status);
/*******************************************************************************
    for(i=0;i<nincl;i++)fprintf(stdout,"%d %g\n", i, thvec[i]);
*******************************************************************************/
/*  check if current value of thetaO is in the range covered by data tables
    if thetaO is "almost" 0 and fits tables are defined for higher values
    then thetaO but from close to 0 we change thetaO to the lowest value
    in the fits tables...(for convenience we put 0.1 equal to 0.0 ...) */
    thetaOO=thetaO;
    if(((thvec[0]-thetaO)<=0.1) && (thetaO<=0.1))thetaOO=(double) thvec[0];
    if((thetaOO<thvec[0]) || (thetaOO>thvec[nincl-1])){
     xs_write("ide: theta_o out of table range",5);ntable=-1;return(1);}
  }
//Move to the extension 'r_vector' and read r_vec values
  ffmrhd(fptr, 1, &hdutype, &status);
  if(ntable_old==-1 || ntable!=ntable_old){
//  read number of r_vector values
    ffgnrw(fptr, &nr, &status);
/*******************************************************************************
    fprintf(stdout, "\nNumber of r_vec values: %d\n", nr);
*******************************************************************************/
/*  We allocate memory for arrays only if we read the tables for the first
    time or if we have changed the dimensions of the arrays (in case the
    tables have changed) */
    if(ntable_old==-1 || nr!=nr_old){
/*    Firstly we have to free allocated memory for these arrays if we have
      already allocated it and the dimensions of the arrays have changed */
//    Free memory from tmp arrays...
      if((ntable_old!=-1) && (nr!=nr_old))free(r_vec);
//    Allocate memory for r_vec...
      r_vec=(double *) malloc((size_t) nr*sizeof(double*));
      if(!r_vec){xs_write("failed to allocate memory for tmp arrays...",5);
        ntable=-1;return(1);}
    }
//  Read the data in the 'r_vector' table
    nelems=nr;
//  ffgcv reads the VALUES from the first column.
    ffgcv(fptr, TDOUBLE, colnum, frow, felem, nelems, &float_nulval, r_vec,
          &anynul, &status);
/*******************************************************************************
    for(i=0;i<nr;i++)fprintf(stdout, "%d %g\n", i, r_vec[i]);
*******************************************************************************/
  }
//Move to the extension 'phi_vector' and read phi_vec values
  ffmrhd(fptr, 1, &hdutype, &status);
  if(ntable_old==-1 || ntable!=ntable_old){
//  read number of phi_vector values
    ffgnrw(fptr, &nph, &status);
/*******************************************************************************
    fprintf(stdout, "\nNumber of phi_vec values: %d\n", nph);
*******************************************************************************/
/*  We allocate memory for arrays only if we read the tables for the first
    time or if we have changed the dimensions of the arrays (in case the
    tables have changed) */
    if(ntable_old==-1 || nph!=nph_old){
/*    Firstly we have to free allocated memory for these arrays if we have
      already allocated it and the dimensions of the arrays have changed */
//    Free memory from tmp arrays...
      if((ntable_old!=-1) && (nph!=nph_old))free(phi_vec);
//    Allocate memory for phi_vec...
      phi_vec=(float *) malloc((size_t) nph*sizeof(float*));
      if(!phi_vec){xs_write("failed to allocate memory for tmp arrays...",5);
        ntable=-1;return(1);}
    }
//  Read the data in the 'phi_vector' table
    nelems=nph;
//  ffgcv reads the VALUES from the first column.
    ffgcv(fptr, TFLOAT, colnum, frow, felem, nelems, &float_nulval, phi_vec,
          &anynul, &status);
/*******************************************************************************
    for(i=0;i<nph;i++)fprintf(stdout, "%d %g\n", i, phi_vec[i]);
*******************************************************************************/
  }
/*Let's read the transfer functions...
  We need to decide which of the functions are to be read:
  we read all present in the fits file if we read the tables for the first
  time or if the tables have changed (ntable has changed)
  - we read delay only if nt>1 */
  if(ntable_old==-1 || ntable!=ntable_old){
    if(column[0]==1)read_alpha=TRUE;
    if(column[0]==0){
      xs_write("ide: alpha is not present in the fits tables,",5);
      xs_write("     we will use alpha = cos(phi)*sin(theta_o))...",5);}
    if(column[1]==1)read_beta=TRUE;
    if(column[1]==0){
      xs_write("ide: beta is not present in the fits tables,",5);
      xs_write("     we will use beta = cos(theta_o)...",5);}
    if(column[2]==1)read_pr=TRUE;
    if(column[2]==0){
      xs_write("ide: p^r is not present in the fits tables,",5);
      xs_write("     we will use flat space p^r instead...",5);}
    if(nt>1 && column[3]==1)read_del=TRUE;
    if(nt>1 && column[3]==0){
      xs_write("ide: delay is not present in the fits tables,",5);
      xs_write("     we will use delay=-r*sin(phi)*sin(theta_o)...",5);}
    if(column[4]==1)read_lens=TRUE;
    if(column[4]==0){
      xs_write("ide: lensing is not present in the fits tables,",5);
      xs_write("     we will use lensing = 1.",5);}
  }else{
    read_alpha=FALSE;
    read_beta=FALSE;
    read_pr=FALSE;
    read_del=FALSE;
    read_lens=FALSE;
  }
/*We allocate memory for arrays only if we read the tables for the first
  time or if we have changed the dimensions of the arrays (in case the
  tables have changed) */
  if(ntable_old==-1 || nhorizon!=nhorizon_old || nincl!=nincl_old || 
     nr!=nr_old || nph!=nph_old){
/*  Firstly we have to free allocated memory for these arrays if we have
    already allocated it and the dimensions of the arrays have changed */
    if((ntable_old!=-1) && (nhorizon!=nhorizon_old || nincl!=nincl_old || 
       nr!=nr_old || nph!=nph_old)){
//    Free memory from tmp arrays...
      if(column_old[0]==1)free(alpha0);
      if(column_old[1]==1)free(beta0);
      if(column_old[2]==1)free(pr0);
      if(nt_old>1 && column_old[3]==1){free(delay);nt_old=1;}
      if(column_old[4]==1)free(lensing);
    }
//  Allocate memory for alpha0, beta0, pr, delay, lensing
    if(column[0]==1){
      alpha0=(float *) malloc((size_t) nr*nph*nhorizon*nincl*sizeof(float*));
      if(!alpha0){xs_write("failed to allocate memory for tmp arrays...",5);
        ntable=-1;return(1);}
    }     
    if(column[1]==1){
      beta0=(float *) malloc((size_t) nr*nph*nhorizon*nincl*sizeof(float*));
      if(!beta0){xs_write("failed to allocate memory for tmp arrays...",5);
        ntable=-1;return(1);}
    }
    if(column[2]==1){
      pr0=(float *) malloc((size_t) nr*nph*nhorizon*nincl*sizeof(float*));
      if(!pr0){xs_write("failed to allocate memory for tmp arrays...",5);
        ntable=-1;return(1);}
    }
    if(nt>1 && column[3]==1){
      delay=(float *) malloc((size_t) nr*nph*nhorizon*nincl*sizeof(float*));
      if(!delay){xs_write("failed to allocate memory for tmp arrays...",5);
        ntable=-1;return(1);}
    }
    if(column[4]==1){
      lensing=(float *) malloc((size_t) nr*nph*nhorizon*nincl*sizeof(float*));
      if(!lensing){xs_write("failed to allocate memory for tmp arrays...",5);
        ntable=-1;return(1);}
    }
  }else{
/*  if we have changed the tables and also we have changed nt from nt>1
    to nt=1 then we have to free memory for delay...*/
    if(ntable!=ntable_old && nt==1 && nt_old>1 && column_old[3]==1){
      free(delay);nt_old=1;}
/*  if we have already read the tables and we haven't changed them and nt=1
    at that time and we have changed nt to be > 1, we need to read delay
    because it hasn't been read yet (we also need to allocate memory for it) */
    if(nt>1 && nt_old==1 && column[3]==1){
      read_del=TRUE;
//    Allocate memory for delay...
      delay=(float *) malloc((size_t) nr*nph*nhorizon*nincl*sizeof(float*));
      if(!delay){xs_write("failed to allocate memory for tmp arrays...",5);
        ntable=-1;return(1);}
    }
    if(nt>1 && nt_old==1 && column[3]==0){
      xs_write("ide: delay is not present in the fits tables,",5);
      xs_write("     we will use delay=-r*sin(phi)*sin(theta_o)...",5);}
  }
// Let's read alpha0, beta0, pr0, delay and lensing...
  for(ihorizon=0;ihorizon<nhorizon;ihorizon++){
    for(iincl=0;iincl<nincl;iincl++){
      ffmrhd(fptr, 1, &hdutype, &status);
/*    to read the file only once we have to read in blocks (all columns
      from the extension are put to buffer together)
      let's find out how many rows are going to be read into the buffer */
      ffgrsz(fptr, &nrow, &status);
      nelems=nrow*nph;
      for(irow=0;irow<nr;irow+=nrow){
//      the last block to read may be smaller:
        if((nr-irow)<nrow)nelems=(nr-irow)*nph;
        icol=0;
        if(column[0]==1)icol++;
        if(read_alpha)
          ffgcv(fptr, TFLOAT, icol, irow+1, 1, nelems, &float_nulval, 
                &(alpha0[irow*nph+nr*nph*ihorizon+nr*nph*nhorizon*iincl]),
                &anynul, &status);
        if(column[1]==1)icol++;
        if(read_beta)          
          ffgcv(fptr, TFLOAT, icol, irow+1, 1, nelems, &float_nulval, 
                &(beta0[irow*nph+nr*nph*ihorizon+nr*nph*nhorizon*iincl]),
                &anynul, &status);
        if(column[2]==1)icol++;
        if(read_pr)          
          ffgcv(fptr, TFLOAT, icol, irow+1, 1, nelems, &float_nulval, 
                &(pr0[irow*nph+nr*nph*ihorizon+nr*nph*nhorizon*iincl]),
                &anynul, &status);
        if(column[3]==1)icol++;        
        if(read_del)
          ffgcv(fptr, TFLOAT, icol, irow+1, 1, nelems, &float_nulval, 
                &(delay[irow*nph+nr*nph*ihorizon+nr*nph*nhorizon*iincl]),
                &anynul, &status);
        if(column[4]==1)icol++;
        if(read_lens)
          ffgcv(fptr, TFLOAT, icol, irow+1, 1, nelems, &float_nulval, 
                &(lensing[irow*nph+nr*nph*ihorizon+nr*nph*nhorizon*iincl]),
                &anynul, &status);
      }
    }
  }
/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   Let us find NULL values in transfer functions
   SET nt=2 IN CALLING ROUTINE
  ------------------------------------------------------------ */
/*
  for(ihorizon=0;ihorizon<nhorizon;ihorizon++){
    for(iincl=0;iincl<nincl;iincl++){
      for(ir=0;ir<nr;ir++){
        for(iph=0;iph<nph;iph++){
          if(isnan(alpha0[iph+ir*nph+nr*nph*ihorizon+nr*nph*nhorizon*iincl]))
            fprintf(stdout, "alpha %14.6g %14.6g %d %d\n",
                    rmvec[ihorizon], thvec[iincl], ir, iph);
          if(isnan(beta0[iph+ir*nph+nr*nph*ihorizon+nr*nph*nhorizon*iincl]))
            fprintf(stdout, "beta %14.6g %14.6g %d %d\n",
                    rmvec[ihorizon], thvec[iincl], ir, iph);
          if(isnan(pr0[iph+ir*nph+nr*nph*ihorizon+nr*nph*nhorizon*iincl]))
            fprintf(stdout, "pr %14.6g %14.6g %d %d\n",
                    rmvec[ihorizon], thvec[iincl], ir, iph);
          if(isnan(delay[iph+ir*nph+nr*nph*ihorizon+nr*nph*nhorizon*iincl]))
            fprintf(stdout, "delay %14.6g %14.6g %d %d\n",
                    rmvec[ihorizon], thvec[iincl], ir, iph);
          if(isnan(lensing[iph+ir*nph+nr*nph*ihorizon+nr*nph*nhorizon*iincl]))
            fprintf(stdout, "lensing %14.6g %14.6g %d %d\n",
                    rmvec[ihorizon], thvec[iincl], ir, iph);
        }
      }
    }
  }
    stop */
//  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
/*******************************************************************************
  outfile1=fopen("alpha.txt","w");
  for(ir=0;ir<nr;ir++){
    for(iph=0;iph<nph;iph++)fprintf(outfile,"%14.6g ",
                                  alpha0[iph+nph*ir+nr*nph*9+nr*nph*nhorizon*9);
    fprintf(outfile1,"\n");
  }
  fclose(outfile1);
*******************************************************************************/
//The FITS file must always be closed before exiting the program.
  ffclos(fptr, &status);
  if(status){ffrprt(stderr, status);ntable=-1;return(1);}
/* Let's create array dc() for interpolated transfer functions:
  We allocate memory for arrays only if we read the tables for the first time or
  if we have changed the dimensions of the arrays (in case the tables have
  changed) */
  if(ntable_old==-1 || nr!=nr_old || nph!=nph_old){
/*  Firstly we have to free allocated memory for these arrays if we have
    already allocated it and the dimensions of the arrays have changed */
    if((ntable_old!=-1) && (nr!=nr_old || nph!=nph_old))free(dc);
//  Allocate memory for dc...
    dc=(double *) malloc((size_t) nr*nph*11*sizeof(double*));
    if(!dc){xs_write("failed to allocate memory for tmp arrays...",5);
      ntable=-1;return(1);}
  }
  for(i=0;i<5;i++)column_old[i]=column[i];
  nhorizon_old=nhorizon;
  nincl_old=nincl;
  nr_old=nr;
  nph_old=nph;
  ntable_old=ntable;
  if(nt>1)nt_old=nt;
  xs_write("     ...initializing finished",5);
}
//end of reading table fits file................................................
/* check if current value of thetaO is in the range covered by data tables
   if thetaO is "almost" 0 and fits tables are defined for higher values
   then thetaO but from close to 0 we change thetaO to the lowest value
   in the fits tables...(for convenience we put 0.1 equal to 0.0 ...) */
thetaOO=thetaO;
if(((thvec[0]-thetaO)<=0.1) && (thetaO<=0.1))thetaOO=(double) thvec[0];
//given am and thetaO, look up the stored tables that are needed now:
imin=0;
imax=nhorizon;
iam0=nhorizon/2;
while((imax-imin)>1){
  if(fabs(am)<=amvec[iam0-1])imin=iam0;
  else imax=iam0;
  iam0=(imin+imax)/2;
}
if(imax==nhorizon && fabs(am)<amvec[nhorizon-1])iam0=nhorizon;
imin=0;
imax=nincl;
ith0=nincl/2;
while((imax-imin)>1){
  if(thetaOO>=thvec[ith0-1])imin=ith0;
  else imax=ith0;
  ith0=(imin+imax)/2;
}
if(imax==nincl && thetaOO>thvec[nincl-1])ith0=nincl;
if(ith0==0)ith0=1;
//if there are tables only for one value of r_horizon
if(nhorizon==1)iam0=1;
/* of all tables, four sets are needed now;
1) a_low/th0_low, 2) a_low/th0_high, 3) a_high/th0_low, 4) a_high/th0_high:
interpolate data by bilinear interpolation to desired values of am and thetaO;
result is in dc(): */
if(nhorizon>1){
  ttmp=(r_plus-rmvec[iam0-1])/(rmvec[iam0]-rmvec[iam0-1]);
  ttmp1=1.-ttmp;
}
utmp=(thetaOO-thvec[ith0-1])/(thvec[ith0]-thvec[ith0-1]);
utmp1=1.-utmp;
for(iph=0;iph<nph;iph++){
  for(ir=1;ir<nr;ir++){      
// alpha
    if(column[0]==1){
      y1=(double) alpha0[iph+nph*ir+nr*nph*(iam0-1)+nr*nph*nhorizon*(ith0-1)];
      y4=(double) alpha0[iph+nph*ir+nr*nph*(iam0-1)+nr*nph*nhorizon*ith0];
      if(nhorizon==1)alpha=utmp1*y1+utmp*y4;
      else{
       y2=(double) alpha0[iph+nph*ir+nr*nph*iam0+nr*nph*nhorizon*(ith0-1)];
       y3=(double) alpha0[iph+nph*ir+nr*nph*iam0+nr*nph*nhorizon*ith0];
       alpha=utmp1*(ttmp1*y1+ttmp*y2)+utmp*(ttmp*y3+ttmp1*y4);
      }
      if(am<0.)alpha=-alpha;
    }
// beta
    if(column[1]==1){
      y1=(double) beta0[iph+nph*ir+nr*nph*(iam0-1)+nr*nph*nhorizon*(ith0-1)];
      y4=(double) beta0[iph+nph*ir+nr*nph*(iam0-1)+nr*nph*nhorizon*ith0];
      if(nhorizon==1)beta=utmp1*y1+utmp*y4;
      else{
        y2=(double) beta0[iph+nph*ir+nr*nph*iam0+nr*nph*nhorizon*(ith0-1)];
        y3=(double) beta0[iph+nph*ir+nr*nph*iam0+nr*nph*nhorizon*ith0];
        beta=utmp1*(ttmp1*y1+ttmp*y2)+utmp*(ttmp*y3+ttmp1*y4);
      }
    }
// pr
    if(column[2]==1){
      y1=(double) pr0[iph+nph*ir+nr*nph*(iam0-1)+nr*nph*nhorizon*(ith0-1)];
      y4=(double) pr0[iph+nph*ir+nr*nph*(iam0-1)+nr*nph*nhorizon*ith0];
      if(nhorizon==1)pr=utmp1*y1+utmp*y4;
      else{
        y2=(double) pr0[iph+nph*ir+nr*nph*iam0+nr*nph*nhorizon*(ith0-1)];
        y3=(double) pr0[iph+nph*ir+nr*nph*iam0+nr*nph*nhorizon*ith0];
        pr=utmp1*(ttmp1*y1+ttmp*y2)+utmp*(ttmp*y3+ttmp1*y4);
      }
    }
// delay (only if nt>1 and column[3]=1)
    if(nt>1 && column[3]==1){
      y1=(double) delay[iph+nph*ir+nr*nph*(iam0-1)+nr*nph*nhorizon*(ith0-1)];
      y4=(double) delay[iph+nph*ir+nr*nph*(iam0-1)+nr*nph*nhorizon*ith0];
      if(nhorizon==1)dc[ir+nr*iph]=utmp1*y1+utmp*y4;
      else{
       y2=(double) delay[iph+nph*ir+nr*nph*iam0+nr*nph*nhorizon*(ith0-1)];
       y3=(double) delay[iph+nph*ir+nr*nph*iam0+nr*nph*nhorizon*ith0];
       dc[ir+nr*iph]=utmp1*(ttmp1*y1+ttmp*y2)+utmp*(ttmp*y3+ttmp1*y4);
      }
    }
// lensing
    if(column[4]==1){
      y1=(double) lensing[iph+nph*ir+nr*nph*(iam0-1)+nr*nph*nhorizon*(ith0-1)];
      y4=(double) lensing[iph+nph*ir+nr*nph*(iam0-1)+nr*nph*nhorizon*ith0];
      if(nhorizon==1)dc[ir+nr*iph+nr*nph*3]=utmp1*y1+utmp*y4;
      else{
        y2=(double) lensing[iph+nph*ir+nr*nph*iam0+nr*nph*nhorizon*(ith0-1)];
        y3=(double) lensing[iph+nph*ir+nr*nph*iam0+nr*nph*nhorizon*ith0];
        dc[ir+nr*iph+nr*nph*3]=utmp1*(ttmp1*y1+ttmp*y2)+utmp*(ttmp*y3+ttmp1*y4);
      }
    }

    rnow=r_vec[ir]+r_plus;
    rnow2=rnow*rnow;
//    delta=r_vec[ir]*r_vec[ir]-2.*r_vec[ir]+am2;
    delta=rnow2-2.*rnow+am2;
    l=alpha*sin(thetaOO/180.*PI);
    q2=beta*beta+cos(thetaOO/180.*PI)*cos(thetaOO/180.*PI)*(alpha*alpha-am2);
    p_r=rnow2/delta*pr;
    if(q2<0.)q2=0.;
    p_th=-sqrt(q2);
//   Keplerian velocity U^\mu
    if(rnow>rms){
      sqr=fabs(rnow2-3.*rnow+2.*am*sqrt(rnow));
      Ut=(rnow2+am*sqrt(rnow))/(rnow*sqrt(sqr));
      Ur=0.;
      Uph=1./sqrt(rnow*sqr);
    }else{
      sqr_ms=fabs(rms2-3.*rms+2.*am*sqrt(rms));
      Ems=(rms2-2.*rms+am*sqrt(rms))/(rms*sqrt(sqr_ms));
      Lms=(rms2+am2-2.*am*sqrt(rms))/sqrt(rms*sqr_ms);
      Ut=(((rnow*(rnow2+am2)+2.*am2)*Ems-2.*am*Lms)/rnow/delta);
      Ur=(-sqrt(fabs((rnow*(rnow2+am2)+2.*am2)*Ems*Ems-4.*am*Ems*Lms
         -(rnow-2.)*Lms*Lms-rnow*delta))/rnow/sqrt(rnow));
      Uph=((2.*am*Ems+(rnow-2.)*Lms)/rnow/delta);
    }
//    if(am<0.)Uph=-Uph;
//redshift factor:
    gfac=1./(Ut-p_r*Ur-l*Uph);
//cosine
    cosine=-gfac*p_th/rnow;    
    if(cosine>1.)cosine=1.;

    AA=(rnow2+am2)*(rnow2+am2)-am2*delta;
    pt=(AA-2.*am*rnow*l)/(delta*rnow2);
    pph=((rnow-2.)*l+2.*am)/(delta*rnow);
//change in azimuthal emission angle (phiphoton):
    phiphoton=atan2(gfac*delta/rnow*(-pt*Uph+pph*Ut),gfac*pr-Ur);
//change in polarization angle:
    pth=p_th/rnow2;
    ft=-cosine*(gfac*pt-Ut);
    fr=-cosine*(gfac*pr-Ur);
    fth=-1./rnow-cosine*gfac*pth;
    fph=-cosine*(gfac*pph-Uph);
    kappa1=rnow*(pth*(am*ft-(rnow2+am2)*fph)-(am*pt-(rnow2+am2)*pph)*fth);
    kappa2=rnow*(pr*(-ft+am*fph)+(pt-am*pph)*fr);
    polarization=atan2((alpha-am*sin(thetaOO/180.*PI))*kappa2+beta*kappa1,
                       -(alpha-am*sin(thetaOO/180.*PI))*kappa1+beta*kappa2);
    dc[ir+nr*iph+nr*nph]=gfac;
    dc[ir+nr*iph+nr*nph*2]=cosine;
    dc[ir+nr*iph+nr*nph*4]=polarization;
    dc[ir+nr*iph+nr*nph*5]=phiphoton;
    dc[ir+nr*iph+nr*nph*9]=alpha;
    dc[ir+nr*iph+nr*nph*10]=beta;
  }
}
//derivatives of gfactor
for(iph=0;iph<nph;iph++){
  for(ir=1;ir<nr;ir++){
//    if(column[1]==1){
    if(ir>0){
      dg1=dc[ir+nr*iph+nr*nph]-dc[ir-1+nr*iph+nr*nph];
      if(iph>0){
        dg5=dc[ir-1+nr*iph+nr*nph]-dc[ir-1+nr*(iph-1)+nr*nph];
        if(iph<nph-1)dg6=dc[ir-1+nr*(iph+1)+nr*nph]-dc[ir-1+nr*iph+nr*nph];
        else dg6=dc[ir-1+nr*(nph-1)+nr*nph]-dc[ir-1+nr*(nph-2)+nr*nph];
      }else{
       dg5=dc[ir-1+nr+nr*nph]-dc[ir-1+nr*nph];
       dg6=dg5;
      }
      if(ir<nr-1){
        kr=((r_vec[ir+1]-r_vec[ir])/(r_vec[ir]-r_vec[ir-1]));
        deltar=r_vec[ir+1]-r_vec[ir-1];
        dg2=dc[ir+1+nr*iph+nr*nph]-dc[ir+nr*iph+nr*nph];
        if(iph>0){
          dg7=dc[ir+1+nr*iph+nr*nph]-dc[ir+1+nr*(iph-1)+nr*nph];
          if(iph<nph-1)dg8=dc[ir+1+nr*(iph+1)+nr*nph]-dc[ir+1+nr*iph+nr*nph];
          else dg8=dc[ir+1+nr*(nph-1)+nr*nph]-dc[ir+1+nr*(nph-2)+nr*nph];
        }else{
          dg7=dc[ir+1+nr+nr*nph]-dc[ir+1+nr*nph];
          dg8=dg7;
        }
      }else{
//      ir=nr-1
        kr=1.;
        deltar=2.*(r_vec[nr-1]-r_vec[nr-2]);
        dg2=dc[nr-1+nr*iph+nr*nph]-dc[nr-2+nr*iph+nr*nph];
        if(iph>0){
          dg7=dc[nr-1+nr*iph+nr*nph]-dc[nr-1+nr*(iph-1)+nr*nph];
          if(iph<nph-1)dg8=dc[nr-1+nr*(iph+1)+nr*nph]-dc[nr-1+nr*iph+nr*nph];
          else dg8=dc[nr-1+nr*(nph-1)+nr*nph]-dc[nr-1+nr*(nph-2)+nr*nph];
        }else{
          dg7=dc[nr-1+nr+nr*nph]-dc[nr-1+nr*nph];
          dg8=dg7;
        }
      }
    }else{
//    ir=1;
      kr=1.;
      deltar=2.*(r_vec[1]-r_vec[0]);
      dg1=dc[1+nr*iph+nr*nph]-dc[nr*iph+nr*nph];
      dg2=dg1;
      if(iph>0){
        dg5=dc[nr*iph+nr*nph]-dc[nr*(iph-1)+nr*nph];
        dg7=dc[1+nr*iph+nr*nph]-dc[1+nr*(iph-1)+nr*nph];
        if(iph<nph-1){
          dg6=dc[nr*(iph+1)+nr*nph]-dc[nr*iph+nr*nph];
          dg8=dc[1+nr*(iph+1)+nr*nph]-dc[1+nr*iph+nr*nph];
        }else{
          dg6=dc[nr*(nph-1)+nr*nph]-dc[nr*(nph-2)+nr*nph];
          dg8=dc[1+nr*(nph-1)+nr*nph]-dc[1+nr*(nph-2)+nr*nph];
        }
      }else{
        dg5=dc[nr+nr*nph]-dc[nr*nph];
        dg6=dg5;
        dg7=dc[1+nr+nr*nph]-dc[1+nr*nph];
        dg8=dg7;
      }
    }
    if(iph>0){
      dg3=dc[ir+nr*iph+nr*nph]-dc[ir+nr*(iph-1)+nr*nph];
      if(iph<nph-1){
        kphi=(phi_vec[iph+1]-phi_vec[iph])/(phi_vec[iph]-phi_vec[iph-1]);
        deltaphi=phi_vec[iph+1]-phi_vec[iph-1];
        dg4=dc[ir+nr*(iph+1)+nr*nph]-dc[ir+nr*iph+nr*nph];
      }else{
        kphi=1.;
        deltaphi=2.*(phi_vec[nph-1]-phi_vec[nph-2]);
        dg4=dc[ir+nr*(nph-1)+nr*nph]-dc[ir+nr*(nph-2)+nr*nph];
      }
    }else{
      kphi=1.;
      deltaphi=2.*(phi_vec[1]-phi_vec[0]);
      dg3=dc[ir+nr+nr*nph]-dc[ir+nr*nph];
      dg4=dg3;
    }
//the first derivative with respect to r
    dc[ir+nr*iph+nr*nph*6]=(kr*dg1+dg2/kr)/deltar;
//the first derivative with respect to phi
    dc[ir+nr*iph+nr*nph*7]=(kphi*dg3+dg4/kphi)/deltaphi;
//the second derivative with respect to r and phi
    dc[ir+nr*iph+nr*nph*8]=(kr*kphi*(dg3-dg5)+kr*(dg4-dg6)/kphi+
                           kphi*(dg7-dg3)/kr+(dg8-dg4)/kr/kphi)/deltar/deltaphi;
//    }
  }
}
/*******************************************************************************
//Let's write the functions to the file for TESTing purposes...
outfile1=fopen("test1.dat","w");
for(iph=0;iph<nph;iph++){
  for(ir=1;ir<nr;ir++)
   fprintf(outfile1,"%14.6g %14.6g %14.6g %14.6g %14.6g %14.6g %14.6g %14.6g\n",
     r_vec[ir]+r_plus, phi_vec[ir], dc[ir+nr*iph], dc[ir+nr*iph+nr*nph],
     dc[ir+nr*iph+nr*nph*2], dc[ir+nr*iph+nr*nph*3], dc[ir+nr*iph+nr*nph*4],
     dc[ir+nr*iph+nr*nph*5]);
}
 fclose(outfile1);
 stop;
*******************************************************************************/
return(0);
}
/*******************************************************************************
*******************************************************************************/
void element(int iloop, const int ne_loc, const int ne, const int nt, 
             double *fce, double *qce, double *uce, double *vce, 
             void (*emissivity)()){

void bcuint(double y[], double y1[], double y2[], double y12[], double x1l,
            double x1u, double x2l, double x2u, double x1, double x2, 
            double *ansy);

double *ear_loc,far_loc[(ne_loc+1)*nt];
double qar_loc[ne_loc*nt], uar_loc[ne_loc*nt], var_loc[ne_loc*nt];
double intfunc1, intfunc2;
double delay, gfac, cosine, transf, alpha, beta;
double phiphoton, polarization, cos2pol, sin2pol;
double rnow, pnow, rnow_rplus, pnow_kerr;
double g[4], g1[4], g2[4], g12[4];
double dfc, drad, dphi, radratio, radrat, facnorm, fact, dener, ddener;
double xlow, xhigh, dshift;
double ttmp, ttmp1, utmp, utmp1, y1, y2, y3, y4, dr2;
double rnow2, delta, ll, q2, Rsgn, pr, p_r, p_th, sqr, Ut, Ur, Uph, sqr_ms, Ems, 
       Lms, AA, pt, pph;//, pth, ft, fr, fth, fph, kappa1, kappa2;
int    imin, imax, ir0, iphi0, iiloop, ir, iphi, it, k, l, iimax, ii;
int    ilow, ihigh, ishiftt;
      
// loop over disc
drad=(rout-rin)/nrad;
dphi=(phi1-phi0)/nphi;
if(rdivision==1)radratio=pow((rout+r_plus)/(rin+r_plus),1./nrad);
iloop=(iloop-1)*nloops+1;
iiloop=iloop;
while(iiloop<(iloop+nloops) && iiloop<=(nrad*nphi)){
  ir=(iiloop-1)/nphi+1;
  iphi=iiloop-((iiloop-1)/nphi)*nphi;
  if(ir>nrad || iphi>nphi)return;
  if(rdivision==0)rnow=rin+r_plus+drad*(ir-0.5);
  if(rdivision==1){
    radrat=pow(radratio,ir-1.);
    rnow=0.5*(rin+r_plus)*(radratio+1.)*radrat;
    drad=(rin+r_plus)*(radratio-1.)*radrat;
  }
  pnow=phi0+dphi*(iphi-0.5);
  facnorm=drad*dphi*rnow;
// given (rnow,pnow), find corresponding indices in r_vec, phi_vec():
// ... the tables are made for rnow-rplus
  rnow_rplus=rnow-r_plus;
//  if(rnow_rplus<0.)rnow_rplus=0.;
//Let's check if rnow is not in the first bin above the horizon where the 
//interpolation does not work due to weird values of gfac, and others at the 
//horizon (NaN, inf, 0.)...
  if((r_vec[0] == 0.) && (rnow_rplus < r_vec[1])) {
    iiloop+=1;
    continue;
  }
  if(rnow_rplus<=r_vec[nr-1]){
    imin=0;
    imax=nr;
    ir0=nr/2;
    while((imax-imin)>1){
      if(rnow_rplus>=r_vec[ir0-1])imin=ir0;
      else imax=ir0;
      ir0=(imin+imax)/2;
    }
    if(imax==nr && rnow_rplus>r_vec[nr-1])ir0=nr;
//  ... we must convert phi from Boyer-Lindquist coordinate to the Kerr coordinate
    if(fabs(am)==1.)pnow_kerr = pnow + am/(rnow-1.);
    else pnow_kerr = pnow - am/(r_plus-r_minus)*log(rnow_rplus/(rnow-r_minus));
//  phi in the data tables is >= 0 and <= 2*PI
//  ...let's modify pnow_kerr so that it lies in this range
    if(am<0.)pnow_kerr=PI-pnow_kerr;
    pnow_kerr = fmod(pnow_kerr,PI2);
    if(pnow_kerr<0.)pnow_kerr+=PI2;
    imin=0;
    imax=nph;
    iphi0=nph/2;
    while((imax-imin)>1){
      if(pnow_kerr>=phi_vec[iphi0-1])imin=iphi0;
      else imax=iphi0;
      iphi0=(imin+imax)/2;
    }
    if(imax==nph && pnow_kerr>phi_vec[nph-1])iphi0=nph;
//  let's interpolate the values in dc() table to desired (rnow,pnow):
    ttmp=(rnow_rplus-r_vec[ir0-1])/(r_vec[ir0]-r_vec[ir0-1]);
    utmp=(pnow_kerr-phi_vec[iphi0-1])/(phi_vec[iphi0]-phi_vec[iphi0-1]);
    ttmp1=1.-ttmp;
    utmp1=1.-utmp;
//  time delay factor (not used for a stationary problem):
    if((nt>1) && (column[3]==1)){
      y1=dc[ir0-1+nr*(iphi0-1)];
      y2=dc[ir0+nr*(iphi0-1)];
      y3=dc[ir0+nr*iphi0];
      y4=dc[ir0-1+nr*iphi0];
      delay=utmp1*(ttmp1*y1+ttmp*y2)+utmp*(ttmp*y3+ttmp1*y4);
    }
    if((nt>1) && (column[3]==0))delay=-rnow*sin(pnow)*sin(thetaO/180.*PI);
//  redshift factor:
//  if(column[1]==1){
/*   //  bilinear interpolation
      y1=dc[ir0-1+nr*(iphi0-1)+nr*nph];
      y2=dc[ir0+nr*(iphi0-1)+nr*nph];
      y3=dc[ir0+nr*iphi0+nr*nph];
      y4=dc[ir0-1+nr*iphi0+nr*nph];
      gfac=utmp1*(ttmp1*y1+ttmp*y2)+utmp*(ttmp*y3+ttmp1*y4)
*/  
//    bicubic interpolation
      g[0]=dc[ir0-1+nr*(iphi0-1)+nr*nph];
      g[1]=dc[ir0+nr*(iphi0-1)+nr*nph];
      g[2]=dc[ir0+nr*iphi0+nr*nph];
      g[3]=dc[ir0-1+nr*iphi0+nr*nph];
      g1[0]=dc[ir0-1+nr*(iphi0-1)+nr*nph*6];
      g1[1]=dc[ir0+nr*(iphi0-1)+nr*nph*6];
      g1[2]=dc[ir0+nr*iphi0+nr*nph*6];
      g1[3]=dc[ir0-1+nr*iphi0+nr*nph*6];
      g2[0]=dc[ir0-1+nr*(iphi0-1)+nr*nph*7];
      g2[1]=dc[ir0+nr*(iphi0-1)+nr*nph*7];
      g2[2]=dc[ir0+nr*iphi0+nr*nph*7];
      g2[3]=dc[ir0-1+nr*iphi0+nr*nph*7];
      g12[0]=dc[ir0-1+nr*(iphi0-1)+nr*nph*8];
      g12[1]=dc[ir0+nr*(iphi0-1)+nr*nph*8];
      g12[2]=dc[ir0+nr*iphi0+nr*nph*8];
      g12[3]=dc[ir0-1+nr*iphi0+nr*nph*8];
      bcuint(g, g1, g2, g12, r_vec[ir0-1], r_vec[ir0], phi_vec[iphi0-1],
                  phi_vec[iphi0], rnow_rplus, pnow_kerr, &gfac);
//  }else gfac=sqrt(rnow-1.)/(sqrt(rnow)-cos(pnow)*sin(thetaO/180.*PI));
    if(gfac < 0.)gfac=0.;
    gfac=zzshift*gfac;
//  cosine
    y1=dc[ir0-1+nr*(iphi0-1)+nr*nph*2];
    y2=dc[ir0+nr*(iphi0-1)+nr*nph*2];
    y3=dc[ir0+nr*iphi0+nr*nph*2];
    y4=dc[ir0-1+nr*iphi0+nr*nph*2];
    cosine=utmp1*(ttmp1*y1+ttmp*y2)+utmp*(ttmp*y3+ttmp1*y4);
//    cosine=cos(thetaO/180.*PI)*sqrt(rnow-1.)/
//           (sqrt(rnow)-cos(pnow)*sin(thetaO/180.*PI))
//  light lensing amplification:
    if(column[2]==1 && column[4]==1){
      y1=dc[ir0-1+nr*(iphi0-1)+nr*nph*2]*dc[ir0-1+nr*(iphi0-1)+nr*nph*3];
      y2=dc[ir0+nr*(iphi0-1)+nr*nph*2]*dc[ir0+nr*(iphi0-1)+nr*nph*3];
      y3=dc[ir0+nr*iphi0+nr*nph*2]*dc[ir0+nr*iphi0+nr*nph*3];
      y4=dc[ir0-1+nr*iphi0+nr*nph*2]*dc[ir0-1+nr*iphi0+nr*nph*3];
      transf=utmp1*(ttmp1*y1+ttmp*y2)+utmp*(ttmp*y3+ttmp1*y4);
    }else transf=cosine;
/*  
    if(column[4]==1){
      y1=dc[ir0-1+nr*(iphi0-1)+nr*nph*3];
      y2=dc[ir0+nr*(iphi0-1)+nr*nph*3];
      y3=dc[ir0+nr*iphi0+nr*nph*3];
      y4=dc[ir0-1+nr*iphi0+nr*nph*3];
      transf=(utmp1*(ttmp1*y1+ttmp*y2)+utmp*(ttmp*y3+ttmp1*y4))*cosine
    }else transf=cosine;
*/  
//  change in polarization angle:
    if(polar){
      y1=dc[ir0-1+nr*(iphi0-1)+nr*nph*4];
      y2=dc[ir0+nr*(iphi0-1)+nr*nph*4];
      y3=dc[ir0+nr*iphi0+nr*nph*4];
      y4=dc[ir0-1+nr*iphi0+nr*nph*4];
      if(fabs(y2-y1)>PI){
       if(y1<0.)y1+=PI2;
       if(y2<0.)y2+=PI2;
      }
      if(fabs(y4-y3)>PI){
       if(y3<0.)y3+=PI2;
       if(y4<0.)y4+=PI2;
      }
      y1=ttmp1*y1+ttmp*y2;
      y2=ttmp*y3+ttmp1*y4;
      if(fabs(y2-y1)>PI){
       if(y1<0.)y1+=PI2;
       if(y2<0.)y2+=PI2;
      }
      polarization=utmp1*y1+utmp*y2;
      if(polarization>PI)polarization-=PI2;
      cos2pol=cos(2.*polarization);
      sin2pol=sin(2.*polarization);
    }
/*  
    if(polar && column[4]==0){
     polarization=-atan2(cos(thetaO/180.*PI)*sin(pnow),
                         -cos(pnow)+sqrt(rnow)*sin(thetaO/180.*PI));
     cos2pol=cos(2.*polarization);
     sin2pol=sin(2.*polarization);
    }
*/  
//   change in azimuthal emission angle (phiphoton):
    y1=dc[ir0-1+nr*(iphi0-1)+nr*nph*5];
    y2=dc[ir0+nr*(iphi0-1)+nr*nph*5];
    y3=dc[ir0+nr*iphi0+nr*nph*5];
    y4=dc[ir0-1+nr*iphi0+nr*nph*5];
    if(fabs(y2-y1)>PI){
     if(y1<0.)y1+=PI2;
     if(y2<0.)y2+=PI2;
    }
    if(fabs(y4-y3)>PI){
     if(y3<0.)y3+=PI2;
     if(y4<0.)y4+=PI2;
    }
    y1=ttmp1*y1+ttmp*y2;
    y2=ttmp*y3+ttmp1*y4;
    if(fabs(y2-y1)>PI){
     if(y1<0.)y1+=PI2;
     if(y2<0.)y2+=PI2;
    }
    phiphoton=utmp1*y1+utmp*y2;
    if(phiphoton>PI)phiphoton-=PI2;
//        if(polar)phiphoton=0.
//  alpha
    y1=dc[ir0-1+nr*(iphi0-1)+nr*nph*9];
    y2=dc[ir0+nr*(iphi0-1)+nr*nph*9];
    y3=dc[ir0+nr*iphi0+nr*nph*9];
    y4=dc[ir0-1+nr*iphi0+nr*nph*9];
    alpha=utmp1*(ttmp1*y1+ttmp*y2)+utmp*(ttmp*y3+ttmp1*y4);
//  beta
    y1=dc[ir0-1+nr*(iphi0-1)+nr*nph*10];
    y2=dc[ir0+nr*(iphi0-1)+nr*nph*10];
    y3=dc[ir0+nr*iphi0+nr*nph*10];
    y4=dc[ir0-1+nr*iphi0+nr*nph*10];
    beta=utmp1*(ttmp1*y1+ttmp*y2)+utmp*(ttmp*y3+ttmp1*y4);
  }else{//compute everything for geodesics approximated by straight lines
//    alpha=rnow*cos(pnow);
//    beta=-rnow*sin(pnow)*cos(thetaO/180.*PI);
//    delay=-rnow*sin(pnow)*sin(thetaO/180.*PI);
//  ... we must convert phi from Boyer-Lindquist coordinate to the Kerr coordinate
    if(fabs(am)==1.)pnow_kerr = pnow + am/(rnow-1.);
    else pnow_kerr = pnow - am/(r_plus-r_minus)*log(rnow_rplus/(rnow-r_minus));
//  phi in the data tables is >= 0 and <= 2*PI
//  ...let's modify pnow_kerr so that it lies in this range
    if(am<0.)pnow_kerr=PI-pnow_kerr;
    pnow_kerr = fmod(pnow_kerr,PI2);
    if(pnow_kerr<0.)pnow_kerr+=PI2;
    imin=0;
    imax=nph;
    iphi0=nph/2;
    while((imax-imin)>1){
      if(pnow_kerr>=phi_vec[iphi0-1])imin=iphi0;
      else imax=iphi0;
      iphi0=(imin+imax)/2;
    }
    if(imax==nph && pnow_kerr>phi_vec[nph-1])iphi0=nph;
    utmp=(pnow_kerr-phi_vec[iphi0-1])/(phi_vec[iphi0]-phi_vec[iphi0-1]);
    utmp1=1.-utmp;
// alpha
    y1=dc[nr-1+nr*(iphi0-1)+nr*nph*9];
    y2=dc[nr-1+nr*iphi0+nr*nph*9];
    alpha=utmp1*y1+utmp*y2;
    alpha = rnow*cos(pnow) + alpha - (r_vec[nr-1]+r_plus)*cos(pnow);
// beta
    y1=dc[nr-1+nr*(iphi0-1)+nr*nph*10];
    y2=dc[nr-1+nr*iphi0+nr*nph*10];
    beta=utmp1*y1+utmp*y2;
    beta = -rnow*sin(pnow)*cos(thetaO/180.*PI) + beta + 
           (r_vec[nr-1]+r_plus)*sin(pnow)*cos(thetaO/180.*PI);
//  time delay factor (not used for a stationary problem):
    if((nt>1) && (column[3]==1)){
      y1=dc[nr-1+nr*(iphi0-1)];
      y2=dc[nr-1+nr*iphi0];
      delay=utmp1*y1+utmp*y2;
    }
    if((nt>1) && (column[3]==0))delay=-rnow*sin(pnow)*sin(thetaO/180.*PI);
    if(nt>1) delay -= ddel;
    delay = -rnow*sin(pnow)*sin(thetaO/180.*PI) +
            delay + (r_vec[nr-1]+r_plus)*sin(pnow)*sin(thetaO/180.*PI);
    Rsgn=1.;
    if(beta>0.)Rsgn = -1.;
    rnow2=rnow*rnow;
    delta=rnow2-2.*rnow+am2;
    ll=alpha*sin(thetaO/180.*PI);
    q2=beta*beta+cos(thetaO/180.*PI)*cos(thetaO/180.*PI)*(alpha*alpha-am2);
    pr=Rsgn*sqrt((rnow2+am2-am*ll)*(rnow2+am2-am*ll)-delta*((ll-am)*(ll-am)+q2))
       /rnow2;
    p_r=rnow2/delta*pr;
    p_th=-sqrt(fabs(q2));
//  Keplerian velocity U^\mu
    if(rnow>rms){
      sqr=fabs(rnow2-3.*rnow+2.*am*sqrt(rnow));
      Ut=(rnow2+am*sqrt(rnow))/(rnow*sqrt(sqr));
      Ur=0.;
      Uph=1./sqrt(rnow*sqr);
    }else{
      sqr_ms=fabs(rms2-3.*rms+2.*am*sqrt(rms));
      Ems=(rms2-2.*rms+am*sqrt(rms))/(rms*sqrt(sqr_ms));
      Lms=(rms2+am2-2.*am*sqrt(rms))/sqrt(rms*sqr_ms);
      Ut=(((rnow*(rnow2+am2)+2.*am2)*Ems-2.*am*Lms)/rnow/delta);
      Ur=(-sqrt(fabs((rnow*(rnow2+am2)+2.*am2)*Ems*Ems-4.*am*Ems*Lms
         -(rnow-2.)*Lms*Lms-rnow*delta))/rnow/sqrt(rnow));
      Uph=((2.*am*Ems+(rnow-2.)*Lms)/rnow/delta);
    }
//  if(am<0.)Uph=-Uph;
//redshift factor:
    gfac=1./(Ut-p_r*Ur-ll*Uph);
//cosine
    cosine=-gfac*p_th/rnow;

    AA=(rnow2+am2)*(rnow2+am2)-am2*delta;
    pt=(AA-2.*am*rnow*ll)/(delta*rnow2);
    pph=((rnow-2.)*ll+2.*am)/(delta*rnow);
//change in azimuthal emission angle (phiphoton):
    phiphoton=atan2(gfac*delta/rnow*(-pt*Uph+pph*Ut),gfac*pr-Ur);
/*
//change in polarization angle:
//THIS DOES NOT WORK PROPERLY
    pth=p_th/rnow2;
    ft=-cosine*(gfac*pt-Ut);
    fr=-cosine*(gfac*pr-Ur);
    fth=-1./rnow-cosine*gfac*pth;
    fph=-cosine*(gfac*pph-Uph);
    kappa1=rnow*(pth*(am*ft-(rnow2+am2)*fph)-(am*pt-(rnow2+am2)*pph)*fth);
    kappa2=rnow*(pr*(-ft+am*fph)+(pt-am*pph)*fr);
    polarization=atan2((alpha-am*sin(thetaO/180.*PI))*kappa2+beta*kappa1,
                       -(alpha-am*sin(thetaO/180.*PI))*kappa1+beta*kappa2);
*/
    gfac=zzshift*gfac;
    transf=cosine;
  }
//  fprintf(stdout,"%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n", 
//          rnow, alpha, beta, delay, gfac, cosine, transf, phiphoton, polarization);
//If we are at position that is obscured go to next step in iiloop
  if(rcloud2>0.){
    dr2=(alpha-alphac)*(alpha-alphac)+(beta-betac)*(beta-betac);
    if(dr2<=rcloud2 && rcloud > 0.){iiloop+=1;continue;}
    else if(dr2>rcloud2 && rcloud < 0.){iiloop+=1;continue;}
  }
/*******************************************************************************
  Let's write the functions to the file for TESTing purposes...
  fprintf(outfile,"%14.6g %14.6g %14.6g %14.6g %14.6g %14.6g %14.6g %14.6g\n",
          rnow,pnow_kerr,delay,gfac,cosine,transf,polarization,phiphoton)
*******************************************************************************/

/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 in the next line there is gfac^2 because we do not deal with
 intensities but with photon density flux
 - there should be actually one more gfac here because we later on forget
   one gfac when rebinning (i.e. we are dealing with integrated flux actually)
 - there should be actually one gfac less here because we integrate in
   coordinates r,phi (dS=r*dr*dphi) and flux is set in local frame...
   (dS_loc/dS=gfac^{-1}, dS_loc is local area)
 => the last two remarks cancel each other...  ;-)
 => see also PhD thesis (Dovciak 2004) and eqs. 2.1-2.11
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
  if(!observed_flux) fact=facnorm*transf*gfac*gfac;
  else fact=facnorm*transf*gfac;

//let's call emissivity subroutine so that we fill far_loc (qar_loc, uar_loc,
//var_loc) array with local photon flux density
  emissivity(&ear_loc, ne_loc, nt, far_loc, qar_loc, uar_loc, var_loc,
             rnow, pnow, cosine, phiphoton, alpha, beta, delay, gfac);
  if(observed_flux){
//  fill fc(:) array (observed spectrum) using the values in far_loc(:)
//  array (local spectrum):
    for(k = 0; k < ne; k++){
/*    before we can add contributions of local photon flux density we have to
      decide into which 'time' element of far(ne,time) it falls */
//    here we loop over local time...
      for(it=0;it<nt;it++){
        if(nt!=1 && far_loc[ne_loc+(ne_loc+1)*it]==0.)continue;
        if(nt==1)ishiftt=0;
        else{
          ishiftt = it + (int) trunc((delay-ddel)/dt);
          dshift = ((delay-ddel)/dt) - (int) trunc((delay-ddel)/dt);
          if(dshift>0.)ishiftt+=1;
          else dshift+=1;
/*        what to do, if shifted time falls outside time span covered by
          the array: */
          if(periodic==1){
/*          the case of a periodic spectrum (e.g. an orbiting spot); in this
            case emissivity subroutine must be constructed in such a way that
            it covers a complete period (or a multiple of it): */
            while(ishiftt<0)ishiftt+=nt;
            while(ishiftt>=nt)ishiftt-=nt;
          }
        }
/*      contribution from the current point in the disc is added to final
        spectrum:
        the case of non-periodic spectrum (e.g. a decaying flare):
        do nothing if spectral points fall outside the array (time)
        boundaries: */
        if((ishiftt<0) || (ishiftt>=nt))continue;
/*      if we are calculating the Stokes parameters as well
        we have to integrate all 4 functions... */
        iimax=1;
        if(polar)iimax=4;
        for(ii=1;ii<=iimax;ii++){
          if(ii==1){
            intfunc1=far_loc[k+(ne_loc+1)*it];
          }else if(ii==2){
            intfunc1=qar_loc[k+ne_loc*it]*cos2pol-
                     uar_loc[k+ne_loc*it]*sin2pol;
          }else if(ii==3){
            intfunc1=qar_loc[k+ne_loc*it]*sin2pol+
                     uar_loc[k+ne_loc*it]*cos2pol;
          }else if(ii==4){
              intfunc1=var_loc[k+ne_loc*it];
          }
          intfunc1 *= fact;
          if(nt==1){
            if(ii==1)fce[k]+=intfunc1;
            else if(ii==2)qce[k]+=intfunc1;
            else if(ii==3)uce[k]+=intfunc1;
            else if(ii==4)vce[k]+=intfunc1;
          }else{
            if(ii==1)fce[k+ne*ishiftt]+=intfunc1*dshift;
            else if(ii==2)qce[k+ne*ishiftt]+=intfunc1*dshift;
            else if(ii==3)uce[k+ne*ishiftt]+=intfunc1*dshift;
            else if(ii==4)vce[k+ne*ishiftt]+=intfunc1*dshift;
            if(ishiftt>0){
              if(ii==1)fce[k+ne*(ishiftt-1)]+=intfunc1*(1-dshift);
              else if(ii==2)qce[k+ne*(ishiftt-1)]+=intfunc1*(1-dshift);
              else if(ii==3)uce[k+ne*(ishiftt-1)]+=intfunc1*(1-dshift);
              else if(ii==4)vce[k+ne*(ishiftt-1)]+=intfunc1*(1-dshift);
            }else{
              if(periodic==1){
                if(ii==1)fce[k+ne*(nt-1)]+=intfunc1*(1-dshift);
                else if(ii==2)qce[k+ne*(nt-1)]+=intfunc1*(1-dshift);
                else if(ii==3)uce[k+ne*(nt-1)]+=intfunc1*(1-dshift);
                else if(ii==4)vce[k+ne*(nt-1)]+=intfunc1*(1-dshift);
              }
            }
          }
        }//next function (far, qar, uar, var)
      }//next it
    }//next energy bin
  }else{
    if(edivision==0)dener=(ear_loc[ne_loc-1]-ear_loc[0])/(ne_loc-1);
    else if(edivision==1)
      dener=(log10(ear_loc[ne_loc-1])-
             log10(ear_loc[0]))/(ne_loc-1);
//  fill fc(:) array (observed spectrum) using the values in far_loc(:)
//  array (local spectrum):
    for(k = 0; k < ne; k++){
/*    let's find lower and upper index of local energy; energy intervals with 
      lower edge between this lower and upper value are shifted and fall into 
      the energy bin eard(k)-eard(k+1) (at least partially) we use real 
      variables xlow, xhigh for these indices so that they don't overflow for 
      very small dener */
      if(edivision==0){
        xlow=(eard[k]/gfac-ear_loc[0])/dener;
        xhigh=(eard[k+1]/gfac-ear_loc[0])/dener;
      }else if(edivision==1){
        xlow=(log10(eard[k]/gfac)-log10(ear_loc[0]))/dener;
        xhigh=(log10(eard[k+1]/gfac)-log10(ear_loc[0]))/dener;
      }
/*    if the lower index is smaller than 0 => the whole first interval of local 
      energies falls into this bin; we set lower index to 1 (real value to 0., 
      we round it to 1 later on) the upper energy must be larger than eard(k) in 
      this case...(otherwise no part of local energy falls into this bin)*/
      if((xlow<0.) && (xhigh>=0.))xlow=0.;
/*    if the upper index is larger than (neloc-1) => the last interval of local 
      energies also falls into this bin; we set upper index to (ne_loc-1) (real 
      value to (ne_loc-2.), we round it to (ne_loc-1) later on), we set it to 
      ne_loc-1 and not to ne_loc because the upper index is index of the lower 
      edge of the last interval of local energies that falls into the bin 
      eard(k)-eard(k+1) the lower energy must be smaller than eard(k) in this 
      case...(otherwise no part of local energy falls into this bin) */
      if((xhigh>=(ne_loc-1.)) && (xlow<(ne_loc-1.)))xhigh=ne_loc-2.;
//    otherwise no local energy falls into this bin:
      if((xlow<0.) || (xhigh>=(ne_loc-1.)))continue;
//    let's change real values of indices to integer values
      ilow=(int) round(xlow+0.5);
      ihigh=(int) round(xhigh+0.5);
/*    before we can add contributions of local photon flux density we have to
      decide into which 'time' element of far(ne,time) it falls */
//    here we loop over local time...
      for(it=0;it<nt;it++){
        if(far_loc[ne_loc+(ne_loc+1)*it]==0. && nt!=1)continue;
        if(nt==1)ishiftt=0;
        else{
          ishiftt = it + (int) trunc((delay-ddel)/dt);
          dshift = ((delay-ddel)/dt) - (int) trunc((delay-ddel)/dt);
          if(dshift>0.)ishiftt+=1;
          else dshift+=1;
/*        what to do, if shifted time falls outside time span covered by
          the array: */
          if(periodic==1){
/*          the case of a periodic spectrum (e.g. an orbiting spot); in this
            case emissivity subroutine must be constructed in such a way that
            it covers a complete period (or a multiple of it): */
            while(ishiftt<0)ishiftt+=nt;
            while(ishiftt>=nt)ishiftt-=nt;
          }
        }
/*      contribution from the current point in the disc is added to final
        spectrum:
        the case of non-periodic spectrum (e.g. a decaying flare):
        do nothing if spectral points fall outside the array (time)
        boundaries: */
        if((ishiftt<0) || (ishiftt>=nt))continue;
        l=ilow;
        while(l<=ihigh){
          ddener=ear_loc[l]-ear_loc[l-1];
/*        if we are calculating the Stokes parameters as well
          we have to integrate all 4 functions... 
          the emissivity routine is in f77 and the arrays there are defined as
          far_loc(0:ne_loc,nt),qar_loc(ne_loc,nt),uar_loc(ne_loc,nt),
          var_loc(ne_loc,nt), i.e. indices in f77 are from 0..ne_loc, 1..ne_loc 
          and 1..nt while in C they are from 0..ne_loc, 0..ne_loc-1, 0..nt-1, so 
          we must treat the indices here differently than in f77 emissivity 
          routine */
          iimax=1;
          if(polar)iimax=4;
          for(ii=1;ii<=iimax;ii++){
            if(ii==1){
              intfunc1=far_loc[l-1+(ne_loc+1)*it];
              intfunc2=far_loc[l+(ne_loc+1)*it];
            }else if(ii==2){
              intfunc1=qar_loc[l-1+ne_loc*it]*cos2pol-
                       uar_loc[l-1+ne_loc*it]*sin2pol;
              intfunc2=qar_loc[l+ne_loc*it]*cos2pol-
                       uar_loc[l+ne_loc*it]*sin2pol;
            }else if(ii==3){
              intfunc1=qar_loc[l-1+ne_loc*it]*sin2pol+
                       uar_loc[l-1+ne_loc*it]*cos2pol;
              intfunc2=qar_loc[l+ne_loc*it]*sin2pol+
                       uar_loc[l+ne_loc*it]*cos2pol;
            }else if(ii==4){
                intfunc1=var_loc[l-1+ne_loc*it];
                intfunc2=var_loc[l+ne_loc*it];
            }
/*          the whole interval of local energies ear_loc[l]-ear_loc[l-1]
            falls into the bin eard[k]-eard[k+1] */
            if((ear_loc[l-1]>=eard[k]/gfac) && (ear_loc[l]<=eard[k+1]/gfac))
              dfc=0.5*(intfunc1+intfunc2)*fact*ddener;
/*          only the higher part of the local energy interval falls into
            the bin eard[k]-eard[k+1] => we have to interpolate far_loc at
            eard[k] */
            else if((ear_loc[l-1]<eard[k]/gfac) && (ear_loc[l]<=eard[k+1]/gfac))
              dfc=fact*(ear_loc[l]-eard[k]/gfac)*
                  (intfunc2-0.5*(intfunc2-intfunc1)*
                  (ear_loc[l]-eard[k]/gfac)/ddener);
/*          only the lower part of the local energy interval falls into the bin 
            eard[k]-eard[k+1] => we have to interpolate far_loc at eard[k+1] */
            else if((ear_loc[l-1]>=eard[k]/gfac) && (ear_loc[l]>eard[k+1]/gfac))
              dfc=fact*(eard[k+1]/gfac-ear_loc[l-1])*
                  (intfunc1+0.5*(intfunc2-intfunc1)*
                  (eard[k+1]/gfac-ear_loc[l-1])/ddener);
/*          only the middle part of the local energy interval falls into
            the bin eard[k]-eard[k+1] => we have to interpolate far_loc at
            eard[k] and eard[k+1] as well */
            else if((ear_loc[l-1]<eard[k]/gfac) && (ear_loc[l]>eard[k+1]/gfac))
              dfc=fact*(eard[k+1]-eard[k])/gfac*
                  (intfunc1+0.5*(intfunc2-intfunc1)*
                  ((eard[k]+eard[k+1])/gfac-2.0*ear_loc[l-1])/ddener);
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    in the next line there should in fact be dfc*gfac because we have to
    calculate the addition of local flux with respect to current energy bin
    - but we accounted for this earlier - see the definition of fact
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
            if(nt==1){
              if(ii==1)fce[k]+=dfc;
              else if(ii==2)qce[k]+=dfc;
              else if(ii==3)uce[k]+=dfc;
              else if(ii==4)vce[k]+=dfc;
            }else{
              if(ii==1)fce[k+ne*ishiftt]+=dfc*dshift;
              else if(ii==2)qce[k+ne*ishiftt]+=dfc*dshift;
              else if(ii==3)uce[k+ne*ishiftt]+=dfc*dshift;
              else if(ii==4)vce[k+ne*ishiftt]+=dfc*dshift;
              if(ishiftt>0){
                if(ii==1)fce[k+ne*(ishiftt-1)]+=dfc*(1-dshift);
                else if(ii==2)qce[k+ne*(ishiftt-1)]+=dfc*(1-dshift);
                else if(ii==3)uce[k+ne*(ishiftt-1)]+=dfc*(1-dshift);
                else if(ii==4)vce[k+ne*(ishiftt-1)]+=dfc*(1-dshift);
              }else{
                if(periodic==1){
                  if(ii==1)fce[k+ne*(nt-1)]+=dfc*(1-dshift);
                  else if(ii==2)qce[k+ne*(nt-1)]+=dfc*(1-dshift);
                  else if(ii==3)uce[k+ne*(nt-1)]+=dfc*(1-dshift);
                  else if(ii==4)vce[k+ne*(nt-1)]+=dfc*(1-dshift);
                }
              }
            }
          }//next function (far, qar, uar, var)
          l+=1;
        }//next energy bin
      }//next it
    }//next energy bin
  }
  iiloop+=1;
}//next iiloop
return;
}
/*******************************************************************************
*******************************************************************************/
void bcucof(double y[], double y1[], double y2[], double y12[], double d1, 
            double d2, double c[]){
  
static int wt[16][16]=
  {{1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
   {0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0},
   {-3,0,0,3,0,0,0,0,-2,0,0,-1,0,0,0,0},
   {2,0,0,-2,0,0,0,0,1,0,0,1,0,0,0,0},
   {0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0},
   {0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0},
   {0,0,0,0,-3,0,0,3,0,0,0,0,-2,0,0,-1},
   {0,0,0,0,2,0,0,-2,0,0,0,0,1,0,0,1},
   {-3,3,0,0,-2,-1,0,0,0,0,0,0,0,0,0,0},
   {0,0,0,0,0,0,0,0,-3,3,0,0,-2,-1,0,0},
   {9,-9,9,-9,6,3,-3,-6,6,-6,-3,3,4,2,1,2},
   {-6,6,-6,6,-4,-2,2,4,-3,3,3,-3,-2,-1,-1,-2},
   {2,-2,0,0,1,1,0,0,0,0,0,0,0,0,0,0},
   {0,0,0,0,0,0,0,0,2,-2,0,0,1,1,0,0},
   {-6,6,-6,6,-3,-3,3,3,-4,4,2,-2,-2,-2,-1,-1},
   {4,-4,4,-4,2,2,-2,-2,2,-2,-2,2,1,1,1,1}};
int l,k,j,i;
double xx,d1d2,cl[16],x[16];

d1d2=d1*d2;
for (i=0;i<4;i++) {
  x[i]=y[i];
  x[i+4]=y1[i]*d1;
  x[i+8]=y2[i]*d2;
  x[i+12]=y12[i]*d1d2;
}
for (i=0;i<16;i++) {
  xx=0.0;
  for (k=0;k<16;k++) xx += wt[i][k]*x[k];
  cl[i]=xx;
}
l=0;
for (i=0;i<4;i++)
  for (j=0;j<4;j++) c[i+4*j]=cl[l++];
return;
}
/*******************************************************************************
*******************************************************************************/
void bcuint(double y[], double y1[], double y2[], double y12[], double x1l,
            double x1u, double x2l, double x2u, double x1, double x2, 
            double *ansy){
  
void bcucof(double y[], double y1[], double y2[], double y12[], double d1,
  double d2, double c[]);
int i;
double t,u,d1,d2,c[16];

d1=x1u-x1l;
d2=x2u-x2l;
bcucof(y,y1,y2,y12,d1,d2,c);
if (x1u == x1l || x2u == x2l){
  xs_write("ide: bad input in routine bcuint",5);exit(1);}
t=(x1-x1l)/d1;
u=(x2-x2l)/d2;
*ansy=0.0;
for (i=3;i>=0;i--) {
  *ansy=t*(*ansy)+((c[i+12]*u+c[i+8])*u+c[i+4])*u+c[i];
}
return;
}
/*******************************************************************************
*******************************************************************************/
