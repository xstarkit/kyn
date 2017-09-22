Table of contents
=================

  * [Description](#description)
  * [Installation](#installation)
    * [Required files](#required-files)
    * [Usage in XSPEC](#usage-in-xspec)
    * [Usage outside of XSPEC](#usage-outside-of-xspec)
  * [KYN models](#kyn-models)
    * [KYNrline](#kynrline)
    * [KYNrlpli](#kynrlpli)
    * [KYNconv](#kynconv)
    * [KYNclp](#kynclp)
    * [KYNlpcr](#kynlpcr)
    * [KYNrefionx](#kynrefionx)
    * [KYNxillver](#kynxillver)
    * [KYNhrefl](#kynhrefl)
    * [KYNphebb](#kynphebb)
    * [KYNbb](#kynbb)


Description
===========

KYN is a set of the XSPEC models for emission from the black hole accretion disc
with the following assumptions and features:

  * space-time around black hole is described by Kerr metric,
  * accretion disc is Keplerian, geometrically thin and optically thick,
  * disc area below ISCO may be chosen to emit radiation, in which case material
    there is assumed to be freely falling and has the same energy
    and angular momentum as the matter which is orbiting at the ISCO,
  * disc may be non-axisymmetric - only part of the disc may be emitting 
    (sections in radius and azimuth),
  * obscuration by circular cloud is possible,
  * full relativistic ray-tracing code in vacuum was used for photon paths to 
    compute the tables of transfer functions used in the models,
  * the pre-calculated tables contain impact parameters - &alpha; and &beta; 
    coordinates thus in principle non-Keplerian (but still geometrically thin) 
    discs are possible - however, the definitions of the disc velocity would 
    have to be changed inside the code,
  * the computations are parallelized with threads.


**Models included:**

  - *Relativistic fluorescent line models:*

      * [KYNrline](#kynrline) - relativistic line with broken power-law 
                                    radial emissivity,
      * [KYNrlpli](#kynrlpli) -  relativistic line in lamp-post geometry.

  - *Relativistic convolution models:*

      * [KYNconv](#kynconv) - relativistic convolution model with broken 
                                  power-law radial emissivity,
      * [KYNclp](#kynclp) - relativistic convolution model in lamp-post 
                                geometry.

  - *Relativistic reflection models:*

      * [KYNlpcr](#kynlpcr) - relativistic reflection model in lamp-post 
                                 geometry for neutral disc (local emissivity 
                                 computed by NOAR),
      * [KYNrefionx](#kynrefionx) - relativistic reflection model in 
                                       lamp-post geometry for ionised disc 
                                       (local emissivity is given by REFLIONX),
      * [KYNxillver](#kynxillver) - relativistic reflection model in 
                                       lamp-post geometry for ionised disc 
                                       (local emissivity is given by XILLVER),
      * [KYNhrefl](#kynhrefl) - relativistic reflection model  with broken 
                                   power-law emissivity for continuum based on 
                                   HREFL(POWERLAW).

  - *Thermal radiation models:*

      * [KYNphebb](#kynphebb) - relativistic thermal radiation with radial 
                                   power-law temperature profile,
      * [KYNbb](#kynbb) - relativistic thermal radiation with Novikov-Thorne 
                              temperature profile (without self-irradiation and 
                              with colour correction factor).

The KYN package is based upon its first version presented in 
Dovciak M., Karas V. & Yaqoob T. 2004, ApJS, 153, 205-221.

Installation
============

Required files
--------------

* Source files in the main repository directory.

* KY tables: [KBHlamp_q.fits](https://owncloud.asu.cas.cz/index.php/s/upMLegCGecuik03) 
  (also [here](http://www.astro.cas.cz/dovciak/pub/KY/KBHlamp_q.fits)) 
  and [KBHtables80.fits](https://owncloud.asu.cas.cz/index.php/s/WP8aLN168MJgcB9) 
  (also [here](http://www.astro.cas.cz/dovciak/pub/KY/KBHtables80.fits)).

* Some models need FITS tables for local re-processing of photons in the 
  accretion disc:

  * tables computed with Monte Carlo code NOAR (Dumont, A.-M., Abrassart, A., & 
    Collin, S. 2000, A&A, 357, 823):

     - [fluorescent_line.fits](https://owncloud.asu.cas.cz/index.php/s/Za0YaAvk2wsj13L) 
       (also [here](http://www.astro.cas.cz/dovciak/pub/KY/fluorescent_line.fits)), 
     - [reflspectra.fits](https://owncloud.asu.cas.cz/index.php/s/svvZLBq2GogNsju) 
       (also [here](http://www.astro.cas.cz/dovciak/pub/KY/reflspectra.fits)),

  * [REFLION(X)](https://heasarc.gsfc.nasa.gov/xanadu/xspec/models/reflion.html) 
    tables (Ross & Fabian 2005, MNRAS, 358, 211) - unpack gzipped files: 

     - [reflion.mod](https://heasarc.gsfc.nasa.gov/xanadu/xspec/models/reflion.mod.gz) (old),
     - [reflionx.mod](https://heasarc.gsfc.nasa.gov/xanadu/xspec/models/reflionx.mod.gz),

     or in case the links are not available or if the tables there are updated and 
     their format/structure has changed:

     - [reflion.mod](https://owncloud.asu.cas.cz/index.php/s/6CWcb0o5Ssjehju) 
       (or [here](http://www.astro.cas.cz/dovciak/pub/KY-external/reflion.mod)) (old),
     - [reflionx.mod](https://owncloud.asu.cas.cz/index.php/s/Q6biiTPM1QBMtiT)
       (or [here](http://www.astro.cas.cz/dovciak/pub/KY-external/reflionx.mod)),

  * [XILLVER](https://hea-www.cfa.harvard.edu/%7Ejavier/xillver/) tables
    (Garcia & Kallman 2010, ApJ, 718, 695 and Garcia et al. 2013, ApJ, 768, 2): 

     - [xillver.fits](https://hea-www.cfa.harvard.edu/%7Ejavier/xillver/tables/xillver.fits),
     - [xillver-a.fits](https://hea-www.cfa.harvard.edu/%7Ejavier/xillver/tables/xillver-a.fits),
     - [xillver-Ec.fits](https://hea-www.cfa.harvard.edu/%7Ejavier/xillver/tables/xillver-Ec.fits),
     - [xillver-a-Ec.fits](https://hea-www.cfa.harvard.edu/%7Ejavier/xillver/tables/xillver-a-Ec.fits),
     - [xillver-a-Ec2.fits](https://hea-www.cfa.harvard.edu/%7Ejavier/xillver/tables/xillver-a-Ec2.fits),
     - [xillver-a-Ec3.fits](https://hea-www.cfa.harvard.edu/%7Ejavier/xillver/tables/xillver-a-Ec3.fits),
     - [xillver-a-Ec4.fits](https://hea-www.cfa.harvard.edu/%7Ejavier/xillver/tables/xillver-a-Ec4.fits),
     - [xillver-a-Ec5.fits](https://hea-www.cfa.harvard.edu/%7Ejavier/xillver/tables/xillver-a-Ec5.fits).


Usage in XSPEC
--------------

The code is compiled inside XSPEC with the following command (assuming all the 
source files and FITS tables are in the directory /path/to/KYN):

* `initpackage kyn lmodel.dat /path/to/KYN`

To use the KYN models inside XSPEC, first the package needs to be loaded 
and directory with KYN set:

* `lmod kyn /path/to/KYN`
* `xset KYDIR /path/to/KYN`

Then any model from KYN package may be used, e.g.:
* `mo kynrline`

_Note_: 
In case of segmentation fault, one may need to increase the stack size, e.g. 
with the command `ulimit -s unlimited` or `ulimit -s 65532`.


Usage outside of XSPEC
----------------------

* One also needs the Makefile and libxspec library included in the directory 
  'other'.

* The library to work with FITS files (libcfitsio.so) is needed, thus one needs 
  to define the name of the library and path to it in the provided Makefile.

* The model parameters have to be changed inside the source file:

  - _energy_ in the following lines:

                  #define NE     30
                  #define E_MIN  0.3
                  #define E_MAX  80.

         and later choose if linear or exponential energy binning should be used:

                  for (ie = 0; ie <= NE; ie++) {
                  //  ear[ie] = E_MIN + ie * (E_MAX-E_MIN) / NE;
                    ear[ie] = E_MIN * pow(E_MAX/E_MIN, ((double) ie) / NE);
                  }

  - _all basic parameters_  of the models (physical ones as well as those 
    defining resolution grid for computations) are defined in the following 
    lines:

                  param[ 0] = 1.;       // a/M
                  param[ 1] = 30.;      // theta_o
                  param[ 2] = 1.;       // rin
                  param[ 3] = 1.;       // ms
                  param[ 4] = 400.;     // rout
                  param[ 5] = 0.;       // phi0
                  param[ 6] = 360.;     // dphi
                  ...
                  ...
                  ...


* Compile with the make command, e.g.:

    * `make kynrline`

* Run the code, e.g.:

    * `./kynrline`

* The model creates the file with the parameters used and the file with the 
  computed spectrum, e.g.:

    * kynrline.txt
    * kynrline_photar.dat

_Note_: 

* To use KYNLPCR, KYNREFIONX and KYNXILLVER outside of XSPEC one still needs the 
  local installation of XSPEC because its libraries are needed (this might be 
  changed in the future). To compile these models one needs to change 
  definitions of the XSPEC directories inside Makefile.

* In case of segmentation fault, one may need to increase the stack size, e.g. 
  with the command `ulimit -s unlimited` or `ulimit -s 65532`.


KYN models
==========

KYNrline
--------

Relativistically broadened emission line from black hole accretion disc.
Black hole may be rotating (Kerr black hole), accretion disc is assumed to be 
Keplerian, geometrically thin and optically thick. Broken power-law radial 
emissivity and several limb darkening/brightening laws for emission 
directionality are implemented. Only part of the disc may be set to be emitting 
radiation (sections defined in radius and azimuth). Obscuration by circular 
cloud is possible. Full relativistic ray-tracing code in vacuum was used for 
photon paths to compute the tables of transfer functions used in the model. Disc 
area below ISCO may be chosen to emit radiation, in which case material there is 
assumed to be freely falling and has the same energy and angular momentum as the 
matter orbiting at the ISCO. The model is based on its first version presented 
in Dovciak M., Karas V. & Yaqoob T. 2004, ApJS, 153, 205-221.

Definition of the parameters:

  * **par1  ... a/M**
    - black hole angular momentum (-1 &le; a/M &le; 1)
  * **par2  ... theta_o**
    - observer inclination in degrees (0&deg;-pole, 90&deg;-disc)
  * **par3  ... rin**
    - inner edge of non-zero disc emissivity (in GM/c^2 or in r~mso~)
  * **par4  ... ms**
    - switch for inner edge
    - 0: we integrate from inner edge = par3 
    - 1: if the inner edge of the disc is below marginally stable orbit (MSO) 
         then we integrate emission above MSO only
    - 2: we integrate from inner edge given in units of MSO, i.e. inner 
         edge = par3 &times; r~mso~ (the same applies for outer edge)
  * **par5  ... rout**
    - outer edge of non-zero disc emissivity (in GM/c^2 or in r~mso~)
  * **par6  ... phi**
    - lower azimuth of non-zero disc emissivity (degrees)
  * **par7  ... dphi**
    - (phi + dphi) is upper azimuth of non-zero disc emissivity 0&deg; &le; 
      dphi &le; 360&deg;
  * **par8  ... Erest**
    - rest energy of the line (keV)
  * **par9  ... sigma** 
    - width of the line - Gaussian sigma (eV)
  * **par10 ... q_out** 
    - power-law index for radial dependence of emissivity for outer region, 
      scales as r^(-q_out)
  * **par11  ... q_in**
    - power-law index for radial dependence of emissivity for inner region, 
      scales as rb^(q_in-q_out) &times; r^(-q_in)
  * **par12 ... rb**
    - boundary between the region with power-law index q_out and q_in
    - if > 0 then the boundary is in units of MSO, i.e. boundary = rb &times; 
      r~mso~
    - if <= 0 then the boundary is equal to -rb+r_horizon where rb is in GM/c^2
  * **par13 ... jump**
    - ratio of local flux in inner region to local flux in outer region at 
      boundary radius defined by rb
  * **par14 ... limb**
    - limb darkening/brightening law for emission directionality
    -  0: for isotropic emission (flux ~ 1)
    - -1: for Laor's limb darkening (flux ~ 1+2.06&mu;)
    - -2: for Haardt's limb brightening (flux ~ ln (1+1/&mu;))
    - if different from 0, -1 and -2 then the local emisivity is ~ mu^(limb)
  * **par15 ... alpha**
    - position of the cloud centre in GM/c^2 in alpha coordinate (alpha being 
      the impact parameter in &phi;-direction, positive for approaching side 
      of the disc)
  * **par16 ... beta**
    - position of the cloud centre in GM/c^2 in beta coordinate (beta being 
      the impact parameter in &theta;-direction, positive in up direction, 
      i.e. above the disc)
  * **par17 ... rcloud**
    - radius of the obscuring cloud
    - the meaning of cloud is inverted for negative values of rcloud, i.e. 
      only the radiation transmitted through the cloud is computed
  * **par18 ... zshift**
    - overall Doppler shift
  * **par19 ... ntable**
    - table of relativistic transfer functions used in the model
      (defines FITS file with tables), 0 &le; ntable &le; 99, currently the 
      tables with ntable=80 are correct for this model
  * **par20 ... nrad**
    - number of grid points in radius
  * **par21 ... division**
    - type of division in radial integration
    - 0: equidistant radial grid (constant linear step)
    - 1: exponential radial grid (constant logarithmic step)
  * **par22 ... nphi**
    - number of grid points in azimuth
  * **par23 ... smooth**
    - whether to smooth the resulting spectrum 
    - 0: no smoothing
    - 1: simple smoothing
  * **par24 ... Stokes** 
    - definition of output
    - 0: photon number density flux per bin (Stokes parameter I devided by 
         energy)
    - 1: Stokes parameter Q devided by energy
    - 2: Stokes parameter U devided by energy
    - 3: Stokes parameter V devided by energy
    - 4: degree of polarisation
    - 5: linear polarisation angle &psi; = 0.5 atan(U/Q)
    - 6: circular polarisation angle &beta; = 0.5 asin(V/sqrt(Q^(2)+U^(2)+V^(2)))
  * **par25 ... nthreads**
    - number of threads used for computations
  * **par26 ... normtype** 
    - how to normalize the spectra
    - 0: normalization to the total photon flux
    - > 0: normalization to the photon flux at par26 keV
    - -1: the photon flux is not re-normalized,
    - -2: normalization to the maximum of the photon flux
  * **par27 ... norm**
    - depending on par26 it sets either the total photon flux, flux at given 
      energy or maximum of the photon flux

_Note:_
  * KYRH (the black hole horizon), KYRIN (the disc inner edge) and 
    KYRMS (the marginally stable orbit) are added to the XSPEC internal
    switches. Use xset command to show their current values.
  * Accuracy vs. speed trade off depends mainly on nrad and nphi.
  * In this model it is assumed that the local emission is completely linearly 
    polarised in the direction perpendicular to the disc.


KYNrlpli
--------

Relativistically broadened emission line from black hole accretion disc
in the lamp-post geometry. Black hole may be rotating (Kerr black hole), 
accretion disc is assumed to be Keplerian, geometrically thin and optically 
thick. The radial emissivity is given by the disc illumination from a point 
source located at some height on the system axis. The fluorescent iron line 
is modelled by Monte Carlo code NOAR (Dumont, A.-M., Abrassart, A., & Collin, S. 
2000, A&A, 357, 823) that gives slight limb brightening emission
directionality that depends on the primary source power-law photon index. 
In this model only part of the disc may be set to be emitting 
radiation (sections defined in radius and azimuth). Obscuration by circular 
cloud is possible. Full relativistic ray-tracing code in vacuum was used for 
photon paths to compute the tables of transfer functions used in the model. Disc 
area below ISCO may be chosen to emit radiation, in which case material there is 
assumed to be freely falling and has the same energy and angular momentum as the 
matter orbiting at the ISCO. The model is presented in Dovciak, M., Svoboda, J., 
Goosmann, R. W., et al.: 2014, in Proceedings of RAGtime 14-16: Workshops on 
black holes and neutron stars (Silesian University in Opava), [arXiv:1412.8627].

Definition of the parameters:

  * **par1  ... a/M**
    - black hole angular momentum (-1 &le; a/M &le; 1)
  * **par2  ... theta_o**
    - observer inclination in degrees (0&deg;-pole, 90&deg;-disc)
  * **par3  ... rin**
    - inner edge of non-zero disc emissivity (in GM/c^2 or in r~mso~)
  * **par4  ... ms**
    - switch for inner edge
    - 0: we integrate from inner edge = par3 
    - 1: if the inner edge of the disc is below marginally stable orbit (MSO) 
         then we integrate emission above MSO only
    - 2: we integrate from inner edge given in units of MSO, i.e. inner 
         edge = par3 &times; r~mso~ (the same applies for outer edge)
  * **par5  ... rout**
    - outer edge of non-zero disc emissivity (in GM/c^2 or in r~mso~)
  * **par6  ... phi**
    - lower azimuth of non-zero disc emissivity (degrees)
  * **par7  ... dphi**
    - (phi + dphi) is upper azimuth of non-zero disc emissivity 0&deg; &le; 
      dphi &le; 360&deg;
  * **par8  ... height**
    - height on the axis (measured from the center) at which the primary 
      source is located (GM/c^(2))
  * **par9  ... PhoIndex**
    - power-law energy index of the primary flux
  * **par10 ... alpha**
    - position of the cloud centre in GM/c^2 in alpha coordinate (alpha being 
      the impact parameter in &phi;-direction, positive for approaching side 
      of the disc)
  * **par11 ... beta**
    - position of the cloud centre in GM/c^2 in beta coordinate (beta being 
      the impact parameter in &theta;-direction, positive in up direction, 
      i.e. above the disc)
  * **par12 ... rcloud**
    - radius of the obscuring cloud
    - the meaning of cloud is inverted for negative values of rcloud, i.e. 
      only the radiation transmitted through the cloud is computed
  * **par13 ... zshift**
    - overall Doppler shift
  * **par14 ... ntable**
    - table of relativistic transfer functions used in the model
      (defines FITS file with tables), 0 &le; ntable &le; 99, currently the 
      tables with ntable=80 are correct for this model
  * **par15 ... nrad**
    - number of grid points in radius
  * **par16 ... division**
    - type of division in radial integration
    - 0: equidistant radial grid (constant linear step)
    - 1: exponential radial grid (constant logarithmic step)
  * **par17 ... nphi**
    - number of grid points in azimuth
  * **par18 ... smooth**
    - whether to smooth the resulting spectrum 
    - 0: no smoothing
    - 1: simple smoothing
  * **par19 ... Stokes** 
    - definition of output
    - 0: photon number density flux per bin (Stokes parameter I devided by 
         energy)
    - 1: Stokes parameter Q devided by energy
    - 2: Stokes parameter U devided by energy
    - 3: Stokes parameter V devided by energy
    - 4: degree of polarisation
    - 5: linear polarisation angle &psi; = 0.5 atan(U/Q)
    - 6: circular polarisation angle &beta; = 0.5 asin(V/sqrt(Q^(2)+U^(2)+V^(2)))
  * **par20 ... nthreads**
    - number of threads used for computations
  * **par21 ... normtype** 
    - how to normalize the spectra
    - 0: normalization to the total photon flux
    - > 0: normalization to the photon flux at par21 keV
    - -1: the photon flux is not re-normalized,
    - -2: normalization to the maximum of the photon flux
  * **par22 ... norm**
    - depending on par21 it sets either the total photon flux, flux at given 
      energy or maximum of the photon flux

_Note:_
  * KYRH (the black hole horizon), KYRIN (the disc inner edge) and 
    KYRMS (the marginally stable orbit) are added to the XSPEC internal
    switches. Use xset command to show their current values.
  * Accuracy vs. speed trade off depends mainly on nrad and nphi.
  * In this model it is assumed that the local emission is completely linearly 
    polarised in the direction perpendicular to the disc.


KYNconv
-------

Convolution model for relativistic broadening of the emission coming from black 
hole accretion disc. Black hole may be rotating (Kerr black hole), accretion 
disc is assumed to be Keplerian, geometrically thin and optically thick. Broken 
power-law radial emissivity and several limb darkening/brightening laws for 
emission directionality are implemented. Only part of the disc may be set to be 
emitting radiation (sections defined in radius and azimuth). Obscuration by 
circular cloud is possible. Full relativistic ray-tracing code in vacuum was 
used for photon paths to compute the tables of transfer functions used in the 
model. Disc area below ISCO may be chosen to emit radiation, in which case 
material there is assumed to be freely falling and has the same energy and 
angular momentum as the matter orbiting at the ISCO. The model is based on its 
first version presented in Dovciak M., Karas V. & Yaqoob T. 2004, ApJS, 153, 
205-221.

Definition of the parameters:

  * **par1  ... a/M**
    - black hole angular momentum (-1 &le; a/M &le; 1)
  * **par2  ... theta_o**
    - observer inclination in degrees (0&deg;-pole, 90&deg;-disc)
  * **par3  ... rin**
    - inner edge of non-zero disc emissivity (in GM/c^2 or in r~mso~)
  * **par4  ... ms**
    - switch for inner edge
    - 0: we integrate from inner edge = par3 
    - 1: if the inner edge of the disc is below marginally stable orbit (MSO) 
         then we integrate emission above MSO only
    - 2: we integrate from inner edge given in units of MSO, i.e. inner 
         edge = par3 &times; r~mso~ (the same applies for outer edge)
  * **par5  ... rout**
    - outer edge of non-zero disc emissivity (in GM/c^2 or in r~mso~)
  * **par6  ... phi**
    - lower azimuth of non-zero disc emissivity (degrees)
  * **par7  ... dphi**
    - (phi + dphi) is upper azimuth of non-zero disc emissivity 0&deg; &le; 
      dphi &le; 360&deg;
  * **par8  ... q_out** 
    - power-law index for radial dependence of emissivity for outer region, 
      scales as r^(-q_out)
  * **par9  ... q_in**
    - power-law index for radial dependence of emissivity for inner region, 
      scales as rb^(q_in-q_out) &times; r^(-q_in)
  * **par10 ... rb**
    - boundary between the region with power-law index q_out and q_in
    - if > 0 then the boundary is in units of MSO, i.e. boundary = rb &times; 
      r~mso~
    - if <= 0 then the boundary is equal to -rb+r_horizon where rb is in GM/c^2
  * **par11 ... jump**
    - ratio of local flux in inner region to local flux in outer region at 
      boundary radius defined by rb
  * **par12 ... limb**
    - limb darkening/brightening law for emission directionality
    -  0: for isotropic emission (flux ~ 1)
    - -1: for Laor's limb darkening (flux ~ 1+2.06&mu;)
    - -2: for Haardt's limb brightening (flux ~ ln (1+1/&mu;))
    - if different from 0, -1 and -2 then the local emisivity is ~ mu^(limb)
  * **par13 ... alpha**
    - position of the cloud centre in GM/c^2 in alpha coordinate (alpha being 
      the impact parameter in &phi;-direction, positive for approaching side 
      of the disc)
  * **par14 ... beta**
    - position of the cloud centre in GM/c^2 in beta coordinate (beta being 
      the impact parameter in &theta;-direction, positive in up direction, 
      i.e. above the disc)
  * **par15 ... rcloud**
    - radius of the obscuring cloud
    - the meaning of cloud is inverted for negative values of rcloud, i.e. 
      only the radiation transmitted through the cloud is computed
  * **par16 ... zshift**
    - overall Doppler shift
  * **par17 ... ntable**
    - table of relativistic transfer functions used in the model
      (defines FITS file with tables), 0 &le; ntable &le; 99, currently the 
      tables with ntable=80 are correct for this model
  * **par18 ... nrad**
    - number of grid points in radius
  * **par19 ... division**
    - type of division in radial integration
    - 0: equidistant radial grid (constant linear step)
    - 1: exponential radial grid (constant logarithmic step)
  * **par20 ... nphi**
    - number of grid points in azimuth
  * **par21 ... ne_loc**
    - number of grid points in local energy (energy resolution of local flux), 
      the grid is equidistant in logarithmic scale
  * **par22 ... smooth**
    - whether to smooth the resulting spectrum 
    - 0: no smoothing
    - 1: simple smoothing
  * **par23 ... Stokes** 
    - definition of output
    - 0: photon number density flux per bin (Stokes parameter I devided by 
         energy)
    - 1: Stokes parameter Q devided by energy
    - 2: Stokes parameter U devided by energy
    - 3: Stokes parameter V devided by energy
    - 4: degree of polarisation
    - 5: linear polarisation angle &psi; = 0.5 atan(U/Q)
    - 6: circular polarisation angle &beta; = 0.5 asin(V/sqrt(Q^(2)+U^(2)+V^(2)))
  * **par24 ... nthreads**
    - number of threads used for computations
  * **par25 ... normtype** 
    - how to normalize the spectra
    - 0: normalization to the total photon flux
    - > 0: normalization to the photon flux at par25 keV
    - -1: the photon flux is not re-normalized,
    - -2: normalization to the maximum of the photon flux

_Note:_
  * KYRH (the black hole horizon), KYRIN (the disc inner edge) and 
    KYRMS (the marginally stable orbit) are added to the XSPEC internal
    switches. Use xset command to show their current values.
  * There are several restrictions that arise from the fact that we use existing
    XSPEC models for definition of the local flux:
    - only the energy dependence of the photon flux can be defined by local XSPEC
      models,
    - only a certain type of radial dependence of the local photon flux can be
      imposed; here we have chosen to use a broken power-law radial dependence,
    - there is no intrinsic azimuthal dependence of the local photon flux, the
      only azimuthal dependence comes through limb darkening/brightening law
      (emission angle depends on azimuth)
    - local flux can highly depend on the energy resolution, i.e. on the energy
      binning used, if the energy resolution is not high enough. This is because
      the flux is defined in the centre of each bin. A large number of bins is
      needed for highly varying local flux with energy.
  * Accuracy vs. speed trade off depends mainly on nrad, nphi and ne_loc.
  * In this model it is assumed that the local emission is completely linearly 
    polarised in the direction perpendicular to the disc.


KYNclp
------

Convolution model for relativistic broadening of the emission from black hole 
accretion disc in the lamp-post geometry. Black hole may be rotating (Kerr black 
hole), accretion disc is assumed to be Keplerian, geometrically thin and 
optically thick. The radial emissivity is given by the disc illumination from a 
point source located at some height on the system axis. The emission 
directionality is modelled by Monte Carlo code NOAR (Dumont, A.-M., 
Abrassart, A., & Collin, S. 2000, A&A, 357, 823) that gives slight limb 
brightening emission law. In this model only part of the disc may be set to be 
emitting radiation (sections defined in radius and azimuth). Obscuration by 
circular cloud is possible. Full relativistic ray-tracing code in vacuum was 
used for photon paths to compute the tables of transfer functions used in the 
model. Disc area below ISCO may be chosen to emit radiation, in which case 
material there is assumed to be freely falling and has the same energy and 
angular momentum as the matter orbiting at the ISCO. The model is presented in 
Dovciak, M., Svoboda, J., Goosmann, R. W., et al.: 2014, in Proceedings of 
RAGtime 14-16: Workshops on black holes and neutron stars (Silesian University 
in Opava), [arXiv:1412.8627].

Definition of the parameters:

  * **par1  ... a/M**
    - black hole angular momentum (-1 &le; a/M &le; 1)
  * **par2  ... theta_o**
    - observer inclination in degrees (0&deg;-pole, 90&deg;-disc)
  * **par3  ... rin**
    - inner edge of non-zero disc emissivity (in GM/c^2 or in r~mso~)
  * **par4  ... ms**
    - switch for inner edge
    - 0: we integrate from inner edge = par3 
    - 1: if the inner edge of the disc is below marginally stable orbit (MSO) 
         then we integrate emission above MSO only
    - 2: we integrate from inner edge given in units of MSO, i.e. inner 
         edge = par3 &times; r~mso~ (the same applies for outer edge)
  * **par5  ... rout**
    - outer edge of non-zero disc emissivity (in GM/c^2 or in r~mso~)
  * **par6  ... phi**
    - lower azimuth of non-zero disc emissivity (degrees)
  * **par7  ... dphi**
    - (phi + dphi) is upper azimuth of non-zero disc emissivity 0&deg; &le; 
      dphi &le; 360&deg;
  * **par8  ... height**
    - height on the axis (measured from the center) at which the primary 
      source is located (GM/c^(2))
  * **par9  ... PhoIndex**
    - power-law energy index of the primary flux
  * **par10 ... alpha**
    - position of the cloud centre in GM/c^2 in alpha coordinate (alpha being 
      the impact parameter in &phi;-direction, positive for approaching side 
      of the disc)
  * **par11 ... beta**
    - position of the cloud centre in GM/c^2 in beta coordinate (beta being 
      the impact parameter in &theta;-direction, positive in up direction, 
      i.e. above the disc)
  * **par12 ... rcloud**
    - radius of the obscuring cloud
    - the meaning of cloud is inverted for negative values of rcloud, i.e. 
      only the radiation transmitted through the cloud is computed
  * **par13 ... zshift**
    - overall Doppler shift
  * **par14 ... ntable**
    - table of relativistic transfer functions used in the model
      (defines FITS file with tables), 0 &le; ntable &le; 99, currently the 
      tables with ntable=80 are correct for this model
  * **par15 ... nrad**
    - number of grid points in radius
  * **par16 ... division**
    - type of division in radial integration
    - 0: equidistant radial grid (constant linear step)
    - 1: exponential radial grid (constant logarithmic step)
  * **par17 ... nphi**
    - number of grid points in azimuth
  * **par18 ... ne_loc**
    - number of grid points in local energy (energy resolution of local flux), 
      the grid is equidistant in logarithmic scale
  * **par19 ... smooth**
    - whether to smooth the resulting spectrum 
    - 0: no smoothing
    - 1: simple smoothing
  * **par20 ... Stokes** 
    - definition of output
    - 0: photon number density flux per bin (Stokes parameter I devided by 
         energy)
    - 1: Stokes parameter Q devided by energy
    - 2: Stokes parameter U devided by energy
    - 3: Stokes parameter V devided by energy
    - 4: degree of polarisation
    - 5: linear polarisation angle &psi; = 0.5 atan(U/Q)
    - 6: circular polarisation angle &beta; = 0.5 asin(V/sqrt(Q^(2)+U^(2)+V^(2)))
  * **par21 ... nthreads**
    - number of threads used for computations
  * **par22 ... normtype** 
    - how to normalize the spectra
    - 0: normalization to the total photon flux
    - > 0: normalization to the photon flux at par22 keV
    - -1: the photon flux is not re-normalized,
    - -2: normalization to the maximum of the photon flux

_Note:_
  * KYRH (the black hole horizon), KYRIN (the disc inner edge) and 
    KYRMS (the marginally stable orbit) are added to the XSPEC internal
    switches. Use xset command to show their current values.
  * There are several restrictions that arise from the fact that we use existing
    XSPEC models for definition of the local flux:
    - only the energy dependence of the photon flux can be defined by local 
      XSPEC models,
    - only a certain type of radial dependence of the local photon flux can be
      imposed; here we have chosen to use radial emissivity profile given by an 
      illumination by a point source on the system axis,
    - there is no intrinsic azimuthal dependence of the local photon flux, the
      only azimuthal dependence comes through limb darkening/brightening law
      (emission angle depends on azimuth)
    - local flux can highly depend on the energy resolution, i.e. on the energy
      binning used, if the energy resolution is not high enough. This is because
      the flux is defined in the centre of each bin. A large number of bins is
      needed for highly varying local flux with energy.
  * Accuracy vs. speed trade off depends mainly on nrad, nphi and ne_loc.
  * In this model it is assumed that the local emission is completely linearly 
    polarised in the direction perpendicular to the disc.


KYNlpcr
-------
 
Neutral reflection spectrum from black hole accretion disc in the lamp-post 
geometry. Black hole may be rotating (Kerr black hole), accretion disc is 
assumed to be Keplerian, geometrically thin and optically thick. The radial 
emissivity is determined by the disc illumination from a point source located at 
some height on the system axis. The primary source is isotropic and emits 
power-law radiation. Re-processing in the neutral disc is modelled by 
Monte Carlo code NOAR (Dumont, A.-M., Abrassart, A., & Collin, S. 
2000, A&A, 357, 823) that gives slight limb brightening emission
directionality. In this model only part of the disc may be set to be emitting 
radiation (sections defined in radius and azimuth). Obscuration by circular 
cloud is possible. Full relativistic ray-tracing code in vacuum was used for 
photon paths to compute the tables of transfer functions used in the model. Disc 
area below ISCO may be chosen to emit radiation, in which case material there is 
assumed to be freely falling and has the same energy and angular momentum as the 
matter orbiting at the ISCO. The model is based on KY package of models first 
presented in Dovciak M., Karas V. & Yaqoob T. 2004, ApJS, 153, 205-221 and 
later in Dovciak, M., Svoboda, J., Goosmann, R. W., et al. 2014, in Proceedings 
of RAGtime 14-16: Workshops on black holes and neutron stars (Silesian 
University in Opava), [arXiv:1412.8627].

Definition of the parameters:

  * **par1  ... a/M**
    - black hole angular momentum (-1 &le; a/M &le; 1)
  * **par2  ... theta_o**
    - observer inclination in degrees (0&deg;-pole, 90&deg;-disc)
  * **par3  ... rin**
    - inner edge of non-zero disc emissivity (in GM/c^2 or in r~mso~)
  * **par4  ... ms**
    - switch for inner edge
    - 0: we integrate from inner edge = par3 
    - 1: if the inner edge of the disc is below marginally stable orbit (MSO) 
         then we integrate emission above MSO only
    - 2: we integrate from inner edge given in units of MSO, i.e. inner 
         edge = par3 &times; r~mso~ (the same applies for outer edge)
  * **par5  ... rout**
    - outer edge of non-zero disc emissivity (in GM/c^2 or in r~mso~)
  * **par6  ... phi**
    - lower azimuth of non-zero disc emissivity (degrees)
  * **par7  ... dphi**
    - (phi + dphi) is upper azimuth of non-zero disc emissivity 0&deg; &le; 
      dphi &le; 360&deg;
  * **par8  ... M/M8**
       - black hole mass in units of 10^8 solar masses
  * **par9  ... height**
    - height on the axis (measured from the center) at which the primary 
      source is located (GM/c^(2))
  * **par10 ... PhoIndex**
    - power-law energy index of the primary flux
  * **par11 ... L/L~Edd~**
    - dE/dt, the observed (if positive) or the intrinsic local (if negative) 
      primary isotropic flux in the X-ray energy range 2-10keV in units of 
      L~Edd~
  * **par12 ... Np:Nr**
    - ratio of the primary to the reflected normalization
    - 1: self-consistent model for isotropic primary source
    - 0: only reflection, primary source is hidden
    - if positive then L/L~Edd~ (par11) means the luminosity towards the 
      observer
    - if negative then L/L~Edd~ (par11) means the luminosity towards the disc
  * **par13 ... line**
    - whether to include lines and/or reflection continuum in the spectra
    -  0: only continuum
    -  1: K&alpha; Fe line with continuum
    -  2: K&alpha; and K&beta; lines with continuum
    -  3: all lines computed by NOAR for neutral disc
    - -1: only K&alpha; Fe line without the reflection continuum
    - -2: K&alpha; and K&beta; lines without the reflection continuum
    - -3: all lines computed by NOAR for neutral disc without the reflection continuum
  * **par14 ... E_cut**
    - cut-off energy
  * **par15 ... alpha**
    - position of the cloud centre in GM/c^2 in alpha coordinate (alpha being 
      the impact parameter in &phi;-direction, positive for approaching side 
      of the disc)
  * **par16 ... beta**
    - position of the cloud centre in GM/c^2 in beta coordinate (beta being 
      the impact parameter in &theta;-direction, positive in up direction, 
      i.e. above the disc)
  * **par17 ... rcloud**
    - radius of the obscuring cloud
    - the meaning of cloud is inverted for negative values of rcloud, i.e. 
      only the radiation transmitted through the cloud is computed
  * **par18 ... zshift**
    - > 0: the overall Doppler shift z; the distance of the source is computed 
           from the Hubble law, D = zc/H, with the Hubble constant 
           H = 70 km/s/Mpc; this distance is used to compute the correct 
           normalisation of the spectrum, thus norm parameter (the normalisation 
           of the model) should be frozen to unity
    - = 0: the overall Doppler shift is set to z = 0; the distance to the source 
           and thus the correct normalisation is defined by the normalisation 
           parameter
    - < 0: the negative value of the overall Doppler shift z used only to 
           compute the distance of the source from the Hubble law; this distance 
           is used to compute the correct normalisation of the spectrum, thus 
           the normalisation parameter should be frozen to unity; the overall 
           Doppler shift is then set to z = 0 (i.e. there is no shift of 
           spectrum with energy)
  * **par19 ... ntable**
    - table of relativistic transfer functions used in the model
      (defines FITS file with tables), 0 &le; ntable &le; 99, currently the 
      tables with ntable=80 are correct for this model
  * **par20 ... nrad**
    - number of grid points in radius
  * **par21 ... division**
    - type of division in radial integration
    - 0: equidistant radial grid (constant linear step)
    - 1: exponential radial grid (constant logarithmic step)
  * **par22 ... nphi**
    - number of grid points in azimuth
  * **par23 ... smooth**
    - whether to smooth the resulting spectrum 
    - 0: no smoothing
    - 1: simple smoothing
  * **par24 ... Stokes** 
    - definition of output
    - 0: photon number density flux per bin (Stokes parameter I devided by 
         energy)
    - 1: Stokes parameter Q devided by energy
    - 2: Stokes parameter U devided by energy
    - 3: Stokes parameter V devided by energy
    - 4: degree of polarisation
    - 5: linear polarisation angle &psi; = 0.5 atan(U/Q)
    - 6: circular polarisation angle &beta; = 0.5 asin(V/sqrt(Q^(2)+U^(2)+V^(2)))
  * **par25 ... nthreads**
    - number of threads used for computations
  * **par26 ... norm**
    - if the overall Doppler shift zshift = 0, then norm = 1/D^2 where D is the 
      distance to the source in Mpc; **in all other cases this parameter should 
      be frozen to 1!**

_Note:_
  * KYRH (the black hole horizon), KYRIN (the disc inner edge), 
    KYRMS (the marginally stable orbit),
    KYLXLAMP (intrinsic luminosity of the primary source in 2-10keV) and 
    KYREFL (ratio of the reflection photon flux to the primary flux at the 
    observer) are added to the XSPEC internal
    switches. XSPEC xset command shows their current values.
  * Accuracy vs. speed trade off depends mainly on nrad and nphi.
  * This model includes a physical model of polarisation based on Rayleigh 
    scattering in single scattering approximation.


KYNrefionx
----------

Ionised reflection spectrum from black hole accretion disc in the lamp-post 
geometry. Black hole may be rotating (Kerr black hole), accretion disc is 
assumed to be Keplerian, geometrically thin and optically thick. The radial 
emissivity is given by the disc illumination from a point source located at some 
height on the system axis. The primary source is isotropic and emits 
power-law radiation. The disc ionisation state changes with radius and depends 
on the illumination patern as well as the radial density profile of the disc. 
Re-processing in the ionised disc is taken from
REFLION(X) tables, see Ross & Fabian 2005, MNRAS, 358, 211. Several limb 
darkening/brightening laws for emission directionality are implemented.
In this model only part of the disc may be set to be emitting 
radiation (sections defined in radius and azimuth). Obscuration by circular 
cloud is possible. Full relativistic ray-tracing code in vacuum was used for 
photon paths to compute the tables of transfer functions used in the model. Disc 
area below ISCO may be chosen to emit radiation, in which case material there is 
assumed to be freely falling and has the same energy and angular momentum as the 
matter orbiting at the ISCO. The model is based on KY package of models first 
presented in Dovciak M., Karas V. & Yaqoob T. 2004, ApJS, 153, 205-221 and 
later in Dovciak, M., Svoboda, J., Goosmann, R. W., et al.: 2014, in Proceedings 
of RAGtime 14-16: Workshops on black holes and neutron stars (Silesian 
University in Opava), [arXiv:1412.8627].

Definition of the parameters:

  * **par1  ... a/M**
    - black hole angular momentum (-1 &le; a/M &le; 1)
  * **par2  ... theta_o**
    - observer inclination in degrees (0&deg;-pole, 90&deg;-disc)
  * **par3  ... rin**
    - inner edge of non-zero disc emissivity (in GM/c^2 or in r~mso~)
  * **par4  ... ms**
    - switch for inner edge
    - 0: we integrate from inner edge = par3 
    - 1: if the inner edge of the disc is below marginally stable orbit (MSO) 
         then we integrate emission above MSO only
    - 2: we integrate from inner edge given in units of MSO, i.e. inner 
         edge = par3 &times; r~mso~ (the same applies for outer edge)
  * **par5  ... rout**
    - outer edge of non-zero disc emissivity (in GM/c^2 or in r~mso~)
  * **par6  ... phi**
    - lower azimuth of non-zero disc emissivity (degrees)
  * **par7  ... dphi**
    - (phi + dphi) is upper azimuth of non-zero disc emissivity 0&deg; &le; 
      dphi &le; 360&deg;
  * **par8  ... M/M8**
       - black hole mass in units of 10^8 solar masses
  * **par9  ... height**
    - height on the axis (measured from the center) at which the primary 
      source is located (GM/c^(2))
  * **par10 ... PhoIndex**
    - power-law energy index of the primary flux
  * **par11 ... L/L~Edd~**
    - dE/dt, the observed (if positive) or the intrinsic local (if negative) 
      primary isotropic flux in the X-ray energy range 2-10keV in units of 
      L~Edd~
  * **par12 ... Np:Nr**
    - ratio of the primary to the reflected normalization
    - 1: self-consistent model for isotropic primary source
    - 0: only reflection, primary source is hidden
    - if positive then L/L~Edd~ (par11) means the luminosity towards the 
      observer
    - if negative then L/L~Edd~ (par11) means the luminosity towards the disc
  * **par13 ... density/ionisation**
    - density profile normalization in 10^15 cm^(-3) if positive
    - ionisation profile normalisation if it is negative
    - this parameter cannot be zero
  * **par14 ... den_prof/ion_prof**
    - radial power-law density profile if par13 is positive
    - radial ionisation profile if par13 is negative
    - the radial profiles in both cases are given by abs(par13) &times; r^(par14)
  * **par15 ... abun**
    - Fe abundance (in solar abundance)
  * **par16 ... alpha**
    - position of the cloud centre in GM/c^2 in alpha coordinate (alpha being 
      the impact parameter in &phi;-direction, positive for approaching side 
      of the disc)
  * **par17 ... beta**
    - position of the cloud centre in GM/c^2 in beta coordinate (beta being 
      the impact parameter in &theta;-direction, positive in up direction, 
      i.e. above the disc)
  * **par18 ... rcloud**
    - radius of the obscuring cloud
    - the meaning of cloud is inverted for negative values of rcloud, i.e. 
      only the radiation transmitted through the cloud is computed
  * **par19 ... zshift**
    - > 0: the overall Doppler shift z; the distance of the source is computed 
           from the Hubble law, D = zc/H, with the Hubble constant 
           H = 70 km/s/Mpc; this distance is used to compute the correct 
           normalisation of the spectrum, thus norm parameter (the normalisation 
           of the model) should be frozen to unity
    - = 0: the overall Doppler shift is set to z = 0; the distance to the source 
           and thus the correct normalisation is defined by the normalisation 
           parameter
    - < 0: the negative value of the overall Doppler shift z used only to 
           compute the distance of the source from the Hubble law; this distance 
           is used to compute the correct normalisation of the spectrum, thus 
           the normalisation parameter should be frozen to unity; the overall 
           Doppler shift is then set to z = 0 (i.e. there is no shift of 
           spectrum with energy)
  * **par20 ... limb**
    - 0: for isotropic emission (flux ~ 1)
    - 1: for Laor's limb darkening (flux ~ 1+2.06&mu;)
    - 2: for Haardt's limb brightening (flux ~ ln (1+1/&mu;))
  * **par21 ... tab**
    - which REFLION table to use
    - 1: REFLION (the old one, lower cut-off energy at 1eV, not good for 
         PhoIndex > 2)
    - 2: REFLIONX (the newer one, lower cut-off energy at 100eV)
  * **par22 ... sw**
    - switch for the way how to compute the refl. spectra
    - 1: use the computed ionisation parameter, &xi;, for the interpolation 
         in REFLION, i.e. use proper total incident intensity with the 
         shifted cut-offs
    - 2: use the ionisation parameter, &xi;, correspondent to the computed 
         normalization of the incident flux, i.e. do not shift the cut-offs 
         when computing the total incident intensity
  * **par23 ... ntable**
    - table of relativistic transfer functions used in the model
      (defines FITS file with tables), 0 &le; ntable &le; 99, currently the 
      tables with ntable=80 are correct for this model
  * **par24 ... nrad**
    - number of grid points in radius
  * **par25 ... division**
    - type of division in radial integration
    - 0: equidistant radial grid (constant linear step)
    - 1: exponential radial grid (constant logarithmic step)
  * **par26 ... nphi**
    - number of grid points in azimuth
  * **par27 ... smooth**
    - whether to smooth the resulting spectrum 
    - 0: no smoothing
    - 1: simple smoothing
  * **par28 ... nthreads**
    - number of threads used for computations
  * **par29 ... norm**
    - if the overall Doppler shift zshift = 0, then norm = 1/D^2 where D is the 
      distance to the source in Mpc; **in all other cases this parameter should 
      be frozen to 1!**

_Note:_
  * KYRH (the black hole horizon), KYRIN (the disc inner edge), 
    KYRMS (the marginally stable orbit),
    KYLXLAMP (intrinsic luminosity of the primary source in 2-10keV) and 
    KYREFL (ratio of the reflection photon flux to the primary flux at the 
    observer) are added to the XSPEC internal
    switches. XSPEC xset command shows their current values.
  * Accuracy vs. speed trade off depends mainly on nrad and nphi.


KYNxillver
----------

Ionised reflection spectrum from black hole accretion disc in the lamp-post 
geometry. Black hole may be rotating (Kerr black hole), accretion disc is 
assumed to be Keplerian, geometrically thin and optically thick. The radial 
emissivity is given by the disc illumination from a point source located at some 
height on the system axis. The primary source is isotropic and emits 
power-law radiation. The disc ionisation state changes with radius and depends 
on the illumination patern as well as the radial density profile of the disc. 
Re-processing in the ionised disc is taken from
XILLVER tables, see Garcia & Kallman 2010, ApJ, 718, 695 and Garcia et al. 2013, 
ApJ, 768, 2. In this model only part of the disc may be set to be emitting 
radiation (sections defined in radius and azimuth). Obscuration by circular 
cloud is possible. Full relativistic ray-tracing code in vacuum was used for 
photon paths to compute the tables of transfer functions used in the model. Disc 
area below ISCO may be chosen to emit radiation, in which case material there is 
assumed to be freely falling and has the same energy and angular momentum as the 
matter orbiting at the ISCO. The model is based on KY package of models first 
presented in Dovciak M., Karas V. & Yaqoob T. 2004, ApJS, 153, 205-221 and 
later in Dovciak, M., Svoboda, J., Goosmann, R. W., et al.: 2014, in Proceedings 
of RAGtime 14-16: Workshops on black holes and neutron stars (Silesian 
University in Opava), [arXiv:1412.8627].

Definition of the parameters:

  * **par1  ... a/M**
    - black hole angular momentum (-1 &le; a/M &le; 1)
  * **par2  ... theta_o**
    - observer inclination in degrees (0&deg;-pole, 90&deg;-disc)
  * **par3  ... rin**
    - inner edge of non-zero disc emissivity (in GM/c^2 or in r~mso~)
  * **par4  ... ms**
    - switch for inner edge
    - 0: we integrate from inner edge = par3 
    - 1: if the inner edge of the disc is below marginally stable orbit (MSO) 
         then we integrate emission above MSO only
    - 2: we integrate from inner edge given in units of MSO, i.e. inner 
         edge = par3 &times; r~mso~ (the same applies for outer edge)
  * **par5  ... rout**
    - outer edge of non-zero disc emissivity (in GM/c^2 or in r~mso~)
  * **par6  ... phi**
    - lower azimuth of non-zero disc emissivity (degrees)
  * **par7  ... dphi**
    - (phi + dphi) is upper azimuth of non-zero disc emissivity 0&deg; &le; 
      dphi &le; 360&deg;
  * **par8  ... M/M8**
       - black hole mass in units of 10^8 solar masses
  * **par9  ... height**
    - height on the axis (measured from the center) at which the primary 
      source is located (GM/c^(2))
  * **par10 ... PhoIndex**
    - power-law energy index of the primary flux
  * **par11 ... L/L~Edd~**
    - dE/dt, the observed (if positive) or the intrinsic local (if negative) 
      primary isotropic flux in the X-ray energy range 2-10keV in units of 
      L~Edd~
  * **par12 ... Np:Nr**
    - ratio of the primary to the reflected normalization
    - 1: self-consistent model for isotropic primary source
    - 0: only reflection, primary source is hidden
    - if positive then L/L~Edd~ (par11) means the luminosity towards the 
      observer
    - if negative then L/L~Edd~ (par11) means the luminosity towards the disc
  * **par13 ... density/ionisation**
    - density profile normalization in 10^15 cm^(-3) if positive
    - ionisation profile normalisation if it is negative
    - this parameter cannot be zero
  * **par14 ... den_prof/ion_prof**
    - radial power-law density profile if par13 is positive
    - radial ionisation profile if par13 is negative
    - the radial profiles in both cases are given by abs(par13) &times; r^(par14)
  * **par15 ... abun**
    - Fe abundance (in solar abundance)
  * **par16 ... E_cut**
    - the observed (if positive) or intrinsic local at the source (if negative) 
      cut-off energy of the primary X-ray radiation
  * **par17 ... alpha**
    - position of the cloud centre in GM/c^2 in alpha coordinate (alpha being 
      the impact parameter in &phi;-direction, positive for approaching side 
      of the disc)
  * **par18 ... beta**
    - position of the cloud centre in GM/c^2 in beta coordinate (beta being 
      the impact parameter in &theta;-direction, positive in up direction, 
      i.e. above the disc)
  * **par19 ... rcloud**
    - radius of the obscuring cloud
    - the meaning of cloud is inverted for negative values of rcloud, i.e. 
      only the radiation transmitted through the cloud is computed
  * **par20 ... zshift**
    - > 0: the overall Doppler shift z; the distance of the source is computed 
           from the Hubble law, D = zc/H, with the Hubble constant 
           H = 70 km/s/Mpc; this distance is used to compute the correct 
           normalisation of the spectrum, thus norm parameter (the normalisation 
           of the model) should be frozen to unity
    - = 0: the overall Doppler shift is set to z = 0; the distance to the source 
           and thus the correct normalisation is defined by the normalisation 
           parameter
    - < 0: the negative value of the overall Doppler shift z used only to 
           compute the distance of the source from the Hubble law; this distance 
           is used to compute the correct normalisation of the spectrum, thus 
           the normalisation parameter should be frozen to unity; the overall 
           Doppler shift is then set to z = 0 (i.e. there is no shift of 
           spectrum with energy)
  * **par21 ... limb**
    - only used for angle averaged XILLVER tables
    - 0: for isotropic emission (flux ~ 1)
    - 1: for Laor's limb darkening (flux ~ 1+2.06&mu;)
    - 2: for Haardt's limb brightening (flux ~ ln (1+1/&mu;))
  * **par22 ... tab**
    - which XILLVER table to use
    - 1: xillver.fits, angle averaged with cut-off energy at 
         300 keV
    - 2: xillver-a.fits, angle dependent with cut-off energy at
         300 keV
    - 3: xillver-Ec.fits, angle averaged with free cut-off energy 
    - 4: xillver-a-Ec.fits, angle dependent with free cut-off 
         energy 
    - 5: xillver-a-Ec2.fits, angle dependent with free cut-off
         energy 
    - 6: xillver-a-Ec3.fits, angle dependent with free cut-off
         energy 
    - 7: xillver-a-Ec4.fits, angle dependent with free cut-off
         energy 
    - 8: xillver-a-Ec5.fits, angle dependent with free cut-off
         energy 
  * **par23 ... ntable**
    - table of relativistic transfer functions used in the model
      (defines FITS file with tables), 0 &le; ntable &le; 99, currently the 
      tables with ntable=80 are correct for this model
  * **par24 ... nrad**
    - number of grid points in radius
  * **par25 ... division**
    - type of division in radial integration
    - 0: equidistant radial grid (constant linear step)
    - 1: exponential radial grid (constant logarithmic step)
  * **par26 ... nphi**
    - number of grid points in azimuth
  * **par27 ... smooth**
    - whether to smooth the resulting spectrum 
    - 0: no smoothing
    - 1: simple smoothing
  * **par28 ... nthreads**
    - number of threads used for computations
  * **par29 ... norm**
    - if the overall Doppler shift zshift = 0, then norm = 1/D^2 where D is the 
      distance to the source in Mpc; **in all other cases this parameter should 
      be frozen to 1!**

_Note:_
  * KYRH (the black hole horizon), KYRIN (the disc inner edge), 
    KYRMS (the marginally stable orbit),
    KYLXLAMP (intrinsic luminosity of the primary source in 2-10keV), 
    KYREFL (ratio of the reflection photon flux to the primary flux at the 
    observer), KYXIIN (ionisation parameter at the inner disc edge) and 
    KYXIOUT (ionisation parameter at the outer disc edge) are added to the XSPEC 
    internal switches. XSPEC xset command shows their current values.
  * Accuracy vs. speed trade off depends mainly on nrad and nphi.


KYNhrefl
--------

Compton reflection model for a black hole accretion disc. Black hole may be 
rotating (Kerr black hole), accretion disc is assumed to be Keplerian, 
geometrically thin and optically thick.
This model is based on an existing multiplicative HREFL model in combination
with the POWERLAW model. Local emission is the same as in HREFL*POWERLAW with
the parameters thetamin = 0 and thetamax = 90 (i.e. it is assumed that the
disc is illuminated from all directions isotropically) and with a broken
power-law radial dependence added. This model can be interpreted as a
Compton-reflection model for which the source of primary irradiation is near
above the disc, in contrast to the lamp-post scheme with the source on the
axis. The approximations for Compton reflection used in HREFL (and therefore
also in this model) are valid below ~15keV (in the disc rest frame).
Only part of the disc may be set to be emitting 
radiation (sections defined in radius and azimuth). Obscuration by circular 
cloud is possible. Full relativistic ray-tracing code in vacuum was used for 
photon paths to compute the tables of transfer functions used in the model. Disc 
area below ISCO may be chosen to emit radiation, in which case material there is 
assumed to be freely falling and has the same energy and angular momentum as the 
matter orbiting at the ISCO. The model is based on its first version presented 
in Dovciak M., Karas V. & Yaqoob T. 2004, ApJS, 153, 205-221. 

Definition of the parameters:

  * **par1  ... a/M**
    - black hole angular momentum (-1 &le; a/M &le; 1)
  * **par2  ... theta_o**
    - observer inclination in degrees (0&deg;-pole, 90&deg;-disc)
  * **par3  ... rin**
    - inner edge of non-zero disc emissivity (in GM/c^2 or in r~mso~)
  * **par4  ... ms**
    - switch for inner edge
    - 0: we integrate from inner edge = par3 
    - 1: if the inner edge of the disc is below marginally stable orbit (MSO) 
         then we integrate emission above MSO only
    - 2: we integrate from inner edge given in units of MSO, i.e. inner 
         edge = par3 &times; r~mso~ (the same applies for outer edge)
  * **par5  ... rout**
    - outer edge of non-zero disc emissivity (in GM/c^2 or in r~mso~)
  * **par6  ... phi**
    - lower azimuth of non-zero disc emissivity (degrees)
  * **par7  ... dphi**
    - (phi + dphi) is upper azimuth of non-zero disc emissivity 0&deg; &le; 
      dphi &le; 360&deg;
  * **par8  ... PhoIndex**
    - power-law energy index of the primary flux
  * **par9  ... q_out** 
    - power-law index for radial dependence of emissivity for outer region, 
      scales as r^(-q_out)
  * **par10 ... q_in**
    - power-law index for radial dependence of emissivity for inner region, 
      scales as rb^(q_in-q_out) &times; r^(-q_in)
  * **par11 ... rb**
    - boundary between the region with power-law index q_out and q_in
    - if > 0 then the boundary is in units of MSO, i.e. boundary = rb &times; 
      r~mso~
    - if <= 0 then the boundary is equal to -rb+r_horizon where rb is in GM/c^2
  * **par12 ... jump**
    - ratio of local flux in inner region to local flux in outer region at 
      boundary radius defined by rb
  * **par13 ... Feabun**
    - iron abundance relative to Solar
  * **par14 ... FeKedge**
    - iron K-edge energy
  * **par15 ... Escfrac**
    - normalization of the original powerlaw emission
  * **par16 ... covfac**
    - normalization of the reflected emission
  * **par17 ... alpha**
    - position of the cloud centre in GM/c^2 in alpha coordinate (alpha being 
      the impact parameter in &phi;-direction, positive for approaching side 
      of the disc)
  * **par18 ... beta**
    - position of the cloud centre in GM/c^2 in beta coordinate (beta being 
      the impact parameter in &theta;-direction, positive in up direction, 
      i.e. above the disc)
  * **par19 ... rcloud**
    - radius of the obscuring cloud
    - the meaning of cloud is inverted for negative values of rcloud, i.e. 
      only the radiation transmitted through the cloud is computed
  * **par20 ... zshift**
    - overall Doppler shift
  * **par21 ... ntable**
    - table of relativistic transfer functions used in the model
      (defines FITS file with tables), 0 &le; ntable &le; 99, currently the 
      tables with ntable=80 are correct for this model
  * **par22 ... nrad**
    - number of grid points in radius
  * **par23 ... division**
    - type of division in radial integration
    - 0: equidistant radial grid (constant linear step)
    - 1: exponential radial grid (constant logarithmic step)
  * **par24 ... nphi**
    - number of grid points in azimuth
  * **par25 ... smooth**
    - whether to smooth the resulting spectrum 
    - 0: no smoothing
    - 1: simple smoothing
  * **par26 ... Stokes** 
    - definition of output
    - 0: photon number density flux per bin (Stokes parameter I devided by 
         energy)
    - 1: Stokes parameter Q devided by energy
    - 2: Stokes parameter U devided by energy
    - 3: Stokes parameter V devided by energy
    - 4: degree of polarisation
    - 5: linear polarisation angle &psi; = 0.5 atan(U/Q)
    - 6: circular polarisation angle &beta; = 0.5 asin(V/sqrt(Q^(2)+U^(2)+V^(2)))
  * **par27 ... nthreads**
    - number of threads used for computations
  * **par28 ... norm**

_Note:_
  * KYRH (the black hole horizon), KYRIN (the disc inner edge) and 
    KYRMS (the marginally stable orbit) are added to the XSPEC internal
    switches. Use xset command to show their current values.
  * Accuracy vs. speed trade off depends mainly on nrad and nphi.
  * In this model it is assumed that the local emission is completely linearly 
    polarised in the direction perpendicular to the disc.


KYNphebb
--------

Phenomenological thermal emission from a black hole accretion disc. Black hole 
may be rotating (Kerr black hole), accretion disc is assumed to be 
Keplerian, geometrically thin and optically thick. The radial dependence of the 
disc black-body temperature is a simple powerlaw. The local flux is defined as
flux ~ E^2/(exp(E/kT)-1), where T=Tin*(r/rin)^(-BBindex).
Only part of the disc may be set to be emitting radiation (sections defined in 
radius and azimuth). Obscuration by circular cloud is possible. Full 
relativistic ray-tracing code in vacuum was used for photon paths to compute the 
tables of transfer functions used in the model. Disc area below ISCO may be 
chosen to emit radiation, in which case material there is assumed to be freely 
falling and has the same energy and angular momentum as the matter orbiting at 
the ISCO. The model is based on KY package of models first 
presented in Dovciak M., Karas V. & Yaqoob T. 2004, ApJS, 153, 205-221.

Definition of the parameters:

  * **par1  ... a/M**
    - black hole angular momentum (-1 &le; a/M &le; 1)
  * **par2  ... theta_o**
    - observer inclination in degrees (0&deg;-pole, 90&deg;-disc)
  * **par3  ... rin**
    - inner edge of non-zero disc emissivity (in GM/c^2 or in r~mso~)
  * **par4  ... ms**
    - switch for inner edge
    - 0: we integrate from inner edge = par3 
    - 1: if the inner edge of the disc is below marginally stable orbit (MSO) 
         then we integrate emission above MSO only
    - 2: we integrate from inner edge given in units of MSO, i.e. inner 
         edge = par3 &times; r~mso~ (the same applies for outer edge)
  * **par5  ... rout**
    - outer edge of non-zero disc emissivity (in GM/c^2 or in r~mso~)
  * **par6  ... phi**
    - lower azimuth of non-zero disc emissivity (degrees)
  * **par7  ... dphi**
    - (phi + dphi) is upper azimuth of non-zero disc emissivity 0&deg; &le; 
      dphi &le; 360&deg;
  * **par8  ... Tin**
    - temperature in keV at the inner edge of the disc
  * **par9  ... BBindex**
    - radial power-law index for radial dependence of the black-body temperature
  * **par10 ... alpha**
    - position of the cloud centre in GM/c^2 in alpha coordinate (alpha being 
      the impact parameter in &phi;-direction, positive for approaching side 
      of the disc)
  * **par11 ... beta**
    - position of the cloud centre in GM/c^2 in beta coordinate (beta being 
      the impact parameter in &theta;-direction, positive in up direction, 
      i.e. above the disc)
  * **par12 ... rcloud**
    - radius of the obscuring cloud
    - the meaning of cloud is inverted for negative values of rcloud, i.e. 
      only the radiation transmitted through the cloud is computed
  * **par13 ... zshift**
    - overall Doppler shift
  * **par14 ... ntable**
    - table of relativistic transfer functions used in the model
      (defines FITS file with tables), 0 &le; ntable &le; 99, currently the 
      tables with ntable=80 are correct for this model
  * **par15 ... nrad**
    - number of grid points in radius
  * **par16 ... division**
    - type of division in radial integration
    - 0: equidistant radial grid (constant linear step)
    - 1: exponential radial grid (constant logarithmic step)
  * **par17 ... nphi**
    - number of grid points in azimuth
  * **par18 ... smooth**
    - whether to smooth the resulting spectrum 
    - 0: no smoothing
    - 1: simple smoothing
  * **par19 ... nthreads**
    - number of threads used for computations
  * **par20 ... norm**
    - equals to 1/D^2 where D is a source distance in 10kpc

_Note:_
  * KYRH (the black hole horizon), KYRIN (the disc inner edge) and 
    KYRMS (the marginally stable orbit) are added to the XSPEC internal
    switches. Use xset command to show their current values.
  * Accuracy vs. speed trade off depends mainly on nrad and nphi.


KYNbb
-----

Thermal emission from a black hole accretion disc. Black hole may be rotating 
(Kerr black hole), accretion disc is assumed to be Keplerian, geometrically thin 
and optically thick. The local flux is a black body radiation with the 
Novikov-Thorne radial temperature profile.
Only part of the disc may be set to be emitting radiation (sections defined in 
radius and azimuth). Obscuration by circular cloud is possible. Full 
relativistic ray-tracing code in vacuum was used for photon paths to compute the 
tables of transfer functions used in the model. Disc area below ISCO may be 
chosen to emit radiation, in which case material there is assumed to be freely 
falling and has the same energy and angular momentum as the matter orbiting at 
the ISCO. The model is based on KY package of models first 
presented in Dovciak M., Karas V. & Yaqoob T. 2004, ApJS, 153, 205-221.

Definition of the parameters:

  * **par1  ... a/M**
    - black hole angular momentum (-1 &le; a/M &le; 1)
  * **par2  ... theta_o**
    - observer inclination in degrees (0&deg;-pole, 90&deg;-disc)
  * **par3  ... rin**
    - inner edge of non-zero disc emissivity (in GM/c^2 or in r~mso~)
  * **par4  ... ms**
    - switch for inner edge
    - 0: we integrate from inner edge = par3 
    - 1: if the inner edge of the disc is below marginally stable orbit (MSO) 
         then we integrate emission above MSO only
    - 2: we integrate from inner edge given in units of MSO, i.e. inner 
         edge = par3 &times; r~mso~ (the same applies for outer edge)
  * **par5  ... rout**
    - outer edge of non-zero disc emissivity (in GM/c^2 or in r~mso~)
  * **par6  ... phi**
    - lower azimuth of non-zero disc emissivity (degrees)
  * **par7  ... dphi**
    - (phi + dphi) is upper azimuth of non-zero disc emissivity 0&deg; &le; 
      dphi &le; 360&deg;
  * **par8  ... BHmass**
    - the black hole mass in units of Solar mass
  * **par9  ... arate**
    - accretion rate in units of Solar mass per Julian year (365.25days)
  * **par10 ... f_col**
    - spectral hardening factor
  * **par11 ... alpha**
    - position of the cloud centre in GM/c^2 in alpha coordinate (alpha being 
      the impact parameter in &phi;-direction, positive for approaching side 
      of the disc)
  * **par12 ... beta**
    - position of the cloud centre in GM/c^2 in beta coordinate (beta being 
      the impact parameter in &theta;-direction, positive in up direction, 
      i.e. above the disc)
  * **par13 ... rcloud**
    - radius of the obscuring cloud
    - the meaning of cloud is inverted for negative values of rcloud, i.e. 
      only the radiation transmitted through the cloud is computed
  * **par14 ... zshift**
    - overall Doppler shift
  * **par15 ... ntable**
    - table of relativistic transfer functions used in the model
      (defines FITS file with tables), 0 &le; ntable &le; 99, currently the 
      tables with ntable=80 are correct for this model
  * **par16 ... nrad**
    - number of grid points in radius
  * **par17 ... division**
    - type of division in radial integration
    - 0: equidistant radial grid (constant linear step)
    - 1: exponential radial grid (constant logarithmic step)
  * **par18 ... nphi**
    - number of grid points in azimuth
  * **par19 ... smooth**
    - whether to smooth the resulting spectrum 
    - 0: no smoothing
    - 1: simple smoothing
  * **par20 ... Stokes** 
    - definition of output
    - 0: photon number density flux per bin (Stokes parameter I devided by 
         energy)
    - 1: Stokes parameter Q devided by energy
    - 2: Stokes parameter U devided by energy
    - 3: Stokes parameter V devided by energy
    - 4: degree of polarisation
    - 5: linear polarisation angle &psi; = 0.5 atan(U/Q)
    - 6: circular polarisation angle &beta; = 0.5 asin(V/sqrt(Q^(2)+U^(2)+V^(2)))
  * **par21 ... nthreads**
    - number of threads used for computations
  * **par22 ... norm**
    - equals to 1/D^2 where D is a source distance in 10kpc

_Note:_
  * KYRH (the black hole horizon), KYRIN (the disc inner edge) and 
    KYRMS (the marginally stable orbit) are added to the XSPEC internal
    switches. Use xset command to show their current values.
  * Accuracy vs. speed trade off depends mainly on nrad and nphi.
  * In this model it is assumed that the local emission is completely linearly 
    polarised in the direction perpendicular to the disc.
