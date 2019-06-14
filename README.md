# 1D Spectral Simulator for VLT-MOONS
Produces a 1D (optimally extracted) MOONS spectrum for a set of observing conditions, starting from a provided template. The sky is modeled using the ESO skycalc model or the user can provide the Sky template. The code produces a SNR calculation and creates a FITS file as output with the resulting simulated spectrum (sky corrected or not). An output plot is also produced with a summary of the simulation performance.

## Getting Started

These instructions will get you a copy of the simulator up and running on your local machine.

### Prerequisites

The simulator code is written in Python 3. The following dependencies need to be installed for the simulator to run properly. This can be done, for example, using pip:

```
pip install PyAstronomy numpy astropy sys scipy matplotlib math warnings
```

## Running the simulator

This includes all the necessary files and subdirectories needed to run the simulator:

```
moons_sim.py 
Inst_setup/ which includes all the throughput information for Instrument, detectors, and telescope.
Skymodel/ which inscludes all the ESO sky model templates for different airmass values at average Paranal conditions.
Example/ containing a FITS file with an example input spectrum to run the simulator (see command line example below).
```

The simulator can be run from the command line:

```
python moons_sim.py arg1 arg2 arg3 ...
```

Input arguments must be provided in the following specific order:

```
1: Template name (and local path)
2: Output fits spectrum name (and local path)
3: LR or HR mode for MOONS
4: RI, YJ, or H band for MOONS
5: Magnitude of the source (AB)
6: DIT (s)
7: NDIT
8: Seeing
9: Airmass
10: Wavelength at which perform atmospheric diffraction correction (in um)
11: Sky residual percentage. If -1 then sky subtraction is NOT performed.
12: ESO_skycalc for ESO skycalc model included in the bundle or Name (and local path) to use your own model
13: value=1 to perform telluric correction; -1 to switch off telluric correction
```

### Example execution

To run a LR H-band simulation of an observation of 3600sec (in 200sx18 DITxNDIT) at 0.8 seeing and airmass 1.2, correcting atmospheric difraction at 1.5um, starting from the example template in Example/input_stellar_template.fits, using the ESO sky model, applying sky subtraction with 2% sky residual, and telluric correction, to produce output file named output_spec.fits

```
moons_sim Example/input_stellar_template.fits output_spec.fits LR H 18 200 18 0.8 1.2 1.5 2.0 eso_skycalc 1
```

Same simulation as above, without applying sky subtraction or telluric correction:

```
moons_sim Example/input_stellar_template.fits output_spec.fits LR H 18 200 18 0.8 1.2 1.5 -1 eso_skycalc -1
```

## Input templates

The simulator requires input templates for the science source input (and for the sky, when user-defined) in specific units and expects the corresponding header keywords.

### Input Source templates

The package contains an example of the input source spectrum in the Example/ folder. 

```
Single extension fits file with
hdu[0] = source flux in erg/s/cm2/A
```

The mandatory header keywords (and an example value) are:
```
hdu.header['CODE'] = ('Cigale', 'Simulation code')
hdu.header['MTYPE'] = ('Galaxy', 'Star, Galaxy, or Quasar')
hdu.header['MNAME'] = ('M11_NoConvolved_model', 'Model name')
hdu.header['CRVAL'] = CRVAL
hdu.header['CRDEL'] = CRDEL
hdu.header['CRPIX'] = 1
hdu.header['R'] = 20000 Resolving power 
hdu.header['Sampling'] = 2. #Pixels per element of resolution
hdu.header['TUNIT1'] = 'A'
hdu.header['TUNIT2'] = 'erg/s/cm2/A'
```

### Input Sky template (OPTIONAL)

If the sky spectrum is provided by the user (so the input argument is different from "ESO_skycalc" and corresponds to the path to an user defined fits-file), then this file should be:

```
Two-extension fits file with:
hdu[0]=SKY emission in erg/s/cm2/A/arcsec2
hdu[1]=Fraction of atmopheric transmission
```

The mandatory header keywords (and an example value) in that FITS-file are:
```
hdu.header['CRVAL'] = CRVAL
hdu.header['CRDEL'] = CRDEL
hdu.header['CRPIX'] = 1
hdu.header['R'] = 100000 Resolving power 
hdu.header['Sampling'] = 2. #Pixels per element of resolution
hdu.header['TUNIT1'] = 'A'
hdu.header['TUNIT2'] = 'erg/s/cm2/A/'
```


## Authors

* **Oscar A. Gonzalez** - *UKATC* -

## Acknowledgments

* Makes use of the Cerro Paranal Advanced Sky Model, which was developed in particular to be used in the ESO Exposure Time Calculators, by a team of astronomers at the Institute for Astro- and Particle Physics at the University of Innsbruck, as part of an Austrian in-kind contribution to ESO. References: Noll et al. (2012, A&A 543, A92) and Jones et al. (2013, A&A 560, A91).
