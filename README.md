# 1D Spectral Simulator for VLT-MOONS
Produces a 1D (optimally extracted) MOONS spectrum for a set of observing conditions, starting from a provided template. The sky is modeled using the ESO skycalc model or the user can provide the Sky template. The code produces a SNR calculation and creates a FITS file as output with the resulting simulated spectrum (sky corrected or not). An output plot is also produced with a summary of the simulation performance.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine.

### Prerequisites

The following Python packages need to be installed for the simulator to run. This can be done for example using pip:

```
pip install PyAstronomy numpy astropy sys scipy matplotlib math warnings
```

### Installing

Download this package. This includes all the necessary files and subdirectories needed to run the simulator:

```
moons_sim.py 
Inst_setup/ which includes all the throughput information for Instrument, detectors, and telescope.
Skymodel/ which inscludes all the ESO sky model templates for different airmass values at average Paranal conditions.
```

## Running the simulator

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

To run a LR H-band simulation of an observation of 3600sec (in 200sx18 DITxNDIT) at 0.8 seeing and airmass 1.2, correcting atmospheric difraction at 1.5um, starting from a template named stellar_template.fits, using the ESO sky model, applying sky subtraction with 2% sky residual, and telluric correction, to produce output file named output_spec.fits

```
moons_sim stellar_template.fits output_spec.fits LR H 18 200 18 0.8 1.2 1.5 2.0 eso_skycalc 1
```

Same simulation as above, without applying sky subtraction or telluric correction:

```
moons_sim stellar_template.fits output_spec.fits LR H 18 200 18 0.8 1.2 1.5 -1 eso_skycalc -1
```

## Authors

* **Oscar A. Gonzalez** - *UKATC* -

## Acknowledgments

* Makes use of the Cerro Paranal Advanced Sky Model, which was developed in particular to be used in the ESO Exposure Time Calculators, by a team of astronomers at the Institute for Astro- and Particle Physics at the University of Innsbruck, as part of an Austrian in-kind contribution to ESO. References: Noll et al. (2012, A&A 543, A92) and Jones et al. (2013, A&A 560, A91).
