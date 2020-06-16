<img src="logo.png" width="20%">

# Paint the Atmospheres of Rotating Stars (PARS)

This software quickly computes magnitudes and spectra of rotating stellar models. The models incorporate **a Roche model of mass distribution** (in which all mass is at the center of the star), **solid body rotation**, and **collinearity of effective gravity and energy flux**.

## Getting Started

The following instructions describe how to install PARS on macOS 10.15.4. These may have to be modified for other operating systems.

### Prerequisites

* git
* python 3

### Installing

Go to the directory where you want to install PARS, clone it into that directory, and go to the software's top directory.
```
cd <directory name>
git clone https://github.com/mlipatov/paint_atmospheres
cd paint_atmospheres
```

Install PARS.
```
pip install .
```

Place a file with limb darkening information from atmosphere models in the data directory.
```
cp ~/im01k2.pck ./data
```

Print the command line syntax for the executables.
```
calc_limbdark -h
calc_star -h
calc_spectra -h
```

## Deployment

The following instructions describe how to access the functionality of PARS in several different ways.

### Compute a Spectrum

Compute fits of intensity versus surface inclination from the limb darkening information.
```
calc_limbdark data/im01k2.pck data/limbdark_m01.pkl 0.1 0.4
```

Perform the inclination-independent computations for a stellar model.
```
calc_star 'data/limbdark_m01.pkl' 'data/vega.pkl' 0.6151 40.346 2.165 2.815 2.3694e19 100
```

Perform the inclination-dependent computations.
```
calc_spectra 'data/vega.pkl' 'data/vega/' -i 0.088418
```

Look at the resulting spectrum
```
cat data/vega/*.txt | more
```

### Compute magnitudes

Create files with filter transmission curves, such as [](other_file.md)

### Break down into end to end tests

Explain what these tests test and why

```
Give an example
```

### And coding style tests

Explain what these tests test and why

```
Give an example
```


Features:
	a closed-form expression for the azimuthal integral,
	a high-order approximation of the longitudinal integral,
	a precise calculation of surface effective temperature.

Executables:
	calc_limbdark computes fits of intensity w.r.t surface inclination
	calc_star performs inclination-independent computations
	calc_spectra computes spectra for multiple inclinations
	<command> -h provides usage instructions

Scripts:
	Located in paint_atmospheres/pa/usr
	Create Figures 3-10 in [LB] and other figures
	Each contains instructions in the comments at the top of the file 

## Authors

* **Timothy D. Brandt**
* [Mikhail Lipatov](https://github.com/mlipatov/)

References: 
	Lipatov M & Brandt TD, Submitted to ApJ [LB]
	Espinosa Lara F & Rieutord M (2011), A&A, 533, A43
	Castelli F & Kurucz RL (2004), arXiv:astro-ph/0405087 
	http://mathworld.wolfram.com/CubicFormula.html
	Wikipedia:Bilinear interpolation:Algorithm
	Wikipedia:Newton's method
	Press WH et al, Numerical Recipes, 3rd ed. (2007) 
