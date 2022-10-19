# Super Capt'n General

Super Capt'n General: fork of [Capt'n General](https://github.com/aaronvincent/captngen), designed to handle supernovae boosting dark matter via elastic scattering to above galactic escape velocities.

If you use this code, you can cite [arXiv:2210.09448](https://arxiv.org/abs/2210.09448) where it was first deployed.

Run `make` to compile the library on its own.

If you want the test file as well, run `make supertest.x` which compiles the test script in `supermain.f90`.
The test main file shows how to work with Super Capt'n.

The program can be run to get a velocity spectrum using the `runsuper.f90` main file. To do this, edit the parameters you desire at the top of the file, then build using the `make getSpectrum.x` command, and finally run the new `getSpectrum.x` executable. The created output file includes a helpful tip in the header on how to read it into a python numpy array alongside a record of all parameters used in the run.

## Main File Structure

The structure of the `runsuper.f90` main file is as follows:

1. First the option to change the output file name and dark matter velocity range/resolution.
1. Next, the option to change dark matter parameters such as the mass, spin, and couplings. Note that these values are input as $\text{GeV}^{-2}$, not $\text{cm}^2$.
You can use $c_i^\tau \approx \sqrt{4\pi\sigma}\,/\,(\mu\hbar c)$ to get a rough translation from cross section to coupling.
Also, any combination of couplings can be run at once, including the isoscalar and isovector ones.
1. Finally, the option to change the supernova parameters, including the local dark matter density, the ejecta mass, the local Interstellar Medium (ISM) density, the distance to the supernova, the total energy released by the supernova, the age of the supernova, the supernova type (either 1 for a TypeIa, 2 for TypeII), the pre-supernova stellar wind, and the pre-supernova stellar mass loss rate.
The values we used can be found in [arXiv:2210.09448](https://arxiv.org/abs/2210.09448).
1. The rest of the script calls and runs Super Capt'n General.
Feel free to modify it to suit your needs.

## Modifying `supercap.f90`

To facilitate modification of the code that performs the velocity calculation, it is summarised here.

### `supercaptn_init`:

The first function call is made to `supercaptn_init`.
This function initialises some data for the actual calculations.
It loads in all input values and converts to internal _natural units_.
It also loads the choice of nova into a shared variable that can be read later.
It then loads all the constants of the _W functions_ ($W_{m}^{\tau\tau^\prime}(y)$) (see [arXiv:1501.03729](https://arxiv.org/abs/1501.03729)) into an array to be quickly read out during calculation.
It then creates a zeroed out $14\times 2$ array to house the couplings.
Finally it creates an array of parameters that contains the conversion between $y$ and $q^2$, one for each isotope in the calculation (again [arXiv:1501.03729](https://arxiv.org/abs/1501.03729)).

### `populate_array_super`:

The next call is to `populate_array_super`, which can populate the array of coupling parameters.
This converts from the coupling number and isospin to the internal array indices and does error checking on the numbers chosen.

### `supercaptn`:

Finally `supercaptn` is called, which does the calculation to find the number of dark matter particles scattered for a given dark matter mass and dark matter velocity at Earth.
It first loads in the dark matter mass, velocity, and spin.
It then finds the time the dark matter scattered, and breaks out if the dark matter would have scattered before the explosion.
Next it zeros out an array of prefactors that correspond to factors that appear as $A_{i,q_\text{pow},w_\text{pow}} q^{2(q_\text{pow}-1)} w^{2(w_\text{pow}-1)}$, where $q_\text{pow}$ runs $1\rightarrow 11$ and $w_\text{pow}$ runs $1\rightarrow 2$.
Then the isotopes, _R_ and _W_ function types, isospin indices, _W_ function terms, and _R_ function terms are iterated over, constructing the prefactors.
Once the prefactors are computed, they are iterated over to calculate the scattering (see [arXiv:1501.03729](https://arxiv.org/abs/1501.03729) for details of the cross section, where here we have factorised the sum into powers of $q^2$ and $w^2$).
Final factors are multiplied to get the correct result, and it is converted to output in $[(\text{cm} \ \text{s}^{-1})^{-1} \cdot (\text{s})^{-1} \cdot (\text{cm}^2)^{-1}]$

### `supermod`:

This is a module that contains useful shared parameters and functions used in the calculation.
It has parameters like physical constants, isotope information, and initialisation of shared parameters and arrays.
It contains the calculation for the shockwave expansion radius and velocity.
It also contains generalised form factor integrals (GFFIs), which are from [arXiv:1504.04378](https://arxiv.org/abs/1504.04378).

## TODO:

- In future: change $V_\text{shock}\left( t \right)$, or include a varying $V_\text{shock}\left( x, t \right)$ with position in the shockwave
- Double check the uses of $V_{\rm Shock}$ and $R_{\rm Shock}$ in the interpolation loop to see if it should be the maximum value instead
- Double check that the opacity condition still behaves in the variable vshock version
