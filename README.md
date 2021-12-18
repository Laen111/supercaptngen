# supercaptngen

Super Capt'n General: fork of [Capt'n General](https://github.com/aaronvincent/captngen), designed to handle supernovae boosting dark matter via scattering to above galactic escape velocities.

If you use this code, you can cite **_NOT USED YET_** where it was first deployed.

Run `make` to compile the library on its own.

If you want the test file as well, run `make supertest.x` which compiles the test script in *supermain.f90*.
The test main file shows how to work with Super Capt'n.

### TODO:

- Figure out why the $c_{5}^{0}$, $c_{7}^{0}$, $c_{8}^{0}$, $c_{13}^{0}$, and $c_{14}^{0}$ couplings have negative values (something to do with $(v^\perp_T)^2 = w^2 - \left(\frac{q}{2\mu_T}\right)^2$)
- Refine user interface to Super Capt'n
- Double check the units (especially on the energy of the supernova [ergs -> GeV?] and the $V_\text{shock}$ and $R_\text{shock}$ functions)
