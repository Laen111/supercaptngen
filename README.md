# supercaptngen

Super Capt'n General: fork of [Capt'n General](https://github.com/aaronvincent/captngen), designed to handle supernovae boosting dark matter via scattering to above galactic escape velocities.

If you use this code, you can cite **_NOT USED YET_** where it was first deployed.

Run `make` to compile the library on its own.

If you want the test file as well, run `make supertest.x` which compiles the test script in _supermain.f90_.
The test main file shows how to work with Super Capt'n.

The program can be run to get a velocity spectrum using the `runsuper.f90` main file. To do this, edit the parameters you desire at the top of the file, then build using the `make getSpectrum.x` command, and finally run the new `getSpectrum.x` executable. The create output file includes a helpful tip in the header on how to read it into a python numpy array alongside a record of all parameters used in the run.

## TODO

- In future: change $V_\text{shock}\left( t \right)$, or include a varying $V_\text{shock}\left( x, t \right)$ with position in the shockwave
- Double check the uses of $V_{\rm Shock}$ and $R_{\rm Shock}$ in the interpolation loop to see if it should be the maximum value instead
