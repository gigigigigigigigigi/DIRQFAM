# Spherical test
We select a generic nucleus <sup>84</sup>Zr, constrain its
shape to spherical by constraining its quadrupole
moment to zero value and perform QFAM calculations
for fixed J and variable 0 <= K <= J.


Due to the Wigner-Eckart theorem, spherical nuclei
should exhibit the strength function response invariant to multipolarity
K for fixed J. In this spherical test, we verify
this assertion numerically.


The demonstrated agreement within 7 most significant digits in
strength response function is satisfying, considering all the
numerical integrations being performed. One can obtain even higher
level of agreement if the induced Coulomb interaction is ignored since
it requires numerical integration of highly oscillating function.



# Note
In this calculation, we set Broyden's iteration tolerance (variable "tol" in iter_fam.f)
to a slightly lower value in order to achieve better
agreement. Default value in the original <code>DIRQFAM</code> code is 1.e-5,
but in this test we lowered it to 1.e-8.


