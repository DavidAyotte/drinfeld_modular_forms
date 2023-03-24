# Drinfeld Modular Forms

This SageMath package provides an implementation for computing with Drinfeld modular forms for the full modular group.

## Installation

This package has been tested on SageMath version 9.8 and higher. It is
not guaranteed to work on previous versions.

### Install from PyPI

The easiest way to install this package is via PyPI. You simply have to run SageMath first and then type the following command

`sage: pip install drinfeld-modular-forms`

### Install from source code

You can also install this package by cloning the source code from the [Github repo](https://github.com/DavidAyotte/drinfeld_modular_forms).

Next, you have to run `make install` inside the project's folder. You can also run the following command:

`sage -pip install --upgrade --no-index -v .`

If there is any changes to the current repo, you will then simply need to pull the changes and run the above command again.

## Usage

After running SageMath, you can import the functionalities of this package by typing the following command:

`sage: from drinfeld_modular_forms import *`

## Documentation

The documentation is available at this address:

https://davidayotte.github.io/drinfeld_modular_forms

## Examples

One may create the ring of Drinfeld modular forms:

```
    sage: from drinfeld_modular_forms import DrinfeldModularFormsRing
    sage: A = GF(3)['T']; K = Frac(A); T = K.gen()
    sage: M = DrinfeldModularFormsRing(K, 2)
    sage: M.ngens()  # number of generators
    2
```

The elements of this ring are viewed as multivariate polynomials in a choice of generators for the ring. The current implemented generators are the coefficient forms of a universal Drinfeld module over the Drinfeld period domain (see theorem 17.5 in \[1\]). In the computation below, the forms `g1` and `g2` corresponds to the weight `q - 1` Eisenstein series and the Drinfeld modular discriminant of weight `q^2 - 1` respectively.
```
    sage: M.inject_variables()
    Defining g1, g2
    sage: F = (g1 + g2)*g1; F
    g1*g2 + g1^2
```
Note that elements formed with polynomial relations `g1` and `g2` may not be homogeneous in the weight and may not define a Drinfeld modular form. We will call elements of this ring *graded Drinfeld modular forms*.

In the case of rank 2, one can compute the expansion at infinity of any graded form:

```
    sage: g1.expansion()
    1 + ((2*T^3+T)*t^2) + O(t^7)
    sage: g2.exansion()
    t^2 + 2*t^6 + O(t^8)
    sage: ((g1 + g2)*g2).expansion()
    1 + ((T^3+2*T+1)*t^2) + ((T^6+T^4+2*T^3+T^2+T)*t^4) + 2*t^6 + O(t^7)
```
This is achieved via the `A`-expansion theory developed by López-Petrov in \[3\] and \[4\]. We note that the returned expansion is a lazy power series. This means that it will compute on demands any coefficient up to any precision:
```
    sage: g2[600]  # 600-th coefficient
    T^297 + 2*T^279 + T^273 + T^271 + T^261 + 2*T^253 + T^249 + 2*T^243 + 2*T^171 + T^163 + T^153 + 2*T^147 + 2*T^145 + T^139 + T^135 + T^129 + 2*T^123 + 2*T^121 + T^117 + T^115 + T^111 + 2*T^109 + T^105 + 2*T^99 + 2*T^97 + T^93 + T^91 + T^87 + 2*T^85 + T^81 + 2*T^75 + T^69 + T^67 + T^63 + 2*T^61 + 2*T^51 + 2*T^45 + T^43 + T^39 + T^29 + T^27 + 2*T^21 + T^19 + T^13 + 2*T^11 + T^9 + T^7 + 2*T^3 + 2*T
```

In rank 2, it is also possible to compute the normalized Eisenstein series of weight `q^k - 1` (see (6.9) in \[2\]):

```
    sage: from drinfeld_modular_forms import DrinfeldModularFormsRing
    sage: q = 3
    sage: A = GF(q)['T']; K = Frac(A); T = K.gen()
    sage: M = DrinfeldModularFormsRing(K, 2)
    sage: M.eisenstein_series(q^3 - 1)  # weight q^3 - 1
    g1^13 + (-T^9 + T)*g1*g2^3
```

## Notes

This package is based on the intial implementation of Alex Petrov.

Drinfeld modules are currently being implemented in SageMath. See https://github.com/sagemath/sage/pull/350263. As of March 2023, this PR is merged in the current latest development version of SageMath.

## Further Developments

* Add Hecke operators computations.
* Add general Goss polynomials

## References

* \[1\] Basson D., Breuer F., Pink R., Drinfeld modular forms of arbitrary rank, Part III: Examples, https://arxiv.org/abs/1805.12339
* \[2\] Gekeler, E.-U., On the coefficients of Drinfelʹd modular forms. Invent. Math. 93 (1988), no. 3, 667–700
* \[3\] López, B. A non-standard Fourier expansion for the Drinfeld discriminant function. Arch. Math. 95, 143–150 (2010). https://doi.org/10.1007/s00013-010-0148-7
* \[4\] Petrov A., A-expansions of Drinfeld modular forms. J. Number Theory 133 (2013), no. 7, 2247–2266
