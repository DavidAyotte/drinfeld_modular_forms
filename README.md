# Drinfeld Modular Forms

This SageMath package provides an implementation for computing A-expansion of Modular Forms

To install this package, first clone this repository and then run the `make`
command (inside the repository main folder). Next, to use it in your SageMath
session you just need to type:

`sage: from drinfeld_modules import *`.

This package is still in development and some parts of the code is
based on the initial implementation of Alex Petrov located here:

`petrov/AlexPetrov-original-code-drinfeld-modular-forms.sage`.

## Examples:

One may create the ring of Drinfeld modular forms:

```sage
    sage: from drinfeld_modular_forms.all import *
    sage: A = GF(3)['T']; K = Frac(A); T = K.gen()
    sage: M = DrinfeldModularFormsRing(K, 2); M
    Ring of Drinfeld modular forms of rank 2 over Fraction Field of Univariate Polynomial Ring in T over Finite Field of size 3
    sage: M.ngens()  # number of generators
    2
```

The elements of this ring are viewed as multivariate polynomials in a choice of generators for the ring. The default generators are the coefficient forms of a universal Drinfeld module over the Drinfeld period domain.

```sage
    sage: g0, g1 = M.gens()
    sage: F = (g0 + g1)*g0; F
    g0*g1 + g0^2
```
Note that such elements are not necessarily modular forms as they can have mixed weight components. We will call elements of this ring *graded Drinfeld modular forms*.

In the case of rank 2, one can compute the expansion at infinity of any graded form:

```sage
    sage: g0.t_expansion()
    1/(T^3 + 2*T) + t^2 + O(t^12)
    sage: g1.t_exansion()
    t^2 + 2*t^6 + (T^3 + 2*T)*t^8 + O(t^12)
    sage: ((g0 + g1)*g0).t_expansion()
    1/(T^6 + T^4 + T^2) + 2*t^4 + (2/(T^3 + 2*T))*t^6 + (T^3 + 2*T)*t^10 + O(t^12)
```

It is also possible to compute the normalized Eisenstein serie of weight `q^k - 1`:

```sage
    sage: from drinfeld_modular_forms.all import *
    sage: A = GF(3)['T']; K = Frac(A); T = K.gen()
    sage: M = DrinfeldModularFormsRing(K, 2)
    sage: M.weighted_eisenstein_serie(3)  # weight q^3 - 1
    g0^13 + (-T^9 + T)*g0*g1^3
```

## Note:

Drinfeld modules are currently being implemented in SageMath. See https://trac.sagemath.org/ticket/33713.

## Further Developments:

Add Hecke operators computations.
