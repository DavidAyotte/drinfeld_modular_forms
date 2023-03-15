# Drinfeld Modular Forms

This SageMath package provides an implementation for computing with Drinfeld modular forms for `GL_r(A)`.

To install this package, first clone this repository and then run the following command (inside the project folder):

`sage -pip install --upgrade --no-index -v .`

If there is any changes to the current repo, you will then simply need to pull the changes and run the above command again.

Next, to use its functionalities in your SageMath session you just need to type:

`sage: from drinfeld_modular_forms.all import *`.

## Examples:

One may create the ring of Drinfeld modular forms:

```
    sage: from drinfeld_modular_forms.all import *
    sage: A = GF(3)['T']; K = Frac(A); T = K.gen()
    sage: M = DrinfeldModularFormsRing(K, 2)
    sage: M.ngens()  # number of generators
    2
```

The elements of this ring are viewed as multivariate polynomials in a choice of generators for the ring. The default generators are the coefficient forms of a universal Drinfeld module over the Drinfeld period domain (see theorem 17.5 in \[1\]).

```
    sage: g0, g1 = M.gens()
    sage: F = (g0 + g1)*g0; F
    g0*g1 + g0^2
```
Note that such elements are not necessarily modular forms as they can have mixed weight components. We will call elements of this ring *graded Drinfeld modular forms*.

In the case of rank 2, one can compute the expansion at infinity of any graded form:

```
    sage: g0.expansion()
    1/(T^3 + 2*T) + t^2 + O(t^12)
    sage: g1.exansion()
    t^2 + 2*t^6 + (T^3 + 2*T)*t^8 + O(t^12)
    sage: ((g0 + g1)*g0).expansion()
    1/(T^6 + T^4 + T^2) + 2*t^4 + (2/(T^3 + 2*T))*t^6 + (T^3 + 2*T)*t^10 + O(t^12)
```
This is achieved via the `A`-expansion theory developed by Petrov in \[3\].

It is also possible to compute the normalized Eisenstein serie of weight `q^k - 1` (see (6.9) in \[2\]):

```
    sage: from drinfeld_modular_forms.all import *
    sage: A = GF(3)['T']; K = Frac(A); T = K.gen()
    sage: M = DrinfeldModularFormsRing(K, 2)
    sage: M.weighted_eisenstein_serie(3)  # weight q^3 - 1
    g0^13 + (-T^9 + T)*g0*g1^3
```

## Note:

Drinfeld modules are currently being implemented in SageMath. See https://github.com/sagemath/sage/pull/350263.
Whenever this ticket gets officially merged in SageMath, the functionalities on Drinfeld modules implemented in this project will be deprecated.

This package is still in development and some parts of the code is
based on the initial implementation of Alex Petrov located here:

`petrov/AlexPetrov-original-code-drinfeld-modular-forms.sage`.

## Further Developments:

* Add Hecke operators computations.
* Speed up the code for `t`-expansions.
* Implement Sturm bound of a Drinfeld modular forms space of fixed weight.

## References:

* \[1\] Basson D., Breuer F., Pink R., Drinfeld modular forms of arbitrary rank, Part III: Examples, https://arxiv.org/abs/1805.12339
* \[2\] Gekeler, E.-U., On the coefficients of Drinfelʹd modular forms. Invent. Math. 93 (1988), no. 3, 667–700
* \[3\] Petrov A., A-expansions of Drinfeld modular forms. J. Number Theory 133 (2013), no. 7, 2247–2266
