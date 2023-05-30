=================================
Background material and notations
=================================

In this section, we define some notations and give some definitions.
We will that assume the reader possesses some knowledge of Drinfeld
modules and their analytic theory. A classical reference is chapter 3
and 4 of [Gos1998]_.

.. RUBRIC:: Function field setting

Let `q` be a prime power and let `A` be the ring of functions of
`\mathbb{P}^1/\mathbb{F}_q` which are regular outside `\infty`. This
ring is the polynomial ring `\mathbb{F}_q[T]`. We denote by
`K := \mathbb{F}_q(T)` its fraction field and let
`K_{\infty} := \mathbb{F}_q((1/T))`. We define `\mathbb{C}_{\infty}` to
be the completion of an algebraic closure of `K_{\infty}`. Lastly, We
denote by `\tau : x\mapsto x^q` the `q`-Frobenius.

.. RUBRIC:: Drinfeld period domain and group action

The *Drinfeld period domain of rank* `r > 1` is
defined by

.. MATH::

    \Omega^r(\mathbb{C}_{\infty}) :=
    \mathbb{P}^{r-1}(\mathbb{C}_{\infty})
    \setminus \{K_{\infty}\text{-rational hyperplanes}\}.

This space is a rigid analytic space
and plays the role of the complex upper half plane. We identify the elements
of `\Omega^r(\mathbb{C}_{\infty})` with the set of column vectors
`(w_1,\ldots, w_{r-1}, w_{r})^{\mathrm{T}}` in `\mathbb{C}_{\infty}^r`
such that the `w_i` are `K_{\infty}`-linearly independant and
`w_r = \xi`, a nonzero constant in `\mathbb{C}_{\infty}`.

We define a left action of `\mathrm{GL}_r(K_{\infty})` on
`\Omega^r(\mathbb{C}_{\infty})` by setting

.. MATH::

    \gamma(w) := j(\gamma, w)^{-1}\gamma w

where `j(\gamma, w) := \xi^{-1} \cdot (\text{last entry of }\gamma w)`.

.. RUBRIC:: Universal Drinfeld module over `\Omega^r(\mathbb{C}_{\infty})` and modular forms

For any `w = (w_1, \ldots, w_{r-1}, \xi)` in
`\Omega^r(\mathbb{C}_{\infty})` we have a corresponding free `A`-lattice
of rank `r`:

.. MATH::

    \Lambda_w := Aw_1 \oplus \cdots \oplus Aw_{r-1} \oplus A\xi.

By analytic uniformization, there exists a corresponding Drinfeld module

.. MATH::

    \phi_w : T \mapsto T + g_1(w)\tau + \cdots
    + g_{r - 1}(w)\tau^{r-1} + g_{r}(w)\tau^{r}

where the coefficients
`g_i : \Omega^r(\mathbb{C}_{\infty}) \rightarrow \mathbb{C}_{\infty}`
are rigid analytic functions satisfying the invariance property:

.. MATH::

    g_i(\gamma(w)) = j(\gamma, w)^{1 - q^i} g_i(w),
    ~\forall \gamma\in \mathrm{GL}_r(A)

where `g_{r}(w)` never vanishes. Moreover, these coefficients `g_i`
admits an expansion at infinity, analogous to the classical theory. In
the rank two case, this expansion at infinity is of the form

.. MATH::

    g_i(w) = \sum_{i = 0}^{\infty} a_n(g_i)t(w)^i

where `t(w) := e(w)^{-1}` and `e` is the exponential of the Carlitz
module. Any rigid analytic function
`f : \Omega^r(\mathbb{C}_{\infty}) \rightarrow \mathbb{C}_{\infty}`
that satisfies the invariance property and the expansion at infinity
is called a Drinfeld modular form of weight `k`. The forms `g_i` are
called the *coefficients forms*.

This Drinfeld module is the universal
Drinfeld `\mathbb{F}_q[T]`-module over `\Omega^r(\mathbb{C}_{\infty})`.
The coefficients
`g_i : \Omega^r(\mathbb{C}_{\infty}) \rightarrow \mathbb{C}_{\infty}`
are rigid analytic function which satisfies a modular invariance
properties under the action of the group `\mathrm{GL}_r(A)`. These
function are examples of *modular forms*.

**Definition.**

A *Drinfeld modular form* of rank `r`, weight `k`, type `m` for
`\mathrm{GL}_r(A)` is a rigid analytic function
`f:\Omega^r(\mathbb{C}_{\infty}) \rightarrow \mathbb{C}_{\infty}`
such that

* `f(\gamma(w)) = \mathrm{det}(\gamma)^m j(\gamma, w)^k f(w)` for all
  `\gamma` in `\mathrm{GL}_r(A)` and
  `w\in \Omega^r(\mathbb{C}_{\infty})`;

* `f` is holomorphic at infinity.

The second condition is similar to the classical case. In the rank two
situation, this expansion is simply given by
`f = \sum_{n\geq 0} a_n(f) t^n` where `a_n(f)\in \mathbb{C}_{\infty}`.

The reader is refered to
part I of [BRP2018]_ for more information about the analytic theory of
Drinfeld modular form of arbitrary rank.

.. RUBRIC:: Ring of Drinfeld modular forms

Letting `M_k^{r, m}(\mathrm{GL}_r(A))` denote the space of rank `r`,
weight `k` and type `m` (`k\geq 0` and `m\in (q - 1)\mathbb{Z}`)
Drinfeld modular forms, we define

.. MATH::

    M^{r, 0}(\mathrm{GL}_r(A)) :=
    \bigoplus_{k\in ZZ} M_k^{r}(\mathrm{GL}_r(A))

to be the graded ring of all Drinfeld modular forms. By theorem 17.5 in
part III of [BRP2018]_, we have

.. MATH::

    M^{r, 0}(\mathrm{GL}_r(A))
    = \mathbb{C}_{\infty}[g_1,\ldots, g_{r-1}, g_{r}].

Furthermore, in the rank two case, we also have

.. MATH::

    M^r(\mathrm{GL}_r(A)) :=
    \bigoplus_{k, m} M_k^{r, m}(\mathrm{GL}_r(A))
    = \mathbb{C}_{\infty}[g_1, h]

where `h` is a modular form of weight `q+1`, type 1 with expansion
`t + O(t^2)` (see (5.13) of [Gek1988]_).

.. RUBRIC:: Rank two examples

- Drinfeld Eisenstein series

For `k \equiv 0` modulo `q - 1`. The *Drinfeld Eisenstein series of
weight* `k` and rank 2 is defined by

.. MATH::

    E_{k}(w) :=
    \sum_{\substack{ (c, c)\in A^{2} \\ (c, c) \neq (0, 0) }}
    \frac{1}{(cw + d)^k}.

This series is absolutely and uniformly convergent and admits the
expansion

.. MATH::

    E_k(w) = \tilde{\pi}^k\delta_k
    - \tilde{\pi}^k \sum_{\substack{a\in A \\a\text{ monic}}} G_k(t(aw))

where `\tilde{\pi}` is the Carlitz period (analogue of `\pi`) and `G_k` is
the `k`-th Goss polynomial and `\delta_k \in K` is some constant
depending on `k`. See section 6 of [Gek1988]_ for the proof of this
fact. We will denote by

.. MATH::

    g_k := \tilde{\pi}^{q^k - 1}\delta_{q^k - 1} E_{q^k - 1}

the *normalized* Eisenstein series. For `k = 1,\ldots r-1`, these forms
corresponds to the coefficients forms defined above.

- Modular discriminant

The *modular discriminant*
`\Delta : \Omega^2(\mathbb{C}_{\infty}) \rightarrow \mathbb{C}` is the
leading coefficient form of the rank 2 universal Drinfeld module over
`\Omega^2(\mathbb{C}_{\infty})`:

.. MATH::

    \phi^w : T \mapsto T + g_1(w)\tau + \Delta(w)\tau^2.

By the work of López in [Lop2010]_, the discriminant function admits an
expansion of the form

.. MATH::

    -\tilde{\pi}^{1 - q^2}\Delta(w)
    = \sum_{\substack{a\in A\\a \text{ monic}}}
    a^{q(q-1)} t(aw)^{q-1}.

- Petrov `A`-expansions

We say that a Drinfeld modular forms of weight `k` admits a
*Petrov expansion* or an `A`-*expansion* if there exists an integer `n`
and elements `c_{a}(f)\in \mathbb{C}_{\infty}` such that

.. MATH::

    f =
    \sum_{\substack{a\in \mathbb{F}_q[T] \\ a\text{ monic}}}
    c_a(f)G_n(t(az)).

In [Pet2013]_, Petrov showed that

.. MATH::

    f_{k, n} :=
    \sum_{\substack{a\in \mathbb{F}_q[T] \\ a\text{ monic}}}
    a^{k - n}G_n(t(az))

defines an infinite family of Drinfeld modular forms of weight `k`
provided that `k - 2n \equiv 0` modulo `q - 1` and
`n \leq p^{v_p(k - n)}`. See theorem 1.3 of loc. cit. for more details.

.. RUBRIC:: References

.. [BRP2018] Basson D., Breuer F., Pink R., Drinfeld modular forms of
             arbitrary rank:
             Part I: `arxiv:1805.12335 <https://arxiv.org/abs/1805.12335>`_,
             Part II: `arxiv:1805.12337 <https://arxiv.org/abs/1805.12337>`_,
             Part III: `arxiv:1805.12339 <https://arxiv.org/abs/1805.12339>`_,
             2018.

.. [Gek1988] Gekeler, EU. On the coefficients of Drinfeld modular forms.
             Invent Math 93, 667-700 (1988).
             `doi.org/10.1007/BF01410204 <https://doi.org/10.1007/BF01410204>`_

.. [Gos1998] Goss D. Basic structures of function field arithmetic.
             Springer, 1998.
             `doi.org/10.1007/978-3-642-61480-4 <https://doi.org/10.1007/978-3-642-61480-4>`_

.. [Lop2010] López, B. A non-standard Fourier expansion for the Drinfeld
             discriminant function. Arch. Math. 95, 143-150 (2010).
             `doi.org/10.1007/s00013-010-0148-7 <https://doi.org/10.1007/s00013-010-0148-7>`_

.. [Pet2013] Petrov A., A-expansions of Drinfeld modular forms,
             Journal of Number Theory, Volume 133, Issue 7, 2013,
             `doi.org/10.1016/j.jnt.2012.12.012 <https://doi.org/10.1016/j.jnt.2012.12.012>`_
