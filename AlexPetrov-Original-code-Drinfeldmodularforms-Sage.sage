runr"""
	AUTHOR: Alexander Petrov

       You should pick/define 'q' before loading the following routines in SAGE.
       The routines are based on:
         [Gek] On the coefficients of Drinfeld modular forms, Invent. Math. Volume 93, Issue 3, pp 667-700, 1988.
                   http://link.springer.com/article/10.1007%2FBF01410204
         [Tha]  Function Field Arithmetic, World Scientific, ISBN: 978-981-238-839-1.
                   http://www.worldscientific.com/worldscibooks/10.1142/5535
"""
q=3
F.<prim_root_mod_q> = FiniteField(q)
P.<T> = PolynomialRing(F)
P1.<u> = PowerSeriesRing(FractionField(P))
P2.<X> = PolynomialRing(FractionField(P))

print("----->")
print("The following objects have been initiated:")
print("--- 'F' a finite field of 'q' elements with generator 'prim_root_mod_q';")
print("--- 'P' a univariate polynomial ring in 'T' over 'F'; ")
print("--- 'P1' a power series ring in 'u' over 'FractionField(P)'.")
print("<-----")
print("To see a list of available functions type 'AvailableFunctions()'")
print("To see a descriptions of a specific function 'Fun_Name' ")
print("type 'Fun_Name?'")
print(" ")

def AvailableFunctions():
    print("Br")
    print("D")
    print("L")
    print("rho")
    print("fpol")
    print("PrimesOfDegree")
    print("Fstar")
    print("AA")
    print("GossPolynomial")
    print("Normalize")
    print("ta")
    print("Aexpansion")
    print("hh")
    print("gg")
    print("Ee")
    print("GossPolynomialList")
    print("HeckeAction")
    print("ListOfPowersOf")
    print("BasisFor_M_k_m")
    print("DimensionOf_M_k_m")
    print("BasisRepresentation")
    print("HeckeMatrix")
    print("FactorWithSingular")
    print("Eigenforms")
    print ("hdif_sub1")
    print ("hdif_sub2")
    print ("Hdiff")

def Br(n):
    r"""
    Returns '[n] = T^(q^n) - T'

    INPUT:  - A positive integer 'n'.

    OUTPUT: - A polynomial in 'Fq[T]'.

    EXAMPLES: sage: q = 3
    	      sage: runfile DMFv1.sage
	      ----->
	      The following objects have been initiated:
	      --- 'F' a finite field of 'q' elements with generator 'prim_root_mod_q';
	      --- 'P' a univariate polynomial ring in 'T' over 'F';
	      --- 'P1' a power series ring in 'u' over 'FractionField(P)'.
	      <-----
	      To see a list of available functions type 'AvailableFunctions()'
	      To see a descriptions of a specific function 'Fun_Name'
	      type 'Fun_Name?'

	      sage: Br(0)
	      ---------------------------------------------------------------------------
	      NotImplementedError                       Traceback (most recent call last)
	      <ipython-input-3-05cfc50ae8ec> in <module>()
	      ----> 1 Br(Integer(0))

	      <string> in Br(n)

	      NotImplementedError:  Input must be positive.
	      sage: Br(1)
	      T^3 + 2*T
	      sage: Br(5)
	      T^243 + 2*T
    """
    if n <= 0:
       raise NotImplementedError(" Input must be positive. ")
    return T^(q^n) - T

def D(n):
    r"""
    Returns the product of all monic polynomials in 'Fq[T]' of degree 'n'. See [Tha] for reference. Uses 'Br'.

    INPUT: - A non-negative integer 'n'.

    OUTPUT: - A polynomial in 'Fq[T]'

    EXAMPLES: sage: q = 3
    	      sage: runfile DMFv1.sage
	      ----->
	      The following objects have been initiated:
	      --- 'F' a finite field of 'q' elements with generator 'prim_root_mod_q';
	      --- 'P' a univariate polynomial ring in 'T' over 'F';
	      --- 'P1' a power series ring in 'u' over 'FractionField(P)'.
	      <-----
	      To see a list of available functions type 'AvailableFunctions()'
	      To see a descriptions of a specific function 'Fun_Name'
	      type 'Fun_Name?'

	       sage: D(1)
	       T^3 + 2*T
	       sage: D(2)
	       T^18 + 2*T^12 + 2*T^10 + T^4
	       sage: factor(D(2))
	       T^4 * (T + 1)^4 * (T + 2)^4 * (T^2 + 1) * (T^2 + T + 2) * (T^2 + 2*T + 2)
    """
    if n < 0:
       raise NotImplementedError(" Input must be non-negative. ")
    if n == 0:
       return 1 + 0*T
    ans = 1
    for j in range(1, n+1):
    	ans = ans*Br(j)^(q^(n-j))
    return ans

def L(n):
    r"""
    Returns the least common multiple (monic) of all monic polynomials of degree 'n'. See [Tha] for reference. Uses 'Br'.

    INPUT: - A non-negative integer 'n'.

    OUTPUT: - A polynomial in 'Fq[T]'.

    EXAMPLES: sage: q = 5

    	      sage: runfile DMFv1.sage

	      ----->

	      The following objects have been initiated:

	      --- 'F' a finite field of 'q' elements with generator 'prim_root_mod_q';

	      --- 'P' a univariate polynomial ring in 'T' over 'F';

	      --- 'P1' a power series ring in 'u' over 'FractionField(P)'.

	      <-----

	      To see a list of available functions type 'AvailableFunctions()'

	      To see a descriptions of a specific function 'Fun_Name'

	      type 'Fun_Name?'


	      sage: L(2)

	      T^30 + 4*T^26 + T^25 + 4*T^6 + T^5 + T^2 + 3*T + 1

	      sage: L(3)

	      T^155 + 4*T^151 + T^150 + 4*T^131 + T^130 + T^127 + 3*T^126 + T^125 + 4*T^31 + T^30 + T^27 + 3*T^26 + T^25 + T^7 + 3*T^6 + T^5 + 4*T^3 + 3*T^2 + 2*T + 1
	      sage:
    """
    if n < 0:
       raise NotImplementedError(" Input must be non-negative. ")
    if n == 0:
       return 1 + 0*T
    ans = 1
    for i in range(1, n+1):
    	ans = ans*Br(i)
    return ans

def rho(poly):
    r"""
    Returns the value of the Carlitz module for 'poly', where 'poly' is a polynomial in 'Fq[T]'. The Carlitz module value is returned as a polynomial in 'Fq[T, u]'. See [Tha] for reference. Uses 'Br'.

    INPUT: - 'poly' a polynomial in 'Fq[T]'.

    OUTPUT: - An element of 'Fq[T][[u]]', which is really an element of 'Fq[T, u]'.

    EXAMPLES: sage: q = 5
    	      sage: runfile DMFv1.sage
	      ----->
	      The following objects have been initiated:
	      --- 'F' a finite field of 'q' elements with generator 'prim_root_mod_q';
	      --- 'P' a univariate polynomial ring in 'T' over 'F';
	      --- 'P1' a power series ring in 'u' over 'FractionField(P)'.
	      <-----
	      To see a list of available functions type 'AvailableFunctions()'
	      To see a descriptions of a specific function 'Fun_Name'
	      type 'Fun_Name?'

	      sage: rho(T)
	      T*u + u^5
	      sage: rho(2*T^2  + 2) == 2*rho(T^2) + rho(2)
	      True
	      sage: rho(T)(rho(T)) == rho(T^2)
	      True
	      sage:
    """
    poly_used = P(poly)
    ans = poly_used*u
    aux = (poly_used^q - poly_used)/Br(1)
    d = poly_used.degree()
    for i in range(1, d+1):
    	ans = ans + aux*u^(q^i)
	aux = (aux^q - aux)/Br(i + 1)
    return ans

def fpol(poly):
    r"""
    Returns the inverse cyclotomic polynomial of a 'poly'. See [Gek] for reference. Uses 'rho'.

    INPUT: - 'poly' a polynomial in 'Fq[T]'.

    OUTPUT: - An element of 'Fq[T][[u]]', which is really in 'Fq[T, u].

    EXAMPLES: sage: q
    	      3
	      sage: fpol(T)
	      1 + T*u^2
	      sage: rho(T)
	      T*u + u^3
	      sage: fpol(T^5 + T^3 + 1) == fpol(T^5) + u^(q^5 - q^3)*fpol(T^3 + 1)
	      True
	      sage:
    """
    poly_used = P(poly)
    d = poly_used.degree()
    if d == -1:
       return 0*u^0
    return (rho(poly)(1/u))*u^(q^d)

def PrimesOfDegree(d):
    r"""
    Returns a list of the monic primes in 'Fq[T]' of degree 'd'.

    INPUT: - A positive integer 'd'.

    OUTPUT: - A list of elements in 'Fq[T]'.

    EXAMPLES: sage: q
    	      3
	      sage: PrimesOfDegree(1)
	      [T, T + 1, T + 2]
	      sage: PrimesOfDegree(2)
	      [T^2 + 1, T^2 + T + 2, T^2 + 2*T + 2]
	      sage: PrimesOfDegree(3)
	      [T^3 + 2*T + 1,
	       T^3 + 2*T + 2,
	       T^3 + T^2 + 2,
	       T^3 + T^2 + T + 2,
	       T^3 + T^2 + 2*T + 1,
	       T^3 + 2*T^2 + 1,
	       T^3 + 2*T^2 + T + 1,
	       T^3 + 2*T^2 + 2*T + 2]
	      sage:
    """
    if d<= 0:
       raise NotImplementedError(" Input must be a positive integer. ")
    Ans = []
    for pp in P.monics(of_degree = d):
    	if factor(pp)[0][0] == pp:
	   Ans = Ans + [pp]
    return Ans

def Fstar():
    r"""
    Returns a list of all the non-zero elements of 'F'.

    INPUT: - None.

    OUTPUT: - A list of elements in 'F'.

    EXAMPLES: sage: q = 4
    	      sage: load DMFv1.sage
	      ----->
	      The following objects have been initiated:
	      --- 'F' a finite field of 'q' elements with generator 'prim_root_mod_q';
	      --- 'P' a univariate polynomial ring in 'T' over 'F';
	      --- 'P1' a power series ring in 'u' over 'FractionField(P)'.
	      <-----
	      To see a list of available functions type 'AvailableFunctions()'
	      To see a descriptions of a specific function 'Fun_Name'
	      type 'Fun_Name?'

	      sage: Fstar()
	      [prim_root_mod_q, prim_root_mod_q + 1, 1]
	      sage: q = 5
	      sage: load DMFv1.sage
	      ----->
	      The following objects have been initiated:
	      --- 'F' a finite field of 'q' elements with generator 'prim_root_mod_q';
	      --- 'P' a univariate polynomial ring in 'T' over 'F';
	      --- 'P1' a power series ring in 'u' over 'FractionField(P)'.
	      <-----
	      To see a list of available functions type 'AvailableFunctions()'
	      To see a descriptions of a specific function 'Fun_Name'
	      type 'Fun_Name?'

	      sage: Fstar()
	      [1, 2, 3, 4]
	      sage:
    """
    Ans = []
    for theta in F:
    	if theta != 0:
	   Ans = Ans + [theta]
    return Ans

def AA(d):
    r"""
    Returns a list of all the polynomials in 'Fq[T]' of degree strictly less than 'd'. Uses 'Fstar'.

    INPUT: - A non-negative integer 'd'.

    OUTPUT: - A list of elements in 'Fq[T]'.

    EXAMPLES: sage: q
    	      3
	      sage: AA(0)
	      [0]
	      sage: AA(1)
	      [0, 1, 2]
	      sage: AA(2)
	      [0, 1, 2, T, 2*T, T + 1, 2*T + 2, T + 2, 2*T + 1]
	      sage:
    """
    if d < 0:
       raise NotImplementedError(" Input must be a non-negative integer. ")
    Ans = [P(0)]
    for j in range(0, d):
    	for pp in P.monics(of_degree = j):
	    for theta in Fstar():
	    	Ans = Ans + [P(theta)*pp]
    return Ans

def GossPolynomial(n):
    r"""
    Returns the 'n'-th Goss polynomial for the lattice '\pi Fq[T]'. Uses 'D'.

    INPUT: - 'n' a non-negative integer.

    OUTPUT: - An element of 'Fq[T][[u]]', which is really a polynomial in 'Fq[T, u]'.

    EXAMPLES: sage: q
    	      3
	      sage: GossPolynomial(2)
	      u^2
	      sage: GossPolynomial(4)
	      (1/(T^3 + 2*T))*u^2 + u^4
	      sage: GossPolynomial(5)
	      (2/(T^3 + 2*T))*u^3 + u^5
	      sage
    """
    if n == 0:
       return 0*u^0
    if n <= q-1:
       return u^n
    if (n % q) == 0:
       return GossPolynomial(n/q)^q
    ans = 0
    for j in range(0, floor(N(log(n, q))) + 1):
    	ans = ans + u*GossPolynomial(n-q^j)*(1/D(j))
    return ans

def Normalize(powerseries):
    r"""
    Returns the input 'powerseries' normalized so that its least-degree, non-zero coefficient is 1.

    INPUT: - 'powerseries' an element of 'Fq[T][[u]]'.

    OUTPUT: - An element of 'Fq[T][[u]]'.

    EXAMPLES: sage: q
    	      3
              sage: GossPolynomial(4)
	      (1/(T^3 + 2*T))*u^2 + u^4
	      sage: Normalize(GossPolynomial(4))
	      u^2 + (T^3 + 2*T)*u^4
	      sage: Normalize(GossPolynomial(5))
	      u^3 + (2*T^3 + T)*u^5
	      sage: GossPolynomial(5)
	      (2/(T^3 + 2*T))*u^3 + u^5
	      sage:
    """
    return (1/powerseries.coefficients()[0])*powerseries

def ta(poly, preCision):
    r"""
    Returns 't_(poly)' truncated up to 'u^n'. See [Gek] for reference. Uses 'fpol' and PARI/GP.

    INPUT: - 'poly' an element of 'Fq[T]'; 'preCision' a positive integer.

    OUTPUT: - A power series in 'Fq[T][[u]]'.

    EXAMPLES: sage: q
    	      3
	      sage: ta(T^2 + T + 1, 20)
	      u^9 + (2*T^3 + 2*T + 2)*u^15 + (2*T^2 + 2*T + 2)*u^17 + (T^6 + 2*T^4 + 2*T^3 + T^2 + 2*T + 1)*u^21 + (2*T^5 + 2*T^4 + T^3 + T^2 + T + 2)*u^23 + (T^4 + 2*T^3 + 2*T + 1)*u^25 + (2*T^9 + 2*T^3 + 2)*u^27
	      sage: ta(T^2 + T + 1, 20)*fpol(T^2 + T + 1)
	      u^9 + (2*T^12 + 2*T^10 + 2*T^9 + 2*T^4 + 2*T^3 + 2*T)*u^33 + (2*T^11 + 2*T^10 + 2*T^9 + 2*T^5 + 2*T^4 + 2*T^3 + 2*T^2 + 2*T + 2)*u^35
	      sage: ta(T^3 + T + 1, 20)*fpol(T^3 + T + 1) + O(u^20)
	      O(u^20)
	      sage: ta(T^3 + T + 1, 20)*fpol(T^3 + T + 1) + O(u^40)
	      u^27 + O(u^40)
	      sage: % time ta(T^2 + 1, 1000);
	      CPU times: user 2.27 s, sys: 0.18 s, total: 2.45 s
	      Wall time: 2.48 s
	      sage:
    """
    poly_used = P(poly)
    d = poly_used.degree()
    if d == -1:
       raise NotImplementedError(" The polynomial must be non-zero. ")
    return u^(q^d)*(1/(fpol(poly) + O(u^preCision)))

def Aexpansion(k, n, deg, preCision):
    r"""
    Returns the A-expansion '\sum a^(k-n) G_n (t_a)' by using polynomials of degree less than or equal to 'deg' and up to precision 'preCision'. Uses 'ta', 'Normalize', 'GossPolynomial'.

    INPUT: - 'k', 'n', 'deg', 'preCision' positive integers.

    OUTPUT: - A power series in 'Fq[T][[u]]'.

    EXAMPLES: sage: q
    	      3
	      sage: Aexpansion(q+1, 1, 3, 50) + O(u^11)
	      u + u^5 + (2*T^3 + T)*u^7 + u^9 + O(u^11)
	      sage: %time Aexpansion(q+1, 1, 3, 50) + O(u^11)
	      CPU times: user 1.93 s, sys: 0.02 s, total: 1.96 s
	      Wall time: 1.94 s
	      u + u^5 + (2*T^3 + T)*u^7 + u^9 + O(u^11)
	      sage: %time Aexpansion(q+1, 1, 3, 250) + O(u^11)
	      Warning: Cannot ensure the required precision. Choose 'deg', so that 'preCision <= q^(deg + 1)'.
	      CPU times: user 7.19 s, sys: 0.16 s, total: 7.36 s
	      Wall time: 7.31 s
	      u + u^5 + (2*T^3 + T)*u^7 + u^9 + O(u^11)
	      sage: %time Aexpansion(q+1, 1, 4, 250) + O(u^11)
	      Warning: Cannot ensure the required precision. Choose 'deg', so that 'preCision <= q^(deg + 1)'.
	      CPU times: user 46.09 s, sys: 0.31 s, total: 46.40 s
	      Wall time: 46.29 s
	      u + u^5 + (2*T^3 + T)*u^7 + u^9 + O(u^11)
	      sage: %time Aexpansion(q+1, 1, 5, 250) + O(u^11)
	      CPU times: user 451.42 s, sys: 1.11 s, total: 452.53 s
	      Wall time: 452.43 s
	      u + u^5 + (2*T^3 + T)*u^7 + u^9 + O(u^11)
	      sage:
    """
    if preCision >= q^(deg + 1):
       print("Warning: Cannot ensure the required precision. Choose 'deg', so that 'preCision <= q^(deg + 1)'.")
    characteristic = F.characteristic()
    if n > characteristic^((k-n).valuation(characteristic)):
       print("Warning: The computed A-expansion is not a modular form.")
    goss_polynomial = Normalize(GossPolynomial(n))
    ans = 0
    for j in range(0, deg + 1):
    	for p in P.monics(of_degree = j):
	    ans = ans + p^(k-n)*goss_polynomial(ta(p, preCision))
    return ans + O(u^preCision)

def hh(deg,  preCision):
    r"""
    Returns the expansion of 'h' up to precision 'preCision' by using polynomials of degree less than or equal to 'deg'. Uses 'Aexpansion'.

    INPUT: - 'deg', 'preCision' positive integers.

    OUTPUT: - A power series in 'Fq[T][[u]]'.

    EXAMPLES: sage: q
    	      5
	      sage: h = hh(3, 100)
	      sage: for j in range(1, q+2):
	      ....:     h^j - Aexpansion((q+1)*j, j, 3, 100) + O(u^100)
	      ....:
	      O(u^100)
	      O(u^100)
	      O(u^100)
	      O(u^100)
	      O(u^100)
	      Warning: The computed A-expansion is not a modular form.
	      4*u^2 + (4*T^5 + T + 1)*u^6 + 3*u^18 + (T^25 + T^5 + 3*T + 1)*u^22 + (4*T^5 + T)*u^26 + 2*u^34 + (2*T^25 + 3*T + 1)*u^38 + (2*T^30 + 3*T^26 + 3*T^10 + 2*T^6 + 3*T^5 + 2*T)*u^42 + (T^35 + 3*T^31 + T^27 + 4*T^11 + T^10 + 2*T^7 + 3*T^6 + 4*T^3 + T^2)*u^46 + u^50 + (3*T^25 + T^5 + T + 1)*u^54 + (2*T^30 + 3*T^26 + 3*T^6 + 2*T^5 + 2*T^2 + 3*T)*u^58 + (2*T^35 + T^31 + 2*T^27 + 3*T^15 + 4*T^11 + 3*T^10 + 3*T^7 + 4*T^6 + 3*T^2)*u^62 + (3*T^40 + T^36 + 4*T^32 + 2*T^28 + 2*T^16 + 4*T^15 + 4*T^12 + 3*T^11 + T^8 + 2*T^7 + 3*T^4 + T^3)*u^66 + (4*T^25 + 4*T^5 + 2*T + 1)*u^70 + (3*T^10 + 4*T^6 + T^5 + 3*T^2 + 4*T)*u^74 + (T^35 + 3*T^31 + T^27 + 2*T^15 + 3*T^11 + T^10 + 3*T^7 + 3*T^6 + 2*T^3 + T^2)*u^78 + (2*T^40 + 4*T^36 + T^32 + 3*T^28 + T^20 + 4*T^16 + T^15 + 2*T^12 + 2*T^11 + 3*T^7 + 3*T^4 + 4*T^3 + 4)*u^82 + (3*T^45 + 3*T^41 + 3*T^37 + 3*T^33 + 3*T^29 + 2*T^21 + T^20 + 2*T^17 + T^16 + 2*T^13 + T^12 + 2*T^9 + T^8 + T^4 + 2*T + 2)*u^86 + (T^30 + 4*T^26 + 4*T^10 + T^6)*u^90 + (T^35 + 3*T^31 + T^27 + 4*T^15 + 2*T^11 + 4*T^7)*u^94 + (T^40 + 2*T^36 + 3*T^32 + 4*T^28 + 4*T^20 + 3*T^16 + 2*T^12 + T^8 + 3)*u^98 + O(u^100)
	      sage:
    """
    return Aexpansion(q+1, 1, deg, preCision)

def Ee(deg,  preCision):
    r"""
    Returns the expansion of 'E' (Gekeler's false Eisenstein series) up to precision 'preCision' by using polynomials of degree less than or equal to 'deg'. Uses 'Aexpansion'.

    INPUT: - 'deg', 'preCision' positive integers.

    OUTPUT: - A power series in 'Fq[T][[u]]'.

    EXAMPLES: sage: q
    	      4
	      sage: E = Ee(2, 100)
	      sage: for j in range(1, q+2):
	      ....:     E^j - Aexpansion(2*j, j, 2, 100) + O(u^100)
              ....:
              Warning: Cannot ensure the required precision. Choose 'deg', so that 'preCision <= q^(deg + 1)'.
              O(u^100)
              Warning: Cannot ensure the required precision. Choose 'deg', so that 'preCision <= q^(deg + 1)'.
              O(u^100)
              Warning: Cannot ensure the required precision. Choose 'deg', so that 'preCision <= q^(deg + 1)'.
              Warning: The computed A-expansion is not a modular form.
              O(u^100)
             Warning: Cannot ensure the required precision. Choose 'deg', so that 'preCision <= q^(deg + 1)'.
             O(u^100)
             Warning: Cannot ensure the required precision. Choose 'deg', so that 'preCision <= q^(deg + 1)'.
             Warning: The computed A-expansion is not a modular form.
             u^2 + (T^4 + T + 1)*u^5 + (T^4 + T + 1)*u^14 + u^20 + (T^4 + T + 1)*u^23 + (T^4 + T)*u^26 + (T^4 + T + 1)*u^32 + (T^12 + T^9 + T^8 + T^6 + T^3 + T^2 + 1)*u^38 + (T^8 + T^4 + T^2 + T)*u^44 + (T^12 + T^9 + T^8 + T^6 + T^3 + T^2)*u^47 + (T^12 + T^9 + T^6 + T^3)*u^50 + u^56 + (T^20 + T^17 + T^16 + T^5 + T^4 + T^2)*u^62 + (T^16 + T^4)*u^68 + (T^20 + T^17 + T^16 + T^8 + T^5 + T^4)*u^71 + (T^20 + T^17 + T^8 + T^5 + 1)*u^74 + (T^4 + T + 1)*u^77 + (T^20 + T^17 + T^16 + T^8 + T^5 + T^4)*u^80 + (T^28 + T^25 + T^24 + T^22 + T^19 + T^18 + T^16 + T^13 + T^12 + T^10 + T^7 + T^6 + T^4 + T + 1)*u^86 + (T^24 + T^20 + T^18 + T^17 + T^12 + T^8 + T^6 + T^5 + 1)*u^92 + (T^28 + T^25 + T^24 + T^22 + T^19 + T^18 + T^16 + T^13 + T^12 + T^10 + T^7 + T^6 + T^4 + T + 1)*u^95 + (T^28 + T^25 + T^22 + T^19 + T^16 + T^13 + T^10 + T^7 + T^4 + T)*u^98 + O(u^100)
             sage: E
             u + u^10 + u^19 + (T^4 + T)*u^22 + u^28 + (T^8 + T^2)*u^34 + u^37 + (T^4 + T)*u^40 + (T^8 + T^2)*u^43 + (T^12 + T^9 + T^6 + T^3 + 1)*u^46 + u^55 + (T^16 + T)*u^58 + u^64 + (T^16 + T^4)*u^67 + (T^20 + T^17 + T^5 + T^2)*u^70 + u^73 + (T^16 + T)*u^76 + (T^8 + T^2)*u^79 + (T^24 + T^18 + T^9 + T^3 + 1)*u^82 + (T^16 + T^4)*u^85 + (T^20 + T^17 + T^8 + T^5)*u^88 + (T^24 + T^18 + T^12 + T^6 + 1)*u^91 + (T^28 + T^25 + T^22 + T^19 + T^16 + T^13 + T^10 + T^7 + T^4 + T)*u^94 + O(u^100)
    """
    return Aexpansion(2, 1, deg, preCision)


def gg(deg, preCision):
    r"""
    Returns the expansion of 'g' up to precision 'preCision' by using polynomials of degree less than or equal to 'deg'. Uses 'ta', 'Br'.

    INPUT: - 'def', 'preCision' positive integers.

    OUTPUT: A power series in 'Fq[T][[u]].

    EXAMPLES: sage: q
    	      4
	      sage: g = gg(3, 50)
	      sage: g
	      1 + (T^4 + T)*u^3 + (T^4 + T)*u^39 + (T^4 + T)*u^48 + O(u^50)
	      sage: gg(3, 50) - gg(2, 50)
	      O(u^50)
	      sage:
    """
    if preCision >= q^(deg + 1):
       print("Warning: Cannot ensure the required precision. Choose 'deg', so that 'preCision <= q^(deg + 1)'.")
    ans = 1
    for j in range(0, deg + 1):
    	for p in P.monics(of_degree = j):
	    ans = ans - Br(1)*ta(p, preCision)^(q-1) + O(u^preCision)
    return ans

def GossPolynomialList(pp, n, previous_list = []):
    r"""
    Returns a list of Goss polynomials for the lattice of 'pp'-torsion points. The list has size 'n' + 'len(previous_list)', where 'previous_list' is an optional input. Uses 'rho'.

    INPUT: - 'pp' a monic prime in 'Fq[T]'; 'n' a positive integer; 'previous_list' is a list of elements in 'Fq[T][[u]]' (optional input).

    OUTPUT: - A list of elements in 'Fq[T][[u]]'.

    EXAMPLES: sage: q
    	      3
	      sage: %time list1 = GossPolynomialList(T^2 + 1, 100)
	      CPU times: user 1.81 s, sys: 0.01 s, total: 1.82 s
	      Wall time: 1.82 s
	      sage: %time list2 = GossPolynomialList(T^2 + 1, 100)
	      CPU times: user 1.62 s, sys: 0.01 s, total: 1.63 s
	      Wall time: 1.63 s
	      sage: %time list2 = GossPolynomialList(T^2 + 1, 200)
	      CPU times: user 12.96 s, sys: 0.04 s, total: 12.99 s
	      Wall time: 13.01 s
	      sage: %time list3 = GossPolynomialList(T^2 + 1, 400)
	      CPU times: user 95.34 s, sys: 0.19 s, total: 95.53 s
	      Wall time: 95.62 s
	      sage: %time list3 == GossPolynomialList(T^2 + 1, 300, list1)
	      CPU times: user 82.08 s, sys: 0.13 s, total: 82.21 s
	      Wall time: 82.25 s
	      True
	      sage: %time list3 == GossPolynomialList(T^2 + 1, 200, list2)
	      CPU times: user 62.62 s, sys: 0.04 s, total: 62.66 s
	      Wall time: 62.65 s
	      True
	      sage:
    """
    if previous_list == []:
       previous_list = [0, u]
    Ans = previous_list
    expoNential= (1/pp)*rho(pp)
    for j in range(0, n):
    	place_in_list = len(previous_list) + j
	if (place_in_list % q) == 0:
	   next_goss_polynomial = Ans[Integer(place_in_list/q)]^q
	if (place_in_list % q) != 0:
	   next_goss_polynomial = 0
	   for i in range(0, floor(N(log(place_in_list, q)))+1):
	       next_goss_polynomial = next_goss_polynomial + u*Ans[place_in_list-q^i]*expoNential[q^i]
        Ans = Ans + [next_goss_polynomial]
    return Ans

def HeckeAction(pp, ff, k, m, preCision, Goss_polynomial_list = []):
    r"""
    Returns 'T_pp ff' , where 'ff' is a Drinfeld modular form of weight 'k', up to precision 'O(u^preCision)'. The algorithm is an implementation of formula (7.3) on page 685 of [Gek]. Uses 'GossPolynomialList', 'ta'.

    INPUT: - 'pp' a monic prime in 'Fq[T]'; 'ff' an element of 'Fq[T][[u]]'; 'k', 'm', 'preCision' positive integers; 'Goss_polynomial_list' is a list of elements in 'Fq[T][[u]]' (optional input).

    OUTPUT: - A power series in 'Fq[T][[u]]'.

    EXAMPLES: sage: q
    	      3
	      sage: h = hh(3, 80)
	      sage: HeckeAction(T^2+1, h, 4, 1, 20) - (T^2+1)*h + O(u^20)
	      ---------------------------------------------------------------------------
	      NotImplementedError                       Traceback (most recent call last)
	      <ipython-input-37-1916a727bed2> in <module>()
	      ----> 1 HeckeAction(T**Integer(2)+Integer(1), h, Integer(4), Integer(20)) - (T**Integer(2)+Integer(1))*h + O(u**Integer(20))

	      <string> in HeckeAction(pp, ff, k, preCision, Goss_polynomial_list)

	      NotImplementedError:  Insufficient precision of 'ff', compute the expansion of 'ff' up to precision 181.
	      sage: h = hh(4, 182)
	      sage: HeckeAction(T^2+1, h, 4, 1, 20) - (T^2+1)*h + O(u^20)
	      O(u^20)
	      sage: g = gg(4, 200)
	      sage: HeckeAction(T^2+1, g, 2, 0, 20) - (T^2+1)^2*g + O(u^20)
	      O(u^20)
	      sage: %time f_14_5 = Aexpansion(14, 5, 4, 200)
	      CPU times: user 652.34 s, sys: 1.12 s, total: 653.46 s
	      Wall time: 654.10 s
	      sage: %time HeckeAction(T^2+1, f_14_5, 14, 1,  20) - (T^2 + 1)^5*f_14_5
	      CPU times: user 7.29 s, sys: 0.02 s, total: 7.30 s
	      Wall time: 7.30 s
	      O(u^20)
	      sage:
    """
    d = pp.degree()
    bound_first_sum = ceil(preCision/q^d) + 1
    bound_second_sum = preCision*q^d + 1
    if bound_second_sum > ff.degree():
       raise NotImplementedError(" Insufficient precision of 'ff', compute the expansion of 'ff' up to precision " + str(bound_second_sum) + ".")
    if Goss_polynomial_list == []:
       Goss_polynomial_list = GossPolynomialList(pp, min(bound_second_sum, ff.degree()))
    ans = 0*u^0
    for j in range(0, bound_first_sum + 1):
    	if ff[j] != 0:
            #Inserting a P! below around ta(?,?) to force sage to recognize that it's in the power series ring, and not the fraction field of it!
	   ans = ans + pp^k*ff[j]*P1((ta(pp, preCision) + O(u^preCision))^j) + O(u^preCision)
    int_med = ans
    for j in range(0, bound_second_sum+1):
    	if ff[j] != 0:
	   ans = ans + ff[j]*(Goss_polynomial_list[j] + O(u^preCision))(pp*u) + O(u^preCision)
    return ans

def UpAction(pp, ff, k, m, preCision, Goss_polynomial_list = []):
    r"""
    Returns 'U_pp ff' , where 'ff' is a Drinfeld modular form of weight 'k', up to precision 'O(u^preCision)'. The algorithm is an implementation of formula (7.3) on page 685 of [Gek]. Uses 'GossPolynomialList', 'ta'.

    INPUT: - 'pp' a monic prime in 'Fq[T]'; 'ff' an element of 'Fq[T][[u]]'; 'k', 'm', 'preCision' positive integers; 'Goss_polynomial_list' is a list of elements in 'Fq[T][[u]]' (optional input).

    OUTPUT: - A power series in 'Fq[T][[u]]'.
    """
    d = pp.degree()
#    bound_first_sum = ceil(preCision/q^d) + 1
    bound_second_sum = preCision*q^d + 1
    if bound_second_sum > ff.degree():
       raise NotImplementedError(" Insufficient precision of 'ff', compute the expansion of 'ff' up to precision " + str(bound_second_sum) + ".")
    if Goss_polynomial_list == []:
       Goss_polynomial_list = GossPolynomialList(pp, min(bound_second_sum, ff.degree()))
    ans = 0*u^0
#    for j in range(0, bound_first_sum + 1):
#    	if ff[j] != 0:
#   ans = ans + pp^k*ff[j]*(ta(pp, preCision) + O(u^preCision))^j + O(u^preCision)
#    int_med = ans
    for j in range(0, bound_second_sum+1):
    	if ff[j] != 0:
	   ans = ans + ff[j]*(Goss_polynomial_list[j] + O(u^preCision))(pp*u) + O(u^preCision)
    return ans

def ListOfPowersOf(ff, n, preCision):
    r"""
    Returns a list of all powers up to the 'n'-th power of the fiven parameter 'ff', each power is computed up to precision 'preCision'.

    INPUT: - 'ff' an element of 'Fq[T][[u]]'; 'n', 'preCision' positive integers.

    OUTPUT: - A list of elements of 'Fq[T][[u]]'.

    EXAMPLES: sage: q
    	      3
	      sage: ff = u + T
	      sage: ListOfPowersOf(ff, 20, 10)[0]
	      1
	      sage: ListOfPowersOf(ff, 20, 10)[1] == ff
	      True
	      sage: ListOfPowersOf(ff, 20, 10)[2] == ff^2
	      True
	      sage: ListOfPowersOf(ff, 20, 10)[3] == ff^3
	      True
	      sage: ListOfPowersOf(ff, 20, 10)[10] == ff^10
	      True
	      sage: ListOfPowersOf(ff, 20, 10)[11] == ff^11
	      True
	      sage:
    """
    ff_new = ff + O(u^preCision)
    Ans = [1*u^0]
    for j in range(0, n):
    	Ans = Ans + [Ans[j]*ff_new + O(u^preCision)]
    return Ans

def BasisFor_M_k_m(k, m, preCision, GG = [], HH = []):
    r"""
    Returns basis in the form of a list for the space 'M_{k, m} (GL_2 (A))'. The basis is computed to precision 'preCision'. Uses 'ListOfPowersOf', 'gg', 'hh'.

    INPUT: - 'k', 'm', 'preCision' positive integers; 'GG', 'HH' list of elements of 'Fq[T][[u]]' (optional inputs).

    OUTPUT: - A list of elements in 'Fq[T][[u]]'.

    EXAMPLES: sage: q
    	      3
	      sage: BasisFor_M_k_m(4, 1, 50)[0] == hh(3, 50)
	      True
	      sage: %time basis1 = BasisFor_M_k_m(22, 1, 200);
	      CPU times: user 216.51 s, sys: 0.59 s, total: 217.09 s
	      Wall time: 217.49 s
	      sage: g = gg(4, 200)
	      sage: h = hh(4, 200)
	      sage: GG = ListOfPowersOf(g, 20, 200)
	      sage: HH = ListOfPowersOf(h, 20, 200)
	      sage: %time basis2 = BasisFor_M_k_m(22, 1, 200, GG, HH);
	      CPU times: user 2.47 s, sys: 0.00 s, total: 2.48 s
	      Wall time: 2.48 s
	      sage: basis1 == basis2
	      True
	      sage:
    """
    auxilary_precision = ceil(log(preCision, q))-1
    biggest_power_of_g = floor(k/(q-1))
    biggest_power_of_h = floor((k-m*(q+1))/(q^2 - 1))
    if GG == []:
       GG = ListOfPowersOf(gg(auxilary_precision, preCision), biggest_power_of_g, preCision)
    if HH == []:
       HH = ListOfPowersOf(hh(auxilary_precision, preCision), biggest_power_of_h*(q-1) + m, preCision)
    Ans = []
    for j in range(0, biggest_power_of_h + 1):
    	hpower = m + j*(q-1)
	gpower = (k-m*(q+1) - j*(q^2 - 1))/(q-1)
	Ans = Ans + [GG[gpower]*HH[hpower] + O(u^preCision)]
    return Ans

def DimensionOf_M_k_m(k, m):
    r"""
    Returns the dimension of the space 'M_{k, m} (GL_2 (A))'.

    INPUT: - 'k', 'm' positive integers.

    OUTPUT: - A non-negative integer.

    EXAMPLES: sage: q
    	      3
	      sage: DimensionOf_M_k_m(4, 0)
	      1
	      sage: DimensionOf_M_k_m(4, 1)
	      1
	      sage: DimensionOf_M_k_m(2, 1)
	      0
	      sage: DimensionOf_M_k_m(3, 1)
	      0
	      sage: DimensionOf_M_k_m(22, 1)
	      3
	      sage: DimensionOf_M_k_m(22, 0)
	      3
	      sage: DimensionOf_M_k_m(24, 0)
	      4
	      sage:
    """
    if ((k - 2*m) % (q-1)) != 0:
       return 0
    biggest_power_of_h = floor((k-m*(q+1))/(q^2 - 1))
    return 1 + biggest_power_of_h


def BasisRepresentation(ff, basis):
    r"""
    Returns a vector representation of 'ff' in the basis 'basis'.

    INPUT: - 'ff' an element of 'Fq[T][[u]]'; 'basis' a list of elements in 'Fq[T][[u]]'. Note that 'basis' has to be computed in the same manner as the computation in 'BasisFor_M_k_m' above.

    OUTPUT: - A vector with elements in 'Fq[T]'.

    EXAMPLES: sage: q
    	      3
	      sage: basis = BasisFor_M_k_m(22, 1, 200)
	      sage: BasisRepresentation(basis[0] + T^2*basis[1], basis)
	      (1, T^2, 0)
	      sage: BasisRepresentation(T*basis[0] + T^2*basis[1] - T^3*basis[2], basis)
	      (T, T^2, 2*T^3)
	      sage: BasisRepresentation(T*basis[0] + T^2*basis[1] - T^3*basis[0], basis)
	      (2*T^3 + T, T^2, 0)
	      sage:
    """
    #Dimensions are flexible here if one needs higher precision in u...
    d = len(basis)
    Matrix_space = MatrixSpace(FractionField(P), d, d)
    Vector_space = VectorSpace(FractionField(P), d)
    m = basis[0].valuation()
    matrix_list = []
    vector_list = []
    for j in range(0, d):
    	for i in range(0, d):
	    matrix_list = matrix_list + [basis[j][m + i*(q-1)]]
    for i in range(0, d):
    	vector_list = vector_list + [ff[m + i*(q-1)]]
    return Matrix_space(matrix_list).solve_left(Vector_space(vector_list))

def HeckeMatrix(pp, k, m, basis = [], Goss_polynomial_list = []):
    r"""
    Returns a matrix representing the action of 'T_pp' on 'M_{k, m} (GL_2(A))' relative to the basis 'basis'. Uses 'GossPolynomialList', 'BasisFor_M_k_m', 'HeckeAction', 'BasisRepresentation', 'DimensionOf_M_k_m'.

    INPUT: - 'pp' a monic prime in 'Fq[T]'; 'k', 'm' positive integers; 'basis', 'Goss_polynomial_list' lists of elements in 'Fq[T][[u]]' (optional inputs).

    OUTPUT: - A matrix with entries in 'Fq(T)'.

    EXAMPLES: sage: q
    	      3
	      sage: %time hecke_matrix1 = HeckeMatrix(T, 22, 1);
	      CPU times: user 1.23 s, sys: 0.00 s, total: 1.23 s
	      Wall time: 1.23 s
	      sage: hecke_matrix1
	      [                           T                            0                            0]
	      [                T^22 + 2*T^4                          T^7                        2*T^4]
	      [                2*T^25 + T^7 2*T^16 + T^14 + T^10 + 2*T^8                2*T^7 + 2*T^5]
	      sage: basis = BasisFor_M_k_m(22, 1, 200)
	      sage: %time hecke_matrix1 = HeckeMatrix(T, 22, 1, basis);
	      CPU times: user 0.17 s, sys: 0.00 s, total: 0.17 s
	      Wall time: 0.17 s
	      sage: hecke_matrix1
	      [                           T                            0                            0]
	      [                T^22 + 2*T^4                          T^7                        2*T^4]
	      [                2*T^25 + T^7 2*T^16 + T^14 + T^10 + 2*T^8                2*T^7 + 2*T^5]
	      sage: %time HeckeMatrix(T, 22, 1);
	      CPU times: user 8.91 s, sys: 0.01 s, total: 8.92 s
	      Wall time: 8.92 s
	      sage: %time HeckeMatrix(T, 42, 1);
	      CPU times: user 11.53 s, sys: 0.01 s, total: 11.54 s
	      Wall time: 11.54 s
	      sage: %time HeckeMatrix(T, 62, 1);
	      CPU times: user 81.45 s, sys: 0.08 s, total: 81.52 s
	      Wall time: 81.53 s
	      sage: %time HeckeMatrix(T^2 + 1, 42, 1);
	      CPU times: user 126.75 s, sys: 0.13 s, total: 126.88 s
	      sage:
    """
    d = pp.degree()
    dim = DimensionOf_M_k_m(k, m)
    required_precision = (m + dim*(q-1)+2)*q^d + q^3 + (q-1)
    if basis == []:
       basis = BasisFor_M_k_m(k, m, required_precision)
    if Goss_polynomial_list == []:
       Goss_polynomial_list = GossPolynomialList(pp, required_precision)
    dim = len(basis)
    Matrix_space = MatrixSpace(FractionField(P), dim, dim)
    matrix_list = []
    for i in range(0, dim):
    	matrix_list = matrix_list + [BasisRepresentation(HeckeAction(pp, basis[i], k, m,  m + dim*(q-1) + 1, Goss_polynomial_list), basis)]
    return Matrix_space(matrix_list)

def UpMatrix(pp, k, m, basis = [], Goss_polynomial_list = []):
    r"""
    Returns a matrix representing the action of 'U_pp' on 'M_{k, m} (GL_2(A))' relative to the basis 'basis'. Uses 'GossPolynomialList', 'BasisFor_M_k_m', 'HeckeAction', 'BasisRepresentation', 'DimensionOf_M_k_m'.

    INPUT: - 'pp' a monic prime in 'Fq[T]'; 'k', 'm' positive integers; 'basis', 'Goss_polynomial_list' lists of elements in 'Fq[T][[u]]' (optional inputs).

    OUTPUT: - A matrix with entries in 'Fq(T)'.

    EXAMPLES: sage: q
    	      3
	      sage: %time hecke_matrix1 = HeckeMatrix(T, 22, 1);
	      CPU times: user 1.23 s, sys: 0.00 s, total: 1.23 s
	      Wall time: 1.23 s
	      sage: hecke_matrix1
	      [                           T                            0                            0]
	      [                T^22 + 2*T^4                          T^7                        2*T^4]
	      [                2*T^25 + T^7 2*T^16 + T^14 + T^10 + 2*T^8                2*T^7 + 2*T^5]
	      sage: basis = BasisFor_M_k_m(22, 1, 200)
	      sage: %time hecke_matrix1 = HeckeMatrix(T, 22, 1, basis);
	      CPU times: user 0.17 s, sys: 0.00 s, total: 0.17 s
	      Wall time: 0.17 s
	      sage: hecke_matrix1
	      [                           T                            0                            0]
	      [                T^22 + 2*T^4                          T^7                        2*T^4]
	      [                2*T^25 + T^7 2*T^16 + T^14 + T^10 + 2*T^8                2*T^7 + 2*T^5]
	      sage: %time HeckeMatrix(T, 22, 1);
	      CPU times: user 8.91 s, sys: 0.01 s, total: 8.92 s
	      Wall time: 8.92 s
	      sage: %time HeckeMatrix(T, 42, 1);
	      CPU times: user 11.53 s, sys: 0.01 s, total: 11.54 s
	      Wall time: 11.54 s
	      sage: %time HeckeMatrix(T, 62, 1);
	      CPU times: user 81.45 s, sys: 0.08 s, total: 81.52 s
	      Wall time: 81.53 s
	      sage: %time HeckeMatrix(T^2 + 1, 42, 1);
	      CPU times: user 126.75 s, sys: 0.13 s, total: 126.88 s
	      sage:
    """
    d = pp.degree()
#    dim = DimensionOf_M_k_m(k, m)
    dim = len(basis)
    required_precision = (m + dim*(q-1)+2)*q^d + q^3 + (q-1)
    if basis == []:
       basis = BasisFor_M_k_m(k, m, required_precision)
    if Goss_polynomial_list == []:
       Goss_polynomial_list = GossPolynomialList(pp, required_precision)
#    dim = len(basis)
    Matrix_space = MatrixSpace(FractionField(P), dim, dim)
    matrix_list = []
    for i in range(0, dim*(q-1)):
    	matrix_list = matrix_list + [BasisRepresentation(UpAction(pp, basis[i], k, m,  m + dim*(q-1) + 1, Goss_polynomial_list), basis)]
    return Matrix_space(matrix_list)

def FactorWithSingular(poly):
    r"""
    Returns the factorization of a polynomial in 'Fq(T)[X]'. Uses SINGULAR. Currently only works when q is prime.

    INPUT: - 'poly' a polynomial in 'Fq(T)[X]'.

    OUTPUT: - A generic list.

    EXAMPLES: sage: q
    	      3
	      sage: hecke_matrix = HeckeMatrix(T, 22, 1);
	      sage: FactorWithSingular(hecke_matrix.charpoly())
	      [(1, 1), (T + 2*u, 1), (T^20 + 2*T^18 + 2*T^12 + 2*T^5*u + 2*u^2, 1)]
	      sage: FactorWithSingular(HeckeMatrix(T^2 + 1, 22, 1).charpoly())
	      [(2, 1),
 	      (T^2 + 1 + 2*u, 1),
 	      (T^40 + T^36 + 2*T^34 + 2*T^32 + 2*T^28 + 2*T^26 + 2*T^24 + 2*T^22 + 2*T^18 + 2*T^16 + 2*T^14 + 2*T^10 + 2*T^8 + 2*T^6 + T^4 + 1 + (T^20 + 2*T^16 + T^14 + 2*T^12 + T^4 + 1)*u + u^2,
  	      1)]
	      sage:
    """
    if q != F.characteristic():
    	raise NotImplementedError(" Currently this routine only works when 'q' is prime. ")
    SingRing = singular.ring(P.characteristic(), '(X, T)', 'dp')
    poly_sing = poly(X)
    factorization_list = poly_sing._singular_().factorize()
    Ans = []
    ans = poly(X)
    for j in range(1, len(factorization_list[1]) + 1):
    	fac = P2(factorization_list[1][j].sage_polystring())
	expo = factorization_list[2][j]
	if (ans % (fac^expo)) == 0:
	   ans = P2(ans/(fac^expo))
	   Ans = Ans + [(fac, expo)]
    if ans.degree() == 0:
       return Ans
    if ans.degree() != 0:
       print 'special'
       return Ans + [(ans, 1, 'special')]

def Eigenforms(maTrix, eigenValue):
    r"""
    Returns a list of eigenvectors for 'maTrix' with eigenvalue 'eigenValue'.

    INPUT: - 'maTrix' a square matrix with entries in 'Fq[T]'; 'eigenValue' an element of 'Fq[T]'.

    OUTPUT: - A generic list.

    EXAMPLES: sage: q
    	      3
	      sage: hecke_matrix = HeckeMatrix(T, 24, 1)
	      sage: FactorWithSingular(hecke_matrix.charpoly())
	      [(2, 1), (T + 2*u, 1), (T^7 + 2*u, 1), (T^3 + 2*u, 1)]
	      sage: Eigenforms(hecke_matrix, T)
	      Vector space of degree 3 and dimension 1 over Fraction Field of Univariate Polynomial Ring in T over Finite Field of size 3
	      Basis matrix:
	      [                                          1 2*T^21 + T^19 + 2*T^13 + T^11 + T^9 + 2*T^5                           T^18 + T^12 + T^6]
	      sage: Eigenforms(hecke_matrix, T^7)
	      Vector space of degree 3 and dimension 1 over Fraction Field of Univariate Polynomial Ring in T over Finite Field of size 3
	      Basis matrix:
	      [        0         1 T^3 + 2*T]
	      sage:
    """
    dim = len(maTrix.rows())
    return kernel((maTrix - eigenValue*MatrixSpace(FractionField(P), dim, dim).identity_matrix()).transpose())

def hdif_sub1(r, i):
    r"""
    Returns a list of powers 'i_1, ..., i_r' such that 'i = q^{i_1} + .... + q^{i_r}'. Used by 'Hdiff' below.

    INPUT: - 'r' and 'i' two positive integers.

    OUTPUT: - A generic list.

    EXAMPLES: sage: q
			3
			sage: hdif_sub1(2, 4)
			[[0, 1], [1, 0]]
			sage: hdif_sub1(1, 4)
			[]
			sage: hdif_sub1(4, 4)
			[[0, 0, 0, 0]]
			sage: hdif_sub1(4, 20)
			[[0, 0, 2, 2],
 			[0, 2, 0, 2],
 			[0, 2, 2, 0],
 			[2, 0, 0, 2],
			[2, 0, 2, 0],
			[2, 2, 0, 0]]
			sage:
    """
    if i == 0:
       return []
    if r == 1:
       if log(i, F.characteristic()) in ZZ:
       	  return [[log(i, F.characteristic())]]
       return []
    ans = []
    for j in range(0, floor(log(i, F.characteristic()))+1):
        for l in hdif_sub1(r-1, i-q^j):
	    ans = ans + [[j] + l]
    return ans

def hdif_sub2(List):
    r"""
    Returns an element of 'Fq(T)', which is used in 'Hdiff' below.

    INPUT: - 'List' a list of integers.

    OUTPUT: - an element of 'Fq(T)'
    """
    if List == []:
       return 0
    ans = 0
    for j in List:
    	term = 1
    	for l in j:
    	    term = term*(1/D(l))
        ans = ans + term
    return ans

def Hdiff(nn, ff, preCision = []):
    r"""
    Returns the 'nn'-th hyper derivative of 'ff', where 'ff' is a u-expansion up to precision 'preCision'.

    INPUT: -'nn' a positive integer; 'ff' an element of 'Fq(T)[[u]]'; 'precision' a positive integer that is an optional parameter.

    OUTPUT: An element of 'Fq(T)[[u]]'.

    EXAMPLES: sage: q = 3

    sage: E = Ee(3, 100)

    Warning: Cannot ensure the required precision. Choose 'deg', so that 'preCision <= q^(deg + 1)'.

    sage: Hdiff(0, E, 100) - E + O(u^100)

    sage: Hdiff(1, E, 100) - E^2 + O(u^100)

    sage: Hdiff(2, E, 100) - E^3*(T^36 + 2*T^30 + 2*T^18 + T^12 + 2)*u^99 + O(u^100)
    """
    if nn == 0:
    	return ff
    ans = [0,0]
    if preCision == []:
       preCision = ff.degree() + 1
    for j in range(2, min(ff.degree(), preCision)):
    	coef_j = 0
	for r in range(1, j):
	    coef_j = coef_j + (-1)^(nn+r)*binomial(j-1, r)*hdif_sub2(hdif_sub1(r, nn))*ff[j-r]
    	ans = ans + [coef_j]
    return P1(ans) + O(u^preCision)
