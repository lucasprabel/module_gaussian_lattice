#ifndef __COMMON_H__
#define __COMMON_H__

#include <inttypes.h>
#include <complex.h>
#include <math.h>

//	-------------------
//	DEFINING PARAMETERS
//	-------------------

/*// q = 268432897, d = 4, NIST I
#define PARAM_Q 268432897 // modulus q
#define PARAM_K 28 // size of q
#define PARAM_N 256 // degree of polynomials
#define PARAM_D 4 // rank of the module
#define PARAM_SIGMA 4.95 // Gaussian parameter (generation of the trapdoor)
#define PARAM_ZETA (57676.0) // Gaussian parameter (presampling)
#define PARAM_T 12 // Tailcut

#define MONT 40944
#define MONT2 65813754
#define QINV 664532479
#define BARRETT_MULT 16
#define BARRETT_SHIFT 32
*/
/*
#define PARAM_Q 1073738753 // modulus q
#define PARAM_K 30 // size of q
#define PARAM_N 256 // degree of polynomials

#define PARAM_D 4 // rank of the module

#define PARAM_SIGMA 7.00 // Gaussian parameter (generation of the trapdoor)

#define PARAM_ALPHA (54.35) // Gaussian parameter (sampling perturbations)
#define PARAM_ZETA (83832.0) // Gaussian parameter (presampling)

#define PARAM_T 12 // Tailcut
*/


/*// q = 7681, n = 256, d = 1
#define PARAM_Q 7681 // modulus q
#define PARAM_K 13 // size of q
#define PARAM_N 256 // degree of polynomials
#define PARAM_D 1 // rank of the module
#define PARAM_SIGMA 7.00 // Gaussian parameter (generation of the trapdoor)
#define PARAM_ALPHA (48.34) // Gaussian parameter (sampling perturbations)
#define PARAM_ZETA (83832.0) // Gaussian parameter (presampling)
#define PARAM_T 12 // Tailcut
*/

/*// q = 13313, n = 512, d = 1
#define PARAM_Q 13313 // modulus q
#define PARAM_K 14 // size of q
#define PARAM_N 512 // degree of polynomials
#define PARAM_D 1 // rank of the module
#define PARAM_SIGMA 7.00 // Gaussian parameter (generation of the trapdoor)
#define PARAM_ALPHA (48.34) // Gaussian parameter (sampling perturbations)
#define PARAM_ZETA (83832.0) // Gaussian parameter (presampling)
#define PARAM_T 12 // Tailcut
*/

/*// q = 12289, n = 1024, d = 1
#define PARAM_Q 12289 // modulus q
#define PARAM_K 14 // size of q
#define PARAM_N 1024 // degree of polynomials
#define PARAM_D 1 // rank of the module
#define PARAM_SIGMA 7.00 // Gaussian parameter (generation of the trapdoor)
#define PARAM_ALPHA (48.34) // Gaussian parameter (sampling perturbations)
#define PARAM_ZETA (83832.0) // Gaussian parameter (presampling)
#define PARAM_T 12 // Tailcut
*/

// q = 1073707009, n = 1024, d = 1
#define PARAM_Q 1073707009 // modulus q
#define PARAM_K 30 // size of q
#define PARAM_N 1024 // degree of polynomials
#define PARAM_D 1 // rank of the module
#define PARAM_SIGMA 7.00 // Gaussian parameter (generation of the trapdoor)
#define PARAM_ALPHA (48.34) // Gaussian parameter (sampling perturbations)
#define PARAM_ZETA (83832.0) // Gaussian parameter (presampling)
#define PARAM_T 12 // Tailcut


/*// q = 1073707009, n = 512, d = 2
#define PARAM_Q 1073707009 // modulus q
#define PARAM_K 30 // size of q
#define PARAM_N 512 // degree of polynomials
#define PARAM_D 2 // rank of the module
#define PARAM_SIGMA 7.00 // Gaussian parameter (generation of the trapdoor)
#define PARAM_ALPHA (48.34) // Gaussian parameter (sampling perturbations)
#define PARAM_ZETA (83832.0) // Gaussian parameter (presampling)
#define PARAM_T 12 // Tailcut
*/

/*// q = 1073707009, n = 256, d = 4
#define PARAM_Q 1073707009 // modulus q
#define PARAM_K 30 // size of q
#define PARAM_N 256 // degree of polynomials
#define PARAM_D 4 // rank of the module
#define PARAM_SIGMA 7.00 // Gaussian parameter (generation of the trapdoor)
#define PARAM_ALPHA (48.34) // Gaussian parameter (sampling perturbations)
#define PARAM_ZETA (83832.0) // Gaussian parameter (presampling)
#define PARAM_T 12 // Tailcut
*/

/*// q = 1073738753, d = 2
#define PARAM_Q 1073738753 // modulus q
#define PARAM_K 30 // size of q
#define PARAM_N 256 // degree of polynomials
#define PARAM_D 4 // rank of the module
#define PARAM_SIGMA 7.00 // Gaussian parameter (generation of the trapdoor)
#define PARAM_ALPHA (48.34) // Gaussian parameter (sampling perturbations)
#define PARAM_ZETA (83832.0) // Gaussian parameter (presampling)
#define PARAM_T 12 // Tailcut
*/

/*// q = 1073738753, d = 4
#define PARAM_Q 1073738753 // modulus q
#define PARAM_K 30 // size of q
#define PARAM_N 256 // degree of polynomials
#define PARAM_D 4 // rank of the module
#define PARAM_SIGMA 7.00 // Gaussian parameter (generation of the trapdoor)
#define PARAM_ALPHA (48.34) // Gaussian parameter (sampling perturbations)
#define PARAM_ZETA (83832.0) // Gaussian parameter (presampling)
#define PARAM_T 12 // Tailcut
*/


/*// q = 1073738753, d = 5
#define PARAM_Q 1073738753 // modulus q
#define PARAM_K 30 // size of q
#define PARAM_N 256 // degree of polynomials
#define PARAM_D 5 // rank of the module
#define PARAM_SIGMA 5.55 // Gaussian parameter (generation of the trapdoor)
#define PARAM_ALPHA (54.35) // Gaussian parameter (sampling perturbations)
#define PARAM_ZETA (83290.0) // Gaussian parameter (presampling)
#define PARAM_T 14 // Tailcut
*/

/*// q = 1073738753, d = 6
#define PARAM_Q 1073738753 // modulus q
#define PARAM_K 30 // size of q
#define PARAM_N 256 // degree of polynomials
#define PARAM_D 6 // rank of the module
#define PARAM_SIGMA 6.15 // Gaussian parameter (generation of the trapdoor)
#define PARAM_ALPHA (60.50) // Gaussian parameter (sampling perturbations)
#define PARAM_ZETA (112522.0) // Gaussian parameter (presampling)
#define PARAM_T 15 // Tailcut
*/

/*// q = 1073738753, d = 8
#define PARAM_Q 1073738753 // modulus q
#define PARAM_K 30 // size of q
#define PARAM_N 256 // degree of polynomials
#define PARAM_D 8 // rank of the module
#define PARAM_SIGMA 6.15 // Gaussian parameter (generation of the trapdoor)
#define PARAM_ALPHA (60.50) // Gaussian parameter (sampling perturbations)
#define PARAM_ZETA (500000.0) // Gaussian parameter (presampling)
#define PARAM_T 15 // Tailcut
*/

/*// q = 1073738753, d = 10
#define PARAM_Q 1073738753 // modulus q
#define PARAM_K 30 // size of q
#define PARAM_N 256 // degree of polynomials
#define PARAM_D 10 // rank of the module
#define PARAM_SIGMA 6.15 // Gaussian parameter (generation of the trapdoor)
#define PARAM_ALPHA (60.50) // Gaussian parameter (sampling perturbations)
#define PARAM_ZETA (500000.0) // Gaussian parameter (presampling)
#define PARAM_T 15 // Tailcut
*/

#if (PARAM_Q == 1073738753)

	#define MONT 12284
	#define MONT2 150896656
	#define QINV 2138043391
	#define BARRETT_MULT 4
	#define BARRETT_SHIFT 32

#elif (PARAM_Q == 1073707009)

	#define MONT 139260
	#define MONT2 66621438
	#define QINV 4156520447
	#define BARRETT_MULT 4
	#define BARRETT_SHIFT 32

#elif (PARAM_Q == 13313)

	#define MONT 7114
	#define MONT2 6283
	#define QINV 896545791
	#define BARRETT_MULT 322614
	#define BARRETT_SHIFT 32

#elif (PARAM_Q == 12289)

	#define MONT 10952
	#define MONT2 5664
	#define QINV 4143984639
	#define BARRETT_MULT 349496
	#define BARRETT_SHIFT 32

#elif (PARAM_Q == 7681)

	#define MONT 5569
	#define MONT2 5564
	#define QINV 1954291199
	#define BARRETT_MULT 559167
	#define BARRETT_SHIFT 32

#endif


/*// q = 2147483137
#define PARAM_Q 2147483137 // modulus q
#define PARAM_K 31 // size of q
#define PARAM_N 256 // degree of polynomials
#define PARAM_D 3 // rank of the module
#define PARAM_SIGMA 4.5 // Gaussian parameter (generation of the trapdoor)
#define PARAM_ZETA (4441.0) // Gaussian parameter (presampling)
#define PARAM_T 12 // Tailcut

#define MONT 1022
#define MONT2 1044484
#define QINV 2013003263
#define BARRETT_MULT 2
#define BARRETT_SHIFT 32
#define DOUBLE_BARRETT_MULT 8
#define DOUBLE_BARRETT_SHIFT 34
*/

/*// q = 1073738753
#define PARAM_Q 1073738753 // modulus q
#define PARAM_K 30 // size of q
#define PARAM_N 256 // degree of polynomials
#define PARAM_D 3 // rank of the module
#define PARAM_SIGMA 4.5 // Gaussian parameter (generation of the trapdoor)
#define PARAM_ZETA (4441.0) // Gaussian parameter (presampling)
#define PARAM_T 12 // Tailcut

#define MONT 12284
#define MONT2 150896656
#define QINV 2138043391
#define BARRETT_MULT 4
#define BARRETT_SHIFT 32
#define DOUBLE_BARRETT_MULT 16
#define DOUBLE_BARRETT_SHIFT 34
*/

#define NTT_ALWAYS_REDUCE
#define ALWAYS_CSUBQ


/*// q = 8380417
#define PARAM_Q 8380417 // modulus q
#define PARAM_K 23 // size of q
#define PARAM_N 256 // degree of polynomials
#define PARAM_R 8 // number of irreducible factors of x^n + 1 in F_q[x]
#define PARAM_D 3 // rank of the module
#define PARAM_SIGMA 4.5 // Gaussian parameter (generation of the trapdoor)
#define PARAM_ZETA (4441.0) // Gaussian parameter (presampling)
#define PARAM_T 12 // Tailcut

#define MONT 4193792
#define MONT2 2365951
#define QINV 4236238847
#define BARRETT_MULT 512
#define BARRETT_SHIFT 32

#define NTT_ALWAYS_REDUCE
*/

/*// q = 32257
#define PARAM_Q 32257 // modulus q
#define PARAM_K 15 // size of q
#define PARAM_N 256 // degree of polynomials
#define PARAM_D 3 // rank of the module
#define PARAM_SIGMA 4.5 // Gaussian parameter (generation of the trapdoor)
#define PARAM_ZETA 9129.0 // Gaussian parameter (presampling)
#define PARAM_T 12 // Tailcut

#define MONT 12260 // 2^32 % Q
#define MONT2 22237 // 2^64 % Q
#define QINV 3120332287 // -q^(-1) mod 2^32
*/

#define NEEDED_XOF_BYTES (2 * PARAM_N * PARAM_D * PARAM_K)
#define NEEDED_AES_FOR_SALT 3
#define SALT_BYTES (16 * NEEDED_AES_FOR_SALT)


/*// Pauline parameter set, 81-bit classical security
#define PARAM_Q 134218289 // modulus q
#define PARAM_K 28 // size of q
#define PARAM_N 256 // degree of polynomials
#define PARAM_R 8 // number of irreducible factors of x^n + 1 in F_q[x]
#define PARAM_D 6 // rank of the module
#define PARAM_SIGMA 4.5 // Gaussian parameter (generation of the trapdoor)
#define PARAM_ZETA 9129.0 // Gaussian parameter (presampling)
#define PARAM_T 12 // Tailcut
*/

/*// Pauline parameter set, 97-bit classical security
#define PARAM_Q 134218289 // modulus q
#define PARAM_K 28 // size of q
#define PARAM_N 256 // degree of polynomials
#define PARAM_R 7 // number of irreducible factors of x^n + 1 in F_q[x]
#define PARAM_D 6 // rank of the module
#define PARAM_SIGMA 4.5 // Gaussian parameter (generation of the trapdoor)
#define PARAM_ZETA 12165.0 // Gaussian parameter (presampling)
#define PARAM_T 12 // Tailcut
*/

//	------------------
//	DERIVED PARAMETERS
//	------------------

#define PARAM_M (PARAM_D * (PARAM_K + 2)) // parameter m ( = d(k+2) for the computational instantiation)
//#define PARAM_ALPHA (10.0 * PARAM_SIGMA) // Gaussian parameter (sampling perturbations)

/*// Parameters for n = 256, q = 8388593
#define PARAM_Q 8388593 // modulus q
#define PARAM_K 23 // size of q
#define PARAM_N 256 // degree of polynomials
#define PARAM_R 8 // number of irreducible factors of x^n + 1 in F_q[x]
#define SMALL_DEGREE (PARAM_N / PARAM_R) // degree of irreducible factors of x^n + 1 in F_q[x]
*/

/*
// Parameters for n = 8, q = 8000053
#define PARAM_Q 8000053 // modulus q
#define PARAM_K 23 // size of q
#define PARAM_N 8 // degree of polynomials
#define PARAM_R 2 // number of irreducible factors of x^n + 1 in F_q[x]
#define SMALL_DEGREE (PARAM_N / PARAM_R) // degree of irreducible factors of x^n + 1 in F_q[x]
*/

/*// Parameters for n = 8, q = 13
#define PARAM_Q 13 // modulus q
#define PARAM_K 4 // size of q
#define PARAM_N 8 // degree of polynomials
#define PARAM_R 2 // number of irreducible factors of x^n + 1 in F_q[x]
#define SMALL_DEGREE (PARAM_N / PARAM_R) // degree of irreducible factors of x^n + 1 in F_q[x]
*/

/*// Parameters for n = 2, q = 13
#define PARAM_Q 13 // modulus q
#define PARAM_K 4 // size of q
#define PARAM_N 2 // degree of polynomials
#define PARAM_R 2 // number of irreducible factors of x^n + 1 in F_q[x]
#define SMALL_DEGREE (PARAM_N / PARAM_R) // degree of irreducible factors of x^n + 1 in F_q[x]
*/

/*// Parameters for n = 256, q = 17
#define PARAM_Q 17 // modulus q
#define PARAM_K 5 // size of q
#define PARAM_N 256 // degree of polynomials
#define PARAM_R 8 // number of irreducible factors of x^n + 1 in F_q[x]
#define SMALL_DEGREE (PARAM_N / PARAM_R) // degree of irreducible factors of x^n + 1 in F_q[x]
*/

/*// Parameters for n = 256, q = 7681
#define PARAM_Q 7681 // modulus q
#define PARAM_K 13 // size of q
#define PARAM_N 256 // degree of polynomials
#define PARAM_R 256 // number of irreducible factors of x^n + 1 in F_q[x]
#define SMALL_DEGREE (PARAM_N / PARAM_R) // degree of irreducible factors of x^n + 1 in F_q[x]
*/

/*// Parameters for n = 256, q = 1032193
#define PARAM_Q 1032193 // modulus q
#define PARAM_K 20 // size of q
#define PARAM_N 256 // degree of polynomials
#define PARAM_R 256 // number of irreducible factors of x^n + 1 in F_q[x]
#define SMALL_DEGREE (PARAM_N / PARAM_R) // degree of irreducible factors of x^n + 1 in F_q[x]
*/


// LOG_N = base 2 log of n, defines the depth of the complex factorisation tree
#if (PARAM_N == 1024)
	#define LOG_N 10
#elif (PARAM_N == 512)
	#define LOG_N 9
#elif (PARAM_N == 256)
	#define LOG_N 8
#elif (PARAM_N == 8)
	#define LOG_N 3
#else
	#error "LOG_N is not defined for this value of PARAM_N"
#endif

#if ((1 << LOG_N) != PARAM_N)
	#error "LOG_N is not base 2 log of PARAM_N"
#endif

#if ((1 << PARAM_K) < PARAM_Q) || ((1 << (PARAM_K - 1)) >= PARAM_Q)
	#error "PARAM_K is not base 2 log of PARAM_Q"
#endif

#if (PARAM_K < 2)
	#error "PARAM_K < 2, not enough space in cplx_p_coeffs (sample_perturb)"
	#endif

#define GET_BIT(x, i) (((x) >> (i)) & 1)
#define Q_BIT(i) GET_BIT(PARAM_Q, i)

//#define DOUBLE_ZERO (((double_scalar) PARAM_Q) << (8*sizeof(double_scalar) - PARAM_K - 1)) // a double_scalar that is 0 mod q, such that there will be no underflow nor overflow during computations
#define DOUBLE_ZERO (((double_scalar) PARAM_Q) * PARAM_Q)  // a double_scalar that is 0 mod q, such that there will be no underflow nor overflow during computations
#define SIGNED_DOUBLE_ZERO (((double_scalar) PARAM_Q) << (8*sizeof(double_scalar) - PARAM_K - 2))

typedef uint32_t scalar;
typedef uint64_t double_scalar;
typedef int32_t signed_scalar;
typedef int64_t signed_double_scalar;

typedef double complex cplx;
typedef cplx *cplx_poly;
typedef cplx *cplx_poly_matrix;

#if (PARAM_Q > 0x7FFFFFFF)
	#error "Modulus PARAM_Q does not fit in 31 bits"
#endif
/*
#if (LOG_N - LOG_R + 2 * PARAM_K) >= 63
	#error "(n/r)*q^2 >= 2^63, there might be overflow in polynomial multiplication"
#endif
*/
// Polynomial in standard (coefficients) representation
typedef scalar *poly;

// Non-normalized polynomial with double_scalar coefficients
typedef double_scalar *double_poly;

// Polynomial with signed coefficients (used in Karatsuba)
typedef signed_scalar *signed_poly;

// Polynomial with signed_double_scalar coefficients (used in Karatsuba)
typedef signed_double_scalar *signed_double_poly;

// Matrix of polynomials (single pointer on scalars, everything is accessed through macros)
typedef scalar *poly_matrix;

// Matrix of signed polynomials, used when sampling since Gaussian values are signed
typedef signed_scalar *signed_poly_matrix;

#define poly_matrix_element(M, nb_col, i, j) (&M[(PARAM_N)*(((i)*(nb_col)) + (j))])

#define triangular_poly_matrix_element(M, i, j) (&M[(PARAM_N) * ((i)*((i)+1)/2 + (j))])

extern cplx cplx_roots_of_unity[2*PARAM_N - 1];

//#define cplx_root_of_unity_of_order(k, i) (cplx_roots_of_unity[(2 * (PARAM_N) / (k) * (i))])

#define cplx_cyclotomic_factorisation_tree(i, j) (cplx_roots_of_unity[(1 << (i)) + (j) - 1])

/*
	Stuff for Gaussian sampling
*/
typedef long double RR_t;
typedef double real;
#define LDRMX ((RR_t) RAND_MAX)
#define LOG_2 0.6931471805599453094172321214581765680755001343602552541206800094933936219696947156058633269964186875420014810205706857336855202357581305570326707516L
#define SIGMA_1 0.84932180028801904272150283410288961971514109378435394286159953238339383120795466719298223538163406787061691601172910413284884326532697308797136114023L //sqrt(1/(2*log(2)))

#endif
