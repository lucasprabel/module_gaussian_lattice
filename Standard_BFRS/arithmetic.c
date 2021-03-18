#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
#include <string.h>
#include <stdbool.h>

//#include "nmod_poly.h"

#include "arithmetic.h"
#include "random.h"

// Code from Dilithium

#if (PARAM_Q == 1073740609)

static const scalar zetas[PARAM_R] = {0, 1020052844, 944408404, 948631413, 403573592, 274911572, 146915528, 959538471, 635909418, 1056965504, 860001759, 282512986, 768341009, 19646585, 893480053, 15712167, 314121943, 221587804, 221703000, 1060116065, 1012438086, 815929587, 871007317, 645855111, 74611346, 572586120, 80209151, 398251412, 900326132, 631316064, 559573890, 272344229};

static const scalar zetas_inv[PARAM_R] = {801396380, 514166719, 442424545, 173414477, 675489197, 993531458, 501154489, 999129263, 427885498, 202733292, 257811022, 61302523, 13624544, 852037609, 852152805, 759618666, 1058028442, 180260556, 1054094024, 305399600, 791227623, 213738850, 16775105, 437831191, 114202138, 926825081, 798829037, 670167017, 125109196, 129332205, 53687765, 0};

#elif (PARAM_Q == 1073739937)

static const scalar zetas[PARAM_R] = {0, 1002457797, 907562409, 246308725, 741357731, 738775482, 564424612, 584648601, 513104813, 938628629, 557351748, 209762326, 773848353, 180017277, 906952126, 747433592};

static const scalar zetas_inv[PARAM_R] = {326306345, 166787811, 893722660, 299891584, 863977611, 516388189, 135111308, 560635124, 489091336, 509315325, 334964455, 332382206, 827431212, 166177528, 71282140, 0};

#endif


/*************************************************
* Name:        ntt
*
* Description: Forward NTT, in-place. No modular reduction is performed after
*              additions or subtractions. Hence output coefficients can be up
*              to 16*Q larger than the coefficients of the input polynomial.
*              Output vector is in bitreversed order.
*
* Arguments:   - uint32_t p[N]: input/output coefficient array
**************************************************/
/*void ntt(uint32_t p[PARAM_N])
	{
	unsigned int len, start, j, k;
	scalar zeta, t;
	
	#ifdef ALWAYS_CSUBQ
	scalar t1;
	#endif

	k = 1;
	for(len = PARAM_N / 2; len > SMALL_DEGREE - 1; len >>= 1)
		{
		for(start = 0; start < PARAM_N; start = j + len)
			{
			zeta = zetas[k++];
			for(j = start; j < start + len; ++j)
				{
				t = montgomery_reduce((double_scalar)zeta * p[j + len]);
				#ifdef ALWAYS_CSUBQ
				t1 = csubq(t);
				p[j + len] = csubq(p[j] + PARAM_Q - t1);
				p[j] = csubq(p[j] + t1);
				#elif defined NTT_ALWAYS_REDUCE
				p[j + len] = barrett_reduce(p[j] + (scalar) 2*PARAM_Q - t);
				p[j] = barrett_reduce(p[j] + t);
				#else
				p[j + len] = p[j] + 2*PARAM_Q - t;
				p[j] = p[j] + t;
				#endif
				}
			}
		}
	}*/

/*************************************************
* Name:        invntt_frominvmont
*
* Description: Inverse NTT and multiplication by Montgomery factor 2^32.
*              In-place. No modular reductions after additions or
*              subtractions. Input coefficient need to be smaller than 2*Q.
*              Output coefficient are smaller than 2*Q.
*
* Arguments:   - uint32_t p[N]: input/output coefficient array
**************************************************/
/*void invntt_frominvmont(uint32_t p[PARAM_N])
	{
	unsigned int start, len, j, k;
	scalar t, zeta;
	const scalar f = (((double_scalar)MONT*MONT % PARAM_Q) * (PARAM_Q-1) % PARAM_Q) * ((PARAM_Q-1) >> 8) % PARAM_Q;
	
	#ifdef ALWAYS_CSUBQ
	scalar t1;
	#endif

	k = 0;
	for(len = SMALL_DEGREE; len < PARAM_N; len <<= 1)
		{
		for(start = 0; start < PARAM_N; start = j + len)
			{
			zeta = zetas_inv[k++];
			for(j = start; j < start + len; ++j)
				{
				t = p[j];
				#ifdef ALWAYS_CSUBQ
				p[j] = csubq(t + p[j + len]);
				p[j + len] = csubq(t + PARAM_Q - p[j + len]);
				p[j + len] = montgomery_reduce((double_scalar)zeta * p[j + len]);
				p[j + len] = csubq(p[j + len]);
				#elif defined NTT_ALWAYS_REDUCE
				p[j] = barrett_reduce(t + p[j + len]);
				p[j + len] = t + 2 * PARAM_Q - p[j + len];
				p[j + len] = montgomery_reduce((double_scalar) zeta * p[j + len]);
				#else
				p[j] = t + p[j + len];
				p[j + len] = t + 256*PARAM_Q - p[j + len];
				p[j + len] = montgomery_reduce((double_scalar)zeta * p[j + len]);
				#endif
				}
			}
		}
	
	for(j = 0; j < PARAM_N; ++j)
		{
		p[j] = montgomery_reduce((double_scalar)f * p[j]);
		#ifdef ALWAYS_CSUBQ
		p[j] = csubq(p[j]);
		#endif
		}
	}*/
/*
void divide_by_2pow32(poly f, int deg)
	{
	for(int i = 0 ; i < deg + 1 ; ++i)
		{
		f[i] = montgomery_reduce(f[i]);
		}
	}

void multiply_by_2pow32(poly f, int deg)
	{
	for(int i = 0 ; i < deg + 1 ; ++i)
		{
		f[i] = montgomery_reduce((double_scalar) MONT2 * f[i]);
		}
	}

scalar barrett_reduce(double_scalar a)
	{
	scalar u = ((double_scalar) a * BARRETT_MULT) >> BARRETT_SHIFT;
	return a - u * PARAM_Q;
	}
*/
/*************************************************
* Name:        csubq
*
* Description: Subtract Q if input coefficient is bigger than Q.
*
* Arguments:   - uint32_t: finite field element a
*
* Returns r.
**************************************************/
/*scalar csubq(scalar a)
	{
	a -= PARAM_Q;
	a += ((signed_scalar)a >> 31) & PARAM_Q;
	return a;
	}*/

/*************************************************
* Name:        montgomery_reduce
*
* Description: For finite field element a with 0 <= a <= Q*2^32,
*              compute r \equiv a*2^{-32} (mod Q) such that 0 <= r < 2*Q.
*
* Arguments:   - uint64_t: finite field element a
*
* Returns r.
**************************************************/
/*uint32_t montgomery_reduce(uint64_t a)
	{
	uint64_t t;

	t = a * QINV;
	t &= (1ULL << 32) - 1;
	t *= PARAM_Q;
	t = a + t;
	t >>= 32;
	return t;
	}

void ntt_poly_mul_invmontgomery(poly h, const poly f, const poly g)
	{
	double_scalar double_prod_coeffs[2 * PARAM_R];
	double_poly double_prod = double_prod_coeffs;
	
	for(int i = 0 ; i < PARAM_R ; ++i)
		{
		poly f_i = crt_poly_component(f, SMALL_DEGREE, i);
		poly g_i = crt_poly_component(g, SMALL_DEGREE, i);
		poly h_i = crt_poly_component(h, SMALL_DEGREE, i);
		
		// prod <- f_i * g_i (non reduced)
		mul_poly_schoolbook(double_prod, f_i, g_i, SMALL_DEGREE - 1, SMALL_DEGREE - 1);
		
		// h_i <- prod (reduced mod X^(n/r) + c_i mod q)
		//barrett_upper_half_double_poly(double_prod, SMALL_DEGREE);
		//modulo_double_poly(double_crt_f_i, deg, c);
		//freeze_double_poly(crt_f_i, double_crt_f_i, deg-1);
		}
	}

void matrix_ntt(poly_matrix A, int l1, int l2)
	{
	for(int i = 0 ; i < l1 * l2 ; ++i)
		{
		poly A_ij = poly_matrix_element(A, 1, i, 0);
		
		ntt(A_ij);
		}
	}

void matrix_invntt(poly_matrix A, int l1, int l2)
	{
	for(int i = 0 ; i < l1 * l2 ; ++i)
		{
		poly A_ij = poly_matrix_element(A, 1, i, 0);
		
		invntt_frominvmont(A_ij);
		}
	}*/

/*
	Reduce double_f mod q and mod X^(n/r) + c
*/
void reduce_ntt_double_poly(poly f, double_poly double_f, scalar c)
	{
	
	}

/*
	Reduce integer input mod q in a naive way
*/
scalar reduce_naive(scalar x)
	{
	return (x % PARAM_Q);
	}

signed_scalar reduce_signed_double_naive(signed_double_scalar x)
	{
	return (x % PARAM_Q);
	}

scalar reduce_signed_double_to_positive_naive(signed_double_scalar x)
	{
	signed_scalar xmod = x % PARAM_Q;
	if(xmod < 0)
		{
		return xmod + PARAM_Q;
		}
	else
		{
		return xmod;
		}
	}

scalar reduce_signed_naive(signed_scalar x)
	{
	signed_scalar xmod = x % PARAM_Q;
	if(xmod < 0)
		{
		return xmod + PARAM_Q;
		}
	else
		{
		return xmod;
		}
	}
	
/*
	Reduce integer input mod q using the structure of q
*/
scalar reduce_sparse(scalar x)
	{
	return 0;
	}

/*
	Reduce integer input mod q
*/
scalar reduce_double_naive(double_scalar x)
	{
	return (scalar) (x % PARAM_Q);
	}

scalar reduce_double_sparse(double_scalar x)
	{
	return 0;
	}

scalar reduce_double_montgomery(double_scalar x)
	{
	return 0;
	}

/*
	Computes x_inv such that :
		x_inv * x = 1 mod q
		|x_inv| < q
*/
signed_scalar signed_scalar_inverse(signed_scalar x)
	{
	// Extended Euclid with x and n as inputs
	// starting with r0 = q, r1 = x costs one fewer loop iteration
	signed_scalar u0 = 1, u1 = 0, r0 = x, r1 = PARAM_Q, quo, t1, t2;
	
	while(r1 != 0)
		{
		quo = r0 / r1;
		t1 = u0 - quo * u1;
		t2 = r0 - quo * r1;
		u0 = u1;
		r0 = r1;
		u1 = t1;
		r1 = t2;
		}
	
	if(r0 == 1)
		{
		return u0;
		}
	else
		{
		return -u0;
		}
	}

/*
	Compute quo and rem such that :
		a = b * quo + rem
		deg(rem) < deg(b)
	where a is of degree <= deg_a and b is of degree deg_b
*/
void polynomial_division(signed_poly quo, signed_poly rem, signed_poly a, signed_poly b, int deg_a, int deg_b)
	{
	signed_double_poly r;
	signed_double_scalar r_coeffs[deg_a+1];
	int deg_r, deg_t;
	signed_scalar lc_b, lc_b_inv, lc_t;
	
	r = r_coeffs;
	
	// q <- 0
	zero_poly((poly) quo, deg_a - deg_b);
	// r <- a
	for(int i=0 ; i < deg_a+1 ; ++i)
		{
		r[i] = SIGNED_DOUBLE_ZERO + a[i];
		}
	
	// Compute degree and leading coefficient of b
	lc_b = b[deg_b];
	lc_b_inv = signed_scalar_inverse(lc_b);
	
	
	for(deg_r = deg_a ; deg_r >= deg_b ; deg_r--)
		{
		
		signed_scalar lc_r = reduce_signed_double_naive(r[deg_r]);
		
		// t <- lc(r) / lc(b) * x^(deg(r) - deg(b))
		deg_t = deg_r - deg_b;
		lc_t = reduce_signed_double_naive(((signed_double_scalar) lc_r) * lc_b_inv);
		
		
		// q <- q + t
		quo[deg_t] = lc_t;
		
		// r <- r - t*b
		for(int i=0 ; i < deg_b+1 ; ++i)
			{
			r[deg_t + i] -= ((double_scalar) lc_t) * b[i];
			}
	
		freeze_signed_double_poly_to_positive((poly) rem, r, deg_r);
		}
	
	freeze_signed_double_poly(rem, r, deg_r);
	

	}

/*
	Invert f of degree < deg modulo x^deg + c
		f is assumed to be reduced mod q
		f is assumed to be invertible mod x^deg + c
*/
void invert_poly(poly f_inv, poly f, int deg, scalar c)
	{
	signed_poly r0, r1, v0, v1, quo, tmp;
	signed_double_poly p;
	signed_scalar r0_coeffs[deg+1], r1_coeffs[deg+1], v0_coeffs[deg], v1_coeffs[deg], quo_coeffs[deg], tmp_coeffs[deg+1], *tmp_ptr;
	signed_double_scalar p_coeffs[deg];
	int deg_r0, deg_r1, deg_quo;
	
	r0 = r0_coeffs;
	r1 = r1_coeffs;
	v0 = v0_coeffs; // store v0 in f_inv ? so we can do one fewer allocation
	v1 = v1_coeffs;
	quo = quo_coeffs;
	tmp = tmp_coeffs;
	p = p_coeffs;
	
	// Compute degree of f (= r1 soon)
	for(deg_r1 = deg-1 ; f[deg_r1] == 0 ; --deg_r1);
	
	deg_r0 = deg;
	
	// (r0, r1, v0, v1) <- (x^deg + c, f, 0, 1)
	zero_poly((poly) r0, deg);
	r0[deg] = 1;
	r0[0] = c;
	memcpy(r1, f, deg * sizeof(scalar));
	zero_poly((poly) v0, deg-1);
	zero_poly((poly) v1, deg-1);
	v1[0] = 1;
	
	while(deg_r0 != 0)
		{
		memcpy(tmp, r1, (deg_r1 + 1) * sizeof(scalar));
		polynomial_division(quo, r1, r0, tmp, deg_r0, deg_r1);
		
		deg_quo = deg_r0 - deg_r1;
		
		
		tmp_ptr = r0;
		r0 = tmp;
		tmp = tmp_ptr;
		memcpy(tmp, v1, deg * sizeof(scalar));
		mul_signed_poly_schoolbook(p, quo, v1, deg_r0 - deg_r1, deg - 1 - deg_r0 + deg_r1); // can get a better bound on the degree of v1
		
		for(int i=0 ; i < deg ; ++i)
			{
			v1[i] = reduce_signed_double_naive(((signed_double_scalar) v0[i]) - p[i]);
			}
		tmp_ptr = v0;
		v0 = tmp;
		tmp = tmp_ptr;
		
		// update the degrees
		deg_r0 = deg_r1;
		if(deg_r1 != 0)
			{
			for( ; r1[deg_r1] == 0 ; --deg_r1)
				{
				//printf("r1[%d] = %" PRId32 "\n", deg_r1, r1.coeffs[deg_r1]);
				}
			}
		
		}
	
	// f_inv <- v0 / r0 (r0 is a scalar)
	signed_double_scalar r0_inv = signed_scalar_inverse(r0[0]);
	
	for(int i=0 ; i < deg ; ++i)
		{
		f_inv[i] = reduce_signed_double_to_positive_naive(v0[i] * r0_inv);
		}
	}

void alloc_poly(poly f, int deg)
	{
	f = malloc((deg + 1) * sizeof(scalar));
	}

void alloc_double_poly(double_poly f, int deg)
	{
	f = malloc((deg + 1) * sizeof(double_scalar));
	}

void free_poly(poly f)
	{
	free(f);
	}

void free_double_poly(double_poly f)
	{
	free(f);
	}

void print_poly(poly f, int deg)
	{
	printf("[");
	for(int i=0 ; i < deg ; ++i)
		{
		printf("%" PRIu32 ", ", f[i]);
		}
	printf("%" PRIu32 "]", f[deg]);
	}

void print_double_poly(double_poly f, int deg)
	{
	for(int i=0 ; i < deg + 1 ; ++i)
		{
		printf("%" PRIu64 ", ", f[i]);
		}
	printf("\n");
	}

void print_signed_poly(signed_poly f, int deg)
	{
	printf("[");
	for(int i=0 ; i < deg ; ++i)
		{
		printf("%" PRId32 ", ", f[i]);
		}
	printf("%" PRId32 "]", f[deg]);
	}

void print_signed_double_poly(signed_double_poly f, int deg)
	{
	for(int i=0 ; i < deg + 1 ; ++i)
		{
		printf("%" PRId64 ", ", f[i]);
		}
	printf("\n");
	}

void print_poly_matrix(poly_matrix A, int l1, int l2)
	{
	printf("[");
	for(int i = 0 ; i < l1 - 1 ; ++i)
		{
		printf("[");
		for(int j = 0 ; j < l2 - 1 ; ++j)
			{
			printf("#(%d, %d) :\n", i, j);
			print_poly(poly_matrix_element(A, l2, i, j), PARAM_N - 1);
			printf(",\n");
			}
		printf("#(%d, %d) :\n", i, l2 - 1);
		print_poly(poly_matrix_element(A, l2, i, l2 - 1), PARAM_N - 1);
		printf("\n],\n");
		}
	printf("[");
	for(int j = 0 ; j < l2 - 1 ; ++j)
		{
		printf("#(%d, %d) :\n", l1 - 1, j);
		print_poly(poly_matrix_element(A, l2, l1 - 1, j), PARAM_N - 1);
		printf(",\n");
		}
	printf("#(%d, %d) :\n", l1 - 1, l2 - 1);
	print_poly(poly_matrix_element(A, l2, l1 - 1, l2 - 1), PARAM_N - 1);
	printf("\n]\n]\n");
	}

void print_signed_poly_matrix(signed_poly_matrix A, int l1, int l2)
	{
	printf("[");
	for(int i = 0 ; i < l1 - 1 ; ++i)
		{
		printf("[");
		for(int j = 0 ; j < l2 - 1 ; ++j)
			{
			printf("#(%d, %d) :\n", i, j);
			print_signed_poly(poly_matrix_element(A, l2, i, j), PARAM_N - 1);
			printf(",\n");
			}
		printf("#(%d, %d) :\n", i, l2 - 1);
		print_signed_poly(poly_matrix_element(A, l2, i, l2 - 1), PARAM_N - 1);
		printf("\n],\n");
		}
	printf("[");
	for(int j = 0 ; j < l2 - 1 ; ++j)
		{
		printf("#(%d, %d) :\n", l1 - 1, j);
		print_signed_poly(poly_matrix_element(A, l2, l1 - 1, j), PARAM_N - 1);
		printf(",\n");
		}
	printf("#(%d, %d) :\n", l1 - 1, l2 - 1);
	print_signed_poly(poly_matrix_element(A, l2, l1 - 1, l2 - 1), PARAM_N - 1);
	printf("\n]\n]\n");
	}

void print_cplx_poly(cplx_poly f, int deg)
	{
	printf("[");
	for(int i=0 ; i < deg ; ++i)
		{
		printf("%f %+f * I, ", creal(f[i]), cimag(f[i]));
		}
	printf("%f %+f * I]", creal(f[deg]), cimag(f[deg]));
	}

void print_cplx_poly_matrix(cplx_poly_matrix A, int l1, int l2)
	{
	printf("[");
	for(int i = 0 ; i < l1 - 1 ; ++i)
		{
		printf("[");
		for(int j = 0 ; j < l2 - 1 ; ++j)
			{
			printf("#(%d, %d) :\n", i, j);
			print_cplx_poly(poly_matrix_element(A, l2, i, j), PARAM_N - 1);
			printf(",\n");
			}
		printf("#(%d, %d) :\n", i, l2 - 1);
		print_cplx_poly(poly_matrix_element(A, l2, i, l2 - 1), PARAM_N - 1);
		printf("\n],\n");
		}
	printf("[");
	for(int j = 0 ; j < l2 - 1 ; ++j)
		{
		printf("#(%d, %d) :\n", l1 - 1, j);
		print_cplx_poly(poly_matrix_element(A, l2, l1 - 1, j), PARAM_N - 1);
		printf(",\n");
		}
	printf("#(%d, %d) :\n", l1 - 1, l2 - 1);
	print_cplx_poly(poly_matrix_element(A, l2, l1 - 1, l2 - 1), PARAM_N - 1);
	printf("\n]\n]\n");
	}

void print_triangular_cplx_poly_matrix(cplx_poly_matrix A, int l)
	{
	for(int i = 0 ; i < l ; ++i)
		{
		for(int j = 0 ; j <= i ; ++j)
			{
			printf("#(%d, %d) :\n", i, j);
			print_cplx_poly(triangular_poly_matrix_element(A, i, j), PARAM_N - 1);
			printf("\n");
			}
		}
	}

bool equals_poly(poly f, poly g, int deg)
	{
	for(int i=0 ; i < deg + 1 ; i++)
		{
		if(f[i] != g[i])
			{
			return false;
			}
		}
	return true;
	}


bool is_zero_poly(poly f, int deg)
	{
	for(int i=0 ; i < deg + 1 ; i++)
		{
		if(f[i] != 0)
			{
			return false;
			}
		}
	return true;
	}

/*
	f <- 0
*/
void zero_poly(poly f, int deg)
	{
	memset(f, 0, (deg + 1) * sizeof(scalar));
	}

void zero_double_poly(double_poly f, int deg)
	{
	memset(f, 0, (deg + 1) * sizeof(double_scalar));
	}

void freeze_poly(poly f, int deg)
	{
	for(int i=0 ; i < deg + 1 ; ++i)
		{
		f[i] = reduce_naive(f[i]);
		}
	}

void freeze_double_poly(poly f, double_poly double_f, int deg)
	{
	for(int i=0 ; i < deg + 1; ++i)
		{
		f[i] = reduce_double_naive(double_f[i]);
		}
	}

/*
	Freeze the upper half of a poly of degree < 2*deg
		from the coeff x^deg to the coeff x^(2*deg - 1)
	IN PLACE
*/
void freeze_upper_half_double_poly(double_poly f, int deg)
	{
	for(int i=deg ; i < 2*deg ; ++i)
		{
		f[i] = reduce_double_naive(f[i]);
		}
	}

/*
	Freeze signed_poly f of degree <= deg
*/
void freeze_signed_poly(poly f_out, signed_poly f_in, int deg)
	{
	for(int i=0 ; i < deg + 1 ; ++i)
		{
		f_out[i] = reduce_signed_naive(f_in[i]);
		}
	}

void freeze_signed_double_poly(signed_poly f_out, signed_double_poly f_in, int deg)
	{
	for(int i=0 ; i < deg + 1 ; ++i)
		{
		f_out[i] = reduce_signed_double_naive(f_in[i]);
		}
	}

void freeze_signed_double_poly_to_positive(poly f_out, signed_double_poly f_in, int deg)
	{
	for(int i=0 ; i < deg + 1; ++i)
		{
		f_out[i] = reduce_signed_double_to_positive_naive(f_in[i]);
		}
	}

void random_poly(poly f, int deg)
	{
	for(int i=0 ; i < deg + 1 ; ++i)
		{
		f[i] = uniform_mod_q();
		}
	}

/*
	h <- f + g
		degrees <= deg
*/
void add_poly(poly h, poly f, poly g, int deg)
	{
	for(int i=0 ; i < deg + 1 ; ++i)
		{
		//h[i] = reduce_naive(f[i] + g[i]);
		h[i] = f[i] + g[i];
		}
	}

void add_double_poly(double_poly h, double_poly f, double_poly g, int deg)
	{
	for(int i = 0 ; i < deg + 1 ; ++i)
		{
		h[i] = f[i] + g[i];
		}
	}

void add_signed_poly(signed_poly h, signed_poly f, signed_poly g, int deg)
	{
	for(int i=0 ; i < deg + 1 ; ++i)
		{
		h[i] = f[i] + g[i];
		}
	}

void sub_signed_poly(signed_poly h, signed_poly f, signed_poly g, int deg)
	{
	for(int i=0 ; i < deg + 1 ; ++i)
		{
		h[i] = f[i] - g[i];
		}
	}

/*
	h <- f - g
	pls make sure g[i] < q
*/
void sub_poly(poly h, poly f, poly g, int deg)
	{
	for(int i=0 ; i < deg + 1 ; ++i)
		{
		h[i] = reduce_naive(PARAM_Q + f[i] - g[i]);
		}
	}

/*
	Reduce f of degree < 2*deg modulo x^deg + c
		returns a non-normalized double_poly
*/
void modulo_poly(double_poly f_out, poly f_in, int deg, scalar c)
	{
	for(int i = 0 ; i < deg ; ++i)
		{
		scalar low_coeff = f_in[i];
		double_scalar high_coeff = f_in[deg+i];
		double_scalar new_coeff = DOUBLE_ZERO + low_coeff - (c * high_coeff);
		f_out[i] = new_coeff;
		}
	}

/*
	Reduce f of degree < 2*deg modulo x^deg + c
		returns a non-normalized double_poly
	IN PLACE
*/
void modulo_double_poly(double_poly f, int deg, scalar c)
	{
	for(int i=0 ; i < deg - 1 ; ++i)
		{
		double_scalar low_coeff = f[i];
		double_scalar high_coeff = f[deg+i];
		double_scalar new_coeff = DOUBLE_ZERO + low_coeff - (c * high_coeff);
		f[i] = new_coeff;
		}
	}

void modulo_signed_poly(signed_double_poly f_out, signed_poly f_in, int deg, scalar c)
	{
	for(int i=0 ; i < deg - 1 ; ++i)
		{
		signed_scalar low_coeff = f_in[i];
		signed_double_scalar high_coeff = f_in[deg+i];
		signed_double_scalar new_coeff = DOUBLE_ZERO + low_coeff - (c * high_coeff);
		f_out[i] = new_coeff;
		}
	f_out[deg-1] = f_in[deg-1];
	}

/*
	Reduce f of degree < 2*n modulo x^n + 1
		returns a non-normalized poly
*/
void modulo_cyclotomic_poly(poly f)
	{
	for(int i=0 ; i < PARAM_N ; ++i)
		{
		f[i] += PARAM_Q - f[PARAM_N + i];
		}
	}

void modulo_cyclotomic_signed_poly(signed_poly f)
	{
	for(int i=0 ; i < PARAM_N ; ++i)
		{
		f[i] -= f[PARAM_N + i];
		}
	}

void modulo_cyclotomic_double_poly(double_poly f)
	{
	for(int i=0 ; i < PARAM_N ; ++i)
		{
		f[i] += PARAM_Q - f[PARAM_N + i];
		}
	}

/*
	Invert the CRT map
	i.e. the isomorphism from Z_q[x]/(x^(2*deg) + c) to Z_q[x]/(x^deg + c1) X Z_q[x]/(x^deg + c2)
	u and v are the (scalar) Bézout coefficients such that u*(x^deg + c1) + v*(x^deg + c2) = 1
*/
void invert_crt(poly f, poly f1, poly f2, int deg, scalar c1, scalar c2, scalar u, scalar v)
	{
	scalar uc1 = reduce_double_naive(((double_scalar) u) * c1);
	scalar vc2 = reduce_double_naive(((double_scalar) v) * c2);
	
	for(int i=0 ; i < deg ; ++i)
		{
		double_scalar coef = (((double_scalar)f2[i]) * uc1) + (((double_scalar)f1[i]) * vc2);
		f[i] = reduce_double_naive(coef);
		}
	
	for(int i=deg ; i < 2*deg ; ++i)
		{
		double_scalar coef = (((double_scalar)f2[i-deg]) * u) + (((double_scalar)f1[i-deg]) * v);
		f[i] = reduce_double_naive(coef);
		}
	}

/*
	Convert a polynomial from its coefficients representation to its CRT representation
		- IN PLACE
		- input does not need to be reduced mod q
*/
void crt_representation(poly f, int max_depth)
	{
	for(int depth = 0 ; depth < max_depth ; ++depth)
		{
		int deg = PARAM_N >> (depth + 1);
		for(int breadth = 0 ; breadth < (1 << depth) ; ++breadth)
			{
			// p is the polynomial we are treating, p1 & p2 hold results temporarily
			// p is at tree[depth][breadth]
			// p1 & p2 will be at tree[depth+1][2*breadth] and tree[depth+1][2*breadth+1]
			scalar p1_coeffs[deg], p2_coeffs[deg];
			double_scalar double_p_coeffs[deg];
			poly p = &f[2*deg*breadth], p1 = p1_coeffs, p2 = p2_coeffs;
			double_poly double_p = double_p_coeffs;
			
			// Reduce p into smaller rings
			scalar c1 = cyclotomic_factorisation_tree[depth+1][2*breadth];
			scalar c2 = cyclotomic_factorisation_tree[depth+1][2*breadth+1];
			modulo_poly(double_p, p, deg, c1);
			freeze_double_poly(p1, double_p, deg-1);
			modulo_poly(double_p, p, deg, c2);
			freeze_double_poly(p2, double_p, deg-1);
			
			
			// Store new results in place of p
			memcpy(p, p1, deg * sizeof(scalar));
			memcpy(&p[deg], p2, deg * sizeof(scalar));
			}
		}
	}

/*
	Convert a polynomial from its CRT representation to its coefficients representation
		IN PLACE
*/
void coeffs_representation(poly f, int max_depth)
	{
	for(int depth = max_depth ; depth > 0 ; --depth)
		{
		int deg = PARAM_N >> depth;
		
		for(int breadth = 0 ; breadth < (1 << depth) ; breadth += 2)
			{
			// p1 and p2 are the polynomials we are treating, p holds a result temporarily
			scalar p_coeffs[2*deg];
			poly p1, p2, p;
			p = p_coeffs;
			p1 = &f[deg*breadth];
			p2 = &f[deg*(breadth+1)];
			
			// Combine p1 and p2 into the bigger ring by inverting the CRT map
			scalar c1 = cyclotomic_factorisation_tree[depth][breadth];
			scalar c2 = cyclotomic_factorisation_tree[depth][breadth+1];
			scalar u1 = bezout_coefficients_tree[depth][breadth];
			scalar u2 = bezout_coefficients_tree[depth][breadth+1];
			invert_crt(p, p1, p2, deg, c1, c2, u1, u2);
			
			
			// Store result where we fetched p1 and p2
			memcpy(&f[deg*breadth], p, 2 * deg * sizeof(scalar));
			}
		}
	}

/*
	Multiply crt_f and crt_g in the CRT domain
		no reduction
		returns a non-normalised double_poly
*/
void mul_crt_poly(double_poly crt_h, poly crt_f, poly crt_g, int depth)
	{
	int deg = PARAM_N >> depth;
	
	// Multiply crt_f and crt_g coordinate by coordinate
	for(int i=0 ; i < (1 << depth) ; ++i)
		{
		poly f_i = crt_poly_component(crt_f, deg, i);
		poly g_i = crt_poly_component(crt_g, deg, i);
		double_poly h_i = crt_poly_component(crt_h, 2*deg, i);
		
		
		mul_poly_schoolbook(h_i, f_i, g_i, deg-1, deg-1);

		}
	}

/*
	Takes a non-normalised double_poly in the CRT domain as input
	Outputs its reduction where each component is :
		- reduced mod x^deg + c (for appropriate c)
		- reduced mod q
*/
void reduce_double_crt_poly(poly crt_f, double_poly double_crt_f, int depth)
	{
	int deg = PARAM_N >> depth;
	
	// Reduce double_crt_f coordinate by coordinate
	for(int i=0 ; i < (1 << depth) ; ++i)
		{
		poly crt_f_i = crt_poly_component(crt_f, deg, i);
		double_poly double_crt_f_i = crt_poly_component(double_crt_f, 2*deg, i);
		
		scalar c = cyclotomic_factorisation_tree[depth][i];
		
		freeze_upper_half_double_poly(double_crt_f_i, deg);
		modulo_double_poly(double_crt_f_i, deg, c);
		freeze_double_poly(crt_f_i, double_crt_f_i, deg-1);
		}
	}


/*
	Multiply f and g of degree <= deg into h of degree <= 2*deg
		returns a non-normalized double_poly
*/
void mul_poly_schoolbook(double_poly h, poly f, poly g, int deg_f, int deg_g)
	{
	int deg_h = deg_f + deg_g;
	
	memset(h, 0, (deg_h + 1) * sizeof(double_scalar));
	
	for(int i=0 ; i < deg_f + 1 ; ++i)
		{
		for(int j=0 ; j < deg_g + 1 ; ++j)
			{
			h[i+j] += ((double_scalar) f[i]) * g[j];
			}
		}

	}

/*
	Multiply f and g of degree <= deg into h of degree <= 2*deg
		returns a non-normalized double_signed_poly
*/
void mul_signed_poly_schoolbook(signed_double_poly h, signed_poly f, signed_poly g, int deg_f, int deg_g)
	{
	int deg_h = deg_f + deg_g;
	
	memset(h, 0, (deg_h + 1) * sizeof(signed_double_scalar));
	
	for(int i=0 ; i < deg_f + 1 ; ++i)
		{
		for(int j=0 ; j < deg_g + 1 ; ++j)
			{
			h[i+j] += ((signed_double_scalar) f[i]) * g[j];
			}
		}
	}

#define KARATSUBA_TO_NAIVE_CUTOFF 64
#define LOG_KARATSUBA_TO_NAIVE_CUTOFF 8
/*
	The coefficients of a, b, c can gain 2 bits per level of recursion
	Be careful, it can overflow a simple signed scalar if q * 4^(nb_levels) does not fit in a signed scalar
*/

void mul_signed_poly_karatsuba(signed_poly h, signed_poly f, signed_poly g, int deg)
	{
	// If the degree is too small do the multiplication naively
	if(deg <= KARATSUBA_TO_NAIVE_CUTOFF)
		{
		signed_double_poly double_h;
		signed_double_scalar double_h_coeffs[2*deg];
		
		double_h = double_h_coeffs;
		
		mul_signed_poly_schoolbook(double_h, f, g, deg-1, deg-1);
		
		freeze_signed_double_poly(h, double_h, 2*deg-1);
		return;
		}
	
	//printf("real karatsuba %d\n", deg);
	int small_deg = deg >> 1;
	signed_poly f0, f1, g0, g1, a, b, c;
	
	// Allocate some extra space for the polynomial b
	signed_scalar b_coeffs[deg];
	b = b_coeffs;
	
	// Split f and g into their high and low parts
	// which each have degree < small_deg = deg / 2
	f0 = f;
	f1 = &f[small_deg];
	g0 = g;
	g1 = &g[small_deg];
	
	// Store a and c into h so that we don't need to allocate space for them
	a = h;
	c = &h[deg];
	
	// Multiply recursively
	// i.e. compute a, b, c polynomials of degree < deg - 1
	// such that f * g = a * x^deg + b * x^(deg/2) + c
	add_signed_poly(a, f0, f1, small_deg-1);
	add_signed_poly(c, g0, g1, small_deg-1);

	mul_signed_poly_karatsuba(b, a, c, small_deg-1);
	
	zero_poly((poly) a, deg-1);
	zero_poly((poly) c, deg-1);
	


	mul_signed_poly_karatsuba(a, f0, g0, small_deg); // a = f0 * g0
	mul_signed_poly_karatsuba(c, f1, g1, small_deg); // c = f1 * g1
	sub_signed_poly(b, b, a, deg-2);
	sub_signed_poly(b, b, c, deg-2); // b = (f0 + f1) * (g0 + g1) - a - c
	

	
	// Combine results
	// (a and c are already stored in h so we only need to add b)
	for(int i = 0 ; i < deg-1 ; ++i)
		{
		h[small_deg + i] += b[i];
		}

	}

/*
	Multiply f and g of degree < deg (a power of 2) into h of degree < 2*deg
*/
void mul_poly_karatsuba(poly h, poly f, poly g, int deg)
	{
	// If the degree is too small do the multiplication naively
	if(deg <= KARATSUBA_TO_NAIVE_CUTOFF)
		{
		double_poly double_h;
		double_scalar double_h_coeffs[2*deg];
		double_h = double_h_coeffs;
		
		mul_poly_schoolbook(double_h, f, g, deg-1, deg-1);
		
		freeze_double_poly(h, double_h, 2*deg-1);
		
		return;
		}
	
	//printf("real karatsuba %d\n", deg);
	int small_deg = deg >> 1;
	poly f0, f1, g0, g1, a, b, c;
	
	// Allocate some extra space for the polynomial b
	scalar b_coeffs[deg];
	b = b_coeffs;
	
	// Split f and g into their high and low parts
	// which each have degree < small_deg = deg / 2
	f0 = f;
	f1 = &f[small_deg];
	g0 = g;
	g1 = &g[small_deg];

	// Store a and c into h so that we don't need to allocate space for them
	a = h;
	c = &h[deg];
	
	// Multiply recursively
	// i.e. compute a, b, c polynomials of degree < deg - 1
	// such that f * g = a * x^deg + b * x^(deg/2) + c
	add_poly(a, f0, f1, small_deg-1);
	add_poly(c, g0, g1, small_deg-1);

	mul_poly_karatsuba(b, a, c, small_deg-1);
	
	zero_poly(a, deg-1); // only set one or two coeffs to 0
	zero_poly(c, deg-1);
	


	mul_poly_karatsuba(a, f0, g0, small_deg); // a = f0 * g0
	mul_poly_karatsuba(c, f1, g1, small_deg); // c = f1 * g1
	sub_poly(b, b, a, deg-2);
	sub_poly(b, b, c, deg-2); // b = (f0 + f1) * (g0 + g1) - a - c
	

	
	// the output polynomials a, b, c are not 0 so it doesn't work when called recursively
	

	
	// Combine results
	// (a and c are already stored in h so we only need to add b)
	for(int i = 0 ; i < deg-1 ; ++i)
		{
		h[small_deg + i] = reduce_naive(h[small_deg + i] + b[i]);
		}

	}

/*
	Turn all the polynomials of the matrix A into their CRT representation
*/
void matrix_crt_representation(poly_matrix A, int l1, int l2, int depth)
	{
	for(int i=0 ; i < l1 ; ++i)
		{
		for(int j = 0 ; j < l2 ; ++j)
			{
			crt_representation(poly_matrix_element(A, l2, i, j), depth);
			}
		}
	}

/*
	Turn all the polynomials of the matrix A into their coeffs representation
*/
void matrix_coeffs_representation(poly_matrix A, int l1, int l2, int depth)
	{
	for(int i=0 ; i < l1 ; ++i)
		{
		for(int j = 0 ; j < l2 ; ++j)
			{
			coeffs_representation(poly_matrix_element(A, l2, i, j), depth);
			}
		}
	}

/*
	H <- F * G
		F is a (l1, l2) matrix in the CRT domain
		G is a (l2, l3) matrix in the CRT domain
*/
void mul_crt_poly_matrix(poly_matrix H, poly_matrix F, poly_matrix G, int l1, int l2, int l3, int depth)
	{
	scalar prod_coeffs[PARAM_N];
	poly prod = prod_coeffs;
	
	double_scalar double_prod_coeffs[2 * PARAM_N];
	double_poly double_prod = double_prod_coeffs;
	
	// H <- 0
	memset(H, 0, l1 * l3 * PARAM_N * sizeof(scalar));
	
	// ikj matrix multiplication (line by line)
	for(int i = 0 ; i < l1 ; ++i)
		{
		for(int k = 0 ; k < l2 ; ++k)
			{
			for(int j = 0 ; j < l3 ; ++j)
				{
				// H[i][j] += F[i][k] * G[k][j] (in the CRT domain)
				poly f_ik = poly_matrix_element(F, l2, i, k);
				poly g_kj = poly_matrix_element(G, l3, k, j);
				poly h_ij = poly_matrix_element(H, l3, i, j);
				
				mul_crt_poly(double_prod, f_ik, g_kj, depth);
				reduce_double_crt_poly(prod, double_prod, depth);
				add_poly(h_ij, h_ij, prod, PARAM_N - 1);
				freeze_poly(h_ij, PARAM_N - 1);
				

				}
			}
		}
	}

/*
	H <- F * G
		F is a (l1, l2) matrix in the CRT domain
		G is a (l2, l3) matrix in the CRT domain
*/
void mul_crt_poly_matrix2(poly_matrix H, poly_matrix F, poly_matrix G, int l1, int l2, int l3, int depth)
	{
	// Temporary storage for one (non reduced) line of H
	double_scalar h_line_coeffs[l3 * 2 * PARAM_N], prod_coeffs[2 * PARAM_N];
	double_poly prod = prod_coeffs, h_line = h_line_coeffs;
	
	// ikj matrix multiplication (line by line)
	for(int i = 0 ; i < l1 ; ++i)
		{
		// Set H's i-th line to 0
		memset(h_line, 0, l3 * 2 * PARAM_N * sizeof(double_scalar));
		
		for(int k = 0 ; k < l2 ; ++k)
			{
			for(int j = 0 ; j < l3 ; ++j)
				{
				// H[i][j] += F[i][k] * G[k][j] (in the CRT domain)
				poly f_ik = poly_matrix_element(F, l2, i, k);
				poly g_kj = poly_matrix_element(G, l3, k, j);
				double_poly h_ij = &h_line[j * 2 * PARAM_N];
				
				mul_crt_poly(prod, f_ik, g_kj, depth);
				add_double_poly(h_ij, h_ij, prod, 2*PARAM_N - 1);
				
				}
			}
		
		// H's i-th line is now entirely computed, reduce it and store it into H
		for(int j = 0 ; j < l3 ; ++j)
			{
			poly h_ij = poly_matrix_element(H, l3, i, j);
			double_poly h_line_j = &h_line[j * 2 * PARAM_N];
			//reduce_double_crt_poly(h_ij, h_line, depth);
			reduce_double_crt_poly(h_ij, h_line_j, depth);

			}
		}
	}

/*
	A <- A + B where A and B are (l1, l2) matrices
*/
void add_to_poly_matrix(poly_matrix A, poly_matrix B, int l1, int l2)
	{
	for(int i = 0 ; i < l1 * l2 * PARAM_N ; ++i)
		{
		A[i] += B[i];
		}
	}

/*
	Multiply v (in Z^k) by the scalar gadget vector g in Z^k
*/
void multiply_by_scalar_gadget_vector(scalar *prod, signed_scalar *v)
	{
	signed_double_scalar double_prod = 0;
	for(int j = 0 ; j < PARAM_K ; ++j)
		{
		double_prod += ((signed_double_scalar) v[j] * (1 << j));
		}
	*prod = reduce_signed_double_to_positive_naive(double_prod);
	}

/*
	Multiply v (in R^k) by the ring gadget vector g in R^k
*/
void multiply_by_ring_gadget_vector(poly prod, signed_poly_matrix v)
	{
	signed_double_scalar double_prod_coeffs[PARAM_N];
	signed_double_poly double_prod = double_prod_coeffs;
	
	memset(double_prod, 0, PARAM_N * sizeof(double_scalar));
	
	for(int i = 0 ; i < PARAM_K ; ++i)
		{
		signed_poly v_i = poly_matrix_element(v, 1, i, 0);
		// multiply v's i-th coordinate by 2^i
		for(int j = 0 ; j < PARAM_N ; ++j)
			{
			double_prod[j] += ((signed_double_scalar) v_i[j] * (1 << i));
			}
		}
	
	freeze_signed_double_poly_to_positive(prod, double_prod, PARAM_N - 1);
	}

/*
	Multiply v (in R^dk) by the module gadget matrix G in R^{d,dk}
*/
void multiply_by_module_gadget_matrix(poly_matrix prod, signed_poly_matrix v)
	{
	for(int i = 0 ; i < PARAM_D ; ++i)
		{
		signed_poly_matrix v_i = poly_matrix_element(v, 1, i*PARAM_K, 0);
		poly prod_i = poly_matrix_element(prod, 1, i, 0);
		
		multiply_by_ring_gadget_vector(prod_i, v_i);
		}
	}

/*
	y <- A * x
		A is in R^{d,m} in the CRT domain (the first d columns are the identity matrix and are not stored)
		x is in R^m in the CRT domain
		y is in R^d in the CRT domain
*/
void multiply_by_A(poly_matrix y, poly_matrix A, poly_matrix x)
	{
	// y <- A[:,d:] * x[d:]
	poly_matrix x_d = poly_matrix_element(x, 1, PARAM_D, 0);
	
	mul_crt_poly_matrix(y, A, x_d, PARAM_D, PARAM_M - PARAM_D, 1, LOG_R);
	
	// y <- y + x[:d]
	for(int i = 0 ; i < PARAM_N * PARAM_D ; ++i)
		{
		y[i] += x[i];
		}
	
	freeze_poly(y, PARAM_N * PARAM_D - 1);
	}

/*
	y <- TI * x, where TI is the vertical concatenation of T and I_dk
		y is in R^m in the CRT domain
		T is in R^{2d,dk} in the CRT domain
		x is in R^dk in th CRT domain
*/
void multiply_by_TI(poly_matrix y, poly_matrix T, poly_matrix x)
	{
	// y[:2d] <- T * x
	mul_crt_poly_matrix(y, T, x, 2 * PARAM_D, PARAM_D * PARAM_K, 1, LOG_R);
	
	// y[2d:] <- x
	poly y_2d = poly_matrix_element(y, 1, 2 * PARAM_D, 0);
	memcpy(y_2d, x, PARAM_N * PARAM_D * PARAM_K * sizeof(scalar));
	}

/*
	y <- T * x
		y is in R^2d in the CRT domain
		T is in R^{2d,dk} in the CRT domain
		x is in R^dk in th CRT domain
*/
void multiply_by_T(poly_matrix y, poly_matrix T, poly_matrix x)
	{
	mul_crt_poly_matrix(y, T, x, 2 * PARAM_D, PARAM_D * PARAM_K, 1, LOG_R);
	}
/*
	A <- A_m
		where A_m = A + [0 | H_m * G] describes the SIS instance associated to m
		H_m = I_d * m
		m is a polynomial of degree < n / r
		A is in the CRT domain
		m is in the normal domain
*/
void construct_A_m(poly_matrix A, scalar *m)
	{
	// Compute m * g^T in the CRT domain, where g^T = (1, 2, ..., 2^(k-1)) is a vector of k polynomials
	scalar mgT_coeffs[PARAM_K * PARAM_N];
	poly_matrix mgT = mgT_coeffs;
	
	for(int i = 0 ; i < PARAM_K ; ++i)
		{
		// mgT[i] = 2^i * m in the normal domain
		poly mgT_i = poly_matrix_element(mgT, 1, i, 0);
		
		for(int j = 0 ; j < SMALL_DEGREE ; ++j)
			{
			mgT_i[j] = reduce_double_naive( ((double_scalar) 1 << i) * m[j]);
			}
		
		// Put mgT[i] into the CRT domain (easy since it is of degree < n/r)
		for(int j = 1 ; j < PARAM_R ; ++j)
			{
			poly mgT_ij = &mgT_i[j * SMALL_DEGREE];
			
			memcpy(mgT_ij, mgT_i, SMALL_DEGREE * sizeof(scalar));
			}
		}
	

	
	// Add m * g^T to A
	for(int i = 0 ; i < PARAM_D ; ++i)
		{
		// A[i, d+k*i : d+k*(i+1)] <- A[i, d+k*i : d+k*(i+1)] + m * g^T
		poly A_idpluski = poly_matrix_element(A, PARAM_M - PARAM_D, i, PARAM_D + i * PARAM_K);
		

		
		add_poly(A_idpluski, A_idpluski, mgT, PARAM_K * PARAM_N - 1);
		}
	
	// Reduce A's coefficients mod q
	freeze_poly(A, PARAM_N * PARAM_D * (PARAM_M - PARAM_D) - 1);
	}

void deconstruct_A_m(poly_matrix A, scalar *m)
	{
	// minus_m <- -m
	scalar minus_m[SMALL_DEGREE];
	for(int i = 0 ; i < SMALL_DEGREE ; ++i)
		{
		minus_m[i] = PARAM_Q - m[i];
		}
	
	// Deconstructing A_m is equivalent to constructing A_{-m}
	construct_A_m(A, minus_m);
	}

/*
	Returns the squared euclidian norm of v
		norm_squared = v_0^2 + ... + v_{l*n-1}^2
		where v_i is in [-q/2, q/2]
		v is given with coefficients in [0, q-1]
*/
double_scalar norm_squared(poly_matrix v, int l)
	{
	double_scalar n_s = 0;
	
	for(int i = 0 ; i < l * PARAM_N ; ++i)
		{
		scalar plusorminus_v_i = (v[i] > PARAM_Q / 2) ? (PARAM_Q - v[i]) : (v[i]);
		n_s += ((double_scalar) plusorminus_v_i) * plusorminus_v_i;
		}
	
	return n_s;
	}

/*
	f <- 0
*/
void zero_cplx_poly(cplx_poly f, int deg)
	{
	for(int i = 0 ; i < deg + 1 ; ++i)
		{
		f[i] = 0;
		}
	}

/*
	h <- f + g
*/
void add_cplx_poly(cplx_poly h, cplx_poly f, cplx_poly g, int deg)
	{
	for(int i = 0 ; i < deg + 1 ; ++i)
		{
		h[i] = f[i] + g[i];
		}
	}

/*
	h <- h + f * g
*/
void fma_cplx_poly(cplx_poly h, cplx_poly f, cplx_poly g, int deg)
	{
	for(int i = 0 ; i < deg + 1 ; ++i)
		{
		h[i] += f[i] * g[i];
		}
	}

/*
	h <- f - g
*/
void sub_cplx_poly(cplx_poly h, cplx_poly f, cplx_poly g, int deg)
	{
	for(int i = 0 ; i < deg + 1 ; ++i)
		{
		h[i] = f[i] - g[i];
		}
	}

/*
	h <- h - f * g
*/
void fmsub_cplx_poly(cplx_poly h, cplx_poly f, cplx_poly g, int deg)
	{
	for(int i = 0 ; i < deg + 1 ; ++i)
		{
		h[i] -= f[i] * g[i];
		}
	}

/*
	h <- f^(-1)
*/
void inv_cplx_poly(cplx_poly h, cplx_poly f, int deg)
	{
	for(int i = 0 ; i < deg + 1 ; ++i)
		{
		h[i] = 1 / f[i];
		}
	}

/*
	h <- f * g
		inputs and output are in the FFT domain
*/
void mul_cplx_poly(cplx_poly h, cplx_poly f, cplx_poly g, int deg)
	{
	for(int i = 0 ; i < deg + 1 ; ++i)
		{
		h[i] = f[i] * g[i];
		}
	}

/*
	h <- f / g
		in the complex CRT domain
*/
void div_cplx_poly(cplx_poly h, cplx_poly f, cplx_poly g, int deg)
	{
	for(int i = 0 ; i < deg + 1 ; ++i)
		{
		h[i] = f[i] / g[i];
		}
	}

/*
	h <- f * g^T
		inputs and output are in the FFT domain
*/
void mul_cplx_poly_by_transpose(cplx_poly h, cplx_poly f, cplx_poly g, int deg)
	{
	for(int i = 0 ; i < deg + 1 ; ++i)
		{
		h[i] = f[i] * conj(g[i]);
		}
	}


/*
	h <- f * f^T
		inputs and output are in the FFT domain
*/
void mul_cplx_poly_by_own_transpose(cplx_poly h, cplx_poly f, int deg)
	{
	for(int i = 0 ; i < deg + 1 ; ++i)
		{
		h[i] = f[i] * conj(f[i]);
		}
	}

/*
	h <- h + f * g^T
		in the FFT domain
*/
void cplx_fma_transpose(cplx_poly h, cplx_poly f, cplx_poly g, int deg)
	{
	for(int i = 0 ; i < deg + 1 ; ++i)
		{
		h[i] += f[i] * conj(g[i]);
		}
	}

/*
	h <- h - f * g^T
		in the complex CRT domain
*/
void fmsub_transpose_cplx_poly(cplx_poly h, cplx_poly f, cplx_poly g, int deg)
	{
	for(int i = 0 ; i < deg + 1 ; ++i)
		{
		h[i] -= f[i] * conj(g[i]);
		}
	}

/*
	f <- f^T
		in the complex CRT domain
*/
void transpose_cplx_poly(cplx_poly f, int deg)
	{
	for(int i = 0 ; i < deg + 1 ; ++i)
		{
		f[i] = conj(f[i]);
		}
	}

/*
	Compute f0, f1 such that f(x) = f0(x^2) + x*f1(x^2)
		f is of degree < deg, f0 and f1 of degree < deg/2
		in place
		in the complex CRT domain
*/
void stride(cplx_poly f, int depth)
	{
	int deg = PARAM_N >> depth;
	cplx *f0 = f, f1[deg/2];
	
	for(int i = 0 ; i < deg / 2 ; ++i)
		{
		// w_i is such that f[2*i] is f(w_i) and f[2*i+1] is f(-w_i)
		cplx w_i = - cplx_cyclotomic_factorisation_tree(LOG_N - depth, 2*i);
		f1[i] = conj(w_i) * (f[2*i] - f[2*i+1]) / 2;
		f0[i] = (f[2*i] + f[2*i+1]) / 2;
		}

	
	memcpy(&f[deg/2], f1, deg / 2 * sizeof(cplx));
	}

/*
	Compute f such that f(x) = f0(x^2) + x*f1(x^2)
		f is of degree < deg, f0 and f1 of degree < deg/2
		in place
		in the complex CRT domain
*/
void inverse_stride(cplx_poly f, int depth)
	{
	int deg = PARAM_N >> depth;
	cplx f0[deg/2], *f1 = &f[deg/2];
	
	memcpy(f0, f, deg / 2 * sizeof(cplx));
	
	for(int i = 0 ; i < deg / 2 ; ++i)
		{
		// w_i is such that f0[i] is f0(w_i^2) and f1[i] is f1(-w_i^2)
		cplx w_i = - cplx_cyclotomic_factorisation_tree(LOG_N - depth, 2*i);
		f[2*i] = f0[i] + w_i * f1[i];
		f[2*i+1] = f0[i] - w_i * f1[i];
		}
	}

/*
	permute p = (p_0, ... p_(s-1)) to (p_0, p_2, ..., p_(s-2), p_1, p_3, ... p_(s-1))
*/
void scalar_stride(signed_scalar *p, int size)
	{
	signed_scalar aux[size / 2];
	
	memcpy(aux, p, size / 2 * sizeof(signed_scalar));
	for(int i = 0 ; i < size / 2 ; ++i)
		{
		p[2 * i] = aux[i];
		p[2*i + 1] = p[size / 2 + i];
		}
	}

/*
	Reduce f mod X^(deg/2) + c and mod X^(deg/2) - c
		in place
*/
void reduce_mod_sparse_cplx_poly(cplx_poly f, int deg, cplx c)
	{
	for(int i = 0 ; i < deg / 2 ; ++i)
		{
		cplx f_i = f[i];
		f[i] = f_i - c * f[deg/2 + i];
		f[deg/2+i] = f_i + c * f[deg/2 + i];
		}
	}

/*
	Compute the complex CRT respresentation of f, i.e. its evaluation at all the complex primitive (2n)-th roots of unity
		in place
*/
void cplx_crt_representation(cplx_poly f)
	{
	
	for(int depth = 0 ; (1 << depth) < PARAM_N ; depth++)
		{
		int deg = PARAM_N >> depth;
		for(int i = 0 ; i < (1 << depth) ; ++i)
			{
			// Reduce f_i mod X^(deg/2) +- w_i
			cplx_poly f_i = &f[i*deg];
			cplx w_i = cplx_cyclotomic_factorisation_tree(depth+1, 2*i);
			reduce_mod_sparse_cplx_poly(f_i, deg, w_i);
			

			}

		}
	}

void matrix_cplx_crt_representation(cplx_poly_matrix M, int l1, int l2)
	{
	for(int i = 0 ; i < l1 * l2 ; ++i)
		{
		cplx_poly M_i = poly_matrix_element(M, 1, i, 0);
		
		cplx_crt_representation(M_i);
		}
	}
/*
	Construct the Schur complement of (zeta^2 - alpha^2)I in \Sigma_p, and a new center c
		sch_comp <- zeta^2 * I - (alpha^(-2) - zeta^(-2))^(-1) * T * T^T
		c <- - alpha^2 / (zeta^2 - alpha^2) * T * p
		T is in the complex CRT domain
		the result is stored in a lower triangular matrix in the complex CRT domain
*/
void construct_T_schur_complement(cplx_poly_matrix sch_comp, cplx_poly_matrix T)
	{
	// sch_comp <- T * T^T
	zero_cplx_poly(sch_comp, PARAM_N * PARAM_D * (2 * PARAM_D + 1) - 1);
	for(int i = 0 ; i < 2 * PARAM_D ; ++i)
		{
		for(int j = 0 ; j < i + 1 ; ++j)
			{
			for(int k = 0 ; k < PARAM_D * PARAM_K ; k++)
				{
				cplx_poly T_ik = poly_matrix_element(T, PARAM_D * PARAM_K, i, k);
				cplx_poly T_jk = poly_matrix_element(T, PARAM_D * PARAM_K, j, k);
				cplx_poly sch_comp_ij = triangular_poly_matrix_element(sch_comp, i, j);
				
				cplx_fma_transpose(sch_comp_ij, T_ik, T_jk, PARAM_N-1);
				

				}
			}
		}
	

	
	// sch_comp = zeta^2 * I - z * sch_comp
	real z = 1 / (1/(PARAM_ALPHA * PARAM_ALPHA) - 1/(PARAM_ZETA * PARAM_ZETA));
	for(int i = 0 ; i < PARAM_N * PARAM_D * (2 * PARAM_D + 1) ; ++i)
		{
		sch_comp[i] = - z * sch_comp[i];
		}
	for(int i = 0 ; i < 2 * PARAM_D ; ++i)
		{
		cplx_poly sch_comp_ii = triangular_poly_matrix_element(sch_comp, i, i);
		for(int j = 0 ; j < PARAM_N ; ++j)
			{
			sch_comp_ii[j] += PARAM_ZETA * PARAM_ZETA;
			}
		}
	}

/*
	Construct the Schur complement of d in M, and a new center c0
		M <- M - d^(-1) * B * B^T
		c0 <- c0 + d^(-1) * B * (x1 - c1), where c = (c0, c1)
		in place, M is replaced by the Schur complement
		M is lower triangular and in the complex CRT domain
		d is in the complex CRT domain
		M's blocks are A, B, B^T, d where A is l*l, B is l*1, B^T is 1*l, d is 1*1
*/
void construct_schur_complement_and_center(cplx_poly_matrix M, cplx_poly_matrix c, cplx_poly d, cplx_poly x1, int l)
	{
	cplx prod_coeffs[PARAM_N], d_inv_coeffs[PARAM_N];
	cplx_poly prod = prod_coeffs, d_inv = d_inv_coeffs;
	
	// d_inv <- d^(-1)
	inv_cplx_poly(d_inv, d, PARAM_N - 1);
	

	
	// c0 <- c0 + d^(-1) * B * (x1 - c1)
	cplx_poly_matrix c0 = c;
	cplx_poly c1 = poly_matrix_element(c, 1, l, 0);
	
	sub_cplx_poly(x1, x1, c1, PARAM_N - 1);
	for(int i = 0 ; i < l ; ++i)
		{
		cplx_poly B_i = triangular_poly_matrix_element(M, l + 1, i);
		cplx_poly x1_i = poly_matrix_element(x1, 1, i, 0);
		
		mul_cplx_poly(c1, B_i, x1_i, PARAM_N - 1);
		fma_cplx_poly(c0, c1, d_inv, PARAM_N - 1);
		}
	
	// A <- A - B * d^-1 * B^T
	for(int i = 0 ; i < l ; ++i)
		{
		for(int j = 0 ; j < i + 1 ; ++j)
			{
			cplx_poly B_i = triangular_poly_matrix_element(M, l + 1, i);
			cplx_poly B_j = triangular_poly_matrix_element(M, l + 1, j);
			cplx_poly A_ij = triangular_poly_matrix_element(M, i, j);
			

			mul_cplx_poly_by_transpose(prod, B_j, B_i, PARAM_N - 1);

			fmsub_cplx_poly(A_ij, d_inv, prod, PARAM_N - 1);

			}
		}
	}

/*
	Construct the Schur complement of d in M
		M <- M - d^(-1) * B * B^T
		in place, M is replaced by the Schur complement
		M is lower triangular and in the complex CRT domain
		M's blocks are A, B, B^T, d where A is l*l, B is l*1, B^T is 1*l, d is 1*1
*/
void construct_schur_complement(cplx_poly_matrix M, int l)
	{
	cplx prod_coeffs[PARAM_N], d_inv_coeffs[PARAM_N];
	cplx_poly prod = prod_coeffs, d_inv = d_inv_coeffs;
	
	// d_inv <- d^(-1)
	cplx_poly d = triangular_poly_matrix_element(M, l, l);
	inv_cplx_poly(d_inv, d, PARAM_N - 1);

	
	// A <- A - B * d^-1 * B^T
	for(int i = 0 ; i < l ; ++i)
		{
		cplx_poly B_i = triangular_poly_matrix_element(M, l, i);
		
		for(int j = 0 ; j < i + 1 ; ++j)
		//for(int j = 0 ; j < i ; ++j)
			{
			cplx_poly B_j = triangular_poly_matrix_element(M, l, j);
			cplx_poly A_ij = triangular_poly_matrix_element(M, i, j);

			
			mul_cplx_poly_by_transpose(prod, B_j, B_i, PARAM_N - 1);
			fmsub_cplx_poly(A_ij, d_inv, prod, PARAM_N - 1);

			}

		}
	}

/*
	Construct a new centre c0 which depends on the former centre c and the newly sampled q1
		c0 = c0 + B * d^(-1) * (q1 - c1), where c = (c0, c1), c0 is of size l and c1 of size 1
		c1 is overwritten
		M is lower triangular, its blocks are A, B, B^T, d where A is l*l, B is l*1, B^T is 1*l, d is 1*1
		M, c and q1 are in the complex CRT domain
*/
void construct_new_center(cplx_poly_matrix c, cplx_poly_matrix M, cplx_poly q1, int l)
	{
	// c1 <- d^(-1) * (q1 - c1)
	cplx_poly c1 = poly_matrix_element(c, 1, l, 0);
	cplx_poly d = triangular_poly_matrix_element(M, l, l);
	
	sub_cplx_poly(c1, q1, c1, PARAM_N - 1);
	div_cplx_poly(c1, c1, d, PARAM_N - 1);
	

	
	// c0 <- c0 + d^(-1) * B * (q1 - c1)
	for(int i = 0 ; i < l ; ++i)
		{
		// c0[i] <- c0[i] + d^(-1) * B[i] * (q1 - c1)
		cplx_poly B_i = triangular_poly_matrix_element(M, l, i);
		cplx_poly c0_i = poly_matrix_element(c, 1, i, 0);
		

		
		// transpose B_i since we actually store B^T and not B
		cplx_fma_transpose(c0_i, c1, B_i, PARAM_N - 1);
		}
	}

/*
	Construct the first center that is not 0 in the perturbation sampling, which depends on the already sampled p
		c <- - alpha^2 / (zeta^2 - alpha^2) * T * p
		c is returned in the complex CRT domain
		T is given in the complex CRT domain
		p is given in the complex CRT domain
*/
void construct_first_center(cplx_poly_matrix c, cplx_poly_matrix T, cplx_poly_matrix p)
	{
	// c <- 0
	zero_cplx_poly(c, PARAM_N * 2 * PARAM_D - 1);
	
	// c <- T * p
	for(int i = 0 ; i < 2 * PARAM_D ; ++i)
		{
		// c_i <- T[i,:] * p
		cplx_poly c_i = poly_matrix_element(c, 1, i, 0);
		
		for(int j = 0 ; j < PARAM_D * PARAM_K ; ++j)
			{
			cplx_poly T_ij = poly_matrix_element(T, PARAM_D * PARAM_K, i, j);
			cplx_poly p_j = poly_matrix_element(p, 1, j, 0);
			
			fma_cplx_poly(c_i, T_ij, p_j, PARAM_N - 1);

			}
		}
	
	// c <- - alpha^2 / (zeta^2 - alpha^2) * c
	real z = - PARAM_ALPHA * PARAM_ALPHA / ((PARAM_ZETA * PARAM_ZETA) - (PARAM_ALPHA * PARAM_ALPHA));
	for(int i = 0 ; i < PARAM_N * 2 * PARAM_D ; ++i)
		{
		c[i] = z * c[i];
		}

	}

/*
	Construct the Schur complements that are necessary to the Gaussian preimage sampling operation + cplx_T
		cplx_T is returned in the complex CRT domain
		sch_comp is returned in the complex CRT domain
		T is given in the normal domain
*/
void construct_complex_private_key(cplx_poly_matrix cplx_T, cplx_poly_matrix sch_comp, poly_matrix T)
	{
	// Put T in the complex CRT domain
	for(int i = 0 ; i < PARAM_N * 2 * PARAM_D * PARAM_D * PARAM_K ; ++i)
		{
		cplx_T[i] = (signed_scalar) T[i];
		}
	

	
	matrix_cplx_crt_representation(cplx_T, 2 * PARAM_D, PARAM_D * PARAM_K);
	
	// Construct the Schur complement of (zeta^2 - alpha^2)I in \Sigma_p
	construct_T_schur_complement(sch_comp, cplx_T);
	

	// Construct all the 2d - 2 remaining Schur complements
	for(int i = 2 * PARAM_D - 1 ; i > 1 ; --i)
		{
		construct_schur_complement(sch_comp, i);
		

		}
	}
/*
	Initialise the trees used during the CRT multiplication :
		- cyclotomic_factorisation_tree which contains the factors of x^n + 1
		- bezout_coefficients_tree which contains the Bézout coefficients used to invert the CRT map
*/
void init_crt_trees(void)
	{
	int index = 0;
	
	for(int i=0 ; i < LOG_R + 1 ; ++i)
		{
		cyclotomic_factorisation_tree[i] = &cyclotomic_factorisation_array[index];
		bezout_coefficients_tree[i] = &bezout_coefficients_array[index];
		index += (1 << i);
		}
	}

cplx cplx_roots_of_unity[2 * PARAM_N - 1];

/*
	
*/
void init_cplx_roots_of_unity(void)
	{
	cplx_cyclotomic_factorisation_tree(0, 0) = 1;
	cplx_cyclotomic_factorisation_tree(1, 0) = I;
	cplx_cyclotomic_factorisation_tree(1, 1) = -I;
	
	for(int i = 2 ; (1 << i) <= PARAM_N ; ++i)
		{
		for(int j = 0 ; j < (1 << (i-1)) ; ++j)
			{
			cplx w = csqrt(-cplx_cyclotomic_factorisation_tree(i-1, j));
			cplx_cyclotomic_factorisation_tree(i, 2*j) = w;
			cplx_cyclotomic_factorisation_tree(i, 2*j+1) = -w;
			}
		}
	

	}
