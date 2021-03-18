#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
#include <string.h>
#include <stdbool.h>

#include "arithmetic.h"
#include "random.h"

// Code from Dilithium

#include "zetas.c"


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
void ntt(uint32_t p[PARAM_N])
	{
	unsigned int len, start, j, k;
	scalar zeta, t;
	
	#ifdef ALWAYS_CSUBQ
	scalar t1;
	#endif

	k = 1;
	for(len = PARAM_N / 2; len > 0; len >>= 1)
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
	}

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
void invntt_frominvmont(uint32_t p[PARAM_N])
	{
	unsigned int start, len, j, k;
	scalar t, zeta;
	const scalar f = (((double_scalar)MONT*MONT % PARAM_Q) * (PARAM_Q-1) % PARAM_Q) * ((PARAM_Q-1) >> LOG_N) % PARAM_Q;
	
	#ifdef ALWAYS_CSUBQ
	#endif

	k = 0;
	for(len = 1; len < PARAM_N; len <<= 1)
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
				p[j + len] = csubq(montgomery_reduce((double_scalar)zeta * p[j + len]));
				#elif defined NTT_ALWAYS_REDUCE
				p[j] = barrett_reduce(t + p[j + len]);
				p[j + len] = t + 2 * PARAM_Q - p[j + len];
				p[j + len] = montgomery_reduce((double_scalar) zeta * p[j + len]);
				#else
				p[j] = t + p[j + len];
				p[j + len] = t + PARAM_N * PARAM_Q - p[j + len];
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
	}

void invntt_without_mont(scalar p[PARAM_N])
	{
	multiply_by_one(p, PARAM_N - 1);
	
	invntt_frominvmont(p);
	}

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

/*************************************************
* Name:        csubq
*
* Description: Subtract Q if input coefficient is bigger than Q.
*
* Arguments:   - uint32_t: finite field element a
*
* Returns r.
**************************************************/
scalar csubq(scalar a)
	{
	a -= PARAM_Q;
	a += ((signed_scalar)a >> 31) & PARAM_Q;
	return a;
	}

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
uint32_t montgomery_reduce(uint64_t a)
	{
	uint64_t t;

	t = a * QINV;
	t &= (1ULL << 32) - 1;
	t *= PARAM_Q;
	t = a + t;
	t >>= 32;
	return t;
	}

/*************************************************
* Name:        poly_pointwise_invmontgomery
*
* Description: Pointwise multiplication of polynomials in NTT domain
*              representation and multiplication of resulting polynomial
*              with 2^{-32}. Output coefficients are less than 2*Q if input
*              coefficient are less than 22*Q.
*
* Arguments:   - poly *c: pointer to output polynomial
*              - const poly *a: pointer to first input polynomial
*              - const poly *b: pointer to second input polynomial
**************************************************/
void poly_pointwise_invmontgomery(poly c, const poly a, const poly b)
{
	unsigned int i;

	for(i = 0; i < PARAM_N; ++i)
		{
		c[i] = montgomery_reduce((uint64_t)a[i] * b[i]);
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
	}

void matrix_invntt_without_mont(poly_matrix A, int l1, int l2)
	{
	for(int i = 0 ; i < l1 * l2 ; ++i)
		{
		poly A_ij = poly_matrix_element(A, 1, i, 0);
		
		invntt_without_mont(A_ij);
		}
	}

void multiply_by_one(poly f, int deg)
	{
	for(int i = 0 ; i < deg + 1 ; ++i)
		{
		f[i] = montgomery_reduce(f[i]);
		}
	}

/*
	A <- B * C
*/
void mul_crt_poly_matrix(poly_matrix C, poly_matrix A, poly_matrix B, int l1, int l2, int l3)
	{
	scalar prod_coeffs[PARAM_N];
	poly prod = prod_coeffs;
	
	
	memset(C, 0, PARAM_N * l1 * l3 * sizeof(scalar));
	
	for(int i = 0 ; i < l1 ; ++i)
		{
		for(int k = 0 ; k < l2 ; ++k)
			{
			poly A_ik = poly_matrix_element(A, l2, i, k);
			
			for(int j = 0 ; j < l3 ; ++j)
				{
				
				poly B_kj = poly_matrix_element(B, l3, k, j);
				poly C_ij = poly_matrix_element(C, l3, i, j);
				
				poly_pointwise_invmontgomery(prod, A_ik, B_kj);
				add_poly(C_ij, C_ij, prod, PARAM_N - 1);
				freeze_poly(C_ij, PARAM_N - 1);
				}
			}
		}
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
	//printf("subtraction won't work\n");
	for(int i=0 ; i < deg + 1 ; ++i)
		{
		h[i] = reduce_naive(PARAM_Q + f[i] - g[i]);
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
	v <- v + g where g is the ring gadget vector
		v is in R^k in the NTT domain
*/
void add_ring_gadget_vector(poly_matrix v)
	{
	for(int i = 0 ; i < PARAM_K ; ++i)
		{
		poly v_i = poly_matrix_element(v, 1, i, 0);
		
		for(int j = 0 ; j < PARAM_N ; ++j)
			{
			v_i[j] = barrett_reduce((double_scalar) v_i[j] + (1 << i));
			}
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
	
	mul_crt_poly_matrix(y, A, x_d, PARAM_D, PARAM_M - PARAM_D, 1);
	
	// y <- y + x[:d]
	divide_by_2pow32(x, PARAM_N * PARAM_D - 1);
	freeze_poly(x, PARAM_N * PARAM_M - 1);
	
	for(int i = 0 ; i < PARAM_N * PARAM_D ; ++i)
		{
		y[i] += x[i];
		}
	
	freeze_poly(y, PARAM_N * PARAM_D - 1);
	
	multiply_by_2pow32(x, PARAM_N * PARAM_D - 1);
	freeze_poly(x, PARAM_N * PARAM_M - 1);
	
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
	mul_crt_poly_matrix(y, T, x, 2 * PARAM_D, PARAM_D * PARAM_K, 1);
	
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
	mul_crt_poly_matrix(y, T, x, 2 * PARAM_D, PARAM_D * PARAM_K, 1);
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
