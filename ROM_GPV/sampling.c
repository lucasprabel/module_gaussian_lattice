#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
#include <math.h>
#include <assert.h>
#include <string.h>

#include "common.h"
#include "random.h"
#include "sampling.h"
#include "arithmetic.h"

#include "cpucycles.h"

/*
	Generates a matrix A in R_q^{d,m} along with its trapdoor T in R^{2d,dk}
		A's first d columns are implicit, since they are the identity I_d
		A and T are both returned in the CRT domain
*/
void TrapGen(poly_matrix A, poly_matrix T)
	{
	scalar A_hat_coeffs[PARAM_D * PARAM_D * PARAM_N], AprimeT_coeffs[PARAM_D * PARAM_D * PARAM_K * PARAM_N];
	poly_matrix A_hat = A_hat_coeffs, AprimeT = AprimeT_coeffs;
	
	// A_hat <- U(R_q^{d,d}) is considered to be in the CRT domain
	random_poly(A_hat, PARAM_N * PARAM_D * PARAM_D);

		
	// T <- D_{R^{2d,dk},sigma}
	SampleR_matrix_centered((signed_poly_matrix) T, 2*PARAM_D, PARAM_D * PARAM_K, PARAM_SIGMA);
	
	matrix_ntt(T, 2*PARAM_D, PARAM_D * PARAM_K);
	
	// AprimeT = A_hat * T2 + T1, where T1 and T2 are the upper and lower half of T
	poly_matrix T1 = T, T2 = poly_matrix_element(T, PARAM_D * PARAM_K, PARAM_D, 0);
	mul_crt_poly_matrix(AprimeT, A_hat, T2, PARAM_D, PARAM_D, PARAM_D * PARAM_K);
	add_to_poly_matrix(AprimeT, T1, PARAM_D, PARAM_D * PARAM_K);
	

	
	// A = (A_hat | -A'T) ( = (I | A_hat | -A'T) implicitly)
	for(int i = 0 ; i < PARAM_D ; ++i)
		{
		poly_matrix A_i0 = poly_matrix_element(A, PARAM_M - PARAM_D, i, 0);
		poly_matrix A_hat_i = poly_matrix_element(A_hat, PARAM_D, i, 0);
		memcpy(A_i0, A_hat_i, PARAM_D * PARAM_N * sizeof(scalar));
		}
	for(int i = 0 ; i < PARAM_D ; ++i)
		{
		poly_matrix A_i1 = poly_matrix_element(A, PARAM_M - PARAM_D, i, PARAM_D);
		poly_matrix AprimeT_i = poly_matrix_element(AprimeT, PARAM_D * PARAM_K, i, 0);
		for(int j = 0 ; j < PARAM_D * PARAM_K * PARAM_N ; ++j)
			{
			A_i1[j] = 2*PARAM_Q - AprimeT_i[j];
			}
		}
	}


//==============================================================================
// Samples from distribution D_{c,sigma}, ie                                              
// Samples an element in Z with probability proportionnal to e^{-(c-x)^2/2*(sigma^2)}    
//==============================================================================
signed int SampleZ(RR_t c, RR_t sigma)
	{
	//printf("SampleZ with sigma = %Lf, c = %Lf\n", sigma, c);
	
	return algorithmF(c, sigma);
	}

/*
	Samples from D_{R,sigma} where R = Z[X] / <X^n+1> is isomorphic to Z^n
*/
void SampleR_centered(signed_poly f, RR_t sigma)
	{
	for(int i=0 ; i < PARAM_N ; ++i)
		{
		f[i] = SampleZ(0, sigma);
		}
	}

/*
	Samples from D_{R^{l1,l2},sigma}
		coefficients are centered on 0
*/
void SampleR_matrix_centered(signed_poly_matrix A, int l1, int l2, RR_t sigma)
	{
	for(int i = 0 ; i < l1 * l2 * PARAM_N ; ++i)
		{
		A[i] = SampleZ(0, sigma);
		}
	}

real d_coeffs[PARAM_K];
real l_coeffs[PARAM_K];
real h_coeffs[PARAM_K];

/*
	Compute coefficients used in D-sampling and G-perturbation-sampling :
		- d_coeffs : the coefficients of D's last column
		- l_coeffs and h_coeffs : the coefficients defining L
*/
void init_D_lattice_coeffs(void)
	{
	// d_i = (d_{i-1} + q_i) / 2, with d_{-1} = 0
	d_coeffs[0] = Q_BIT(0) / 2.0;
	for(int i = 1 ; i < PARAM_K ; ++i)
		{
		d_coeffs[i] = (d_coeffs[i-1] + Q_BIT(i)) / 2;
		}
	
	// l_0^2 = 2(1 + 1/k) + 1, and l_i^2 = 2(1 + 1/(k-i))
	l_coeffs[0] = sqrt(2*(1 + 1.0 / PARAM_K) + 1);
	for(int i = 0 ; i < PARAM_K ; ++i)
		{
		l_coeffs[i] = sqrt(2*(1 + 1.0 / (PARAM_K - i)));
		}
	
	// h_{i+1}^2 = 2(1 - 1/(k-i)), there is no h_0 so we set it to zero
	h_coeffs[0] = 0;
	for(int i = 0 ; i < PARAM_K - 1 ; ++i)
		{
		h_coeffs[i+1] = sqrt(2*(1 - 1.0 / (PARAM_K - i)));
		}
	}

extern unsigned long long timing_sampleZ_G;
extern unsigned long long timing_sampleZ_P;
/*
	Sample z from the D-lattice (included in Z^k) with center c and parameter sigma
*/
void sample_D(signed_scalar *z, real *c, real sigma)
	{
	unsigned long long begin_timing = 0;
	unsigned long long end_timing   = 0;

	real c_d = - c[PARAM_K - 1] / d_coeffs[PARAM_K - 1];
	real c_d_floor = floor(c_d);
	real c_d_frac = c_d - c_d_floor;
	
	begin_timing = cpucycles_start();
	z[PARAM_K - 1] = c_d_floor + SampleZ(c_d_frac, sigma / d_coeffs[PARAM_K - 1]);
	end_timing = cpucycles_stop();
	timing_sampleZ_G += (end_timing - begin_timing);
	
	
	for(int i = 0 ; i < PARAM_K - 1 ; ++i)
		{
		real c_i = z[PARAM_K - 1] * d_coeffs[i] - c[i];
		real c_i_floor = floor(c_i);
		real c_i_frac = c_i - c_i_floor;

		begin_timing = cpucycles_start();
		z[i] = c_i_floor + SampleZ(c_i_frac, sigma);
		end_timing = cpucycles_stop();
		timing_sampleZ_G += (end_timing - begin_timing);
		
		}
	}

/*
	Samples the perturbation vector p so that the distribution of Sample_G is spherical
*/
void sample_G_perturb(real *p, real sigma)
	{
	unsigned long long begin_timing = 0;
	unsigned long long end_timing   = 0;

	real beta = 0, z[PARAM_K+1];
	z[PARAM_K] = 0;
	
	for(int i = 0 ; i < PARAM_K ; ++i)
		{
		real c_i = beta / l_coeffs[i], c_i_floor = floor(c_i), c_i_frac = c_i - c_i_floor;
		real sigma_i = sigma / l_coeffs[i];
		begin_timing = cpucycles_start();
		z[i] = c_i_floor + SampleZ(c_i_frac, sigma_i);
		end_timing = cpucycles_stop();
		timing_sampleZ_G += (end_timing - begin_timing);
		
		beta = -z[i] * h_coeffs[i];
		}
	
	p[0] = 5*z[0] + 2*z[1];
	for(int i = 1 ; i < PARAM_K ; ++i)
		{
		p[i] = 2 * (z[i-1] + 2*z[i] + z[i+1]);
		}
	}

/*
	Sample t from the scalar G-lattice (included in Z^k) with parameter alpha, such that <g,t> = u mod q
*/
void scalar_sample_G(signed_scalar *t, scalar u)
	{
	real sigma = PARAM_ALPHA / 3, c[PARAM_K], p[PARAM_K];
	signed_scalar z[PARAM_K];
	sample_G_perturb(p, sigma);
	
	c[0] = ((real) GET_BIT(u,0) - p[0]) / 2;
	for(int i = 1 ; i < PARAM_K ; ++i)
		{
		c[i] = (c[i-1] + GET_BIT(u,i) - p[i]) / 2;
		}
	
	sample_D(z, c, sigma);
	
	t[0] = 2*z[0] + Q_BIT(0)*z[PARAM_K - 1] + GET_BIT(u,0);
	for(int i = 1 ; i < PARAM_K - 1 ; ++i)
		{
		t[i] = 2*z[i] - z[i-1] + Q_BIT(i)*z[PARAM_K - 1] + GET_BIT(u,i);
		}
	t[PARAM_K - 1] = Q_BIT(PARAM_K - 1)*z[PARAM_K - 1] - z[PARAM_K - 2] + GET_BIT(u,PARAM_K-1);
	}

/*
	Sample t form the ring G-lattice (included in R^k) with parameter alpha, such that <g,t> = u mod q
*/
void ring_sample_G(signed_poly_matrix t, poly u)
	{
	signed_scalar t_T[PARAM_N * PARAM_K];
	
	// sample n times from the scalar G-lattice
	for(int i = 0 ; i < PARAM_N ; ++i)
		{
		signed_poly t_i = &t_T[i*PARAM_K];
		scalar_sample_G(t_i, u[i]);
		
		}
	
	
	// permute the coefficients of t
	// as if t was a (n, k) matrix of scalars and we transposed it
	transpose_scalar_matrix((scalar *) t, (scalar *) t_T, PARAM_N, PARAM_K);
	}

/*
	Transpose A of size (l0, l1) into A_T of size (l1, l0)
*/
void transpose_scalar_matrix(scalar *A_T, scalar *A, int l0, int l1)
	{
	#define mat_A(i,j) A[i*l1 + j]
	#define mat_A_T(i,j) A_T[i*l0 + j]
	for(int i = 0 ; i < l0 ; ++i)
		{
		for(int j = 0 ; j < l1 ; ++j)
			{
			mat_A_T(j,i) = mat_A(i,j);
			}
		}
	#undef mat_A
	#undef mat_A_T
	}

/*
	Sample t from the module G-lattice (included in R^dk) with parameter alpha, such that Gt = u mod q
*/
void module_sample_G(signed_poly_matrix t, poly_matrix u)
	{
	// sample d times from the ring G-lattice
	for(int i = 0 ; i < PARAM_D ; ++i)
		{
		signed_poly_matrix t_i = poly_matrix_element(t, PARAM_D, 0, i*PARAM_K);
		poly u_i = poly_matrix_element(u, 1, i, 0);
		
		ring_sample_G(t_i, u_i);
		}
	}

/*
	Samples q in Z^(2*deg) with center c and a covariance defined by a, b, d (each of degree < deg)
*/
void sample_2z(signed_scalar *q, cplx_poly cplx_q, cplx_poly a, cplx_poly b, cplx_poly d, cplx *c, int depth)
	{
	int deg = PARAM_N >> depth;
	
	
	// Split q = (q0, q1) and c = (c0, c1) in half, and copy c1
	signed_scalar *q0 = q, *q1 = &q[deg];
	cplx *c0 = c, *c1 = &c[deg], *cplx_q0 = cplx_q, *cplx_q1 = &cplx_q[deg];
	
	// Copy c1 since it will be modified by sample_fz and we need it afterwards
	cplx c1_bis[deg];
	
	memcpy(c1_bis, c1, deg * sizeof(cplx));

	
	// b_times_d_inv <- b * d^(-1)
	cplx b_times_d_inv_coeffs[deg];
	cplx_poly b_times_d_inv = b_times_d_inv_coeffs;
	
	div_cplx_poly(b_times_d_inv, b, d, deg-1);
	
	// Sample q1 with covariance d and center c1
	sample_fz(q1, cplx_q1, d, c1_bis, depth);
	
	
	// Compute the new covariance
	// a <- a - b^T * d^(-1) * b
	fmsub_transpose_cplx_poly(a, b_times_d_inv, b, deg-1);
	
	// Compute the new center
	// c0 <- c0 + b * d^(-1) * (q1 - c1)
	sub_cplx_poly(c1, cplx_q1, c1, deg - 1); // we don't need c1 anymore so we can overwrite it with (q1 - c1)
	fma_cplx_poly(c0, b_times_d_inv, c1, deg - 1);
	
	
	// Sample q0 with the new covariance and the new center
	sample_fz(q0, cplx_q0, a, c0, depth);
	
	// Update cplx_q by merging q0 and q1 the FFT way
	inverse_stride(cplx_q, depth - 1);
	}

/*
	Samples p in Z^deg with center c and covariance f, where f is of degree < deg
*/
void sample_fz(signed_scalar *p, cplx_poly cplx_p, cplx_poly f, cplx *c, int depth)
	{
	unsigned long long begin_timing = 0;
	unsigned long long end_timing   = 0;

	int deg = PARAM_N >> depth;
	
	
	if(deg == 1)
		{
		// f is a real polynomial of degree 0, so Re(f[0]) = f[0] = f(-1) = f
		begin_timing = cpucycles_start();
		p[0] = SampleZ(*c, sqrt(creal(f[0])));
		end_timing = cpucycles_stop();
		timing_sampleZ_P += (end_timing - begin_timing);
		
		cplx_p[0] = p[0];
		
		return;
		}
	
	stride(f, depth);
	stride(c, depth);
	cplx_poly f0 = f, f1 = &f[deg/2];

	// copy f0 before calling sample_2z
	cplx f0_bis_coeffs[deg/2];
	cplx_poly f0_bis = f0_bis_coeffs;
	
	memcpy(f0_bis, f0, deg/2 * sizeof(cplx));
	
	sample_2z(p, cplx_p, f0, f1, f0_bis, c, depth + 1);
	scalar_stride(p, deg);
	}

/*
	Sample a perturbation p in R^m with the complementary covariance \Sigma_p defined by T and represented by sch_comp
		p is returned in the normal domain
		T is given in the complex CRT domain
		sch_comp is given in the complex CRT domain
*/
void sample_perturb(signed_poly_matrix p, cplx_poly_matrix T, cplx_poly_matrix sch_comp)
	{
	unsigned long long begin_timing = 0;
	unsigned long long end_timing   = 0;
	
	cplx c_coeffs[PARAM_N * 2 * PARAM_D], cplx_p_coeffs[PARAM_N * PARAM_D * PARAM_K];
	cplx_poly_matrix c = c_coeffs, cplx_p = cplx_p_coeffs;
	
	// First sample dk independant centred polynomials with covariance (zeta^2 - alpha^2)
	signed_poly_matrix p_2d = poly_matrix_element(p, 1, 2 * PARAM_D, 0);
	
	real param = sqrt((PARAM_ZETA * PARAM_ZETA) - (PARAM_ALPHA * PARAM_ALPHA));
	for(int i = 0 ; i < PARAM_N * PARAM_D * PARAM_K ; ++i)
		{
		begin_timing = cpucycles_start();
		p_2d[i] = SampleZ(0, param); 
		end_timing = cpucycles_stop();
		timing_sampleZ_P += (end_timing - begin_timing);
		// add q so that the coefficients are positive
		}
	
	
	// Compute the complex CRT transform of the sampled polynomials
	for(int i = 0 ; i < PARAM_N * PARAM_D * PARAM_K ; ++i)
		{
		cplx_p[i] = p_2d[i];
		}
	
	
	matrix_cplx_crt_representation(cplx_p, PARAM_D * PARAM_K, 1);
	
	// Construct the new center (depends on the dk polynomials sampled before)
	construct_first_center(c, T, cplx_p);
	
	
	// Sample 2d - 2 polynomials iteratively
	for(int i = 2 * PARAM_D - 1 ; i > 1 ; --i)
		{
		
		// Sample p[i] with covariance sch_comp[i,i] and center c[i]
		// (copying sch_comp[i,i] and c[i] since they're going to be modified)
		signed_scalar *p_i = poly_matrix_element(p, 1, i, 0);
		cplx_poly sch_comp_ii = triangular_poly_matrix_element(sch_comp, i, i);
		cplx_poly c_i = poly_matrix_element(c, 1, i, 0);
		
		cplx covariance_coeffs[PARAM_N], center_coeffs[PARAM_N];
		cplx_poly covariance = covariance_coeffs, center = center_coeffs;
		memcpy(covariance, sch_comp_ii, PARAM_N * sizeof(cplx));
		memcpy(center, c_i, PARAM_N * sizeof(cplx));
		
		sample_fz(p_i, cplx_p, covariance, center, 0);
		
		// Compute the CRT representation of the newly-sampled polynomial
		
		// Update the center
		construct_new_center(c, sch_comp, cplx_p, i);
		}
	

	// Sample the last 2 polynomials with the specialized sample_2z algorithm (do not forget to copy the covariance first)
	cplx sch_comp_copy_coeffs[PARAM_N * 3];
	cplx_poly_matrix sch_comp_copy = sch_comp_copy_coeffs;
	
	memcpy(sch_comp_copy, sch_comp, PARAM_N * 3 * sizeof(cplx));
	
	cplx_poly sch_comp_00 = triangular_poly_matrix_element(sch_comp_copy, 0, 0);
	cplx_poly sch_comp_01 = triangular_poly_matrix_element(sch_comp_copy, 1, 0);
	transpose_cplx_poly(sch_comp_01, PARAM_N - 1); // transpose it since sch_comp[0, 1] = sch_comp[1, 0]^T
	cplx_poly sch_comp_11 = triangular_poly_matrix_element(sch_comp_copy, 1, 1);
	
	sample_2z((signed_scalar *) p, cplx_p, sch_comp_00, sch_comp_01, sch_comp_11, c, 0);
	}



extern unsigned long long timing_sampleG;
extern unsigned long long timing_samplePerturb;
extern unsigned long long timing_sampleArith;

/*
	Sample x in R^m with parameter zeta such that A * x = u mod q (using the trapdoor T and h_inv the inverse of the tag for A)
		x is in the CRT domain
		A is in the CRT domain
		u is in the CRT domain
		T is in the CRT domain
		cplx_T is in the complex CRT domain
		sch_comp is in the complex CRT domain
		h_inv is in the CRT domain
*/
void sample_pre(poly_matrix x, poly_matrix A, poly_matrix T, cplx_poly_matrix cplx_T, cplx_poly_matrix sch_comp, poly_matrix u)
	{
		unsigned long long begin_timing = 0;
		unsigned long long end_timing   = 0;
	
	// Sample a perturbation p in R^m
	signed_scalar p_coeffs[PARAM_N * PARAM_M];
	signed_poly_matrix p = p_coeffs;
	
	begin_timing = cpucycles_start();
	sample_perturb(p, cplx_T, sch_comp);
	end_timing = cpucycles_stop();
	timing_samplePerturb += (end_timing - begin_timing);
	

	
	// Add q to p's coeffs so that they are positive, and put p in the CRT domain
	begin_timing = cpucycles_start();
	for(int i = 0 ; i < PARAM_N * PARAM_M ; ++i)
		{
		p[i] += PARAM_Q;
		}
	
	freeze_signed_poly((poly) p, p, PARAM_N * PARAM_M - 1); // only for test parameters, useless if zeta << q

	matrix_ntt((poly_matrix) p, PARAM_M, 1);
	freeze_poly((poly_matrix) p, PARAM_N * PARAM_M - 1);
	

	
	// v <- u - A * p (in the CRT domain)
	poly_matrix v = x; // store v at the beginning of x
	
	
	multiply_by_A(v, A, (poly_matrix) p);
	
	
	sub_poly(v, u, v, PARAM_N * PARAM_D - 1);
	freeze_poly(v, PARAM_N * PARAM_D - 1);
	
	
	// Put v back into the normal domain
	matrix_invntt(v, PARAM_D, 1);
	freeze_poly(v, PARAM_N * PARAM_D - 1);
	
	
	// Sample z from the G-lattice with target v
	signed_poly_matrix z = (signed_poly_matrix) poly_matrix_element(x, 1, 2 * PARAM_D, 0); // store z at the end of x to get TI * z for cheaper
	end_timing = cpucycles_stop();
	timing_sampleArith += (end_timing - begin_timing);
	

	begin_timing = cpucycles_start();
	module_sample_G(z, v);
	end_timing = cpucycles_stop();
	timing_sampleG += (end_timing - begin_timing);
	
	
	// Make sure z has positive coefficients and put it in the CRT domain
	begin_timing = cpucycles_start();
	for(int i = 0 ; i < PARAM_N * PARAM_D * PARAM_K ; ++i)
		{
		z[i] += PARAM_Q;
		}
	
	matrix_ntt((poly_matrix) z, PARAM_D * PARAM_K, 1);
	
	
	// x <- p + TI * z (in the CRT domain)
	multiply_by_T(x, T, (poly_matrix) z); // only need to multiply by T since z is already at the end of x
	multiply_by_2pow32(x, PARAM_N * 2 * PARAM_D - 1);
	
	
	add_poly(x, x, (poly) p, PARAM_N * PARAM_M - 1);

	// Reduce x mod q one final time but keep it in the CRT domain
	freeze_poly(x, PARAM_N * PARAM_M - 1);

	end_timing = cpucycles_stop();
	timing_sampleArith += (end_timing - begin_timing);
	}