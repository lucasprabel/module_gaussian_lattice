#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
#include <time.h>
#include <string.h>

#include "common.h"
#include "random.h"
#include "sampling.h"
#include "arithmetic.h"

void test_random_bytes(int n)
	{
	uint8_t r[16];
	
	for(int i = 0 ; i < n ; ++i)
		{
		random_bytes(r);
		for(int j = 0 ; j < 16 ; ++j)
			{
			printf("%d ", r[j]);
			}
		printf("\n");
		}
	}

void time_SampleZ(int n)
	{
	clock_t t0, t1;
	real c = 0.5, sigma = PARAM_SIGMA;
	
	printf("Sampling from the discrete Gaussian over Z with center %f and parameter %f...\n", c, sigma);
	
	t0 = clock();
	for(int i = 0 ; i < n ; ++i)
		{
		SampleZ(c, sigma);
		}
	t1 = clock() - t0;
	
	printf("Done : %ld\n", t1);
	}

void test_SampleZ(int n)
	{
	real c = -10, sigma = PARAM_SIGMA;
	
	printf("Sampling from the discrete Gaussian over Z with center %f and parameter %f...\n", c, sigma);
	
	for(int i = 0 ; i < n ; ++i)
		{
		printf("%d\n", SampleZ(c, sigma));
		}
	}

void test_sample_D(int n)
	{
	signed_scalar z[PARAM_K];
	real c[PARAM_K], sigma = PARAM_ALPHA / 3;
	
	for(int i = 0 ; i < PARAM_K ; ++i)
		{
		c[i] = 0;
		}
	
	for(int i = 0 ; i < n ; ++i)
		{
		sample_D(z, c, sigma);
		print_signed_poly(z, PARAM_K - 1);
		printf("\n");
		}
	}

void test_scalar_sample_G(int n)
	{
	signed_scalar z[PARAM_K];
	scalar u, prod;
	
	for(int i = 0 ; i < n ; ++i)
		{
		random_poly(&u, 0);
		scalar_sample_G(z, u);
		
		//printf("\n");
		print_signed_poly(z, PARAM_K - 1);
		
		multiply_by_scalar_gadget_vector(&prod, z);
		
		printf("(%" PRIu32 ", %" PRIu32 ")\n", u, prod);
		}
	}

void test_ring_sample_G(int n)
	{
	signed_scalar z_coeffs[PARAM_K * PARAM_N];
	scalar prod_coeffs[PARAM_N], u_coeffs[PARAM_N];
	signed_poly_matrix z = z_coeffs;
	poly prod = prod_coeffs, u = u_coeffs;
	
	for(int i = 0 ; i < n ; ++i)
		{
		random_poly(u, PARAM_N - 1);
		// Sample z such that <g,z> = 0 mod q
		ring_sample_G(z, u);
		
		//print_signed_poly_matrix(z, PARAM_K, 1);
		
		// Check that <g,z> = u mod q
		multiply_by_ring_gadget_vector(prod, z);
		for(int j = 0 ; j < PARAM_N ; ++j)
			{
			if(prod[j] != u[j])
				{
				printf("%d %d : %" PRIu32 "\n", i, j, prod[j]);
				}
			}
		}
	}

void test_module_sample_G(int n)
	{
	signed_scalar z_coeffs[PARAM_D * PARAM_K * PARAM_N];
	scalar prod_coeffs[PARAM_D * PARAM_N], u_coeffs[PARAM_D * PARAM_N];
	signed_poly_matrix z = z_coeffs;
	poly_matrix prod = prod_coeffs, u = u_coeffs;
	
	for(int i = 0 ; i < n ; ++i)
		{
		random_poly(u, PARAM_D * PARAM_N - 1);
		// Sample z such that Gz = 0 mod q
		module_sample_G(z, u);
		
		// Check that Gz = u mod q
		multiply_by_module_gadget_matrix(prod, z);
		for(int j = 0 ; j < PARAM_D * PARAM_N ; ++j)
			{
			if(prod[j] != u[j])
				{
				printf("%d %d : %" PRIu32 "\n", i, j, prod[j]);
				}
			}
		}
	}

void time_module_sample_G(int n)
	{
	
	printf("Presampling from the module gadget lattice... (not really)");
	fflush(stdout);
	
	clock_t t0 = clock();
	for(int i = 0 ; i < n ; ++i)
		{
		// Sample z such that Gz = 0 mod q
		//module_sample_G(z);
		}
	clock_t t1 = clock() - t0;
	
	printf("Done : %ld\n", t1);
	}

void test_cplx_crt_representation()
	{
	cplx c[PARAM_N] = {0, 1};
	
	printf("c\n");
	print_cplx_poly(c, PARAM_N - 1);
	printf("\n\n");
	
	cplx_crt_representation(c);
	}

void test_stride()
	{
	cplx c[PARAM_N] = {-18.905300 -5.445592 * I, -2.879636 -18.912596 * I, -19.592393 -13.546650 * I, 4.393183 +10.168026 * I, 0.136162 -20.316123 * I, 3.686480 -0.978678 * I, 1.866518 +14.151483 * I, -23.605831 +26.153469 * I, 26.077725 -2.359032 * I, 0.935550 -42.339880 * I, 16.450492 +23.200703 * I, -30.146933 +12.744050 * I, 16.575329 +26.822235 * I, 7.684506 -30.880119 * I, 24.880183 -0.490263 * I, -30.051783 -11.144942 * I, 5.053629 +23.913111 * I, -14.522386 -8.780045 * I, 26.740038 +2.756713 * I, -22.119335 -11.544306 * I, -1.454303 +26.514699 * I, -22.126534 +10.208652 * I, -0.128176 -6.008698 * I, -0.598240 +14.701610 * I, 17.123404 +2.447041 * I, 1.310250 -1.233796 * I, -15.474274 -11.018091 * I, -18.199295 -18.885344 * I, -7.565432 -2.588454 * I, 6.736067 +32.195708 * I, -18.718270 -8.365083 * I, -0.119658 -9.539082 * I, -0.135211 +5.677274 * I, 4.899913 +35.867994 * I, 0.766370 -7.047868 * I, -3.924297 -5.574359 * I, 15.891594 +24.031707 * I, -16.889394 -9.214467 * I, 30.876032 +2.127784 * I, 4.412871 -11.316174 * I, 9.659418 -11.248520 * I, 11.736900 -1.703379 * I, 0.047987 -3.578209 * I, -21.354478 -8.680834 * I, 17.931173 -15.922143 * I, 11.086809 +10.779151 * I, 13.585724 -4.719505 * I, -20.261057 +3.066047 * I, -15.563322 -11.770299 * I, 1.697781 -3.486217 * I, -3.359731 -17.690137 * I, -5.850095 -7.690241 * I, -17.029591 +0.286818 * I, -9.129248 -2.186703 * I, 7.113253 -29.564135 * I, -29.604096 +3.055580 * I, 15.542087 +14.194458 * I, -14.025619 -3.674872 * I, 2.491711 -35.626212 * I, -11.420998 -8.325095 * I, 11.131402 -1.632941 * I, 0.791437 -22.735221 * I, -1.116591 -12.096086 * I, -16.697887 -7.442130 * I, 13.499630 +3.790266 * I, -0.784560 +1.939879 * I, 26.078747 +13.118937 * I, -16.586872 -10.726357 * I, 3.262905 +17.753371 * I, -9.476851 -10.528257 * I, -9.565060 +35.057456 * I, -7.186776 -6.232244 * I, -5.779794 -24.553142 * I, 10.850690 +19.532908 * I, -3.707857 +16.565053 * I, -9.593560 -12.896251 * I, 3.695477 +3.822116 * I, -21.589690 -18.257808 * I, 9.282997 -16.757499 * I, -12.372171 -0.535008 * I, 21.628090 +13.263022 * I, -0.681887 +14.109435 * I, -2.359898 -26.263666 * I, -1.368716 +0.460977 * I, 15.366485 -27.972409 * I, -15.488563 -16.103426 * I, -10.516728 +24.503821 * I, -29.643079 -31.866245 * I, 8.541624 -10.187008 * I, -18.151082 -9.052437 * I, -30.444488 -0.220254 * I, -12.923441 +11.772260 * I, -21.728899 +5.103852 * I, -11.927816 -18.923924 * I, -12.987039 +2.099621 * I, -5.919494 -4.199156 * I, -10.512516 +9.092107 * I, -9.327002 -11.033699 * I, 1.082346 -10.329277 * I, 0.146523 +14.564256 * I, -0.885394 -14.915647 * I, -19.295058 +39.632039 * I, 16.394718 -8.357002 * I, -14.648808 +10.368257 * I, 22.907285 -11.693653 * I, 19.282105 -5.731400 * I, -1.104540 +8.656422 * I, -10.154189 +0.968388 * I, 4.250382 +13.429783 * I, -5.868074 +3.224903 * I, 40.692294 -4.957691 * I, -18.502260 +5.358976 * I, -2.025901 +12.257041 * I, -20.311446 +12.891419 * I, -4.097861 +13.433358 * I, -24.148133 +6.621441 * I, 3.134287 +18.092785 * I, 15.584336 +10.226329 * I, 6.873673 +27.573397 * I, 1.523667 +16.010889 * I, -17.291153 -10.394464 * I, -15.253666 +5.080843 * I, -10.720671 -20.769515 * I, -29.890443 +4.051944 * I, -33.106716 +1.603919 * I, 14.394574 -2.295278 * I, 0.912826 -33.955949 * I, -0.254709 -12.664393 * I, -18.905300 +5.445592 * I, -2.879636 +18.912596 * I, -19.592393 +13.546650 * I, 4.393183 -10.168026 * I, 0.136162 +20.316123 * I, 3.686480 +0.978678 * I, 1.866518 -14.151483 * I, -23.605831 -26.153469 * I, 26.077725 +2.359032 * I, 0.935550 +42.339880 * I, 16.450492 -23.200703 * I, -30.146933 -12.744050 * I, 16.575329 -26.822235 * I, 7.684506 +30.880119 * I, 24.880183 +0.490263 * I, -30.051783 +11.144942 * I, 5.053629 -23.913111 * I, -14.522386 +8.780045 * I, 26.740038 -2.756713 * I, -22.119335 +11.544306 * I, -1.454303 -26.514699 * I, -22.126534 -10.208652 * I, -0.128176 +6.008698 * I, -0.598240 -14.701610 * I, 17.123404 -2.447041 * I, 1.310250 +1.233796 * I, -15.474274 +11.018091 * I, -18.199295 +18.885344 * I, -7.565432 +2.588454 * I, 6.736067 -32.195708 * I, -18.718270 +8.365083 * I, -0.119658 +9.539082 * I, -0.135211 -5.677274 * I, 4.899913 -35.867994 * I, 0.766370 +7.047868 * I, -3.924297 +5.574359 * I, 15.891594 -24.031707 * I, -16.889394 +9.214467 * I, 30.876032 -2.127784 * I, 4.412871 +11.316174 * I, 9.659418 +11.248520 * I, 11.736900 +1.703379 * I, 0.047987 +3.578209 * I, -21.354478 +8.680834 * I, 17.931173 +15.922143 * I, 11.086809 -10.779151 * I, 13.585724 +4.719505 * I, -20.261057 -3.066047 * I, -15.563322 +11.770299 * I, 1.697781 +3.486217 * I, -3.359731 +17.690137 * I, -5.850095 +7.690241 * I, -17.029591 -0.286818 * I, -9.129248 +2.186703 * I, 7.113253 +29.564135 * I, -29.604096 -3.055580 * I, 15.542087 -14.194458 * I, -14.025619 +3.674872 * I, 2.491711 +35.626212 * I, -11.420998 +8.325095 * I, 11.131402 +1.632941 * I, 0.791437 +22.735221 * I, -1.116591 +12.096086 * I, -16.697887 +7.442130 * I, 13.499630 -3.790266 * I, -0.784560 -1.939879 * I, 26.078747 -13.118937 * I, -16.586872 +10.726357 * I, 3.262905 -17.753371 * I, -9.476851 +10.528257 * I, -9.565060 -35.057456 * I, -7.186776 +6.232244 * I, -5.779794 +24.553142 * I, 10.850690 -19.532908 * I, -3.707857 -16.565053 * I, -9.593560 +12.896251 * I, 3.695477 -3.822116 * I, -21.589690 +18.257808 * I, 9.282997 +16.757499 * I, -12.372171 +0.535008 * I, 21.628090 -13.263022 * I, -0.681887 -14.109435 * I, -2.359898 +26.263666 * I, -1.368716 -0.460977 * I, 15.366485 +27.972409 * I, -15.488563 +16.103426 * I, -10.516728 -24.503821 * I, -29.643079 +31.866245 * I, 8.541624 +10.187008 * I, -18.151082 +9.052437 * I, -30.444488 +0.220254 * I, -12.923441 -11.772260 * I, -21.728899 -5.103852 * I, -11.927816 +18.923924 * I, -12.987039 -2.099621 * I, -5.919494 +4.199156 * I, -10.512516 -9.092107 * I, -9.327002 +11.033699 * I, 1.082346 +10.329277 * I, 0.146523 -14.564256 * I, -0.885394 +14.915647 * I, -19.295058 -39.632039 * I, 16.394718 +8.357002 * I, -14.648808 -10.368257 * I, 22.907285 +11.693653 * I, 19.282105 +5.731400 * I, -1.104540 -8.656422 * I, -10.154189 -0.968388 * I, 4.250382 -13.429783 * I, -5.868074 -3.224903 * I, 40.692294 +4.957691 * I, -18.502260 -5.358976 * I, -2.025901 -12.257041 * I, -20.311446 -12.891419 * I, -4.097861 -13.433358 * I, -24.148133 -6.621441 * I, 3.134287 -18.092785 * I, 15.584336 -10.226329 * I, 6.873673 -27.573397 * I, 1.523667 -16.010889 * I, -17.291153 +10.394464 * I, -15.253666 -5.080843 * I, -10.720671 +20.769515 * I, -29.890443 -4.051944 * I, -33.106716 -1.603919 * I, 14.394574 +2.295278 * I, 0.912826 +33.955949 * I, -0.254709 +12.664393 * I};
	
	printf("c\n");
	print_cplx_poly(c, PARAM_N - 1);
	printf("\n\n");
	
	for(int depth = 0 ; depth < 1 ; ++depth)
		{
		int deg = PARAM_N >> depth;
		
		for(int j = 0 ; j < (1 << depth) ; ++j)
			{
			cplx_poly c_j = &c[j * deg];
			
			stride(c_j, depth);
			}
		
		printf("strided c (%d)\n", depth);
		print_cplx_poly(c, PARAM_N - 1);
		printf("\n\n");
		}
	
	for(int depth = 0 ; depth >= 0 ; --depth)
		{
		int deg = PARAM_N >> depth;
		
		for(int j = 0 ; j < (1 << depth) ; ++j)
			{
			cplx_poly c_j = &c[j * deg];
			
			stride(c_j, depth);
			}
		
		printf("unstrided c (%d)\n", depth);
		print_cplx_poly(c, PARAM_N - 1);
		printf("\n\n");
		}
	}

void test_inverse_stride(void)
	{
	cplx c[PARAM_N] = {-216.0, 6233.0};
	
	printf("c\n");
	print_cplx_poly(c, 1);
	printf("\n\n");
	
	inverse_stride(c, LOG_N - 1);
	
	printf("c\n");
	print_cplx_poly(c, 1);
	printf("\n\n");
	}

void test_trap_gen(void)
	{
	// Generate A and its trapdoor (in the CRT domain)
	scalar A_coeffs[PARAM_N * PARAM_D * (PARAM_M - PARAM_D)], T_coeffs[PARAM_N * 2 * PARAM_D * PARAM_D * PARAM_K];
	poly_matrix A = A_coeffs, T = T_coeffs;
	
	TrapGen(A, T);
	
	printf("crt A\n");
	print_poly_matrix(A, PARAM_D, PARAM_M - PARAM_D);
	printf("crt T\n");
	print_poly_matrix(T, 2 * PARAM_D, PARAM_D * PARAM_K);
	
	// Check that A * TI = 0
	scalar TI_i_coeffs[PARAM_M * PARAM_N], prod_i_coeffs[PARAM_D * PARAM_N];
	poly_matrix TI_i = TI_i_coeffs, prod_i = prod_i_coeffs;
	
	for(int i = 0 ; i < PARAM_D * PARAM_K ; ++i)
		{
		poly_matrix T_i = poly_matrix_element(T, PARAM_D * PARAM_K, 0, i), TI_i_2d = poly_matrix_element(TI_i, 1, 2 * PARAM_D, 0), TI_i_2dplusi = poly_matrix_element(TI_i, 1, 2 * PARAM_D + i, 0);
		
		memcpy(TI_i, T_i, 2 * PARAM_D * PARAM_N * sizeof(scalar));
		memset(TI_i_2d, 0, PARAM_D * PARAM_K * PARAM_N * sizeof(scalar));
		
		for(int j = 0 ; j < PARAM_N ; j += SMALL_DEGREE)
			{
			TI_i_2dplusi[j] = 1;
			}
		
		multiply_by_A(prod_i, A, TI_i);
		
		printf("prod[%d]\n", i);
		print_poly_matrix(prod_i, PARAM_D, 1);
		}
	
	matrix_coeffs_representation(A, PARAM_D, PARAM_M - PARAM_D, LOG_R);
	matrix_coeffs_representation(T, 2 * PARAM_D, PARAM_D * PARAM_K, LOG_R);
	
	printf("A\n");
	print_poly_matrix(A, PARAM_D, PARAM_M - PARAM_D);
	printf("T\n");
	print_poly_matrix(T, 2 * PARAM_D, PARAM_D * PARAM_K);
	}


/*
	Checks whether the l by l triangular matrix M has only positive coefficients on its diagonal
*/
bool is_diagonally_positive(cplx_poly_matrix M, int l)
	{
	bool diag_pos = true;
	
	for(int i = 0 ; i < l ; ++i)
		{
		cplx_poly M_ii = triangular_poly_matrix_element(M, i, i);
		
		for(int j = 0 ; j < PARAM_N ; ++j)
			{
			if((fabs(cimag(M_ii[j])) != 0.0) || (creal(M_ii[j]) <= 0.0))
				{
				diag_pos = false;
				}
			}
		}
	
	return diag_pos;
	}

void test_sch_comp_computations(void)
	{
	scalar *T_coeffs = malloc(PARAM_N * 2 * PARAM_D * PARAM_D * PARAM_K * sizeof(scalar));
	cplx *sch_comp_coeffs = malloc(PARAM_N * PARAM_D * (2 * PARAM_D + 1) * sizeof(cplx)), *cplx_T_coeffs = malloc(PARAM_N * 2 * PARAM_D * PARAM_D * PARAM_K * sizeof(cplx));
	
	poly_matrix T = T_coeffs;
	cplx_poly_matrix sch_comp = sch_comp_coeffs, cplx_T = cplx_T_coeffs;
	
	// Generate T
	printf("Generating gaussian T...\n");
	
	SampleR_matrix_centered((signed_poly_matrix) T, 2*PARAM_D, PARAM_D * PARAM_K, PARAM_SIGMA);
	
	printf("Done\n");
	
	// Compute cplx_T (CRT domain)
	printf("Computing cplx_T in the CRT domain...\n");
	
	for(int i = 0 ; i < PARAM_N * 2 * PARAM_D * PARAM_D * PARAM_K ; ++i)
		{
		cplx_T[i] = (signed_scalar) T[i];
		}
	
	matrix_cplx_crt_representation(cplx_T, 2 * PARAM_D, PARAM_D * PARAM_K);
	
	printf("Done\n");
	
	// Compute the Schur complements and check if they are diagonally positive
	printf("Computing the Schur complements...\n");
	
	construct_T_schur_complement(sch_comp, cplx_T);
	
	printf("first sch_comp done\n");
	
	if(!is_diagonally_positive(sch_comp, 2 * PARAM_D))
		{
		print_triangular_cplx_poly_matrix(sch_comp, 2 * PARAM_D);
		return;
		}
	
	for(int i = 2 * PARAM_D - 1 ; i > 1 ; --i)
		{
		construct_schur_complement(sch_comp, i);
		
		printf("sch_comp %d done\n", i);
		
		if(!is_diagonally_positive(sch_comp, i))
			{
			print_triangular_cplx_poly_matrix(sch_comp, i);
			return;
			}
		}
	}

void is_zeta_big_enough(int n)
	{
	scalar *T_coeffs = malloc(PARAM_N * 2 * PARAM_D * PARAM_D * PARAM_K * sizeof(scalar));
	cplx *sch_comp_coeffs = malloc(PARAM_N * PARAM_D * (2 * PARAM_D + 1) * sizeof(cplx)), *cplx_T_coeffs = malloc(PARAM_N * 2 * PARAM_D * PARAM_D * PARAM_K * sizeof(cplx));
	
	poly_matrix T = T_coeffs;
	cplx_poly_matrix sch_comp = sch_comp_coeffs, cplx_T = cplx_T_coeffs;
	
	scalar p_coeffs[PARAM_M * PARAM_N];
	poly_matrix p = p_coeffs;
	
	bool failed = false;
	
	for(int i = 0 ; i < n ; ++i)
		{
		// Generate (T, cplx_T, sch_comp)
		SampleR_matrix_centered((signed_poly_matrix) T, 2*PARAM_D, PARAM_D * PARAM_K, PARAM_SIGMA);
		construct_complex_private_key(cplx_T, sch_comp, T);
				
		// Try to sample a perturbation
		sample_perturb((signed_poly_matrix) p, cplx_T, sch_comp);
		
		// Check that all Gaussian parameters were positive during the "problematic part" of sample_perturb (i.e. when it's not just calling SampleZ(0, sqrt(zeta^2 - alpha^2)))
		if(!is_zero_poly(p, 2 * PARAM_D * PARAM_N - 1))
			{
			printf("failed @ %d\n", i);
			
			print_poly_matrix(p, 2 * PARAM_D, 1);
			
			failed = true;
			
			break;
			}
		}
	
	if(failed)
		{
		printf("zeta is NOT enough\n");
		}
	else
		{
		printf("did not fail in %d tries\nzeta might be enough\n", n);
		}
	
	free(T);
	free(sch_comp);
	free(cplx_T);
	}

int main(int argc, char **argv)
	{
	if(argc < 2)
		{
		fprintf(stderr, "Usage : %s n\n", argv[0]);
		exit(EXIT_FAILURE);
		}
	
	int n = atoi(argv[1]);
	
	init_crt_trees();
	init_D_lattice_coeffs();
	init_cplx_roots_of_unity();
	
	random_bytes_init();
	
	//time_SampleZ(n);
	
	//test_random_bytes(n);
	//test_SampleZ(n);
	//test_sample_D(n);
	//test_scalar_sample_G(n);
	//test_ring_sample_G(n);
	//test_module_sample_G(n);
	//test_sample_perturb(n);
	//test_trap_gen();
	//test_sample_pre();
	//test_stride();
	//test_inverse_stride();
	//test_cplx_crt_representation();
	//test_sch_comp_computations();
	
	is_zeta_big_enough(n);
	
	//time_module_sample_G(n);
	
	return 0;
	}
