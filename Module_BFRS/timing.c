#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

#include "common.h"
#include "random.h"
#include "sampling.h"
#include "signature.h"
#include "ibe.h"
#include "arithmetic.h"

#include "cpucycles.h"

#define MESSAGE_BYTES 512 // keep it divisible by 16
#define NTESTS 100
#define CPU_CYCLES (1.9 * 1000000000.0)

unsigned long long timing_overhead;
unsigned long long timing_sampleZ_G = 0;
unsigned long long timing_sampleZ_P = 0;
unsigned long long timing_sampleG = 0;
unsigned long long timing_samplePerturb = 0;
unsigned long long timing_sampleArith = 0;

unsigned long long timing_extract = 0;
unsigned long long timing_encrypt = 0;
unsigned long long timing_decrypt = 0;

unsigned long long timing_sampleZ_Setup = 0;
unsigned long long timing_precomp_Setup = 0;
unsigned long long timing_arith_Setup = 0;
unsigned long long timing_setup = 0;

void time_setup(void)
	{
	timing_overhead = cpucycles_overhead();
	unsigned long long begin_timing = 0;
	unsigned long long end_timing = 0;

	
	scalar *A_coeffs = malloc(PARAM_N * PARAM_D * (PARAM_M - PARAM_D) * sizeof(scalar)), *T_coeffs = malloc(PARAM_N * 2 * PARAM_D * PARAM_D * PARAM_K * sizeof(scalar));
	cplx *sch_comp_coeffs = malloc(PARAM_N * PARAM_D * (2 * PARAM_D + 1) * sizeof(cplx)), *cplx_T_coeffs = malloc(PARAM_N * 2 * PARAM_D * PARAM_D * PARAM_K * sizeof(cplx));
	scalar *u_coeffs = malloc(PARAM_D * PARAM_N * sizeof(scalar));
	poly_matrix A = A_coeffs, T = T_coeffs, u = u_coeffs;
	cplx_poly_matrix sch_comp = sch_comp_coeffs, cplx_T = cplx_T_coeffs;



	for(unsigned i = 0; i < NTESTS; ++i) {

		begin_timing = cpucycles_start();
		
		// Generate Keys
		Setup(A, T, cplx_T, sch_comp, u);
		
		end_timing = cpucycles_stop();
		timing_setup += (end_timing - begin_timing);
 	}

 	timing_sampleZ_Setup = timing_sampleZ_Setup/NTESTS - timing_overhead;
 	timing_precomp_Setup = timing_precomp_Setup/NTESTS - timing_overhead;
 	timing_arith_Setup = timing_arith_Setup/NTESTS - timing_overhead;
 	double timing_total = timing_sampleZ_Setup + timing_precomp_Setup + timing_arith_Setup;


 	printf("----------- Setup -----------\n");
 	printf("SampleZ      : %lld cycles (%.2lf ms) (%.2lf %% of Setup)\n", timing_sampleZ_Setup, (timing_sampleZ_Setup*1000)/CPU_CYCLES, (timing_sampleZ_Setup/timing_total)*100.0);
 	printf("Precomp      : %lld cycles (%.2lf ms) (%.2lf %% of Setup)\n", timing_precomp_Setup, (timing_precomp_Setup*1000)/CPU_CYCLES, (timing_precomp_Setup/timing_total)*100.0);
 	printf("Arith        : %lld cycles (%.2lf ms) (%.2lf %% of Setup)\n", timing_arith_Setup, (timing_arith_Setup*1000)/CPU_CYCLES, (timing_arith_Setup/timing_total)*100.0);
 	printf("-----\n");
 	
 	timing_setup = timing_setup/NTESTS - timing_overhead;
 	printf("Total: %lld cycles (%.2lf ms)\n", timing_setup, (timing_setup*1000)/CPU_CYCLES);
 	
 	printf("\n\n");

 	free(A);
	free(T);
	free(sch_comp);
	free(cplx_T);
	free(u);
	}

void time_extract(void)
	{
	timing_overhead = cpucycles_overhead();
	unsigned long long begin_timing = 0;
	unsigned long long end_timing = 0;

	
	scalar *A_coeffs = malloc(PARAM_N * PARAM_D * (PARAM_M - PARAM_D) * sizeof(scalar)), *T_coeffs = malloc(PARAM_N * 2 * PARAM_D * PARAM_D * PARAM_K * sizeof(scalar));
	cplx *sch_comp_coeffs = malloc(PARAM_N * PARAM_D * (2 * PARAM_D + 1) * sizeof(cplx)), *cplx_T_coeffs = malloc(PARAM_N * 2 * PARAM_D * PARAM_D * PARAM_K * sizeof(cplx));
	scalar *u_coeffs = malloc(PARAM_D * PARAM_N * sizeof(scalar));
	poly_matrix A = A_coeffs, T = T_coeffs, u = u_coeffs;
	cplx_poly_matrix sch_comp = sch_comp_coeffs, cplx_T = cplx_T_coeffs;

	// Generate Keys
	Setup(A, T, cplx_T, sch_comp, u);



	for(unsigned i = 0; i < NTESTS; ++i) {
		// Generate an identity
		scalar id_coeffs[PARAM_N] = {0};
		poly id = id_coeffs;
		
		random_poly(id, SMALL_DEGREE - 1);

 		// Compute a signature (nu)
		scalar nu_coeffs[PARAM_N * PARAM_M];
		poly_matrix nu = nu_coeffs;

		begin_timing = cpucycles_start();
		Extract(nu, A, u, T, cplx_T, sch_comp, id);
		end_timing = cpucycles_stop();
		timing_extract += (end_timing - begin_timing);
 	}

 	timing_sampleZ_G = timing_sampleZ_G/NTESTS - timing_overhead;
 	timing_sampleZ_P = timing_sampleZ_P/NTESTS - timing_overhead;
 	unsigned long long timing_sampleZ = timing_sampleZ_G + timing_sampleZ_P;
 	
 	timing_sampleG = timing_sampleG/NTESTS - timing_overhead;
 	timing_samplePerturb = timing_samplePerturb/NTESTS - timing_overhead;
 	timing_sampleArith = timing_sampleArith/NTESTS - timing_overhead;
 	double timing_samplePre = timing_sampleG + timing_samplePerturb + timing_sampleArith;


 	printf("----------- SampleZ -----------\n");
 	printf("SampleZ      : %lld cycles (%.2lf ms) (%.2lf %% of SampleG)\n", timing_sampleZ_G, (timing_sampleZ_G*1000)/CPU_CYCLES, ((double) timing_sampleZ_G/timing_sampleG)*100.0);
 	printf("SampleZ      : %lld cycles (%.2lf ms) (%.2lf %% of SampleP)\n", timing_sampleZ_P, (timing_sampleZ_P*1000)/CPU_CYCLES, ((double) timing_sampleZ_P/timing_samplePerturb)*100.0);
 	printf("SampleZ      : %lld cycles (%.2lf ms) (%.2lf %% of SamplePre)\n", timing_sampleZ, (timing_sampleZ*1000)/CPU_CYCLES, ((double) timing_sampleZ/timing_samplePre)*100.0);


 	printf("\n----------- SamplePre -----------\n");
 	printf("SampleG      : %lld cycles (%.2lf ms) (%.2lf %% of SamplePre)\n", timing_sampleG, (timing_sampleG*1000)/CPU_CYCLES, (timing_sampleG/timing_samplePre)*100.0);
 	printf("G-SampleZ      : %lld cycles (%.2lf ms) (%.2lf %% of SamplePre)\n", timing_sampleZ_G, (timing_sampleZ_G*1000)/CPU_CYCLES, (timing_sampleZ_G/timing_samplePre)*100.0);
 	printf("G-Arith        : %lld cycles (%.2lf ms) (%.2lf %% of SamplePre)\n", timing_sampleG - timing_sampleZ_G, ((timing_sampleG - timing_sampleZ_G)*1000)/CPU_CYCLES, ((timing_sampleG - timing_sampleZ_G)/timing_samplePre)*100.0);
 	printf("-----\n");
 	
 	printf("SamplePerturb: %lld cycles (%.2lf  ms) (%.2lf %% of SamplePre)\n", timing_samplePerturb, (timing_samplePerturb*1000)/CPU_CYCLES, (timing_samplePerturb/timing_samplePre)*100.0);
 	printf("P-SampleZ      : %lld cycles (%.2lf ms) (%.2lf %% of SamplePre)\n", timing_sampleZ_P, (timing_sampleZ_P*1000)/CPU_CYCLES, (timing_sampleZ_P/timing_samplePre)*100.0);
 	printf("P-Arith        : %lld cycles (%.2lf ms) (%.2lf %% of SamplePre)\n", timing_samplePerturb - timing_sampleZ_P, ((timing_samplePerturb - timing_sampleZ_P)*1000)/CPU_CYCLES, ((timing_samplePerturb - timing_sampleZ_P)/timing_samplePre)*100.0);
 	printf("-----\n");
 	
 	printf("Arith        : %lld  cycles (%.2lf  ms) (%.2lf  %% of SamplePre)\n", timing_sampleArith, (timing_sampleArith*1000)/CPU_CYCLES, (timing_sampleArith/timing_samplePre)*100.0);

 	printf("\n----------- Extract -----------\n");
 	timing_extract = timing_extract/NTESTS - timing_overhead;
 	printf("Extract: %lld cycles (%.2lf ms)\n", timing_extract, (timing_extract*1000)/CPU_CYCLES);
 	
 	printf("\n\n");

 	free(A);
	free(T);
	free(sch_comp);
	free(cplx_T);
	free(u);
	}

void time_encrypt(void)
	{
	timing_overhead = cpucycles_overhead();
	unsigned long long begin_timing = 0;
	unsigned long long end_timing = 0;

	
	scalar *A_coeffs = malloc(PARAM_N * PARAM_D * (PARAM_M - PARAM_D) * sizeof(scalar)), *T_coeffs = malloc(PARAM_N * 2 * PARAM_D * PARAM_D * PARAM_K * sizeof(scalar));
	cplx *sch_comp_coeffs = malloc(PARAM_N * PARAM_D * (2 * PARAM_D + 1) * sizeof(cplx)), *cplx_T_coeffs = malloc(PARAM_N * 2 * PARAM_D * PARAM_D * PARAM_K * sizeof(cplx));
	scalar *u_coeffs = malloc(PARAM_D * PARAM_N * sizeof(scalar));
	poly_matrix A = A_coeffs, T = T_coeffs, u = u_coeffs;
	cplx_poly_matrix sch_comp = sch_comp_coeffs, cplx_T = cplx_T_coeffs;

	// Generate Keys
	Setup(A, T, cplx_T, sch_comp, u);
	
	for(unsigned i = 0; i < NTESTS; ++i) {

		// Generate a message
		scalar m_coeffs[PARAM_N] = {0};
		poly m = m_coeffs;
		
		random_poly(m, SMALL_DEGREE - 1);

		// Generate an identity
		scalar id_coeffs[PARAM_N] = {0};
		poly id = id_coeffs;
	
		random_poly(m, SMALL_DEGREE - 1);

 		// Compute b and c
 		scalar *b_coeffs = malloc(PARAM_N * PARAM_M * sizeof(scalar));
 		scalar *c_coeffs = malloc(PARAM_N * sizeof(scalar));

 		poly_matrix b = b_coeffs;
 		poly c = c_coeffs;


		begin_timing = cpucycles_start();
		Encrypt(A, u, id, m, b, c);
		end_timing = cpucycles_stop();
		timing_encrypt += (end_timing - begin_timing);
 	}

 	timing_sampleZ_G = timing_sampleZ_G/NTESTS - timing_overhead;
 	timing_sampleZ_P = timing_sampleZ_P/NTESTS - timing_overhead;
 	unsigned long long timing_sampleZ = timing_sampleZ_G + timing_sampleZ_P;
 	
 	timing_sampleG = timing_sampleG/NTESTS - timing_overhead;
 	timing_samplePerturb = timing_samplePerturb/NTESTS - timing_overhead;
 	timing_sampleArith = timing_sampleArith/NTESTS - timing_overhead;
 	double timing_samplePre = timing_sampleG + timing_samplePerturb + timing_sampleArith;


 	printf("----------- SampleZ -----------\n");
 	printf("SampleZ      : %lld cycles (%.2lf ms) (%.2lf %% of SampleG)\n", timing_sampleZ_G, (timing_sampleZ_G*1000)/CPU_CYCLES, ((double) timing_sampleZ_G/timing_sampleG)*100.0);
 	printf("SampleZ      : %lld cycles (%.2lf ms) (%.2lf %% of SampleP)\n", timing_sampleZ_P, (timing_sampleZ_P*1000)/CPU_CYCLES, ((double) timing_sampleZ_P/timing_samplePerturb)*100.0);
 	printf("SampleZ      : %lld cycles (%.2lf ms) (%.2lf %% of SamplePre)\n", timing_sampleZ, (timing_sampleZ*1000)/CPU_CYCLES, ((double) timing_sampleZ/timing_samplePre)*100.0);


 	printf("\n----------- SamplePre -----------\n");
 	printf("SampleG      : %lld cycles (%.2lf ms) (%.2lf %% of SamplePre)\n", timing_sampleG, (timing_sampleG*1000)/CPU_CYCLES, (timing_sampleG/timing_samplePre)*100.0);
 	printf("G-SampleZ      : %lld cycles (%.2lf ms) (%.2lf %% of SamplePre)\n", timing_sampleZ_G, (timing_sampleZ_G*1000)/CPU_CYCLES, (timing_sampleZ_G/timing_samplePre)*100.0);
 	printf("G-Arith        : %lld cycles (%.2lf ms) (%.2lf %% of SamplePre)\n", timing_sampleG - timing_sampleZ_G, ((timing_sampleG - timing_sampleZ_G)*1000)/CPU_CYCLES, ((timing_sampleG - timing_sampleZ_G)/timing_samplePre)*100.0);
 	printf("-----\n");
 	
 	printf("SamplePerturb: %lld cycles (%.2lf  ms) (%.2lf %% of SamplePre)\n", timing_samplePerturb, (timing_samplePerturb*1000)/CPU_CYCLES, (timing_samplePerturb/timing_samplePre)*100.0);
 	printf("P-SampleZ      : %lld cycles (%.2lf ms) (%.2lf %% of SamplePre)\n", timing_sampleZ_P, (timing_sampleZ_P*1000)/CPU_CYCLES, (timing_sampleZ_P/timing_samplePre)*100.0);
 	printf("P-Arith        : %lld cycles (%.2lf ms) (%.2lf %% of SamplePre)\n", timing_samplePerturb - timing_sampleZ_P, ((timing_samplePerturb - timing_sampleZ_P)*1000)/CPU_CYCLES, ((timing_samplePerturb - timing_sampleZ_P)/timing_samplePre)*100.0);
 	printf("-----\n");
 	
 	printf("Arith        : %lld  cycles (%.2lf  ms) (%.2lf  %% of SamplePre)\n", timing_sampleArith, (timing_sampleArith*1000)/CPU_CYCLES, (timing_sampleArith/timing_samplePre)*100.0);

 	printf("\n----------- Encrypt -----------\n");
 	timing_encrypt = timing_encrypt/NTESTS - timing_overhead;
 	printf("Encrypt: %lld cycles (%.2lf ms)\n", timing_encrypt, (timing_encrypt*1000)/CPU_CYCLES);
 	
 	printf("\n\n");

 	free(A);
	free(T);
	free(sch_comp);
	free(cplx_T);
	free(u);
	}

void time_decrypt(void)
	{
	timing_overhead = cpucycles_overhead();
	unsigned long long begin_timing = 0;
	unsigned long long end_timing = 0;

	
	scalar *A_coeffs = malloc(PARAM_N * PARAM_D * (PARAM_M - PARAM_D) * sizeof(scalar)), *T_coeffs = malloc(PARAM_N * 2 * PARAM_D * PARAM_D * PARAM_K * sizeof(scalar));
	cplx *sch_comp_coeffs = malloc(PARAM_N * PARAM_D * (2 * PARAM_D + 1) * sizeof(cplx)), *cplx_T_coeffs = malloc(PARAM_N * 2 * PARAM_D * PARAM_D * PARAM_K * sizeof(cplx));
	scalar *u_coeffs = malloc(PARAM_D * PARAM_N * sizeof(scalar));
	poly_matrix A = A_coeffs, T = T_coeffs, u = u_coeffs;
	cplx_poly_matrix sch_comp = sch_comp_coeffs, cplx_T = cplx_T_coeffs;

	// Generate Keys
	Setup(A, T, cplx_T, sch_comp, u);



	// Generate a message
	scalar m_coeffs[PARAM_N] = {0};
	poly m = m_coeffs;
		
	random_poly(m, SMALL_DEGREE - 1);

	// Generate an identity
	scalar id_coeffs[PARAM_N] = {0};
	poly id = id_coeffs;
	
	random_poly(id, SMALL_DEGREE - 1);

	// Generate nu = sk_id
	scalar nu_coeffs[PARAM_N * PARAM_M];
	poly_matrix nu = nu_coeffs;


	Extract(nu, A, u, T, cplx_T, sch_comp, id);

 	// Compute b and c
 	scalar *b_coeffs = malloc(PARAM_N * PARAM_M * sizeof(scalar));
 	scalar *c_coeffs = malloc(PARAM_N * sizeof(scalar));

 	poly_matrix b = b_coeffs;
 	poly c = c_coeffs;


	
	Encrypt(A, u, id, m, b, c);


	scalar *M_coeffs = malloc(PARAM_N * sizeof(scalar));
	poly M = M_coeffs;

	for(unsigned i = 0; i < NTESTS; ++i) {


		begin_timing = cpucycles_start();
		
		Decrypt(nu, b, c, M);
		
		end_timing = cpucycles_stop();
		timing_decrypt += (end_timing - begin_timing);
 	}


 	printf("----------- Decrypt -----------\n");
 	
 	timing_decrypt = timing_decrypt/NTESTS - timing_overhead;
 	printf("Total: %lld cycles (%.2lf ms)\n\n\n", timing_decrypt, (timing_decrypt*1000)/CPU_CYCLES);
 	
 	printf("\n\n");

 	free(A);
	free(T);
	free(sch_comp);
	free(cplx_T);
	free(u);
	free(b);
	free(c);
	free(M);
	}



int main(void) {
	
	init_crt_trees();
	init_D_lattice_coeffs();
	init_cplx_roots_of_unity();
	
 	time_setup();
 	
 	time_extract();
 	
 	time_encrypt();

 	time_decrypt();
	

	return 0;
}

