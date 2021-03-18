#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

#include "common.h"
#include "random.h"
#include "sampling.h"
#include "signature.h"
#include "arithmetic.h"
#include "hash.h"

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

unsigned long long timing_sign = 0;

unsigned long long timing_sampleZ_KG = 0;
unsigned long long timing_precomp_KG = 0;
unsigned long long timing_arith_KG = 0;
unsigned long long timing_keygen = 0;

unsigned long long timing_verify = 0;

void time_keygen(void)
	{
	timing_overhead = cpucycles_overhead();
	unsigned long long begin_timing = 0;
	unsigned long long end_timing = 0;

	
	scalar *A_coeffs = malloc(PARAM_N * PARAM_D * (PARAM_M - PARAM_D) * sizeof(scalar)), *T_coeffs = malloc(PARAM_N * 2 * PARAM_D * PARAM_D * PARAM_K * sizeof(scalar));
	cplx *sch_comp_coeffs = malloc(PARAM_N * PARAM_D * (2 * PARAM_D + 1) * sizeof(cplx)), *cplx_T_coeffs = malloc(PARAM_N * 2 * PARAM_D * PARAM_D * PARAM_K * sizeof(cplx));
	poly_matrix A = A_coeffs, T = T_coeffs;
	cplx_poly_matrix sch_comp = sch_comp_coeffs, cplx_T = cplx_T_coeffs;



	for(unsigned i = 0; i < NTESTS; ++i) {

		begin_timing = cpucycles_start();
		
		// Generate Keys
		KeyGen(A, T, cplx_T, sch_comp);
		
		end_timing = cpucycles_stop();
		timing_keygen += (end_timing - begin_timing);
 	}

 	timing_sampleZ_KG = timing_sampleZ_KG/NTESTS - timing_overhead;
 	timing_precomp_KG = timing_precomp_KG/NTESTS - timing_overhead;
 	timing_arith_KG = timing_arith_KG/NTESTS - timing_overhead;
 	double timing_total = timing_sampleZ_KG + timing_precomp_KG + timing_arith_KG;


 	printf("----------- KeyGen -----------\n");
 	printf("SampleZ      : %lld cycles (%.2lf ms) (%.2lf %% of KeyGen)\n", timing_sampleZ_KG, (timing_sampleZ_KG*1000)/CPU_CYCLES, (timing_sampleZ_KG/timing_total)*100.0);
 	printf("Precomp      : %lld cycles (%.2lf ms) (%.2lf %% of KeyGen)\n", timing_precomp_KG, (timing_precomp_KG*1000)/CPU_CYCLES, (timing_precomp_KG/timing_total)*100.0);
 	printf("Arith        : %lld cycles (%.2lf ms) (%.2lf %% of KeyGen)\n", timing_arith_KG, (timing_arith_KG*1000)/CPU_CYCLES, (timing_arith_KG/timing_total)*100.0);
 	printf("-----\n");
 	
 	timing_keygen = timing_keygen/NTESTS - timing_overhead;
 	printf("Total: %lld cycles (%.2lf ms)\n", timing_keygen, (timing_keygen*1000)/CPU_CYCLES);
 	
 	printf("\n\n");

 	free(A);
	free(T);
	free(sch_comp);
	free(cplx_T);
	}

void time_sign(void)
	{
	timing_overhead = cpucycles_overhead();
	unsigned long long begin_timing = 0;
	unsigned long long end_timing = 0;

	
	scalar *A_coeffs = malloc(PARAM_N * PARAM_D * (PARAM_M - PARAM_D) * sizeof(scalar)), *T_coeffs = malloc(PARAM_N * 2 * PARAM_D * PARAM_D * PARAM_K * sizeof(scalar));
	cplx *sch_comp_coeffs = malloc(PARAM_N * PARAM_D * (2 * PARAM_D + 1) * sizeof(cplx)), *cplx_T_coeffs = malloc(PARAM_N * 2 * PARAM_D * PARAM_D * PARAM_K * sizeof(cplx));
	poly_matrix A = A_coeffs, T = T_coeffs;
	cplx_poly_matrix sch_comp = sch_comp_coeffs, cplx_T = cplx_T_coeffs;

	// Generate Keys
	KeyGen(A, T, cplx_T, sch_comp);



	for(unsigned i = 0; i < NTESTS; ++i) {
		// Generate a message
		uint8_t m[MESSAGE_BYTES];
		for(int i = 0 ; i < MESSAGE_BYTES / 16 ; ++i) {
			random_bytes(&m[16 * i]);
		}

 		// Compute a signature (nu, r)
		uint8_t r[SALT_BYTES];
		scalar nu_coeffs[PARAM_N * PARAM_M];
		poly_matrix nu = nu_coeffs;

		begin_timing = cpucycles_start();
		Sign(nu, r, A, T, cplx_T, sch_comp, m, MESSAGE_BYTES);
		end_timing = cpucycles_stop();
		timing_sign += (end_timing - begin_timing);
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

 	printf("\n----------- Signature -----------\n");
 	timing_sign = timing_sign/NTESTS - timing_overhead;
 	printf("Signature: %lld cycles (%.2lf ms)\n", timing_sign, (timing_sign*1000)/CPU_CYCLES);
 	
 	printf("\n\n");

 	free(A);
	free(T);
	free(sch_comp);
	free(cplx_T);
	}

void time_verify(void)
	{
	timing_overhead = cpucycles_overhead();
	unsigned long long begin_timing = 0;
	unsigned long long end_timing = 0;

	
	scalar *A_coeffs = malloc(PARAM_N * PARAM_D * (PARAM_M - PARAM_D) * sizeof(scalar)), *T_coeffs = malloc(PARAM_N * 2 * PARAM_D * PARAM_D * PARAM_K * sizeof(scalar));
	cplx *sch_comp_coeffs = malloc(PARAM_N * PARAM_D * (2 * PARAM_D + 1) * sizeof(cplx)), *cplx_T_coeffs = malloc(PARAM_N * 2 * PARAM_D * PARAM_D * PARAM_K * sizeof(cplx));
	poly_matrix A = A_coeffs, T = T_coeffs;
	cplx_poly_matrix sch_comp = sch_comp_coeffs, cplx_T = cplx_T_coeffs;


	// Generate Keys
	KeyGen(A, T, cplx_T, sch_comp);

	// Generate a message
	uint8_t m[MESSAGE_BYTES];
	for(int i = 0 ; i < MESSAGE_BYTES / 16 ; ++i) {
		random_bytes(&m[16 * i]);
	}
	
	// Compute a signature (nu, r)
	uint8_t r[SALT_BYTES];
	scalar nu_coeffs[PARAM_N * PARAM_M];
	poly_matrix nu = nu_coeffs;
	
	Sign(nu, r, A, T, cplx_T, sch_comp, m, MESSAGE_BYTES);

	for(unsigned i = 0; i < NTESTS; ++i) {

		begin_timing = cpucycles_start();
		
		Verify(nu, r, A, m, MESSAGE_BYTES);
		
		end_timing = cpucycles_stop();
		timing_verify += (end_timing - begin_timing);
 	}


 	printf("----------- Verify -----------\n");
 	
 	timing_verify = timing_verify/NTESTS - timing_overhead;
 	printf("Total: %lld cycles (%.2lf ms)\n", timing_verify, (timing_verify*1000)/CPU_CYCLES);
 	
 	printf("\n\n");

 	free(A);
	free(T);
	free(sch_comp);
	free(cplx_T);
	}

int main(void) {

	init_D_lattice_coeffs();
	init_cplx_roots_of_unity();
	
 	time_keygen();
 	
 	time_sign();
 	
 	time_verify();
	

	return 0;
}

