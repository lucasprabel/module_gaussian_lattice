#include <inttypes.h>

#include "common.h"

void TrapGen(poly_matrix A, poly_matrix T);

//==============================================================================
// Samples from distribution D_{c,sigma}, ie                                              
// Samples an element in Z with probability proportionnal to e^{-(c-x)^2/2*(sigma^2)}    
//==============================================================================
signed int SampleZ(RR_t c, RR_t sigma);

void SampleR_centered(signed_poly f, RR_t sigma);

void SampleR_matrix_centered(signed_poly_matrix A, int l1, int l2, RR_t sigma);

void init_D_lattice_coeffs(void);

void sample_D(signed_scalar *z, real *c, real sigma);

void sample_G_perturb(real *p, real sigma);

void scalar_sample_G(signed_scalar *t, scalar u);

void ring_sample_G(signed_poly_matrix t, poly u);

void transpose_scalar_matrix(scalar *A_T, scalar *A, int l0, int l1);

void module_sample_G(signed_poly_matrix t, poly_matrix u);

void sample_2z(signed_scalar *q, cplx_poly cplx_q, cplx_poly a, cplx_poly b, cplx_poly d, cplx *c, int deg);

void sample_fz(signed_scalar *p, cplx_poly cplx_p, cplx_poly f, cplx *c, int deg);

void sample_perturb(signed_poly_matrix p, cplx_poly_matrix T, cplx_poly_matrix sch_comp);

void sample_pre(poly_matrix x, poly_matrix A, poly_matrix T, cplx_poly_matrix cplx_T, cplx_poly_matrix sch_comp, poly_matrix u);

extern real d_coeffs[PARAM_K];
extern real l_coeffs[PARAM_K];
extern real h_coeffs[PARAM_K];
