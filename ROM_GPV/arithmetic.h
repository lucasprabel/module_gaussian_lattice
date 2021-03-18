#ifndef __ARITHMETIC_H__
#define __ARITHMETIC_H__

#include <inttypes.h>
#include <stdbool.h>


#include "common.h"

void ntt(uint32_t p[PARAM_N]);

void invntt_frominvmont(uint32_t p[PARAM_N]);

void invntt_without_mont(scalar p[PARAM_N]);

scalar barrett_reduce(double_scalar a);

uint32_t csubq(uint32_t a);

uint32_t montgomery_reduce(uint64_t a);

void poly_pointwise_invmontgomery(poly c, const poly a, const poly b);

void matrix_ntt(poly_matrix A, int l1, int l2);

void matrix_invntt(poly_matrix A, int l1, int l2);

void matrix_invntt_without_mont(poly_matrix A, int l1, int l2);

void multiply_by_one(poly f, int deg);

void divide_by_2pow32(poly f, int deg);

void multiply_by_2pow32(poly f, int deg);

void mul_crt_poly_matrix(poly_matrix C, poly_matrix A, poly_matrix B, int l1, int l2, int l3);

scalar reduce_naive(scalar x);

signed_scalar reduce_signed_double_naive(signed_double_scalar x);

scalar reduce_signed_double_to_positive_naive(signed_double_scalar x);

scalar reduce_signed_naive(signed_scalar x);

scalar reduce_sparse(scalar x);

scalar reduce_double_naive(double_scalar x);

scalar reduce_double_sparse(double_scalar x);

scalar reduce_double_montgomery(double_scalar x);

void alloc_poly(poly f, int deg);

void alloc_double_poly(double_poly f, int deg);

void free_poly(poly f);

void free_double_poly(double_poly f);

void print_poly(poly f, int deg);

void print_double_poly(double_poly f, int deg);

void print_signed_poly(signed_poly f, int deg);

void print_signed_double_poly(signed_double_poly f, int deg);

void print_poly_matrix(poly_matrix A, int l1, int l2);

void print_signed_poly_matrix(signed_poly_matrix A, int l1, int l2);

bool equals_poly(poly f, poly g, int deg);

void print_cplx_poly(cplx_poly f, int deg);

void print_cplx_poly_matrix(cplx_poly_matrix A, int l1, int l2);

void print_triangular_cplx_poly_matrix(cplx_poly_matrix A, int l);

void zero_poly(poly f, int deg);

bool is_zero_poly(poly f, int deg);

void freeze_poly(poly f, int deg);

void freeze_double_poly(poly f, double_poly double_f, int deg);

void freeze_upper_half_double_poly(double_poly double_f, int deg);

void freeze_signed_poly(poly f_out, signed_poly f_in, int deg);

void freeze_signed_double_poly(signed_poly f_out, signed_double_poly f_in, int deg);

void freeze_signed_double_poly_to_positive(poly f_out, signed_double_poly f_in, int deg);

void random_poly(poly f, int deg);

void add_poly(poly h, poly f, poly g, int deg);

void add_double_poly(double_poly h, double_poly f, double_poly g, int deg);

void add_signed_poly(signed_poly h, signed_poly f, signed_poly g, int deg);

void sub_signed_poly(signed_poly h, signed_poly f, signed_poly g, int deg);

void sub_poly(poly h, poly f, poly g, int deg);

void add_to_poly_matrix(poly_matrix A, poly_matrix B, int l1, int l2);

void multiply_by_scalar_gadget_vector(scalar *prod, signed_scalar *v);

void multiply_by_ring_gadget_vector(poly prod, signed_poly_matrix v);

void multiply_by_module_gadget_matrix(poly_matrix res, signed_poly_matrix v);

void add_ring_gadget_vector(poly_matrix v);

void multiply_by_A(poly_matrix y, poly_matrix A, poly_matrix x);

void multiply_by_TI(poly_matrix y, poly_matrix T, poly_matrix x);

void multiply_by_T(poly_matrix y, poly_matrix T, poly_matrix x);

double_scalar norm_squared(poly_matrix v, int l);

void zero_cplx_poly(cplx_poly f, int deg);

void add_cplx_poly(cplx_poly h, cplx_poly f, cplx_poly g, int deg);

void fma_cplx_poly(cplx_poly h, cplx_poly f, cplx_poly g, int deg);

void sub_cplx_poly(cplx_poly h, cplx_poly f, cplx_poly g, int deg);

void fmsub_cplx_poly(cplx_poly h, cplx_poly f, cplx_poly g, int deg);

void inv_cplx_poly(cplx_poly h, cplx_poly f, int deg);

void mul_cplx_poly(cplx_poly h, cplx_poly f, cplx_poly g, int deg);

void div_cplx_poly(cplx_poly h, cplx_poly f, cplx_poly g, int deg);

void mul_cplx_poly_by_transpose(cplx_poly h, cplx_poly f, cplx_poly g, int deg);

void mul_cplx_poly_by_own_transpose(cplx_poly h, cplx_poly f, int deg);

void cplx_fma_transpose(cplx_poly h, cplx_poly f, cplx_poly g, int deg);

void fmsub_transpose_cplx_poly(cplx_poly h, cplx_poly f, cplx_poly g, int deg);

void transpose_cplx_poly(cplx_poly f, int deg);

void stride(cplx_poly f, int depth);

void inverse_stride(cplx_poly f, int depth);

void scalar_stride(signed_scalar *p, int size);

void reduce_mod_sparse_cplx_poly(cplx_poly f, int deg, cplx c);

void cplx_crt_representation(cplx_poly f);

void matrix_cplx_crt_representation(cplx_poly_matrix M, int l1, int l2);

void construct_T_schur_complement(cplx_poly_matrix sch_comp, cplx_poly_matrix T);

void construct_schur_complement_and_center(cplx_poly_matrix M, cplx_poly_matrix c, cplx_poly d, cplx_poly x1, int l);

void construct_schur_complement(cplx_poly_matrix M, int l);

void construct_new_center(cplx_poly_matrix c, cplx_poly_matrix M, cplx_poly q1, int l);

void construct_first_center(cplx_poly_matrix c, cplx_poly_matrix T, cplx_poly_matrix p);

void construct_complex_private_key(cplx_poly_matrix cplx_T, cplx_poly_matrix sch_comp, poly_matrix T);

void init_cplx_roots_of_unity(void);

#endif
