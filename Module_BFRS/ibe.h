#include <stdbool.h>

#include "common.h"

void Setup(poly_matrix A, poly_matrix T, cplx_poly_matrix cplx_T, cplx_poly_matrix sch_comp, poly_matrix u);

void Extract(poly_matrix x, poly_matrix A, poly_matrix u, poly_matrix T, cplx_poly_matrix cplx_T, cplx_poly_matrix sch_comp, scalar *id);

void Encrypt(poly_matrix A, poly_matrix u, scalar *id, poly M, poly_matrix b, poly c);

void Decrypt(poly_matrix x, poly_matrix b, poly c, poly M);
