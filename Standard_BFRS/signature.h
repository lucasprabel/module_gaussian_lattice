#include <stdbool.h>

#include "common.h"

void KeyGen(poly_matrix A, poly_matrix T, cplx_poly_matrix cplx_T, cplx_poly_matrix sch_comp);

void Sign(poly_matrix nu, poly_matrix A, poly_matrix T, cplx_poly_matrix cplx_T, cplx_poly_matrix sch_comp, scalar *m);

bool Verify(poly_matrix nu, poly_matrix A, scalar *m);
