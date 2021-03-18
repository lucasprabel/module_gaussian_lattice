#include <stdbool.h>

#include "common.h"

void KeyGen(poly_matrix A, poly_matrix T, cplx_poly_matrix cplx_T, cplx_poly_matrix sch_comp);

void Sign(poly_matrix nu, uint8_t *r, poly_matrix A, poly_matrix T, cplx_poly_matrix cplx_T, cplx_poly_matrix sch_comp, uint8_t *m, int m_len);

bool Verify(poly_matrix nu, uint8_t *r, poly_matrix A, uint8_t *m, int m_len);
