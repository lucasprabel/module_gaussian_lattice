#q = 7681
#q = 8388593
#q = 8000033
#q = 1073741969
#q = 1032193
#q = 17
#q = 134218289
#q = 33554641
#q = 13
#q = 1073740609
#q = 1073739937
#q = 1073740609
q = 1073741441
n = 256
#n = 512
#r = 16
#r = 32
r = 64
log_r = ceil(log(r, 2))
K = ceil(log(q, 2))
d = 2
m = d * (K + 2)

ZZx.<xz> = ZZ[]
GFq = GF(q)
R.<xx> = GF(q)[]
Rq.<x> = R.quotient_ring((xx^(n) + 1))

def my_ext_euclid(f, g):
	(r0, r1, v0, v1) = (g, f, R(0), R(1))
	while r0.degree() != 0:
		(q, r) = r0.quo_rem(r1)
		(v0, v1) = (v1, v0 - q*v1)
		(r0, r1) = (r1, r0 - q*r1)
		print "q\n", q, "\nr\n", r, "\nv0\n", v0, "\nv1\n", v1, "\nr0\n", r0, "\nr1\n", r1, "\n"

def cyclotomic_factorisation_tree(n):
    tree = [[GFq(1)]]
    i = 1
    for i in range(1, n):
        primitive_roots = []
        for r in tree[i-1]:
            primitive_roots += (-r).square_root(extend=False, all=True)
        tree.append(primitive_roots)
    return tree

def bezout_coefficients_tree(c_tree):
    b_tree = [[0]]
    for i in range(1, len(c_tree)):
        b_tree.append(2**i*[0])
        for j in range(2**(i-1)):
            #inv1 = (c_tree[i][2*j] - c_tree[i][2*j+1]).inverse_of_unit()
            inv1 = GFq((ZZ(c_tree[i][2*j] - c_tree[i][2*j+1])).inverse_mod(q))
            b_tree[i][2*j] = inv1
            b_tree[i][2*j+1] = -inv1
    return b_tree

def flatten_tree(tree):
    l = []
    for list in tree:
        l += list
    return l

def ext_euclid(f, nb):
    (r0, r1, v0, v1) = (xx^n + 1, f, R(0), R(1))
    for i in range(nb):
        (quo, rem) = r0.quo_rem(r1)
        (r0, r1) = (r1, rem)
        (v0, v1) = (v1, v0 - quo * v1)
    return (r0, r1, quo, v0, v1)

def crt(f):
	crt_f = [R(f.list()).quo_rem(c_f)[1].list() for c_f in cyclo_factors]
	return [crt_f_i if crt_f_i != [] else [0] for crt_f_i in crt_f]

def my_norm(v):
	coeffs = [ZZ(v_ij) for v_i in v.list() for v_ij in v_i]
	centered_v = vector(ZZ, [v_ij if v_ij < q/2 else (v_ij - q) for v_ij in coeffs])
	return centered_v.norm()

c_tree = cyclotomic_factorisation_tree(log_r + 1)
b_tree = bezout_coefficients_tree(c_tree)

print flatten_tree(c_tree)
print flatten_tree(b_tree)

cyclo_factors = [xx^(n/r) + c for c in c_tree[-1]]

Rx.<xxr> = RLF[]
Rr.<xr> = Rx.quotient_ring((xxr^n + 1))
Cx.<xxc> = CIF[]
Rc.<xc> = Cx.quotient_ring((xxc^n + 1))

if(n == 8):
	cplx_prim_roots = [-e^(I*k*pi/n) for k in [3, 11, 15, 7, 13, 5, 1, 9]]
	cplx_roots = [[1], [I, -I], [-e^(I*k*pi/4) for k in [7, 3, 1, 5]], cplx_prim_roots]

def stride(f):
	return Cx([f[i] for i in range(0, len(f.list()), 2)]), Cx([f[i] for i in range(1, len(f.list()), 2)])

def inv_cplx_crt(f):
	return Cx.lagrange_polynomial(zip(cplx_prim_roots, f))

def is_almost_real(f):
	return max([abs(f_i.imag()) for f_i in f])

def cplx_crt(f):
	return Rc([Cx(f.list())(w) for w in cplx_prim_roots])

def matrix_cplx_crt(A):
	return A.apply_map(cplx_crt)

def extend_matrix(A):
	return block_matrix([[a.matrix().T for a in r] for r in A.rows()])

def poly_matrix(A):
	return matrix(Rr, [[A[i:i+n, j].list() for j in range(0, A.ncols(), n)] for i in range(0, A.nrows(), n)])

stored_A = matrix(Rq)


T = matrix(Rc)

p = matrix()

h_m = Rq()

if stored_A.dimensions() != (d, m-d):
	print "wrong dimensions for A"
	exit

if T.dimensions() != (2*d, d*K):
	print "wrong dimensions for T"
	exit

A = block_matrix([[1, stored_A]])
TI = block_matrix([[T], [1]])

g = matrix(Rq, [2**i for i in range(K)])
G = block_diagonal_matrix([g for i in range(d)])

def construct_A_m(A, h_m):
	zero_matrix = matrix(Rq, d, 2*d, 0)
	return A + block_matrix([[zero_matrix, h_m * G]])

