#q = 1073740609
#q = 1073739937
q = 1073740609
n = 256
r = 32
#r = 16

GFq = GF(q)
root = GFq.zeta(2*r)

def flatten_tree(tree):
    l = []
    for list in tree:
        l += list
    return l

def montgomery_list(l):
	return [l_i * 2^32 for l_i in l]

def brv(a):
	return sum([a_i * 2^i for i, a_i in enumerate(reversed(a.bits()))])

def precompute_zetas(root):
	assert root.multiplicative_order() == 2*r
	
	zetas = []
	
	k = r // 2
	while k > 0:
		zetas = [root^brv(k+i) for i in range(k)] + zetas
		root = zetas[0]^2
		k = k // 2
	
	zetas = [0] + zetas
	zetas_inv = [-zeta for zeta in reversed(zetas)]
	
	return zetas, zetas_inv

def list_order(l):
	return [l_i.multiplicative_order() if l_i.is_unit() else 0 for l_i in l]

mont = mod(2^32, q)
mont2 = mod(2^64, q)
q_inv = - mod(q, 2^32)^(-1)
barrett_mult = floor(2^32 / q)
barrett_shift = 32
double_barrett_mult = floor(2^34 / q)
double_barrett_shift = 34

print "===== common.h =====\n"
print "#define MONT", mont
print "#define MONT2", mont2
print "#define QINV", q_inv
print "#define BARRETT_MULT", barrett_mult
print "#define BARRETT_SHIFT", barrett_shift
print "#define DOUBLE_BARRETT_MULT", double_barrett_mult
print "#define DOUBLE_BARRETT_SHIFT", double_barrett_shift
print "\n"

zetas, zetas_inv = precompute_zetas(root)
montgomery_zetas = montgomery_list(zetas)
montgomery_zetas_inv = montgomery_list(zetas_inv)

m_z_str = "{" + ", ".join(map(str, montgomery_zetas)) + "}"
m_z_i_str = "{" + ", ".join(map(str, montgomery_zetas_inv)) + "}"

print "===== arithmetic.c =====\n"
print "static const scalar zetas[PARAM_R] = " + m_z_str + ";\n"
print "static const scalar zetas_inv[PARAM_R] = " + m_z_i_str + ";\n"
