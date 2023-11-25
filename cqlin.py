import kzg, bls12_381, fields, poly, util, random, ec, time

n=8
p=bls12_381.groupOrder
alpha2=util.nth_root(n*2,p)
alpha=pow(alpha2,2,p)
x=random.randint(0,p-1)

#x^0,x^1,x^2,... -> L_0(x),L_1(x),L_2(x),...
def convert_to_L(xs):
  r = util.dft(xs,alpha,p)
  r = [r[0]] + r[1:][::-1]
  r = [r[i] * util.mult_inv(n,p) for i in range(n)]
  return r

t=time.time()

srs=[bls12_381.G1 * pow(x,i,p) for i in range(n**2 + n)]
srs2=[bls12_381.G2 * pow(x,i,p) for i in range(n**2 + n)]

M = [list(range(i,i+n)) for i in range(n)]
a = list(range(n))
b = util.matmul(M,a)
f = util.inv_dft(a, alpha, p)
g = util.inv_dft(b, alpha, p)
f_com = sum((srs[i] * f[i]  for i in range(n)), start=ec.Point(None,None))
g_com = sum((srs[i] * g[i]  for i in range(n)), start=ec.Point(None,None))

temp = [convert_to_L([srs[i + n*j] for j in range(n)]) for i in range(n)] #L_j(x^n)x^i
U = [convert_to_L([temp[i][j] for i in range(n)]) for j in range(n)] #L_i(x^n)L_j(x)
temp = convert_to_L(srs[n**2-n : n**2])
V = [temp[i] * util.mult_inv(n,p) for i in range(n)]

srs_star = [convert_to_L(srs[n*i:n*i+n]) for i in range(n)]
srs_star = [[srs_star[j][i] for j in range(n)] for i in range(n)]
srs_star = [util.dft(srs_star[i][::-1] + [ec.Point(None,None)]*n, alpha2, p) for i in range(n)]

Ls = [util.inv_dft([0]*i+[1]+[0]*(n-1-i),alpha,p) for i in range(n)]
S = [sum(((U[i][j] * util.mult_inv(pow(alpha,i,p),p) - V[j]) * M[i][j] for j in range(n)), start=ec.Point(None,None)) for i in range(n)]
R = [sum((U[i][j] * M[i][j] for j in range(n)), start=ec.Point(None,None)) for i in range(n)]
C = [[sum([M[i][j] * Ls[i][coeff] for i in range(n)]) %p for coeff in range(n)] for j in range(n)]

temp = [util.dft(C[i]+[0]*n, alpha2, p) for i in range(n)] #1
temp = [[srs_star[i][j]*temp[i][j] for j in range(2*n)] for i in range(n)] #2
temp = [sum([temp[j][i] for j in range(n)], start=ec.Point(None,None)) for i in range(2*n)] #3
temp = util.inv_dft(temp, alpha2, p) #coefficients
temp = util.dft(temp[n:], alpha, p) #evaluate
Q = [temp[i]*(pow(alpha,i,p)*util.mult_inv(n,p)%p) for i in range(n)] #scale by c_i

z = srs2[n**2] - srs2[0]
temp = [convert_to_L([srs2[i + n*j] for j in range(n)]) for i in range(n)] #L_j(x^n)x^i
U2 = [convert_to_L([temp[i][j] for i in range(n)]) for j in range(n)] #L_i(x^n)L_j(x)
m = sum((sum((U2[i][j] * M[i][j] for j in range(n)), start=ec.Point(None,None)) for i in range(n)), start=ec.Point(None,None))
print("Setup time:",time.time()-t);t=time.time()

r = sum((R[i] * a[i] for i in range(n)), start=ec.Point(None,None))
q = sum((Q[i] * a[i] for i in range(n)), start=ec.Point(None,None))
A_com = sum((srs[n*i] * f[i]  for i in range(n)), start=ec.Point(None,None))
print(bls12_381.pairing(A_com,m) == bls12_381.pairing(q,z) * bls12_381.pairing(r,srs2[0]))

s = sum((S[i] * a[i] for i in range(n)), start=ec.Point(None,None))
print(bls12_381.pairing(r - g_com * util.mult_inv(n,p), srs2[0]) == bls12_381.pairing(s, srs2[n]))

p_x = sum((srs[n**2-n + i] * g[i] for i in range(n)), start=ec.Point(None,None))
print(bls12_381.pairing(g_com, srs2[n**2-n]) == bls12_381.pairing(p_x, srs2[0]))

gamma = random.randint(0,p-1)
gamma_n = pow(gamma,n,p)
z = poly.eval(f,gamma_n,p)
h = poly.synth_div([(f[0]-z)%p]+f[1:],gamma_n,p)[0]
pi = sum((srs[i] * h[i]  for i in range(n-1)), start=ec.Point(None,None))
print(bls12_381.pairing(f_com - srs[0]*z + pi*gamma_n,srs2[0]) == bls12_381.pairing(pi,srs2[1]))

temp = sum([[(f[0]-z)%p]] + [[0]*(n-1)+[f[i]] for i in range(1,n)], start=[])
h_1 = poly.div(temp, [p-gamma_n]+[0]*(n-1)+[1], p)[0]
pi_1 = sum((srs[i] * h_1[i]  for i in range(n**2-2*n + 1)), start=ec.Point(None,None))
print(bls12_381.pairing(A_com - srs[0]*z + pi_1*gamma_n,srs2[0]) == bls12_381.pairing(pi_1,srs2[n]))
print("Proof time:",time.time()-t);t=time.time()

