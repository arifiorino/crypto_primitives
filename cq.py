import kzg, bls12_381, ec, fields, poly, util, random

N=32
n=4
p=bls12_381.groupOrder
alpha_2N=util.nth_root(N*2,p)
alpha_N=pow(alpha_2N,2,p)
alpha_n=pow(alpha_N,N//n,p)

table=list(range(N))
f_i=list(range(n))
f=util.inv_dft(f_i,alpha_n,p)

# gen(N, t):
x = random.randint(0,p-1) #1
x_pow = [pow(x,i,p) for i in range(N)] #1
srs = [bls12_381.G1 * y for y in x_pow] #1
srs2 = [bls12_381.G2 * y for y in x_pow] + [bls12_381.G2 * pow(x,N,p)] #1
Z_V_x_2 = bls12_381.G2 * (x*x_pow[-1]-1) #2
T = util.inv_dft(table,alpha_N,p) #3
T_x_2 = sum((srs2[i] * T[i] for i in range(N)), start=ec.Point(None,None)) #3
K_i_x = util.dft(util.toeplitz_mult(T[1:]+[0]*N,x_pow[::-1],alpha_2N,p),alpha_N,p) #4a
Q_i_x = [K_i_x[i] * util.mult_inv(N*pow(alpha_N,i*(N-1),p),p) % p for i in range(N)] #4a
Q_i_x_1 = [bls12_381.G1 * y for y in Q_i_x] #4a
L_i_x = util.inv_dft(x_pow,alpha_N,p) #4b
L_i_x_1 = [bls12_381.G1 * y for y in L_i_x] #4b
L_i_0_x_1 = [L_i_x_1[i] * util.mult_inv(pow(alpha_N,i,p),p) - srs[-1] * util.mult_inv(N,p) for i in range(N)] #4c
f_x_1 = sum((srs[i] * f[i] for i in range(n)), start=ec.Point(None,None)) #commitment to f

# IsInTable(f_x_1, table, srs; f):
# Round 1
m_i = [(i,1) for i in range(n)] #1
m_x_1 = sum((L_i_x_1[i] * y for i,y in m_i), start=ec.Point(None,None)) #2
# Round 2
beta = random.randint(0,p-1) #1
A_i = [(i,y * util.mult_inv(table[i]+beta,p)) for i,y in m_i] #2
A_x_1 = sum((L_i_x_1[i] * y for i,y in A_i), start=ec.Point(None,None)) #3
Q_A_x_1 = sum((Q_i_x_1[i] * y for i,y in A_i), start=ec.Point(None,None)) #4
B_i = [util.mult_inv(f_i[i]+beta,p) for i in range(n)] #5
B = util.inv_dft(B_i, alpha_n, p) #5
B_0 = B[1:]+[0] #6
B_0_x_1 = bls12_381.G1*poly.eval(B_0,x,p) #7
Q_B = poly.mult(B, [(f[0]+beta)%p]+f[1:], p) #8
Q_B = poly.div([(Q_B[0]-1)%p]+Q_B[1:], [p-1]+[0]*(n-1)+[1], p)[0] #8
Q_B_x_1 = sum((srs[i] * Q_B[i] for i in range(n)), start=ec.Point(None,None)) #9
P_x_1 = sum((srs[N-1-(n-2)+i] * B_0[i] for i in range(n-1)), start=ec.Point(None,None)) #10
print(bls12_381.pairing(A_x_1,T_x_2) == bls12_381.pairing(Q_A_x_1,Z_V_x_2) * bls12_381.pairing(m_x_1 - A_x_1 * beta, srs2[0])) #11
print(bls12_381.pairing(B_0_x_1,srs2[N-1-(n-2)]) == bls12_381.pairing(P_x_1,srs2[0])) #12
# Round 3
gamma = random.randint(0,p-1) #1
eta = random.randint(0,p-1) #1
B_0_gamma = poly.eval(B_0,gamma,p) #2
f_gamma = poly.eval(f,gamma,p) #2
A_0 = util.mult_inv(N,p) * sum(y for _,y in A_i) #3
b_0 = N * A_0 * util.mult_inv(n,p) % p #4
Z_H_gamma = pow(gamma,n,p)-1 #5
b_gamma = (B_0_gamma * gamma + b_0) %p #5
Q_b_gamma = (b_gamma * (f_gamma + beta) - 1) * util.mult_inv(Z_H_gamma,p) %p #5
v = (B_0_gamma + eta * f_gamma + pow(eta,2,p) * Q_b_gamma) %p #6a
num = [B_0[i] + eta * f[i] + pow(eta,2,p) * Q_B[i] for i in range(n)] #6b
num[0] -= v #6b
num=[y%p for y in num] #6b
h = poly.synth_div(num, gamma, p)[0] #6b
pi_gamma = sum((srs[i] * h[i] for i in range(n-1)), start=ec.Point(None,None)) #6b
c = B_0_x_1 + f_x_1 * eta + Q_B_x_1 * pow(eta,2,p) #6c
print(bls12_381.pairing(c - bls12_381.G1*v + pi_gamma * gamma, srs2[0]) == bls12_381.pairing(pi_gamma, srs2[1])) #6c
A_0_x = sum((L_i_0_x_1[i] * y for i,y in A_i), start=ec.Point(None,None)) #7a
print(bls12_381.pairing(A_x_1 - bls12_381.G1*A_0, srs2[0]) == bls12_381.pairing(A_0_x, srs2[1])) #7b



