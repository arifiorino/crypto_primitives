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

#t=time.time()
#print(time.time()-t);t=time.time()

srs=[bls12_381.G1 * pow(x,i,p) for i in range(n**2 + n)]

Ls = [util.inv_dft([0]*i+[1]+[0]*(n-1-i),alpha,p) for i in range(n)]
M = [list(range(i,i+n)) for i in range(n)]

temp = [convert_to_L([srs[i + n*j] for j in range(n)]) for i in range(n)] #L_j(x^n)x^i
U = [convert_to_L([temp[i][j] for i in range(n)]) for j in range(n)] #L_i(x^n)L_j(x)
temp = convert_to_L(srs[n**2-n : n**2])
V = [temp[i] * util.mult_inv(n,p) for i in range(n)]

srs_star = [convert_to_L(srs[n*i:n*i+n]) for i in range(n)]
srs_star = [[srs_star[j][i] for j in range(n)] for i in range(n)]
srs_star = [util.dft(srs_star[i][::-1] + [ec.Point(None,None)]*n, alpha2, p) for i in range(n)]

S = [sum(((U[i][j] * util.mult_inv(pow(alpha,i,p),p) - V[j]) * M[i][j] for j in range(n)), start=ec.Point(None,None)) for i in range(n)]
R = [sum((U[i][j] * M[i][j] for j in range(n)), start=ec.Point(None,None)) for i in range(n)]
C = [[sum([M[i][j] * Ls[i][coeff] for i in range(n)]) %p for coeff in range(n)] for j in range(n)]

temp = [util.dft(C[i]+[0]*n, alpha2, p) for i in range(n)] #1
temp = [[srs_star[i][j]*temp[i][j] for j in range(2*n)] for i in range(n)] #2
temp = [sum([temp[j][i] for j in range(n)], start=ec.Point(None,None)) for i in range(2*n)] #3
temp = util.inv_dft(temp, alpha2, p) #coefficients
temp = util.dft(temp[n:], alpha, p) #evaluate
Q = [temp[i]*(pow(alpha,i,p)*util.mult_inv(n,p)%p) for i in range(n)] #scale by c_i

L_x = [poly.eval(Ls[i],x,p) for i in range(n)]
L_x_n = [poly.eval(Ls[i],pow(x,n,p),p) for i in range(n)]
R_x = [sum([M[i][j] * L_x[j] %p for j in range(n)]) %p for i in range(n)]
M_x = bls12_381.G1 * (sum(L_x_n[i] * R_x[i] for i in range(n))%p)
Z_x = pow(x,n*n,p) - 1
#Test Q,R
for i in range(n):
  print(M_x * L_x_n[i] == Q[i] * Z_x + R[i])

for i in range(n):
  print(bls12_381.G1 * ((L_x_n[i] - util.mult_inv(n,p)) * R_x[i] %p) == S[i] * pow(x,n,p))





