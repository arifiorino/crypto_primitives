import kzg, bls12_381, fields, poly, util, random

n=16
p=257
alpha=255
alpha2=240

def calc_Li_x():
  L=[[0]*n for _ in range(n)]
  for j in range(n):
    L[j][j]=1
    L[j]=util.inv_dft(L[j],alpha,p)
  x=69
  print([poly.eval(L[i],x,p) for i in range(n)])
  P=[pow(x,i,p) for i in range(n)]
  print(util.inv_dft(P,alpha,p))

def batched_kzg():
  x=69
  Ti=list(range(1,17))
  T=util.inv_dft(Ti,alpha,p)
  setup=[pow(x,i,p) for i in range(16)]
  h = util.toeplitz_mult(T[1:]+[0]*16,setup[::-1],alpha2,p)
  print(util.dft(h,alpha,p))

  Kis = [poly.div([T[0]-Ti[i]]+T[1:],[p-pow(alpha,i,p),1],p)[0] for i in range(16)]
  print([poly.eval(Ki,x,p) for Ki in Kis])


def cached_quotients():
  x=69
  L=[[0]*n for _ in range(n)]
  for j in range(n):
    L[j][j]=1
    L[j]=util.inv_dft(L[j],alpha,p)
  Ti=list(range(1,17))
  T=util.inv_dft(Ti,alpha,p)
  setup=[pow(x,i,p) for i in range(16)]
  Ki_x = util.dft(util.toeplitz_mult(T[1:]+[0]*16,setup[::-1],alpha2,p),alpha,p)
  Z_V = [p-1]+[0]*(n-1)+[1]
  Z2_V = [0]*(n-1)+[n]
  Qi_x = [Ki_x[i] * util.mult_inv(poly.eval(Z2_V,pow(alpha,i,p),p),p) % p for i in range(n)]
  for i in range(n):
    print(poly.eval(L[i],x,p) * poly.eval(T,x,p) %p,end=' ')
    print((Ti[i]*poly.eval(L[i],x,p) + poly.eval(Z_V,x,p)*Qi_x[i]) %p)

def use_quotients():
  x=68
  L=[[0]*n for _ in range(n)]
  for j in range(n):
    L[j][j]=1
    L[j]=util.inv_dft(L[j],alpha,p)
  Ti=list(range(1,17))
  T=util.inv_dft(Ti,alpha,p)
  setup=[pow(x,i,p) for i in range(16)]
  Ki_x = util.dft(util.toeplitz_mult(T[1:]+[0]*16,setup[::-1],alpha2,p),alpha,p)
  Z_V = [p-1]+[0]*(n-1)+[1]
  Z2_V = [0]*(n-1)+[n]
  Qi_x = [Ki_x[i] * util.mult_inv(poly.eval(Z2_V,pow(alpha,i,p),p),p) % p for i in range(n)]
  Ai=[1]*4+[0]*12
  A=util.inv_dft(Ai,alpha,p)
  Q_x = sum(Ai[i] * Qi_x[i] for i in range(n))
  R_x = sum(Ai[i] * Ti[i] * poly.eval(L[i],x,p)for i in range(n))
  print(poly.eval(A,x,p) * poly.eval(T,x,p) %p, end=' ')
  print((Q_x * poly.eval(Z_V,x,p) + R_x) %p)

