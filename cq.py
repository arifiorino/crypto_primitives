import kzg, bls12_381, fields, poly, util, random

n=16
p=257
alpha=255
alpha2=240

def calc_Li_x():
  L=[[0]*n for _ in range(n)]
  Ln=[[0]*n*n for _ in range(n)]
  for j in range(n):
    L[j][j]=1
    L[j]=util.inv_dft(L[j],alpha,p)
  x=69
  print([poly.eval(L[i],x,p) for i in range(n)])
  P=[pow(x,i,p) for i in range(n)]
  print(util.inv_dft(P,alpha,p))


x=69
Ti=list(range(1,17))
T=util.inv_dft(Ti,alpha,p)
setup=[pow(x,i,p) for i in range(16)]
h = util.toeplitz_mult(T+[0]*15,setup[::-1],alpha2,p)
print(util.dft(h,alpha,p))

print([(poly.eval(T,x,p)-Ti[i])*util.mult_inv((x-pow(alpha,i,p))%p,p)%p for i in range(16)])

Kis = [poly.div([T[0]-Ti[i]]+T[1:], [p-pow(alpha,i,p), 1], p)[0] for i in range(16)]
print([poly.eval(Ki,x,p) for Ki in Kis])




