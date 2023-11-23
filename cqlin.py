import kzg, bls12_381, fields, poly, util, random, ec, time

n=16
p=bls12_381.groupOrder
alpha4=util.nth_root(n*4,p)
alpha2=pow(alpha4,2,p)
alpha=pow(alpha2,2,p)
x=random.randint(0,p-1)

def batched_kzg():
  srs=[bls12_381.G1 * pow(x,i,p) for i in range(n)]
  Pi=list(range(n))
  P=util.inv_dft(Pi,alpha,p)

  h = util.toeplitz_mult(P[1:]+[0]*n,srs[::-1],alpha2,p)
  print(*util.dft(h,alpha,p))

  Kis = [poly.synth_div([P[0]-Pi[i]]+P[1:],pow(alpha,i,p),p)[0] for i in range(n)]
  print(*[bls12_381.G1 * poly.eval(Ki,x,p) for Ki in Kis])

def batched_kzg_An():
  A=list(range(n))
  srs=[bls12_381.G1 * pow(x,n*i,p) * poly.eval(A,x,p) for i in range(n)]
  P=list(range(n))

  h = util.toeplitz_mult(P[1:]+[0]*n,srs[::-1],alpha2,p)
  print(*util.dft(h,alpha,p))

  P_n = []
  for i in range(n):
    P_n+=[P[i]]+[0]*(n-1)
  Kis = [poly.div([P_n[0]-poly.eval(P,pow(alpha,i,p),p)]+P_n[1:],[p-pow(alpha,i,p)]+[0]*(n-1)+[1],p)[0] for i in range(n)]
  print(*[bls12_381.G1 * poly.eval(Ki,x,p) * poly.eval(A,x,p) for Ki in Kis])

def lemma_5_2():
  k=4
  t=time.time()
  srs=[bls12_381.G1 * pow(x,i,p) for i in range(2*n)]
  As = [list(range(i,n+i)) for i in range(k)]

  temp = [util.toeplitz_mult([0]*n + As[j][::-1] + [0]*(2*n-1),srs,alpha4,p)[:n] for j in range(k)]
  srs_star = [util.dft(l[::-1] + [ec.Point(None,None)]*n, alpha2, p) for l in temp]

  Ps = [list(range(i,n+i)) for i in range(k)]
  step_1 = [util.dft(Ps[i]+[0]*n, alpha2, p) for i in range(k)]
  step_2 = [[a*b for a,b in zip(srs_star[i], step_1[i])] for i in range(k)]
  step_3 = [sum([step_2[j][i] for j in range(k)], start=ec.Point(None,None)) for i in range(2*n)]
  h_D_coeff = util.inv_dft(step_3, alpha2, p)
  h_D_a = util.dft(h_D_coeff[n:], alpha, p)

  pi_a = [sum(poly.eval(As[j],x,p) * (poly.eval(Ps[j],x,p) - poly.eval(Ps[j],pow(alpha,i,p),p)) * util.mult_inv(x-pow(alpha,i,p),p)
              for j in range(k)) %p for i in range(n)]
  pi_a = [bls12_381.G1 * pi_a[i] for i in range(n)]
  print(*h_D_a)
  print(*pi_a)



