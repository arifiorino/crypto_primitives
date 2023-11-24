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
  t=time.time()
  srs=[bls12_381.G1 * pow(x,i,p) for i in range(n*2)]
  Ls = [util.inv_dft([0]*i+[1]+[0]*(n-1-i),alpha,p) for i in range(n)]
  Ps = [list(range(i,n+i)) for i in range(n)]

  srs_star = [util.dft(srs[i:i+n],alpha,p) for i in range(n)]
  srs_star = [[srs_star[i][0]] + srs_star[i][1:][::-1] for i in range(n)]
  srs_star = [[srs_star[i][j] * util.mult_inv(n,p) for i in range(n)] for j in range(n)]
  srs_star = [util.dft(srs_star[i][::-1] + [ec.Point(None,None)]*n, alpha2, p) for i in range(n)]

  step_1 = [util.dft(Ps[i]+[0]*n, alpha2, p) for i in range(n)]
  step_2 = [[a*b for a,b in zip(srs_star[i], step_1[i])] for i in range(n)]
  step_3 = [sum([step_2[j][i] for j in range(n)], start=ec.Point(None,None)) for i in range(2*n)]
  h_D_coeff = util.inv_dft(step_3, alpha2, p)
  pi_a = util.dft(h_D_coeff[n:], alpha, p)

  test = [sum(poly.eval(Ls[j],x,p) * (poly.eval(Ps[j],x,p) - poly.eval(Ps[j],pow(alpha,i,p),p)) * util.mult_inv((x-pow(alpha,i,p))%p,p)
              for j in range(n)) %p for i in range(n)]
  test = [bls12_381.G1 * test[i] for i in range(n)]
  for a,b in zip(pi_a, test):
    print(a==b)

def lemma_5_2_moreover():
  srs=[bls12_381.G1 * pow(x,i,p) for i in range(n**2)]
  Ls = [util.inv_dft([0]*i+[1]+[0]*(n-1-i),alpha,p) for i in range(n)]
  Ps = [list(range(i,n+i)) for i in range(n)]

  srs_star = [util.dft(srs[i*n:i*n+n],alpha,p) for i in range(n)]
  srs_star = [[srs_star[i][0]] + srs_star[i][1:][::-1] for i in range(n)]
  srs_star = [[srs_star[i][j] * util.mult_inv(n,p) for i in range(n)] for j in range(n)]
  srs_star = [util.dft(srs_star[i][::-1] + [ec.Point(None,None)]*n, alpha2, p) for i in range(n)]

  step_1 = [util.dft(Ps[i]+[0]*n, alpha2, p) for i in range(n)]
  step_2 = [[a*b for a,b in zip(srs_star[i], step_1[i])] for i in range(n)]
  step_3 = [sum([step_2[j][i] for j in range(n)], start=ec.Point(None,None)) for i in range(2*n)]
  h_D_coeff = util.inv_dft(step_3, alpha2, p)
  pi_a_prime = util.dft(h_D_coeff[n:], alpha, p)

  test = [sum(poly.eval(Ls[j],x,p) * (poly.eval(Ps[j],pow(x,n,p),p) - poly.eval(Ps[j],pow(alpha,i,p),p)) * util.mult_inv((pow(x,n,p)-pow(alpha,i,p))%p,p)
              for j in range(n)) %p for i in range(n)]
  test = [bls12_381.G1 * test[i] for i in range(n)]
  for a,b in zip(pi_a_prime, test):
    print(a==b)


#x^0,x^1,x^2,... ->DFT->  ->reverse/multiply->  L_0(x),L_1(x),L_2(x),...
#x^i,x^{i+1},x^{i+2},... ->DFT->  ->reverse/multiply->  L_0(x)x^i,L_1(x)x^i,L_2(x)x^i,...
#x^{ni},x^{ni+1},x^{ni+2},... ->DFT->  ->reverse/multiply->  L_0(x)x^{ni},L_1(x)x^{ni},L_2(x)x^{ni},...

