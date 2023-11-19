import random, ec

def mult_inv(x,p):
  t,newt=0,1
  r,newr=p,x
  while newr!=0:
    quotient=r//newr
    t,newt=newt,t-quotient*newt
    r,newr=newr,r-quotient*newr
  if r>1:
    raise Exception("not invertible")
  if t<0:
    t+=p
  return t

def nth_root(n, p):
  assert (p-1) % n == 0
  assert n & (n-1) == 0
  while 1:
    x=random.randint(1,p-1)
    g=pow(x,(p-1)//n,p)
    if pow(g,n//2,p) != 1:
      return g

# x is coefficients of f
# eval f(alpha^k)
def dft(y, alpha, prime):
  def dft2(x,N,s):
    if N==1:
      return [x[0]]
    else:
      X= dft2(x, N//2, s*2) + dft2(x[s:], N//2, s*2)
      for k in range(N//2):
        p=X[k]
        q=X[k+N//2] * pow(alpha,s*k,prime)
        X[k]=(p+q)%prime
        X[k+N//2]=(p-q)%prime
    return X
  return dft2(y,len(y),1)

# x is eval f(alpha^k)
# return coefficients of f
def inv_dft(y, alpha, prime):
  def inv_dft2(x,N,s):
    if N==1:
      return [x[0]*N]
    else:
      X= inv_dft2(x, N//2, s*2) + inv_dft2(x[s:], N//2, s*2)
      for k in range(N//2):
        p=X[k]
        q=X[k+N//2] * mult_inv(pow(alpha,s*k,prime),prime)
        X[k]=(p+q)%prime
        X[k+N//2]=(p-q)%prime
    return X
  return [c * mult_inv(len(y),prime) % prime for c in inv_dft2(y,len(y),1)]

matmul=lambda A,x: [sum([x1*x2 for x1,x2 in zip(A[i],x)]) for i in range(len(A))]

def circulant_mult(c,x,alpha,p):
  Lambda=dft(c,alpha,p)
  r=dft(x,alpha,p)
  r=[a*b for a,b in zip(r,Lambda)]
  r=inv_dft(r,alpha,p)
  return r

def test_circulant_mult():
  n=4
  p=257
  alpha=nth_root(n,p)
  c=[2,3,0,1]
  C=[[2,1,0,3],[3,2,1,0],[0,3,2,1],[1,0,3,2]]
  x=[9,2,0,0]
  print(matmul(C,x))
  print(circulant_mult(c,x,alpha,p))

def toeplitz_mult(a,x,alpha,p):
  n=(len(a)+1)//2
  c=a[n-1:]+[0]+a[:n-1]
  if isinstance(x[0],int):
    r=circulant_mult(c,x+[0]*n,alpha,p)
  else:
    r=circulant_mult(c,x+[ec.Point(None,None)]*n,alpha,p)
  return r[:n]

def test_toeplitz_mult():
  n=4
  p=257
  alpha=nth_root(n,p)
  a=[1,2,3]
  T=[[2,1],[3,2]]
  x=[9,2]
  print(matmul(T,x))
  print(toeplitz_mult(a,x,alpha,p))

