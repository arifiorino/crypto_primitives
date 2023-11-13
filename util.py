import random

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
def dft(x, alpha, p):
  n=len(x)
  assert n & (n-1) == 0
  return [sum(x[j]*pow(alpha,j*k,p) % p for j in range(n))%p for k in range(n)]

# x is eval f(alpha^k)
# return coefficients of f
def inv_dft(x, alpha, p):
  n=len(x)
  assert n & (n-1) == 0
  return [mult_inv(n,p) * sum(x[j]*mult_inv(pow(alpha,j*k,p),p) % p for j in range(n))%p for k in range(n)]

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
  r=circulant_mult(c,x+[0]*n,alpha,p)
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

test_circulant_mult()
test_toeplitz_mult()
