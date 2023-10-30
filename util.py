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

