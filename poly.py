import util

#All polynomials are increasing

def eval(f, x, p):
  s,curr=0,1
  for c in f:
    s+=c*curr
    curr*=x
  return s%p

def div(num, den, p):
  quot = []
  divisor = den[-1]
  shiftlen=len(num)-len(den)
  den=[0]*shiftlen+den
  for i in range(shiftlen + 1):
    mult = num[-1] * util.mult_inv(divisor,p)
    quot = [mult] + quot
    if not mult == 0:
      d = [mult * u for u in den]
      num = [u - v for u, v in zip(num, d)]
    num.pop()
    den.pop(0)
  quot = [x%p for x in quot]
  num = [x%p for x in num]
  return quot, num

def mult(a, b, p):
  n=max(len(a),len(b))*2
  alpha=util.nth_root(n,p)
  a+=[0]*(n-len(a))
  b+=[0]*(n-len(b))
  a=util.dft(a,alpha,p)
  b=util.dft(b,alpha,p)
  c=[x*y for x,y in zip(a,b)]
  c=util.inv_dft(c,alpha,p)
  return c


