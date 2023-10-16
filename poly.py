import util

def eval(f, x, p):
  s,curr=0,1
  for c in f:
    s+=c*curr
    curr*=x
  return s%p

def div(num, den, p):
  quot = []
  divisor = den[-1]
  shiftlen=0
  for x in den:
    if x==0:
      shiftlen+=1
    else:
      break
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
  res = [0]*(len(a)+len(b)-1)
  for ai,ac in enumerate(a):
    for bi,bc in enumerate(b):
      res[ai+bi] += ac*bc
  res = [x%p for x in res]
  return res

def lagrange(ys, p):
  k = len(ys)-1
  ls = []
  for j in range(k+1):
    ls.append([1])
    for m in range(k+1):
      if m!=j:
        x=util.mult_inv((j-m)%p, p)
        ls[-1] = mult(ls[-1], [-m * x, x], p)
  L = [0 for _ in range(k+1)]
  for j in range(k+1):
    L = [(prev+ys[j]*l)%p for prev,l in zip(L,ls[j])]
  return L




