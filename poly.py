import bls12_381, fields

def eval(f, x):
  s,curr=0,1
  if isinstance(x, fields.Fq1):
    s=fields.Fq1(0)
    curr=fields.Fq1(1)
  for c in f:
    s+=c*curr
    curr*=x
  return s

def div(num, den):
  quot = []
  divisor = den[-1]
  shiftlen=0
  for x in den:
    if x.x==0:
      shiftlen+=1
    else:
      break
  for i in range(shiftlen + 1):
    mult = num[-1] * divisor.inv()
    quot = [mult] + quot
    if not mult == fields.Fq1(0):
      d = [mult * u for u in den]
      num = [u - v for u, v in zip(num, d)]
    num.pop()
    den.pop(0)
  return quot, num

def mult(a, b):
  res = [fields.Fq1(0)]*(len(a)+len(b)-1)
  for ai,ac in enumerate(a):
    for bi,bc in enumerate(b):
      res[ai+bi] += ac*bc
  return res

def lagrange(ys):
  k = len(ys)-1
  ls = []
  for j in range(k+1):
    ls.append([fields.Fq1(1)])
    for m in range(k+1):
      if m!=j:
        x=(fields.Fq1(j) - fields.Fq1(m)).inv()
        ls[-1] = mult(ls[-1], [fields.Fq1(-m) * x, fields.Fq1(1) * x])
  L = [fields.Fq1(0) for _ in range(k+1)]
  for j in range(k+1):
    L = [prev + ys[j]*l for prev,l in zip(L,ls[j])]
  return L




