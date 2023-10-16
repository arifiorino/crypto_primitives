import bls12_381, fields

def eval(f, x):
  s=0
  curr=1
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



def lagrange(points):
  pass
