import bls12_381, fields, ec, random, poly

def setup(Nmax):
  tau = random.randint(0,bls12_381.groupOrder-1)
  xi = random.randint(0,bls12_381.groupOrder-1)
  return [bls12_381.G1*x for x in [tau**i for i in range(Nmax)]+[xi]] +\
         [bls12_381.G2*x for x in [1,tau,xi]]

def commit(srs, f):
  r = random.randint(0,bls12_381.groupOrder-1)
  C = srs[-4] * r
  for i in range(len(f)):
    C += srs[i]*f[i]
  return C, r

def prover_open(srs, r, f, u):
  v = poly.eval(f, u)
  s = random.randint(0,bls12_381.groupOrder-1)
  f = [fields.Fq1(c) for c in f]
  f[0] -= fields.Fq1(v)
  g = [fields.Fq1(0) for _ in range(len(f)-2)] + [fields.Fq1(-u), fields.Fq1(1)]
  q, rem = poly.div(f, g)
  assert len(rem)==1 and rem[0].x == 0
  pi = srs[-4] * s
  for i in range(len(q)):
    pi += srs[i]*q[i].x
  delta = srs[0]*r - srs[1]*s + srs[0]*(s*u)
  return (v, pi, delta)

def verify_open(srs, C, r, pi, delta, u, v):
  lhs = bls12_381.pairing(C - srs[0] * v, srs[-3])
  rhs = bls12_381.pairing(pi, srs[-2] - srs[-3] * u) * bls12_381.pairing(delta, srs[-1])
  return lhs==rhs

def test():
  f = [1,2,3,4,5,6,7]
  srs=setup(len(f))
  C, r = commit(srs,f)
  u = 567 
  v, pi, delta = prover_open(srs, r, f, u)
  assert(verify_open(srs, C, r, pi, delta, u, v))

