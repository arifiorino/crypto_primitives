import bls12_381, fields, ec, random, poly

def setup(Nmax):
  tau = random.randint(0,bls12_381.groupOrder-1)
  xi = random.randint(0,bls12_381.groupOrder-1)
  srs=[bls12_381.G1]
  for i in range(Nmax-1):
    srs.append(srs[-1] * tau)
  return srs + [bls12_381.G1*xi] + [bls12_381.G2*x for x in [1,tau,xi]]

def commit(srs, f):
  r = random.randint(0,bls12_381.groupOrder-1)
  C = srs[-4] * r
  for i in range(len(f)):
    C += srs[i]*f[i]
  return C, r

def prover_open(srs, r, f, u):
  v = poly.eval(f, u, bls12_381.groupOrder)
  s = random.randint(0,bls12_381.groupOrder-1)
  f[0] = (f[0]-v) % bls12_381.groupOrder
  g = [-u % bls12_381.groupOrder, 1]
  q, rem = poly.div(f, g, bls12_381.groupOrder)
  assert len(rem)==0 or (len(rem)==1 and rem[0] == 0)
  pi = srs[-4] * s
  for i in range(len(q)):
    pi += srs[i]*q[i]
  delta = srs[0]*r - srs[1]*s + srs[0]*(s*u)
  return (v, pi, delta)

def verify_open(srs, C, r, pi, delta, u, v):
  lhs = bls12_381.pairing(C - srs[0] * v, srs[-3])
  rhs = bls12_381.pairing(pi, srs[-2] - srs[-3] * u) * bls12_381.pairing(delta, srs[-1])
  return lhs==rhs

def test():
  #random polynomial of degree 128:
  f = [random.randint(0,bls12_381.groupOrder-1) for _ in range(128)]
  srs=setup(len(f))
  C, r = commit(srs,f)
  u = random.randint(0,bls12_381.groupOrder-1)
  v, pi, delta = prover_open(srs, r, f, u)
  assert(verify_open(srs, C, r, pi, delta, u, v))

