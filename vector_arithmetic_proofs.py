import kzg, bls12_381, fields, poly, util, random

def addvec(f,g,h):
  #Setup
  srs=kzg.setup(len(f))
  Cf, rf = kzg.commit(srs,f)
  Cg, rg = kzg.commit(srs,g)
  Ch, rh = kzg.commit(srs,h)
  #Verifier picks random point:
  u = random.randint(0, bls12_381.groupOrder - 1)
  #Prover opens f,g,h at u:
  vf, pi_f, delta_f = kzg.prover_open(srs, rf, f, u)
  vg, pi_g, delta_g = kzg.prover_open(srs, rg, g, u)
  vh, pi_h, delta_h = kzg.prover_open(srs, rh, h, u)
  #Verifier verifies opens:
  print(kzg.verify_open(srs, Cf, rf, pi_f, delta_f, u, vf))
  print(kzg.verify_open(srs, Cg, rg, pi_g, delta_g, u, vg))
  print(kzg.verify_open(srs, Ch, rh, pi_h, delta_h, u, vh))
  #Verifier verifies f+g=h:
  print((vf+vg) % bls12_381.groupOrder == vh % bls12_381.groupOrder)

def multvec(a,b,c):
  #Setup
  srs=kzg.setup(len(a)*2)
  n,p=len(a),bls12_381.groupOrder
  alpha=util.nth_root(n, p)
  f=util.inv_dft(a, alpha, p)
  g=util.inv_dft(b, alpha, p)
  h=util.inv_dft(c, alpha, p)
  Z=[p-1]+[0]*(n-1)+[1]
  Cf, rf = kzg.commit(srs,f)
  Cg, rg = kzg.commit(srs,g)
  Ch, rh = kzg.commit(srs,h)
  CZ, rZ = kzg.commit(srs,Z)
  #Prover calculates T:
  temp=[(x-y)%p for x,y in zip(poly.mult(f,g,p),h+[0]*len(h))]
  T,_=poly.div(temp,Z,p)
  CT, rT = kzg.commit(srs,T)
  #Verifier picks random point:
  u = random.randint(0, p - 1)
  #Prover opens f,g,h,T,Z at u:
  vf, pi_f, delta_f = kzg.prover_open(srs, rf, f, u)
  vg, pi_g, delta_g = kzg.prover_open(srs, rg, g, u)
  vh, pi_h, delta_h = kzg.prover_open(srs, rh, h, u)
  vZ, pi_Z, delta_Z = kzg.prover_open(srs, rZ, Z, u)
  vT, pi_T, delta_T = kzg.prover_open(srs, rT, T, u)
  #Verifier verifies opens:
  print(kzg.verify_open(srs, Cf, rf, pi_f, delta_f, u, vf))
  print(kzg.verify_open(srs, Cg, rg, pi_g, delta_g, u, vg))
  print(kzg.verify_open(srs, Ch, rh, pi_h, delta_h, u, vh))
  print(kzg.verify_open(srs, CZ, rZ, pi_Z, delta_Z, u, vZ))
  print(kzg.verify_open(srs, CT, rT, pi_T, delta_T, u, vT))
  #Verifier verifies fg-h=TZ:
  print((vf*vg-vh) % p == (vT*vZ) % p)

def test():
  # random vector of length 128:
  a = [random.randint(0,bls12_381.groupOrder-1) for _ in range(128)]
  b = [random.randint(0,bls12_381.groupOrder-1) for _ in range(128)]
  c = [(x+y) % bls12_381.groupOrder for x,y in zip(a,b)]
  addvec(a,b,c)
  c = [(x*y) % bls12_381.groupOrder for x,y in zip(a,b)]
  multvec(a,b,c)

