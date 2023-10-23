import kzg, bls12_381, fields, poly, random

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
  f=poly.lagrange(a, bls12_381.groupOrder)
  g=poly.lagrange(b, bls12_381.groupOrder)
  h=poly.lagrange(c, bls12_381.groupOrder)
  Z=poly.range_poly(len(a),bls12_381.groupOrder)
  Cf, rf = kzg.commit(srs,f)
  Cg, rg = kzg.commit(srs,g)
  Ch, rh = kzg.commit(srs,h)
  CZ, rZ = kzg.commit(srs,Z)
  #Prover calculates T:
  temp=[x-y for x,y in zip(poly.mult(f,g,bls12_381.groupOrder),h+[0]*len(h))]
  T=poly.div(temp,Z,bls12_381.groupOrder)[0]
  CT, rT = kzg.commit(srs,T)
  #Verifier picks random point:
  u = random.randint(0, bls12_381.groupOrder - 1)
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
  print((vf*vg-vh) % bls12_381.groupOrder == (vT*vZ) % bls12_381.groupOrder)

def test():
  # random vector of length 100:
  a = [random.randint(0,bls12_381.groupOrder-1) for _ in range(100)]
  b = [random.randint(0,bls12_381.groupOrder-1) for _ in range(100)]
  c = [(x+y) % bls12_381.groupOrder for x,y in zip(a,b)]
  addvec(a,b,c)
  c = [(x*y) % bls12_381.groupOrder for x,y in zip(a,b)]
  multvec(a,b,c)

