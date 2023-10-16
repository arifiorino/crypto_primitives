import kzg, bls12_381, fields, poly, random

def addvec(f,g,h):
  srs=kzg.setup(len(f))
  Cf, rf = kzg.commit(srs,f)
  Cg, rg = kzg.commit(srs,g)
  Ch, rh = kzg.commit(srs,h)
  u = random.randint(0, bls12_381.groupOrder - 1)
  vf, pi_f, delta_f = kzg.prover_open(srs, rf, f, u)
  vg, pi_g, delta_g = kzg.prover_open(srs, rg, g, u)
  vh, pi_h, delta_h = kzg.prover_open(srs, rh, h, u)
  print(kzg.verify_open(srs, Cf, rf, pi_f, delta_f, u, vf))
  print(kzg.verify_open(srs, Cg, rg, pi_g, delta_g, u, vg))
  print(kzg.verify_open(srs, Ch, rh, pi_h, delta_h, u, vh))
  print((vf+vg) % bls12_381.groupOrder == vh % bls12_381.groupOrder)


def test():
  # random vector of length 100:
  a = [random.randint(0,bls12_381.groupOrder-1) for _ in range(100)]
  b = [random.randint(0,bls12_381.groupOrder-1) for _ in range(100)]
  c = [(x+y) % bls12_381.groupOrder for x,y in zip(a,b)]
  addvec(a,b,c)

