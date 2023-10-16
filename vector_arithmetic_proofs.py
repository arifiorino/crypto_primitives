import kzg, bls12_381, fields, poly, random

def addvec(f,g,h):
  srs=kzg.setup(len(f))
  Cf, rf = kzg.commit(srs,f)
  Cg, rg = kzg.commit(srs,g)
  Ch, rh = kzg.commit(srs,h)
  u = 38141234123
  vf, pi_f, delta_f = kzg.prover_open(srs, rf, f, u)
  vg, pi_g, delta_g = kzg.prover_open(srs, rg, g, u)
  vh, pi_h, delta_h = kzg.prover_open(srs, rh, h, u)
  print(kzg.verify_open(srs, Cf, rf, pi_f, delta_f, u, vf))
  print(kzg.verify_open(srs, Cg, rg, pi_g, delta_g, u, vg))
  print(kzg.verify_open(srs, Ch, rh, pi_h, delta_h, u, vh))
  print(vf+vg == vh)

def mulvec(a,b,c):
  a,b,c=[fields.Fq1(x) for x in a],[fields.Fq1(x) for x in b],[fields.Fq1(x) for x in c]
  f,g,h=poly.lagrange(a),poly.lagrange(b),poly.lagrange(c)
  f,g,h=[x.x for x in f], [x.x for x in g],[x.x for x in h]
  srs=kzg.setup(len(f))
  Cf, rf = kzg.commit(srs,f)
  Cg, rg = kzg.commit(srs,g)
  Ch, rh = kzg.commit(srs,h)
  u = 12362452
  vf, pi_f, delta_f = kzg.prover_open(srs, rf, f, u)
  vg, pi_g, delta_g = kzg.prover_open(srs, rg, g, u)
  vh, pi_h, delta_h = kzg.prover_open(srs, rh, h, u)
  print(kzg.verify_open(srs, Cf, rf, pi_f, delta_f, u, vf))
  print(kzg.verify_open(srs, Cg, rg, pi_g, delta_g, u, vg))
  print(kzg.verify_open(srs, Ch, rh, pi_h, delta_h, u, vh))
  print(vf*vg == vh)

def mulvecscaled(a,b,c,s):
  a,b,c=[fields.Fq1(x) for x in a],[fields.Fq1(x) for x in b],[fields.Fq1(x) for x in c]
  f,g,h=poly.lagrange(a),poly.lagrange(b),poly.lagrange(c)
  f,g,h=[x.x for x in f], [x.x for x in g],[x.x for x in h]
  srs=kzg.setup(len(f))
  Cf, rf = kzg.commit(srs,f)
  Cg, rg = kzg.commit(srs,g)
  Ch, rh = kzg.commit(srs,h)
  u = 12362452
  vf, pi_f, delta_f = kzg.prover_open(srs, rf, f, u)
  vg, pi_g, delta_g = kzg.prover_open(srs, rg, g, u)
  vh, pi_h, delta_h = kzg.prover_open(srs, rh, h, u)
  print(kzg.verify_open(srs, Cf, rf, pi_f, delta_f, u, vf))
  print(kzg.verify_open(srs, Cg, rg, pi_g, delta_g, u, vg))
  print(kzg.verify_open(srs, Ch, rh, pi_h, delta_h, u, vh))
  print(vf*vg//s == vh)
  

addvec([1,2,3,4,2,6,7],
       [2,4,2,2,4,1,1],
       [3,6,5,6,6,7,8])
mulvec([1,2,3,4,5,6,7],
       [1,2,3,4,5,6,7],
       [1,4,9,16,25,36,49])
mulvecscaled([2,4,6,8,10,12,14],
             [1,2,3,4,5,6,7],
             [1,4,9,16,25,36,49],
             2)



