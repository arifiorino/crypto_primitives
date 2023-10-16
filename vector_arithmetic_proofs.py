import kzg

f = [1,2,3,4,5,6,7]
g = [2,4,2,4,4,3,2]
h = [3,6,5,8,9,9,9]

srs=kzg.setup(len(f))

Cf, rf = kzg.commit(srs,f)
Cg, rg = kzg.commit(srs,g)
Ch, rh = kzg.commit(srs,h)

u = 567 

vf, pi_f, delta_f = kzg.prover_open(srs, rf, f, u)
vg, pi_g, delta_g = kzg.prover_open(srs, rg, g, u)
vh, pi_h, delta_h = kzg.prover_open(srs, rh, h, u)

print(kzg.verify_open(srs, Cf, rf, pi_f, delta_f, u, vf))
print(kzg.verify_open(srs, Cg, rg, pi_g, delta_g, u, vg))
print(kzg.verify_open(srs, Ch, rh, pi_h, delta_h, u, vh))
print(vf+vg == vh)

