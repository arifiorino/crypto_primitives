import fields, ec

fields.fieldPrime=0x1A0111EA397FE69A4B1BA7B6434BACD764774B84F38512BF6730D2A0F6B0F6241EABFFFEB153FFFFB9FEFFFFFFFFAAAB
groupOrder=0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001

g1Gen = ec.Point(fields.Fq1(0x17F1D3A73197D7942695638C4FA9AC0FC3688C4F9774B905A14E3A3F171BAC586C55E83FF97A1AEFFB3AF00ADB22C6BB),
                 fields.Fq1(0x08B3F481E3AAA0F1A09E30ED741D8AE4FCF5E095D5D00AF600DB18CB2C04B3EDD03CC744A2888AE40CAA232946C5E7E1))
g2Gen = ec.Point(fields.Fq2(fields.Fq1(0x024AA2B2F08F0A91260805272DC51051C6E47AD4FA403B02B4510B647AE3D1770BAC0326A805BBEFD48056C8C121BDB8),
                            fields.Fq1(0x13E02B6052719F607DACD3A088274F65596BD0D09920B61AB5DA61BBDC7F5049334CF11213945D57E5AC7D055D042B7E)),
                 fields.Fq2(fields.Fq1(0x0CE5D527727D6E118CC9CDC6DA2E351AADFD9BAA8CBDD3A76D429A695160D12C923AC9CC3BACA289E193548608B82801),
                            fields.Fq1(0x0606C4A02EA734CC32ACD2B02BC28B99CB3E287E85A763AF267492AB572E99AB3F370D275CEC1DA1AAA9075FF05F79BE)))

def doubleEval(r, p):
  wideR = r.untwist()
  slope = (wideR.x * wideR.x + wideR.x * wideR.x + wideR.x * wideR.x) * (wideR.y + wideR.y).inv()
  v = wideR.y - slope * wideR.x
  return fields.Fq12.fromInt(p.y.x) - (fields.Fq12.fromInt(p.x.x) * slope) - v

def addEval(r, q, p):
  wideR = r.untwist()
  wideQ = q.untwist()
  if wideR.x == wideQ.x and wideR.y == fields.Fq12.fromInt(0)-wideQ.y:
    return fields.Fq12.fromInt(p.x.x) - wideR.x
  slope = (wideQ.y - wideR.y) * (wideQ.x - wideR.x).inv()
  v = (wideQ.y * wideR.x - wideR.y * wideQ.x) * (wideR.x - wideQ.x).inv()
  return fields.Fq12.fromInt(p.y.x) - (fields.Fq12.fromInt(p.x.x) * slope) - v

def millerLoop(p, q):
  l=[c=='1' for c in bin(0xd201000000010000)[2:]]
  r = q
  result = fields.Fq12.fromInt(1)
  for t in l:
    result *= result * doubleEval(r, p)
    r += r
    if t:
      r += q
      result *= addEval(r+r, q, p)
  return result

def pairing(p, q):
  return millerLoop(p, q) ** ((fields.fieldPrime**12 - 1) // groupOrder)

print(pairing(g1Gen*5, g2Gen*5))
print(pairing(g1Gen*5, g2Gen*3) + pairing(g1Gen*5, g2Gen*2))

