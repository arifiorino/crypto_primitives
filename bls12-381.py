import random

fieldPrime=0x1A0111EA397FE69A4B1BA7B6434BACD764774B84F38512BF6730D2A0F6B0F6241EABFFFEB153FFFFB9FEFFFFFFFFAAAB
groupOrder=0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001

class Fq1():
  def __init__(self,x):
    self.x=x%fieldPrime
  def __add__(self,b):
    return Fq1(self.x+b.x)
  def __sub__(self,b):
    return Fq1(self.x-b.x)
  def __mul__(self,b):
    return Fq1(self.x*b.x)
  def inv(self):
    t,newt=0,1
    r,newr=fieldPrime,self.x
    while newr!=0:
      quotient=r//newr
      t,newt=newt,t-quotient*newt
      r,newr=newr,r-quotient*newr
    if r>1:
      raise Exception("not invertible")
    if t<0:
      t+=fieldPrime
    return Fq1(t)
  def mul_nonres(self):
    return self
  def __str__(self):
    return hex(self.x)
  def __eq__(self,b):
    return self.x==b.x
  @staticmethod
  def fromInt(x):
    return Fq1(x)

class Fq2():
  def __init__(self,x,y):
    self.x=x
    self.y=y
  def __add__(self,b):
    return Fq2(self.x+b.x,self.y+b.y)
  def __sub__(self,b):
    return Fq2(self.x-b.x,self.y-b.y)
  def __mul__(self,b):
    return Fq2(self.x*b.x - self.y*b.y, self.x*b.y + self.y*b.x)
  def inv(self):
    factor=(self.x*self.x + self.y*self.y).inv()
    return Fq2(self.x*factor,(Fq1(0)-self.y)*factor)
  def mul_nonres(self):
    return Fq2(self.x-self.y,self.x+self.y)
  def __str__(self):
    return '%s,%s'%(str(self.x),str(self.y))
  def __eq__(self,b):
    return (self.x,self.y)==(b.x,b.y)
  @staticmethod
  def fromInt(x):
    return Fq2(Fq1.fromInt(x),Fq1.fromInt(0))

class Fq6():
  def __init__(self,x,y,z):
    self.x=x
    self.y=y
    self.z=z
  def __add__(self,b):
    return Fq6(self.x+b.x,self.y+b.y,self.z+b.z)
  def __sub__(self,b):
    return Fq6(self.x-b.x,self.y-b.y,self.z-b.z)
  def __mul__(self,b):
    t0=self.x*b.x
    t1=self.x*b.y + self.y*b.x
    t2=self.x*b.z + self.y*b.y + self.z*b.x
    t3=(self.y*b.z + self.z*b.y).mul_nonres()
    t4=(self.z*b.z).mul_nonres()
    return Fq6(t0+t3, t1+t4, t2)
  def inv(self):
    t0 = self.x*self.x - (self.y*self.z).mul_nonres()
    t1 = (self.z*self.z).mul_nonres() - self.x*self.y
    t2 = self.y*self.y - self.x*self.z
    factor = (self.x*t0 + (self.z*t1).mul_nonres() + (self.y*t2).mul_nonres()).inv()
    return Fq6(t0*factor, t1*factor, t2*factor)
  def mul_nonres(self):
    return Fq6(self.z.mul_nonres(), self.x, self.y)
  def __str__(self):
    return '%s,%s,%s'%(str(self.x),str(self.y),str(self.z))
  def __eq__(self,b):
    return (self.x,self.y,self.z)==(b.x,b.y,b.z)
  @staticmethod
  def fromInt(x):
    return Fq6(Fq2.fromInt(x),Fq2.fromInt(0),Fq2.fromInt(0))

class Fq12():
  def __init__(self,x,y):
    self.x=x
    self.y=y
  def __add__(self,b):
    return Fq12(self.x+b.x,self.y+b.y)
  def __sub__(self,b):
    return Fq12(self.x-b.x,self.y-b.y)
  def __mul__(self,b):
    return Fq12(self.x*b.x + (self.y*b.y).mul_nonres(), self.y*b.x + self.x*b.y)
  def __pow__(self,k):
    res = self
    temp = self
    k-=1
    while k>0:
      if k&1:
        res*=temp
      temp*=temp
      k>>=1
    return res
  def inv(self):
    factor=(self.x*self.x - (self.y*self.y).mul_nonres()).inv()
    zero=Fq6(Fq2(Fq1(0),Fq1(0)),Fq2(Fq1(0),Fq1(0)),Fq2(Fq1(0),Fq1(0)))
    return Fq12(self.x*factor, (zero-self.y)*factor)
  def __str__(self):
    return '%s,%s'%(str(self.x),str(self.y))
  def __eq__(self,b):
    return (self.x,self.y)==(b.x,b.y)
  @staticmethod
  def fromInt(x):
    return Fq12(Fq6.fromInt(x),Fq6.fromInt(0))

class Point():
  def __init__(self,x,y):
    self.x,self.y=x,y
  def __add__(self,b):
    if not self.x:
      return b
    if not b.x:
      return self
    if self==b:
      slope=(b.x*b.x + b.x*b.x + b.x*b.x) * (b.y + b.y).inv() #tangent
    else:
      if self.x==b.x:
        return Point(None,None)
      slope = (b.y-self.y) * (b.x-self.x).inv() #line
    x=slope*slope - self.x - b.x
    return Point(x, slope * (self.x - x) - self.y)
  def __mul__(self,k):
    res = Point(None,None)
    temp = self
    while k>0:
      if k&1:
        res+=temp
      temp+=temp
      k>>=1
    return res
  def check(self):
    t=Fq1(4)
    if isinstance(self.x,Fq2):
      t=Fq2(Fq1(4),Fq1(4))
    return self.y*self.y == self.x*self.x*self.x + t
  def __str__(self):
    return '(%s),(%s)'%(str(self.x),str(self.y))
  def __eq__(self,b):
    return (self.x,self.y)==(b.x,b.y)
  def untwist(self):
    assert isinstance(self.x,Fq2)
    root=Fq6(Fq2.fromInt(0), Fq2.fromInt(1), Fq2.fromInt(0))
    return Point(Fq12(Fq6(self.x, Fq2.fromInt(0), Fq2.fromInt(0)), Fq6.fromInt(0)) * Fq12(root, Fq6.fromInt(0)).inv(),
                 Fq12(Fq6(self.y, Fq2.fromInt(0), Fq2.fromInt(0)), Fq6.fromInt(0)) * Fq12(Fq6.fromInt(0), root).inv())
    


g1Gen = Point(Fq1(0x17F1D3A73197D7942695638C4FA9AC0FC3688C4F9774B905A14E3A3F171BAC586C55E83FF97A1AEFFB3AF00ADB22C6BB),
              Fq1(0x08B3F481E3AAA0F1A09E30ED741D8AE4FCF5E095D5D00AF600DB18CB2C04B3EDD03CC744A2888AE40CAA232946C5E7E1))
g2Gen = Point(Fq2(Fq1(0x024AA2B2F08F0A91260805272DC51051C6E47AD4FA403B02B4510B647AE3D1770BAC0326A805BBEFD48056C8C121BDB8),
                  Fq1(0x13E02B6052719F607DACD3A088274F65596BD0D09920B61AB5DA61BBDC7F5049334CF11213945D57E5AC7D055D042B7E)),
              Fq2(Fq1(0x0CE5D527727D6E118CC9CDC6DA2E351AADFD9BAA8CBDD3A76D429A695160D12C923AC9CC3BACA289E193548608B82801),
                  Fq1(0x0606C4A02EA734CC32ACD2B02BC28B99CB3E287E85A763AF267492AB572E99AB3F370D275CEC1DA1AAA9075FF05F79BE)))

def doubleEval(r, p):
  wideR = r.untwist()
  slope = (wideR.x * wideR.x + wideR.x * wideR.x + wideR.x * wideR.x) * (wideR.y + wideR.y).inv()
  v = wideR.y - slope * wideR.x
  return Fq12.fromInt(p.y.x) - (Fq12.fromInt(p.x.x) * slope) - v

def addEval(r, q, p):
  wideR = r.untwist()
  wideQ = q.untwist()
  if wideR.x == wideQ.x and wideR.y == Fq12.fromInt(0)-wideQ.y:
    return Fq12.fromInt(p.x.x) - wideR.x
  slope = (wideQ.y - wideR.y) * (wideQ.x - wideR.x).inv()
  v = (wideQ.y * wideR.x - wideR.y * wideQ.x) * (wideR.x - wideQ.x).inv()
  return Fq12.fromInt(p.y.x) - (Fq12.fromInt(p.x.x) * slope) - v

def millerLoop(p, q):
  l=[c=='1' for c in bin(0xd201000000010000)[2:]]
  r = q
  result = Fq12.fromInt(1)
  for t in l:
    result *= result * doubleEval(r, p)
    r += r
    if t:
      r += q
      result *= addEval(r+r, q, p)
  return result

def pairing(p, q):
  return millerLoop(p, q) ** ((fieldPrime**12 - 1) // groupOrder)

print(pairing(g1Gen*5, g2Gen*5))
print(pairing(g1Gen*5, g2Gen*3) + pairing(g1Gen*5, g2Gen*2))

