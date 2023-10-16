import random

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
    res = Fq12.fromInt(1)
    temp = self
    while k>0:
      if k&1:
        res*=temp
      temp*=temp
      k>>=1
    return res
  def inv(self):
    factor=(self.x*self.x - (self.y*self.y).mul_nonres()).inv()
    zero=Fq6.fromInt(0)
    return Fq12(self.x*factor, (zero-self.y)*factor)
  def __str__(self):
    return '%s,%s'%(str(self.x),str(self.y))
  def __eq__(self,b):
    return (self.x,self.y)==(b.x,b.y)
  @staticmethod
  def fromInt(x):
    return Fq12(Fq6.fromInt(x),Fq6.fromInt(0))

Fq1Rand=lambda: Fq1(random.randint(0,fieldPrime))
Fq2Rand=lambda:Fq2(Fq1Rand(),Fq1Rand())
Fq6Rand=lambda:Fq6(Fq2Rand(),Fq2Rand(),Fq2Rand())
Fq12Rand=lambda:Fq12(Fq6Rand(),Fq6Rand())

