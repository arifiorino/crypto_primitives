import fields

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
    t=fields.Fq1(4)
    if isinstance(self.x,fields.Fq2):
      t=fields.Fq2(fields.Fq1(4),fields.Fq1(4))
    return self.y*self.y == self.x*self.x*self.x + t
  def __str__(self):
    return '(%s,%s)'%(str(self.x),str(self.y))
  def __eq__(self,b):
    return (self.x,self.y)==(b.x,b.y)
  def untwist(self):
    assert isinstance(self.x,fields.Fq2)
    root=fields.Fq6(fields.Fq2.fromInt(0), fields.Fq2.fromInt(1), fields.Fq2.fromInt(0))
    return Point(fields.Fq12(fields.Fq6(self.x, fields.Fq2.fromInt(0), fields.Fq2.fromInt(0)), fields.Fq6.fromInt(0)) * fields.Fq12(root, fields.Fq6.fromInt(0)).inv(),
                 fields.Fq12(fields.Fq6(self.y, fields.Fq2.fromInt(0), fields.Fq2.fromInt(0)), fields.Fq6.fromInt(0)) * fields.Fq12(fields.Fq6.fromInt(0), root).inv())
    

