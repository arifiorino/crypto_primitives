def mod_inv(a,n): #Using extended euclidian algorithm
  t,newt=0,1
  r,newr=n,a
  while newr!=0:
    quotient=r//newr
    t,newt=newt,t-quotient*newt
    r,newr=newr,r-quotient*newr
  if r>1:
    raise Exception("not invertible")
  if t<0:
    t+=n
  return t

p=0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F
a,b=0,7

verify=lambda x: ((x[1]*x[1])%p)==(x[0]*x[0]*x[0]+a*x[0]+b)%p
def add(P,Q):
  if not P:
    return Q
  if not Q:
    return P
  xp,yp=P
  xq,yq=Q
  if P==Q:
    if yp==0:
      return None
    l = (3*xp*xp+a)*mod_inv((2*yp)%p,p) %p
  else:
    if xp==xq:
      return None
    l = (yq-yp)*mod_inv((xq-xp)%p,p) %p
  xr=(l*l-xp-xq)%p
  return (xr,(l*(xp-xr)-yp) %p)

def mult(P,s):
  res = None
  temp = P
  while s>0:
    if s&1:
      res=add(res,temp)
    temp=add(temp,temp)
    s>>=1
  return res

G=(55066263022277343669578718895168534326250603453777594175500187360389116729240,32670510020758816978083085130507043184471273380659243275938904335757337482424)

secret=0x227dbb8586117d55284e26620bc76534dfbd2394be34cf4a09cb775d593b6f2b
x=mult(G,secret)
print('Private key',hex(secret))
print('Public key',hex(x[0]),hex(x[1]))




