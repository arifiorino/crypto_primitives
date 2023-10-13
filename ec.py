def mod_inv(a,n):
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


# elliptic curves
p,a,b=17,0,7

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
    l = (3*xp*xp+a)*mod_inv(2*yp%p,p) %p
  else:
    if xp==xq:
      return None
    l = (yq-yp)*mod_inv(xq-xp%p,p) %p
  xr=(l*l-xp-xq)%p
  return (xr,(l*(xp-xr)-yp) %p)


u=(15,13)
curr=u
for i in range(20):
  print(i,curr)
  curr=add(curr,u)




