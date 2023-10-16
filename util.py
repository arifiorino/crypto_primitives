def mult_inv(x,p):
  t,newt=0,1
  r,newr=p,x
  while newr!=0:
    quotient=r//newr
    t,newt=newt,t-quotient*newt
    r,newr=newr,r-quotient*newr
  if r>1:
    raise Exception("not invertible")
  if t<0:
    t+=p
  return t

