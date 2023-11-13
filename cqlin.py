import kzg, bls12_381, fields, poly, util, random

n=2
p=257
alpha=4

M=[[1,2],[3,4]]
a=[1,2]
b=[5,11]

L=[[0]*n for _ in range(n)]
Ln=[[0]*n*n for _ in range(n)]
for j in range(n):
  L[j][j]=1
  L[j]=util.inv_dft(L[j],alpha,p)
  for i in range(n):
    Ln[j][i*n]=L[j][i]

print('L',L)
print('Ln',Ln)
R=[[0]*n for _ in range(n)]
for i in range(n):
  for j in range(n):
    R[i]=[(x+M[i][j]*y)%p for x,y in zip(R[i],L[j])]
M=[0]*n*n*n
for i in range(n):
  temp = poly.mult(Ln[i].copy(),R[i].copy(),p)
  M=[(x+y)%p for x,y in zip(M,temp)]
A=[0]*n*n
for k in range(n):
  for i in range(n):
    A[k*n]=(A[k*n]+a[i]*L[i][k]) %p
print('R',R)
print('M',M)
print('A',A)
Y=poly.mult(A.copy(),M.copy(),p)
print('Y',Y)

X=[0]*n*n*n
for i in range(n):
  temp=poly.mult(Ln[i].copy(),R[i].copy(),p)
  temp=[a[i]*temp[j] for j in range(n*n*n)]
  X=[(X[j]+temp[j])%p for j in range(n*n*n)]
print('X',X)
Z=[p-1]+[0]*(n*n-1)+[1]
print('Z',Z)
print('poly.div(Y,Z,p)',poly.div(Y,Z,p))
print('poly.div(X,Z,p)',poly.div(X,Z,p))

