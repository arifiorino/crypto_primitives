import fields, ec

fields.fieldPrime=0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F

G=ec.Point(fields.Fq1(0x79be667ef9dcbbac55a06295ce870b07029bfcdb2dce28d959f2815b16f81798),fields.Fq1(0x483ada7726a3c4655da4fbfc0e1108a8fd17b448a68554199c47d08ffb10d4b8))

secret=0x51897b64e85c3f714bba707e867914295a1377a7463a9dae8ea6a8b914246319
print(f'Bitcoin private key: %x' % secret)
x=G*secret
print('Bitcoin public key: 0%d%x'%(2+x.y.x%2,x.x.x))


