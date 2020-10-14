import math as m 

def leglength(lc,lb,lt,beta):

    Lc = lc/m.cos(m.pi/4)
    Ltot = lt-(lb/2)+Lc
    Ll = m.sin(m.pi/4)*Ltot/(m.sin(m.pi/4+beta))
    return Ll, beta, lb

def CGM(Ll, beta, lb):

    CGM = lb/2 + Ll*m.sin(beta) - Ll*m.cos(beta)

    return CGM

Ll, beta, lb = leglength(50, 260, 1690/2, m.pi*2/9)
CGM = CGM(Ll, beta, lb)
print(CGM, Ll)
