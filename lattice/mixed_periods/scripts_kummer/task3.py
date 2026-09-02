from kummer import *
setup(80)
S3=sqrt(3); S2=sqrt(2)
def forms_k3(D,fam):
    r=mpf(D)**(mpf(1)/3)
    if fam=='M':
        A=log((r+1)/sqrt(r*r-r+1)); T=S3*atan(S3/(2*r-1))
        return {1:(r/D)*(A+T), 2:(r*r/D)*(T-A)}
    else:
        A=log((r-1)/sqrt(r*r+r+1)); T=S3*atan(S3/(2*r+1))
        return {1:(r/D)*(A-T), 2:(r*r/D)*(A+T)}
def forms_k4(D,fam):
    r=mpf(D)**(mpf(1)/4)
    if fam=='M':
        Lg=log((r*r+S2*r+1)/(r*r-S2*r+1))/2; Th=atan(S2*r/(r*r-1))
        return {1:(S2*r/D)*(Lg+Th), 2:(2*r*r/D)*atan(1/(r*r)), 3:(S2*r**3/D)*(Th-Lg)}
    else:
        return {1:(r/D)*(log((r-1)/(r+1))-2*atan(1/r)),
                2:(r*r/D)*log((r*r-1)/(r*r+1)),
                3:(r**3/D)*(log((r-1)/(r+1))+2*atan(1/r))}
print("Task 3: explicit real closed forms for c^{(a)}[1/(1-t)], verified against the general formula")
worst=0
for (k,Ns) in [(3,[9,18,27,54,81,63,126]),(4,[8,16,80,256,512,624])]:
    for N in Ns:
        for fam in ['M','P']:
            h=Host(k,N,fam); D=h.D
            F = forms_k3(D,fam) if k==3 else forms_k4(D,fam)
            for a in range(1,k):
                exact = re(h.cB(a)); got=F[a]
                e=float(abs(exact-got)/abs(exact)); worst=max(worst,e)
                print(f"  k={k} N={N} {fam} D={D} a={a}: {mp.nstr(exact,25)}  vs explicit {mp.nstr(got,25)}  relerr={e:.2e}")
print("worst relative error:",worst)
