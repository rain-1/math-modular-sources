default(parisize,"2G");
N=110;
Lt=vector(N+1);Lt[1]=1;for(n=1,N,Lt[n+1]=lcm(Lt[n],n)); L(n)=Lt[n+1];
intg(s)=my(r=O(x^(N+1)));for(k=0,N-1, r+=polcoeff(s,k)*x^(k+1)/(k+1));r;
exc(s,typ)={my(m=1,v=vector(N-1));for(n=1,N-1,my(dn=denominator(polcoeff(s,n)),T=typ(n));v[n]=dn/gcd(dn,T);m=max(m,v[n]));[m,v];}
report(name,s)={my(e1=exc(s,n->L(n)),e12=exc(s,n->L(n)*L(n\2)),e13=exc(s,n->L(n)*L(n\3)),e2=exc(s,n->L(n)^2));
 print(name,":  max excess over [n]: ",e1[1],"  over [n][n/2]: ",e12[1],"  over [n][n/3]: ",e13[1],"  over [n]^2: ",e2[1]);
 if(e12[1]>1 && e12[1]<10^6, print("     excess seq [n][n/2] (n=1..40): ",vector(40,i,e12[2][i])));
}
host(name,HA)={
 print("=== host ",name);
 my(lg=log(1-x+O(x^(N+2))));
 report("HB=HA*int HA/(1-t)",HA*intg(HA/(1-x)));
 report("HC=HA*int HA log/t",HA*intg(HA*lg/x));
 report("HD=HA*int HA log/(1-t)",HA*intg(HA*lg/(1-x)));
 report("HA*int HA log^2/(1-t)",HA*intg(HA*lg^2/(1-x)));
 report("HA*int HA log^2/t",HA*intg(HA*lg^2/x));
 report("HA*int (HA/t) int(HA/(1-t))/HA",HA*intg(HA/x*intg(HA/(1-x))/HA));
 report("HA*int (HA/(1-t)) int(HA/(1-t))/HA",HA*intg(HA/(1-x)*intg(HA/(1-x))/HA));
}
host("1/sqrt(1-4x)",1/sqrt(1-4*x+O(x^(N+2))));
host("1/sqrt(1-8x)",1/sqrt(1-8*x+O(x^(N+2))));
host("1/sqrt(1-12x)",1/sqrt(1-12*x+O(x^(N+2))));
host("1/sqrt((1-x)(1-9x))",1/sqrt((1-x)*(1-9*x)+O(x^(N+2))));
host("1/sqrt((1-x)(1-25x))",1/sqrt((1-x)*(1-25*x)+O(x^(N+2))));
host("1/sqrt((1+x)(1-3x))",1/sqrt((1+x)*(1-3*x)+O(x^(N+2))));
