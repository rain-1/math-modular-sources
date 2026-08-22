MM=40; default(seriesprecision,MM);
Ap=vector(MM);Ap[1]=1;Ap[2]=5;for(n=1,MM-2,Ap[n+2]=((2*n+1)*(17*n^2+17*n+5)*Ap[n+1]-n^3*Ap[n])/(n+1)^3);
a=vector(MM);a[1]=1;a[2]=10;for(n=1,MM-2,a[n+2]=((136*n^2+68*n+10)*a[n+1]-4*(2*n-1)^2*a[n])/(n+1)^2);
F=sum(n=0,MM-1,Ap[n+1]*(4*u)^n)+O(u^MM); f=sum(n=0,MM-1,a[n+1]*u^n)+O(u^MM);
th(g)=u*deriv(g,u);
L1(g)=(1-136*u+16*u^2)*th(th(g))+(-68*u+16*u^2)*th(g)+(-10*u+4*u^2)*g;
LAp(g)=th(th(th(g)))-4*u*(34*th(th(th(g)))+51*th(th(g))+27*th(g)+5*g)+16*u^2*(th(th(th(g)))+3*th(th(g))+3*th(g)+g);
print("L1 f = 0: ",L1(f)==0,"   LAp F = 0: ",LAp(F)==0,"   f^2 = F: ",f^2-F==0);
\q
