/* Check the claimed second-order recurrence of SOL notes 1 (sec 11) and 2 (sec 5),
   the denominator law den(V_n) | d_{2n}^2, and the characteristic roots. */
read("out_cellular.txt");
Pn(nn) = 3520*nn^6+26752*nn^5+83024*nn^4+134592*nn^3+120196*nn^2+56088*nn+10699;
lead(nn) = (nn+2)^2*(2*nn+3)^2*(20*nn^2+32*nn+13);
trail(nn) = 256*(nn+1)^2*(2*nn+1)^2*(20*nn^2+72*nn+65);
rec(X,nn) = lead(nn)*X[nn+3] - 4*Pn(nn)*X[nn+2] - trail(nn)*X[nn+1];
NM = #Uv-1;
{
print("recurrence residuals (0 = satisfied), n=0..",NM-2);
for(nn=0,NM-2, print(nn,"  U:",rec(Uv,nn),"   V:",rec(Vv,nn)));
print();
print("denominator law: den(V_n) vs d_{2n}^2");
for(nn=0,NM, my(dv=denominator(Vv[nn+1]), dd=lcm(vector(max(1,2*nn),k,k))^2);
   print(nn," den=",dv,"  divides d_{2n}^2 ? ",dd%dv==0, "  quotient=", if(dv!=0,dd/dv)));
print();
print("den(V_n) vs d_{2n}^2 also test d_n^2 and odd-part:");
for(nn=0,NM, my(dv=denominator(Vv[nn+1]), dn=lcm(vector(max(1,nn),k,k))^2);
   print(nn,"  d_n^2 % den = ", if(dv!=0, dn%dv)));
print();
print("char roots of lambda^2-176*lambda-256: ", polroots(x^2-176*x-256));
print("88-40*sqrt(5) = ", 88-40*sqrt(5), "   88+40*sqrt(5) = ", 88+40*sqrt(5));
print("ratio lambda+/|lambda-| = ", (88+40*sqrt(5))/(40*sqrt(5)-88));
print();
print("measured U_{n+1}/U_n and error ratios:");
for(nn=1,NM, print(nn,"  U ratio=",Uv[nn+1]*1./Uv[nn],
   "   err=",(-Vv[nn+1]/Uv[nn+1]-Catalan)));
}
quit();
