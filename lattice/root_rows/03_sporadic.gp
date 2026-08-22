read("lattice/root_rows/lib.gp");
default(realprecision,250);
N=420;

\\ --- six Almkvist-Zudilin rows: (n+1)^3 A_{n+1}=(2n+1)(a n^2+a n+b)A_n - c n^3 A_{n-1}
azrow(a,b,c,N)={my(A=vector(N+1)); A[1]=1; A[2]=b;
 for(n=1,N-1, A[n+2]=((2*n+1)*(a*n^2+a*n+b)*A[n+1]-c*n^3*A[n])/(n+1)^3); A};
AZ=[[17,5,1,"Apery"],[12,4,16,"T"],[10,4,64,"Domb"],[9,3,-27,"AZ(9,3,-27)"],[11,5,125,"AZ(11,5,125)"],[7,3,81,"AZ(7,3,81)"]];

\\ --- Cooper's three
coop(fn,b,N)={my(A=vector(N+1));A[1]=1;A[2]=b;
 for(n=1,N-1,A[n+2]=fn(n,A[n+1],A[n])/(n+1)^3); A};
f10(n,u1,u0)=2*(2*n+1)*(3*n^2+3*n+1)*u1 + 4*n*(16*n^2-1)*u0;
f7(n,u1,u0)=(2*n+1)*(13*n^2+13*n+4)*u1 + 3*n*(9*n^2-1)*u0;
f18(n,u1,u0)=2*(2*n+1)*(7*n^2+7*n+3)*u1 - 12*n*(16*n^2-1)*u0;

ROWS=[];
{for(i=1,#AZ, ROWS=concat(ROWS,[[AZ[i][4], azrow(AZ[i][1],AZ[i][2],AZ[i][3],N), 2]]));}
ROWS=concat(ROWS,[["Cooper s_7", coop(f7,4,N),2],["Cooper s_10",coop(f10,2,N),2],["Cooper s_18",coop(f18,6,N),2]]);

{for(i=1,#ROWS,
  my(nm=ROWS[i][1], A=ROWS[i][2], w=ROWS[i][3]);
  print("\n########## ",nm,"  (w=",w,") ##########");
  print("A_0..A_6 = ", vector(7,j,A[j]));
  print("e_2 = ",emin(A,2,N), "   e_3 = ",emin(A,3,N));
  my(rr=rootrow(A,w,N), lam=rr[1], a=rr[2]);
  print("lambda = ",lam,"   (predicted max(1, 2^(2-e_2)) = ", max(1,2^(2-emin(A,2,N))),")");
  print("a_0..a_5 = ",vector(6,j,a[j]));
  my(chk=1); for(n=1,N, if(denominator(a[n+1])!=1, chk=0;print("NONINT at n=",n);break));
  print("a_n integral to n=",N,": ",chk);
  my(mr=minrec(a,4,4));
  if(type(mr)=="t_INT", print("NO recurrence found (order<=4,deg<=4)"),
    print("minimal recurrence: order=",mr[1]," deg=",mr[2]);
    print("   coeffs: ",mr[3]);
    my(rts=charroots(mr[3]));
    print("   char roots: ",rts);
  );
);}
\q
