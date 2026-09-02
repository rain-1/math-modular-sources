default(parisize,"10G");
read("lib.gp");
N = 3000;
{ CASES = [[1,2],[2,3],[3,3],[5,2],[6,2],[6,3],[7,2],[9,3],[10,2],[11,3],[12,5],[15,3]]; }
XI = vector(#CASES);
{ for(j=1,#CASES,
  my(i=CASES[j][1], p=CASES[j][2], R=ROWS[i], AB=genrow(R,N), A=AB[1], B=AB[2]);
  my(x1=B[N+1]/A[N+1], x0=B[N]/A[N], cau=if(x1==x0, 99999, valuation(x1-x0,p)), vx=if(x1==0,99999,valuation(x1,p)));
  XI[j] = x1;
  print(R[1], " p=", p, " : v_p(xi)=", vx, "  Cauchy v_p(xi_N - xi_{N-1})=", cau, "  xi*p^-v mod p^8 = ", if(x1==0, 0, lift(Mod(x1/p^vx, p^8)))));
}
print("");
print("Conjecture-D ratio tests (v_p of the exact rational combination at N=", N, "):");
print("  p=3  xi_B - xi_C            : ", valuation(XI[2]-XI[3],3));
print("  p=3  xi_B - xi_s18          : ", valuation(XI[2]-XI[12],3));
print("  p=3  4*xi_F - 5*xi_B        : ", valuation(4*XI[6]-5*XI[2],3));
print("  p=3  4*xi_delta - ?         : v_3(xi_delta) = ", valuation(XI[8],3), " (no partner row)");
print("  p=2  3*xi_alpha - 4*xi_eps  : ", valuation(3*XI[7]-4*XI[9],2));
print("  p=2  xi_E                   : v_2 = ", valuation(XI[4],2));
print("  p=2  xi_A (should be 0)     : v_2 = ", valuation(XI[1],2));
print("  p=3  xi_zeta (should be 0)  : v_3 = ", valuation(XI[11],3));
print("  p=2  xi_F                   : v_2 = ", valuation(XI[5],2));
print("  p=5  xi_eta                 : v_5 = ", valuation(XI[10],5));
quit;
