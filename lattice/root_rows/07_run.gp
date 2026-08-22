read("lattice/root_rows/07_scan.gp");
{for(k=1,6,
 my(nms=["beta4_L24","Lchi3_L24","zeta7_L24","zeta5_L16","zeta5_L12p","zeta5_L12m"],
    ws=[3,3,6,4,4,4]);
 my(A=readrow(nms[k]), rr=rootrow(A,ws[k],#A-1), a=rr[2]);
 print("\n==== ",nms[k]," (lambda=",rr[1],") kernel dims (-1 = too few equations) ====");
 scan(a, 8, 14));}
\q
