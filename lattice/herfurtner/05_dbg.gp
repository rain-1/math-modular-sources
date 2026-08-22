read("05_jtest.gp");
{
 my(q,Jh,Jt,W,mat,ker,neq,d,V,U,R);
 q = nome(1,0,0, 9,3,27, NS);
 print("nome(ZagierB) = ", q + O(x^8));
 q = nome(1,0,0, 7,2,-8, NS);
 print("nome(ZagierA) = ", q + O(x^8));
 q = nome(1,0,0, 11,3,-1, NS);
 print("nome(ZagierD) = ", q + O(x^8));
 q = nome(2,1,1, 136,10,4, NS);
 print("nome(Apery-root) = ", q + O(x^8));
}
quit;
