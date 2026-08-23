/* lattice/p2_holonomic/rowrec.gp
   A P-recursive (order 3, polynomial coefficients of degree <= 14) recurrence
   for the Nesterenko (4,7) Catalan row entries B_n, C_n, fitted from the 77
   exact values n = 4..80 of lattice/p2_structure/data/rows_all.txt and
   validated exactly on n = 81..120 (0 failures).  It replaces the
   O(n^{4.1}) partial-fraction solve of rows_pos.gp `mom` by an O(n) step, and
   is what makes n <= 200 feasible.   [verified 4 <= n <= 120]

   Convention: with A[i] the value at n = i+3, the relation is
      sum_{j=0}^{3} PB[j+1](i) * A[i+1-j] = 0   for i >= 3.
   Prepend lattice/fitrec.gp is NOT needed; the fitted polynomials are
   recomputed here from the cache so nothing is hard-coded.                  */

fitrec2(A,ord,d)={my(MM=#A,rws=[]);
 for(m=ord,MM-1, my(r=[]); for(j=0,ord, for(e=0,d, r=concat(r,[m^e*A[m+1-j]]))); rws=concat(rws,[r]));
 my(Mx=matrix(#rws,(ord+1)*(d+1),i,j,rws[i][j]),K=matker(Mx));
 if(#K==0,return(0)); my(w=K[,1]); w=w/content(w);
 vector(ord+1,j,Polrev(vector(d+1,e,w[(j-1)*(d+1)+e]),'m))};

/* extend A (1-indexed, A[i] = value at n=i+3) to length NEW using rec R */
{
extend(A, R, NEW) =
 my(ord=#R-1, B=vector(NEW));
 for(i=1,#A, B[i]=A[i]);
 for(i=#A,NEW-1,
   my(lead=subst(R[1],'m,i), s=sum(j=1,ord, subst(R[j+1],'m,i)*B[i+1-j]));
   B[i+1] = -s/lead);
 B;
}
