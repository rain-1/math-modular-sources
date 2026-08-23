/* Scoring the cyclotomic cellular Catalan row against the CDT/obstruction framework. */
read("out_cellular.txt");
NM = #Uv-1;
{
print("=== singular points of the row (roots of 256x^2+176x-1) ===");
print("  polroots: ", polroots(256*x^2+176*x-1));
print("  t1 = (-11+5*sqrt(5))/32 = ", (-11+5*sqrt(5))/32, "   1/lambda_+ = ", 1/(88+40*sqrt(5)));
print("  t2 = (-11-5*sqrt(5))/32 = ", (-11-5*sqrt(5))/32, "   1/lambda_- = ", 1/(88-40*sqrt(5)));
print("  t1*t2 = ", (-11+5*sqrt(5))/32*(-11-5*sqrt(5))/32, "   (= -1/256)");
print("  16*|t2| = ", 16*abs((-11-5*sqrt(5))/32), "   phi^5 = ", ((1+sqrt(5))/2)^5);
print("  16*|t1| = ", 16*abs((-11+5*sqrt(5))/32), "   phi^-5 = ", ((1+sqrt(5))/2)^-5);
print("  ceiling log(16|t2|) = ", log(16*abs((-11-5*sqrt(5))/32)), "  = 5*log(phi) = ", 5*log((1+sqrt(5))/2));
print();
print("=== denominator growth of the companion V_n ===");
print("  n | (1/n)log den(V_n) | (1/n)log(d_{2n}^2) | ratio | d_{2n}^2/den(V_n)");
for(nn=6,NM, if(nn%4==0 || nn==NM,
  my(dv=denominator(Vv[nn+1]), dd=lcm(vector(2*nn,k,k))^2);
  print("  ",nn,"   ",log(dv*1.)/nn,"   ",log(dd*1.)/nn,"   ",log(dv*1.)/log(dd*1.),"   ",dd/dv)));
print();
print("=== raw Apery score: log|lambda_-| + sigma ===");
print("  log|lambda_-| = ", log(abs(88-40*sqrt(5))));
print("  sigma (den type [1..2n]^2) = 4;  need log|lambda_-| - 0 ... score = ",log(abs(88-40*sqrt(5)))+4);
print("  (positive => no irrationality)");
print();
print("=== entry numbers (width law: entry = (1/w)log16 + log|t2| - k) ===");
print("  w=1, t2=-0.6931, k=sigma=4  ->  entry = ", log(16.)+log(abs((-11-5*sqrt(5))/32))-4);
print("  hypothetical k=2           ->  entry = ", log(16.)+log(abs((-11-5*sqrt(5))/32))-2);
print("  level-8 modular host        ->  entry = ", log(16.)+log(0.25)-2);
}
quit();
