M=300;
sig = sqrt(1-34*t+t^2+O(t^(M+1)));
print("sigma coeffs 0..10: ", vector(11,i,polcoeff(sig,i-1)));
bad=0; {for(i=0,M, if(denominator(polcoeff(sig,i))!=1, bad++));}
print("sigma integral to t^",M,"? bad=",bad);
print("v_2 of sigma coeffs 1..20: ", vector(20,i,valuation(polcoeff(sig,i),2)));
\\ 1-34t+t^2 = (1-17t)^2 - 288 t^2
print("check (1-17t)^2-288t^2 = ", (1-17*t)^2-288*t^2);
quit;
