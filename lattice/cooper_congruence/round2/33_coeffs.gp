\\ 33_coeffs.gp -- rebuilds f4a, f4b coefficient vectors A4A, A4B indexed n=-4..NMAXC
NMAXC = 3000;
BPC = NMAXC + 40;
TH = 1 + 2*sum(n=1, sqrtint(BPC), 'q^(n^2)) + O('q^BPC);
{ E24 = sum(n=1, BPC-1, if(n%2, sigma(n,1), 0)*'q^n) + O('q^BPC); }
{ E4q(t,P) = 1 + 240*sum(n=1, (P-1)\t, sigma(n,3)*'q^(t*n)) + O('q^P); }
{ E6q(t,P) = 1 - 504*sum(n=1, (P-1)\t, sigma(n,5)*'q^(t*n)) + O('q^P); }
E44 = E4q(4,BPC);
E64 = E6q(4,BPC);
DEL4 = (E44^3 - E64^2)/1728;
G0 = TH*(TH^4 - 20*E24);
G1 = TH*E44^2*E64/DEL4;
G2 = G0*E44^3/DEL4;
F4A = 7/8*G0 + 1/768*G1 - 1/768*G2;
F4B = 19/18*G0 - 5/648*G1 - 1/648*G2;
A4A = vector(NMAXC+5, i, polcoeff(F4A, i-5));
A4B = vector(NMAXC+5, i, polcoeff(F4B, i-5));
