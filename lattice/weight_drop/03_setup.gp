/* shared: builds t,F,Phi,Xi q-series for the three Cooper rows to NQ terms */
ee(d) = eta(q^d + O(q^NQ));
E2s(d) = 1 - 24*sum(n=1,(NQ-1)\d, sigma(n)*q^(d*n)) + O(q^NQ);
th(f) = q*deriv(f,q);
{build() =
  XX7 = (ee(1)/ee(7))^4;            T7  = q*XX7/(XX7^2+13*q*XX7+49*q^2);
  F7  = (1/6)*(7*E2s(7)-E2s(1));    PHI7 = F7*th(T7);
  XX10= (ee(2)*ee(5)/(ee(1)*ee(10)))^6; T10 = q*XX10/(XX10-q)^2;
  F10 = (1/12)*(10*E2s(10)+5*E2s(5)-2*E2s(2)-E2s(1)); PHI10 = F10*th(T10);
  XX18= ee(3)^4*ee(6)^4/(ee(1)^2*ee(2)^2*ee(9)^2*ee(18)^2); T18 = q*XX18/(XX18+3*q)^2;
  F18 = (1/4)*(18*E2s(18)-9*E2s(9)-12*E2s(6)+6*E2s(3)+2*E2s(2)-E2s(1)); PHI18 = F18*th(T18);
  DAT = [["s7",7,T7,F7,PHI7,XX7],["s10",10,T10,F10,PHI10,XX10],["s18",18,T18,F18,PHI18,XX18]];}
