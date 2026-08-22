\\ Genera of all Atkin-Lehner quotients of X_0(60):  g(X_0(60)/H) = dim S_2(Gamma_0(60))^H.
mf = mfinit([60,2],1); n = mfdim(mf);
print("g_0(60) = dim S_2 = ", n);
QS = [3,4,5];
{WM = vector(3, i, my(a=mfatkininit(mf,QS[i])); a[2]/a[3]);}
{for(i=1,3, print("W_",QS[i],": involution? ", WM[i]^2==matid(n)));}
{sub = [[],[1],[2],[3],[1,2],[1,3],[2,3],[1,2,3]];
 nm  = ["1","W3","W4","W5","W3,W4 (=+W12)","W3,W5 (=+W15)","W4,W5 (=+W20)","W3,W4,W5 (full)"];
 for(s=1,#sub, my(S=sub[s], Mx);
   Mx = if(#S==0, 0*matid(n), matconcat(vector(#S, i, WM[S[i]]-matid(n))~));
   print("genus X_0(60)/<",nm[s],"> = ", #matker(Mx)));}
\\ also the single Fricke W_60 = W_3 W_4 W_5 and the pair-subgroups <W_5,W_60>
W60 = WM[1]*WM[2]*WM[3];
print("genus X_0(60)/<W60> = ", #matker(W60-matid(n)));
print("genus X_0(60)/<W5,W60> = ", #matker(matconcat([WM[3]-matid(n); W60-matid(n)])));
W12 = WM[1]*WM[2]; W15=WM[1]*WM[3]; W20=WM[2]*WM[3];
{for(i=1,4, my(W=[W12,W15,W20,W60][i], L=[12,15,20,60][i]);
  print("genus X_0(60)/<W",L,"> = ", #matker(W-matid(n))));}
\q
