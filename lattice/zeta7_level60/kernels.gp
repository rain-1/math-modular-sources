\\ Support-depth (= denominator entropy 7/depth) of the three Atkin-Lehner
\\ kernel lines K_3, K_4, K_5 of the 3-dim purified level-60 space.
dl = divisors(60);
P(c,s) = sum(i=1,#dl, c[i]*dl[i]^(-s));
B0=[0,0,2673,-73216,398125,-497664,-3840000,13208832,-14478750,5280000,0,0];
B1=[0,33,0,-39424,276250,-375921,-2900625,9165312,-7796250,0,1670625,0];
B2=[1,0,0,-402688,3062500,-4313088,-33280000,101606400,-79633125,0,0,12960000];
comb(v) = vector(#dl, i, v[1]*B0[i]+v[2]*B1[i]+v[3]*B2[i]);
depth(c) = my(d=1); for(i=1,#dl, if(c[i]!=0, d=dl[i]; break)); d;
{lines = [[-1,-380,33],[14723,-112676,8349],[13,-52,3],[13/3,-52/3,1],[1,0,0],[0,1,0],[0,0,1]];}
{nms = ["K_3","K_4","K_5","Phi60^(5)","B_0","B_1","B_2"];}
{for(i=1,#lines, my(c=comb(lines[i]));
  print(nms[i], "  coords ", lines[i], "  support depth = ", depth(c),
        "   entropy 7/depth = ", 7/depth(c)*1.0,
        "   L(.,7)/zeta(7) = ", -P(c,7)/2));}
print("");
print("c_1 of K_3,K_4,K_5 (=0 iff the line lies in the depth>=2 plane F_2): ", vector(3,i, comb(lines[i])[1]));
print("=> no Atkin-Lehner kernel line meets F_2, so every monodromy-purified source has entropy 7.");
print("Do K_3,K_4,K_5 span the plane? rank = ", matrank(Mat([[-1,-380,33]~,[14723,-112676,8349]~,[13,-52,3]~])));
\q
