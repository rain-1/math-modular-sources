from fractions import Fraction as Fr
exec(open('18_padic2.py').read().split('N = 5000')[0])
N = 5000
A, B = seqs(N)
for a0 in [1,3,5,7,9,11]:
    print("a0 =", a0, ":  (n, v2(b_n/a_n))")
    s = 0; out=[]
    while a0*2**s < N-1:
        n = a0*2**s
        out.append((n, vp(B[n]/A[n], 2)))
        s += 1
    print("   ", out)
