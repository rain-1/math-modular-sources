"""SOL note 2 sec 2: moments I_n = int_0^1 t^n F(t)^2 dt, F=2F1(1/4,3/4;1;t)."""
from mpmath import mp, mpf, hyp2f1, quad, catalan, pi, mpmathify, nstr
mp.dps = 40
Ff = lambda t: hyp2f1(mpf(1)/4, mpf(3)/4, 1, t)
def I(n, tol=None):
    return quad(lambda t: t**n * Ff(t)**2, [0, mpf('0.5'), mpf('0.9'), mpf('0.99'), 1])
Iv = [I(n) for n in range(0,8)]
print("I_n for n=0..7:")
for n,v in enumerate(Iv): print("  I_%d = %s"%(n, nstr(v,30)))
print()
print("S1: I_0 vs 16G/pi^2 :", nstr(Iv[0],30), nstr(16*catalan/pi**2,30), " diff", nstr(Iv[0]-16*catalan/pi**2,8))
print("S4: I_1 vs (8G+4/3)/pi^2:", nstr(Iv[1],30), nstr((8*catalan+mpf(4)/3)/pi**2,30), " diff", nstr(Iv[1]-(8*catalan+mpf(4)/3)/pi**2,8))
print("    1/pi^2 vs (3/4)(I_1-I_0/2):", nstr(1/pi**2,30), nstr(mpf(3)/4*(Iv[1]-Iv[0]/2),30))
print()
print("S3: recurrence 8n^3 I_{n-1} - (2n+1)(8n^2+8n+3) I_n + 2(n+1)(2n+1)(2n+3) I_{n+1} = 8/pi^2")
for n in range(1,7):
    lhs = 8*n**3*Iv[n-1] - (2*n+1)*(8*n**2+8*n+3)*Iv[n] + 2*(n+1)*(2*n+1)*(2*n+3)*Iv[n+1]
    print("  n=%d: lhs=%s   8/pi^2=%s   diff=%s"%(n, nstr(lhs,25), nstr(8/pi**2,25), nstr(lhs-8/pi**2,8)))
