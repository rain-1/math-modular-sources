"""Task 5: how much extra supply would H2/H3 need for a positive margin?"""
import math
from final import best, HOSTS
print("(a) single-layer supply nu (nu=1 is the measured/proved value); measured top block, g=1")
print("| host | nu | best entry@ceil | best entry@contour | best margin@contour | m,a |")
print("|---|---|---|---|---|---|")
for h, k, ceil, BC, D, I in HOSTS:
    for nu in (1, 2, 3, 4, 6, 8, 12):
        r = best(k, ceil, BC, D, I, nu=nu, key='marR')
        e = best(k, ceil, BC, D, I, nu=nu, key='ec')
        print(f"| {h} | {nu} | {e['ec']:+.4f} | {e['er']:+.4f} | {r['marR']:+.3f} | m={r['m']},a={r['a']} |")
print()
print("(b) H3: extra pure atoms at (2,0)  [measured q2 = 4 incl. F_1,F_2,B_4,B_5], g=1")
print("| q2 (atoms at L=2,e=0) | best entry@ceil | best margin@contour | m,a,n2 |")
print("|---|---|---|---|")
h, k, ceil, BC, D, I = HOSTS[2]
for q2 in (3, 4, 6, 8, 12, 16, 24, 32, 48):
    r = best(k, ceil, BC, D, I, q2=q2, key='marR')
    e = best(k, ceil, BC, D, I, q2=q2, key='ec')
    print(f"| {q2} | {e['ec']:+.4f} | {r['marR']:+.3f} | m={r['m']},a={r['a']},n2={r['n2']} |")
print()
print("(c) H2: extra pure atoms at (2,0)=top level -- cannot help u_2, only m")
h, k, ceil, BC, D, I = HOSTS[1]
for qk in (2, 6, 12, 24, 48):
    r = best(k, ceil, BC, D, I, qk=qk, key='marR')
    print(f"|  qk={qk} | best margin {r['marR']:+.3f} at m={r['m']}, a={r['a']} |")
print()
print("(d) sensitivity to the conditional orbit (D derivatives, I integrations), g=1")
for h, k, ceil, BC, D0, I0 in HOSTS:
    for (D, I) in ((D0, I0), (D0, 6), (2*D0, I0), (2*D0, 2*I0)):
        r = best(k, ceil, BC, D, I, key='marR')
        print(f"|  {h:16s} D={D:2d} I={I:2d} | margin {r['marR']:+.3f} at m={r['m']} | entry@ceil {r['ec']:+.4f} |")
