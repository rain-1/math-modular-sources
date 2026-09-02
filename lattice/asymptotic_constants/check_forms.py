from mpmath import mp, mpf, sqrt, pi, gamma, nstr
mp.dps = 50
rows = {
 "alpha": (12, 16, 4, "0.359174244250333123378163967255"),
 "gamma": (6, 17+12*sqrt(2), 17-12*sqrt(2), "0.220043767112643037850689759810"),
 "eps":   (8, 12+8*sqrt(2), 12-8*sqrt(2), "0.257797308917838939757311495404"),
 "zeta":  (9, 9+6*sqrt(3), 9-6*sqrt(3), "0.260201758994440992303217903862"),
 "s7":    (7, 27, -1, "0.233290514929399000627453098124"),
 "s10":   (10, 16, -4, "0.253974543736963879143053219739"),
 "s18":   (18, 16, 12, "0.761923631210891637429159659216"),
}
for k,(N,l1,l2,val) in rows.items():
    K = sqrt(N)/(2*pi**1.5)*sqrt(mpf(l1)/(l1-l2))
    print(k, nstr(K,32), " diff vs table:", nstr(K-mpf(val),3))
closed = {
 "alpha 2/pi^1.5": 2/pi**1.5,
 "gamma (1+sqrt2)^2/2^(9/4)": (1+sqrt(2))**2/(2**mpf(2.25)*pi**1.5),
 "eps sqrt(4+3sqrt2)/2": sqrt(4+3*sqrt(2))/(2*pi**1.5),
 "zeta 3(sqrt2+sqrt6)/8": 3*(sqrt(2)+sqrt(6))/(8*pi**1.5),
 "s7 3sqrt3/4": 3*sqrt(3)/(4*pi**1.5),
 "s10 sqrt2": sqrt(2)/pi**1.5,
 "s18 3sqrt2": 3*sqrt(2)/pi**1.5,
}
for k,v in closed.items(): print(k, nstr(v,32))
# level-12 zeta(5): old vs new closed form
G3 = gamma(mpf(1)/3)**12
old = G3/pi**(mpf(17)/2) * 3**(mpf(8)/3)*(15*sqrt(3)-23)*sqrt((12+7*sqrt(3))/2) / (2**(mpf(37)/3)*(2205-1273*sqrt(3))*(45+26*sqrt(3))**(mpf(1)/3))
new = 9*((3*(45+26*sqrt(3)))/2**17)**(mpf(2)/3)*sqrt((7+4*sqrt(3))/(24*sqrt(3)))*G3/pi**(mpf(17)/2)
print("zeta5 L12 old:", nstr(old,40)); print("zeta5 L12 new:", nstr(new,40)); print("old-new:", nstr(old-new,3))
# level-16 zeta(5)
G4 = gamma(mpf(1)/4)**4
K16 = (10021-7083*sqrt(2))*(2+sqrt(2))*sqrt(4+sqrt(2))/16 * G4/pi**(mpf(9)/2)
print("zeta5 L16:", nstr(K16,40), " table: 2.04997125634010629357510092802")
# level-12 zeta(7)
K7 = 3*sqrt(3)/2*2025*(490730*sqrt(3)-849969)*((3*(45+26*sqrt(3)))/2**17)**(mpf(2)/3)*sqrt(mpf(5)/(2016-1152*sqrt(3)))*G3/pi**(mpf(19)/2)
print("zeta7 L12:", nstr(K7,40), " table: 31.7215208279640907727003835703")
