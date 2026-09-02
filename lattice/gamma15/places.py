# -*- coding: utf-8 -*-
"""The two real places of K = Q(sqrt5) for Zagier's row D on Gamma_1(5).

t1 = phi^-5 = (5 sqrt5 - 11)/2  (near cusp, the fold: the Apery limit lives there)
t2 = -phi^5 = -(11 + 5 sqrt5)/2 (far cusp),  t1 t2 = -1, t1 + t2 = -11, t1 - t2 = 5 sqrt5.
The descent uses s := t2 (the outer singularity, forced: H is regular at t1 and singular
at t2), y = x^2/(x - s), branch point y = 4s, extra point y(t1) = t1^2/(t1 - t2).
Under sigma: sqrt5 -> -sqrt5 one has sigma(t1) = t2 and sigma(t2) = t1, so at the second
real place s is realised as t1 = +0.0901699... and the extra point as t2^2/(t2-t1).
Exact:  Y_v := y_v/s_v  with  Y_1 = -phi^-15/(5 sqrt5) = (682 sqrt5 - 1525)/25,
                              Y_2 = -phi^{15}/(5 sqrt5) = -(682 sqrt5 + 1525)/25.
(phi^15 = 682 + 305 sqrt5; L_15 = 1364, F_15 = 610.)   Y_1 Y_2 = 1/125, Y_2/Y_1 = phi^30."""
import math
s5 = math.sqrt(5)
phi = (1+s5)/2
t1 = phi**-5
t2 = -phi**5
s_1 = t2; y_1 = t1*t1/(t1-t2); Y_1 = y_1/s_1
s_2 = t1; y_2 = t2*t2/(t2-t1); Y_2 = y_2/s_2
PLACES = [('v1', abs(s_1), Y_1), ('v2', abs(s_2), Y_2)]
TAU = 16603/3920
