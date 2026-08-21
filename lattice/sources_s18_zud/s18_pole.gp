default(realprecision,50);
/* x_18(tau) as eta quotient, evaluated numerically */
x18(tau) = my(f=d->eta(d*tau,1)); f(3)^4*f(6)^4/(f(1)^2*f(2)^2*f(9)^2*f(18)^2);
tau0 = 1/2 + I/6;
print("tau0 = ", tau0);
print("x_18(tau0) = ", x18(tau0));
print("q0 = exp(2*Pi*I*tau0) = ", exp(2*Pi*I*tau0));
print("|q0| = ", abs(exp(2*Pi*I*tau0)), "   e^{-Pi/3} = ", exp(-Pi/3));
print("t_18(tau0) = x/(x+3)^2 -> ", my(X=x18(tau0)); X/(X+3)^2);
print("x_18(1/2+I/6) fixed by W_9?  9/x = ", 9/x18(tau0));
/* other AL fixed point x=3 */
print("--- growth of Phi coefficients vs 1/|q0| ---");
print("1/|q0| = ", 1/exp(-Pi/3));
quit;
