default(parisizemax, 6000000000);
read("sc_rows.gp");
{
A4a = r4A(120);
A4b = genseq(R4cf, 6, [1], 120);
print("R4 recurrence matches series cube: ", A4a == A4b);
print("A4 first 8: ", vector(8,i,A4b[i]));
}
