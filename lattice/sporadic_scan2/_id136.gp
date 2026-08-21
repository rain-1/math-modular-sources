default(parisizemax, 4000000000);
read("identlib.gp");
default(realprecision, 370);
x = eval(Str(externstr("cat _lim136.txt")[1]));
print("x = ", x);
print("rational? ", lindep([1,x]));
B = constbasis(370);
print("--- singles ---");
print(tryrel(x, B, 360));
