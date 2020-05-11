function val = optimize_EJ12_symmetric(x,a)
options2 = optimset('TolX',1e-3,'TolFun',1e-3);
fun = @(y)optimize_frequencies_symmetric(y,[x,a(2),a(3),a(4),a(5),a(6)]);
y0 = a(1);
y = fminsearch(fun,y0,options2);

[J]=find_frequencies(y,y,x,a(2),a(2),a(3),a(4),a(5),0);

val=abs(J-a(7));

end
