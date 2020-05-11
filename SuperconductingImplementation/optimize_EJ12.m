function val = optimize_EJ12(x,a)

options2 = optimset('TolX',1e-2,'TolFun',1e-2);

fun = @(y)optimize_frequencies(y,[x,a(3),a(4),a(5),a(6),a(7),a(8),a(9)]);
y0 = [a(1),a(2)];
y = fminsearch(fun,y0,options2);

[J]=find_frequencies(y(1),y(2),x,a(3),a(4),a(5),a(6),a(7),0);

val=abs(J-a(10));

end
