function val = optimize_frequencies(x,a)
[~,Delta1,Delta2]=find_frequencies(x(1),x(2),a(1),a(2),a(3),a(4),a(5),a(6),0);

val=abs(Delta1-a(7))+abs(Delta2-a(8));

end

