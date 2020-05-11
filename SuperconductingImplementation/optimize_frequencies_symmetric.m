function val = optimize_frequencies_symmetric(x,a)
[~,Delta]=find_frequencies(x,x,a(1),a(2),a(2),a(3),a(4),a(5),0);
val=abs(Delta-a(6));

end

