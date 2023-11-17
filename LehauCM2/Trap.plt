set xrange [0:2200]
set ticslevel 0
plot "Trap.dat" u 1:2 w l, x**3*(1-x**2)**(0.5)