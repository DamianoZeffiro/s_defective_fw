function Ay = embedfs(Ay, Ayface, Sfrac, Igyface, fsize, Icurr, totstep, s)
if fsize > 0
Ayface(max(Icurr) + 1:fsize) = totstep * Ayface(max(Icurr) + 1:fsize);
if sum(Ay, 'all') > s + 0.01
         error('ll');
end
Ay(Sfrac(Igyface)) = Ayface; 
if sum(Ay, 'all') > s + 0.01
         warning('ll');
end
end
end