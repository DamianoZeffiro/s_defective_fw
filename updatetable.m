function updatetable(tablename, cnum, nrun, s, cputimes, cliquesizes)
    load(tablename, 'sizetable', 'timestable'); 
    if nrun == 0
    meanC = 0;
    maxC = 0;
    stdC = 0;
    meant = 0;
    stdt = 0;        
    else
    meanC = mean(cliquesizes(cnum, 1:nrun, s));
    maxC = max(cliquesizes(cnum, 1:nrun, s));
    stdC = std(cliquesizes(cnum, 1:nrun, s));
    meant = mean(cputimes(cnum, 1:nrun, s));
    stdt = std(cputimes(cnum, 1:nrun, s));
    end
    sizetable{cnum, s+1} = table([maxC], [meanC], [stdC]);
    timestable{cnum, s+1} = table([meant], [stdt]);    
    save(tablename, 'sizetable', 'timestable');
    
end
