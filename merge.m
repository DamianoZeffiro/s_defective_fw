function Iyfwold2 = merge(Iyfw, Vyfw, Iyfwold, Vyfwold, s)
i = 1;
j = 1;
k = 1;
Iyfwold2 = zeros(1, s);
for k = 1:s
    while i < s && ismember(Iyfw(i), Iyfwold2)
       i = i + 1; 
    end
    while j < s && ismember(Iyfwold(j), Iyfwold2)
        j = j + 1;
    end
    if Vyfw(i) > Vyfwold(j)        
        Iyfwold2(k) = Iyfw(i);
        i = i + 1;
    else
        Iyfwold2(k) = Iyfwold(j);
        j = j + 1;
    end
end