function missingedges = missingedgecount(Q, findx)
missingedges = 0;
for i = 1:size(findx, 1)
    for j = 1:size(findx, 1)
        if (i < j)
            if Q(findx(i), findx(j)) < 0.5
                missingedges = missingedges + 1; 
            end
        end
    end
end
end