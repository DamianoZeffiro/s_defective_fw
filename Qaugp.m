function [g, Iyfwold] = Qaugp(Q, xfind, Iyfwold, sizex, x, xp, s, n, Ac, gammay)
xfindlin = find(Ac(xfind, xfind) > 0);
g = 2 * (Q * x);
Iyfw = zeros(s, 1);
Vyfw = zeros(s, 1);
Vyfwold = zeros(s, 1);
for ij = 1:size(xfindlin, 1)
    ijxf = xfindlin(ij);
    ijxfrem = rem(ijxf - 1, sizex) + 1;
    ijxffix = fix(ijxf / sizex);
    i = ijxfrem;
    j = ijxffix + 1;
%    [i,j] = ind2sub([sizex, sizex], xfindlin(ij));
    if Q(xfind(i), xfind(j)) == 0
        wy = 2 * xp(i) * xp(j);
        if wy > Vyfw(s)
            Vyfw(s) = wy;
            Iyfw(s) = sub2ind([n,n], xfind(i), xfind(j));
            [Vyfw, I] = sort(Vyfw, 'descend');
            Iyfw = Iyfw(I);
        end
    end  
end

for i = 1:s
    if Iyfwold(i) > 0        
       [j,k] = ind2sub([n, n], Iyfwold(i));
       Vyfwold(i) = 2 * x(j) * x(k) + gammay;
    else
        Vyfwold(i) = 0;
    end
end
[Vyfwold, vyO] = sort(Vyfwold, 'descend');
Iyfwold = Iyfwold(vyO);
Iyfwold =  merge(Iyfw, Vyfw, Iyfwold, Vyfwold, s);  
for i = 1:s
    if Iyfwold(i) > 0
        [a1, b1] = ind2sub([n,n], Iyfwold(i));
        g(a1) = g(a1) + 2 * x(b1);
        g(b1) = g(b1) + 2 * x(a1);
    end
end    
end