function [Ay, S1, gyface, Igyface, Ayface, Sfrac, Ayface_s, IAyface_s, Aysums, Aygysums, Icurr, fsize, IAysmax] = buildface(Ay, gy, s, epsc)
S1 = find(Ay > 1 - epsc);
Ay(S1) = 1;
Sfrac = find(Ay > 0 & Ay < 1);
fsize = size(Sfrac, 1);
[gyface, Igyface] = sort(gy(Sfrac));
Ayface = Ay(Sfrac(Igyface));
[Ayface_s, IAyface_s] = sort(Ayface);
Ayface2 = Ayface .* Ayface;
if fsize > 0
    Aysums = partials(flip(Ayface2));
    Aygysums = partials(flip(Ayface.*gyface));
else
    Aysums = [];
    Aygysums = [];
end
Icurr = 1 : (s - size(S1, 1));
IAysmax = 1;
if fsize > 0
    while IAyface_s(IAysmax) <= s - size(S1, 1)
        IAysmax = IAysmax + 1;
        if IAysmax > fsize
            Ayface = [];
            S1 = find(Ay > 0);
            Ay(S1) = 1;
            fsize = 0;
            gyface = [];
            Igyface = [];
            Sfrac = [];
            Aysums = [];
            Aygysums = [];
            Icurr = [];
            break;
        end
    end
end
if size(S1, 1) == s
    fsize = 0;
    gyface = [];
    Igyface = [];
    Sfrac = [];
    Aysums = [];
    Aygysums = [];
    Icurr = [];
end
end