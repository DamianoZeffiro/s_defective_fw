function [x, it, ttot] = FWdc(Q, Ac, x, maxtime, eps, gammay, s, L)
tic;
[~,n] = size(Q);
flagls=0;
% values for the smart computation of the o.f.
it=1;
ssc = 1;
while flagls==0
    % gradient evaluation   
    xfind = find(x > 0);   
    xp = x(xfind);
    xsize = size(xp, 1);  
    if it == 1
        Iyfwold = zeros(1, s);
    end 
    [g, Iyfwold] = Qaugp(Q, xfind, Iyfwold, xsize, x, xp, s, n, Ac, gammay);
    [~,istar]=max(g);
    xstar=zeros(n,1);
    xstar(istar)=1.0;
    continuecycle = 1;
    xbar = x;
    ttot = toc;
    if (ttot > maxtime)
        break;
    end
    % Short Step Chain
    while(continuecycle)        
        indcc=find(x>0e0);
        [~,istaraw] = min(g(indcc));
        xstar2=zeros(n,1);
        indaw=indcc(istaraw);
        xstar2(indaw)=1.0;
        % in face step direction
        dAS = x - xstar2;
        % fw direction
        dFW = xstar - x;
        p1 = g' * dFW;
        findx = find(x > 0.000001);
        if (p1 <= eps)
            if missingedgecount(Q, findx) <= s
                flagls = 1;
                break;
            end
        end
        p2=g'*dAS;
        if (p2<= p1)
            d=dFW;
            gnr = p1;
            alpham=1.0;
            cdir=1;
        else
            d=dAS;
            gnr = p2;
            alpham=x(indaw)/(1-x(indaw));
            cdir=2;
        end
        nd = d'*d;
        gdeltax = g'*(x-xbar);
        ddeltax = d'*(x-xbar);
        xxbar = (x-xbar)'*(x-xbar);
        av = L*nd;
        bv = 2*L*ddeltax - gnr;
        cv = - gdeltax + L * xxbar;
        beta = (- bv + sqrt(bv^2-4*av*cv))/(2*av);
        cvbis = xxbar - gnr^2/(nd*L^2);
        bvbis = 2*ddeltax;
        avbis = nd;
        if cvbis > 0
            beta = 0;
        else
            beta = min([beta,(- bvbis + sqrt(bvbis^2-4*avbis*cvbis))/(2*avbis) ]);
        end
        continuecycle = 0;
        if beta >= alpham && ssc == 1 && cdir == 2
            x = x + alpham * d;
            x(indaw) = 0;
            continuecycle = 1;
        else
            alpha = min([beta, alpham]);
            x = x + alpha * d;
        end
    end
if flagls == 0
    it = it+1;
end
end
ttot = toc;
end