function [x,it, ttot] = FDFW(Q, Ac, x, Ay, maxtime, eps, gammay, s, L)
% control variable for stopping criterion
flagls=0;
% sufficient decrease parameter
gamma = 0.002;
% rounding parameter for fractionary y components
epsc = 10^(-9);
it=1;
IAc = find(Ac == 1);
tic;
while flagls==0     
    if (it==1)
        % gradient and objective computation
        Qx = Q*x;
        Aysim = Ay + Ay';
        Ayx = Aysim * x;
        xx = x*x';
        xxAc = xx.*Ac;
        gx = 2 * Qx + 2 * Ayx;
        gy = 2 * xxAc + gammay * Ay;
        fx = x'*Qx + x'*Ayx + 0.5 * gammay * sum(Ay.^2, 'all');
    end
    % anchor point for the Short Step Chain
    ttot = toc;
    if ttot > maxtime
        break
    end
    xbar = x;
    Aybar = Ay;
    gzbar = dot(gx, x) + sum(gy.*Aybar, 'all'); 
    gz = gzbar;
    % computation of FW vertex
    [~, Ixfw] = max(gx);
    [~, Iyfw] = maxk(gy(IAc), s);
    Iyfw = IAc(Iyfw);
    gvfw = sum(gy(Iyfw), 'all') + gx(Ixfw); % (g, v_fw)
    fwgap = gvfw - gz;  
    findx = find(abs(x) > 0.000001);
    if (fwgap <= eps)
        if missingedgecount(Q, findx) <= s
            break;
        end
    end    
    % preprocessing of current face info to speed up the Short Step Chain 
    [Ay, S1, gyface, Igyface, Ayface, Sfrac, Ayface_s, IAyface_s, Aysums, Aygysums, Icurr, fsize, IAysmax] = buildface(Ay, gy, s, epsc);
    Ayfbar = Ayface;
    szbar0 = 0;
    totstep = 1;
    continuecycle = 1;
    Invertex = 0;
    % Short Step Chain
    while(continuecycle)
        indcc = find(x > 0e0);
        [~,istaraw] = min(gx(indcc));
        istaraw = indcc(istaraw);
        gvaw = sum(gyface(Icurr), 'all') + sum(gy(S1), 'all') + gx(istaraw);
        if gz - gvaw  > (gvfw - gz) 
            % in face direction computation
            gap = gz - gvaw;            
            dir = 2;
            dx = x;
            dx(istaraw) = dx(istaraw) - 1;
            alphaxm = x(istaraw)/(1 - x(istaraw));   
            if fsize > 0
                dycur = Ayface(Icurr) - 1;
                dzzbar = sum(dycur .* (Ayface(Icurr) - Ayfbar(Icurr)), 'all') + totstep * (totstep - 1) * Aysums(fsize - max(Icurr)) + dot(dx, x - xbar);
                nd = norm(dycur)^2 + totstep^2 * Aysums(fsize - max(Icurr)) + norm(dx)^2;                     
                alphaym = min(Ayface(Icurr)./(1 - Ayface(Icurr)));              
                alpha1m = (1 - Ayface_s(IAysmax) * totstep)/(Ayface_s(IAysmax) * totstep);
            else
                nd = norm(dx)^2;
                dzzbar = dot(dx, x - xbar);
                dycur = [];
                alphaym = 1;
                alpha1m = 1;
            end
            alpham = min([alphaxm, alphaym, alpha1m]);       
        else
            % fw direction computation
            gap = gvfw - gz;
            dir = 1;
            dx = - x;
            dx(Ixfw) = 1 + dx(Ixfw);
            if sum(Ayface, 'all') + size(S1, 1) > s + 0.05
                warning('ll');
            end
            Ay = embedfs(Ay, Ayface, Sfrac, Igyface, fsize, Icurr, totstep, s);           
            dy = - Ay;
            dy(Iyfw) = 1 + dy(Iyfw);
            dzzbar = dot(dx, x - xbar) + sum(dy.*(Ay - Aybar), 'all');
            nd = sum(dy.*dy, 'all') + norm(dx)^2;
            alpham = 1;
            continuecycle = 0;
        end
        zzbar = szbar0 + norm(Ayface(Icurr) - Ayfbar(Icurr))^2 + norm(x - xbar)^2;
        if max(Icurr) < fsize
            zzbar = zzbar + (totstep - 1)^2 * Aysums(fsize - max(Icurr));
        end
        gzzbar = gz - gzbar;
        av = L*nd;
        bv = 2*L*dzzbar - gap;
        cv = - gzzbar + L * zzbar;
        if sum(size(bv)) < 2
            warning('l');
        end
        if bv^2 < 0
            warning('w');
        end 
        beta = (- bv + sqrt(bv^2-4*av*cv))/(2*av);    
        cvbis = zzbar - gap^2/(nd*L^2);
        bvbis = 2*dzzbar;
        avbis = nd;
        if cvbis > 0
            beta = 0;
        else
            beta = min([beta,(- bvbis + sqrt(bvbis^2-4*avbis*cvbis))/(2*avbis) ]);
        end
        alpha = min([beta, alpham]);
        x = x + alpha * dx; 
        if dir == 2 
            totstep = totstep * (1 + alpha);
            if alpha == alphaym && fsize > 0
               [~, Ialphaym] = min(Ayface(Icurr));
               Ayface(Icurr) = Ayface(Icurr) + alpha * dycur;
               Ayface(Icurr(Ialphaym)) = 0;
               szbar0 = szbar0 + Ayfbar(Icurr(Ialphaym))^2;
               Icurr(Ialphaym) = max(Icurr) + 1;           
               Ayface(max(Icurr)) = totstep * Ayface(max(Icurr));                  
               while IAyface_s(IAysmax) <= max(Icurr)
                   IAysmax = IAysmax + 1;
                   if IAysmax > fsize                
                       Invertex = 1;
                       break
                   end                  
               end           
            elseif fsize > 0
                Ayface(Icurr) = Ayface(Icurr) + alpha * dycur;
            end       
        else  
            Ay = Ay + alpha * dy; 
        end
        if dir==2   && ~isempty(Icurr) && max(Icurr) < fsize   
            gz = dot(gx, x) + sum(gyface(Icurr) .* Ayface(Icurr), 'all') + totstep * Aygysums(fsize - max(Icurr)) + sum(gy(S1));
        elseif dir==2
            gz = dot(gx, x) + sum(gyface(Icurr) .* Ayface(Icurr), 'all') + sum(gy(S1));
        elseif dir==1
            gz = dot(gx, x) + sum(Ay.*gy, 'all');
        end     
        if beta < alpham || dir == 1 || alpha == alpha1m || alpha < 10^(-20) || Invertex == 1
            continuecycle = 0;
            gzzbar = gz - gzbar;
            if dir == 2                      
                Ay = embedfs(Ay, Ayface, Sfrac, Igyface, fsize, Icurr, totstep, s);          
            end
        end
    end
    % smart computation of new gradient
    Qx = Q*x;
    Aysim = Ay + Ay';
    Ayx = Aysim * x;
    xx2 = (2 * x)*x';
    xx2Ac = xx2.*Ac;
    gxnew = 2 * Qx + 2 * Ayx;
    gynew = xx2Ac + gammay * Ay;
    % new objective value
    fxnew = x'*Qx + x'*Ayx + 0.5 * gammay * sum(Ay.^2, 'all');
    % check sufficient decrease
    if fxnew <= fx + gamma * gzzbar && L < 8
        L = 2*L;
        Ay = Aybar;
        x = xbar;
    else
        fx = fxnew;
        gx = gxnew;
        gy = gynew;
    end 
    it = it+1;
end
ttot = toc;
end