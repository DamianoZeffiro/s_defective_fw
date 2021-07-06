% test the FDFW on a list of graphs from the 2nd DIMACS optimization challenge.
% see butenkoclique in cliquenames.mat for a list of 
% graph names. cnumV must be given as a list of indices
function main(cnumV)
% total number of cliques
numclique = 50;
% time limit for each instance
maxtime = 600;
% max number of runs for each instance, in case time limit is not reached
nrunmax = 100;
% vector with possible number of fake edges
sVector = [1, 2, 3, 4];
% Lipschitz constant starting estimate
L = 1;
% 1 for FDFW, 2 for FWdc
algorithm = 1;
algorithmV = ["FDFW", "FWdc"];
savenum = 1;
% comment if you don't want to reset results
cliquesizes = zeros(numclique, nrunmax, size(sVector, 2));
cputimes = zeros(numclique, nrunmax, size(sVector, 2));
save(strcat("cm", string(savenum), ".mat"), 'cliquesizes', 'cputimes');
% uncomment if you want to load results
% load(strcat("cm", string(savenum), ".mat"), 'cliquesizes', 'cputimes');
tname = strcat(algorithmV(algorithm), "table", string(savenum), ".mat");
% initialize or reset table with results 
inittable(tname, numclique)
for sI = 1:4
    s = sVector(sI);
    for i = 1:size(cnumV, 2)
        % number of completed runs
        nrun = 0;
        % control variable for the time limit
        breakd = 0;
        cnum = cnumV(i);
        %  [adjacency matrix, clique name]
        [A, sN] = clique_init2(cnum);
        n = size(A, 2);        
        % stopping condition: FW gap <= eps and support is a clique
        eps = 2 * 1e-3;
        % regularization coefficients
        gammax = 1;
        gammay = 2/(n^2);
        % tot time for the instance
        Ittot = 0;
        for krun = 1:nrunmax
            %%%%% Gcleneration of the instance: %%%%%
            % seed changes at every run to generate starting point
            rng(krun);
            x0 = rand(n,1);
            x0 = x0/(sum(x0));
            [~,n] = size(A);
            Q = A + 0.5 * gammax * eye(n); 
            filterU = 1 - tril(ones(n));
            Ac = (1 - Q).*filterU;            
            if algorithm == 1
                disp('******************************');
                % initialization of a random y component
                rng(krun);
                Ay = rand(n).*(1 - A);
                Ay = Ay.*filterU;
                Ay = Ay*s/sum(Ay, 'all');
                % [solution found, number of iterations, cputime]
                [x, itf,  cputimes(cnum, krun, sI)] = FDFW(Q, Ac, x0, Ay, maxtime, eps, gammay, s, L);
            else
                % no random y component needed
                [x, itf,  cputimes(cnum, krun, sI)] = FWdc(Q, Ac, x0, maxtime, eps, gammay, s, L);
            end
            cliquesizes(cnum, krun, sI) =  sum((abs(x)>=0.000001));
            % Print results:%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            fprintf(1, 'Missing edges: %d\n',   missingedgecount(A, find(abs(x) >=0.000001)));
            fprintf(1,'Number of non-zero components of x = %d\n',...
                cliquesizes(cnum, krun, sI));
            fprintf(1,'Number of iterations = %d\n',...
                itf);
            fprintf(1,'CPU time = %10.3e\n', cputimes(cnum, krun, sI));
            disp(sN);
            Ittot = Ittot + cputimes(cnum, krun, sI);
            if Ittot > maxtime
                nrun = krun - 1;
                breakd = 1;
                break
            end               
        end
        if breakd == 0
            nrun = krun;
        end
        save(strcat(algorithmV(algorithm), string(savenum), ".mat"), 'cliquesizes', 'cputimes');       
        updatetable(tname, cnum, nrun, s, cputimes, cliquesizes);     
    end
end
load(tname, 'sizetable', 'timestable');
disp(sizetable);
disp(timestable);
end