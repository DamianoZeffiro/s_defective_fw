function [Q, cliquename] = clique_init2(a)
    cliquenames = load('cliquenames.mat');
    cliquename = cliquenames.butenkoclique{a}; 
    E = graphcreator2(cliquename);
    n = max(E(:,1));
    [edges, ~] = size(E);
    Q = zeros(n);   
    for i = 1:edges
        Q(E(i, 1), E(i, 2)) = 1;
        Q(E(i, 2), E(i, 1)) = 1;
    end
    