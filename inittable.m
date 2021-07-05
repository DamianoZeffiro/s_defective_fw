function inittable(tname, totcliques)
    Max = zeros(totcliques, 1);
    Mean = zeros(totcliques, 1); 
    Std = zeros(totcliques, 1);
    s1 = table(Max, Mean, Std);
    s2 = table(Max, Mean, Std);
    s3 = table(Max, Mean, Std);
    s4 = table(Max, Mean, Std);    
    load('cliquenames.mat');
    Graph_names = butenkoclique;
    sizetable = table(Graph_names, s1, s2, s3, s4);
    s1 = table(Mean, Std);
    s2 = table(Mean, Std);
    s3 = table(Mean, Std);
    s4 = table(Mean, Std); 
    timestable = table(Graph_names, s1, s2, s3, s4);
    save(tname, 'sizetable', 'timestable');
end