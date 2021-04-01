function [ U,t,k,complex ] = ComparMethode(N,A,b,C,rho,eps,kmax)
    %%% Pas du tout fini
    
    
    U = zeros(3,N);
    t = zeros(3,1);
    k = zeros(3,1);
    complex = [];
    
    [U(1,:),~,k(1),~,t(1)] = DecompositionPrix(N4,A4,b4,C4,rho,eps,kmax);
    [U(2,:),~,k(2),~,t(2)] = DecompositionQuantites(N4,A4,b4,C4,eps,kmax);
    [U(3,:),~,k(3),~,t(3)] = DecompositionPrediction(N4,A4,b4,C4,eps,kmax);
    
    
    
end

