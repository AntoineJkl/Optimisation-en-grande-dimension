function [ U,t,k,complex ] = ComparMethode(rho,eps,kmax)
    %%% Pas du tout fini
    
    % Comparaison evolution des solutions pour chaque méthode pour N=5
    N=5;
    [A,b,C] = CreateInstance(N);
    t = zeros(3,1);
    k = zeros(3,1);
    complex = [];
    
    [~,~,k(1),~,t(1),U1] = DecompositionPrix(N,A,b,C,rho,eps,kmax,true);
    [~,~,k(2),~,t(2),U2,~] = DecompositionQuantites(N,A,b,C,rho,eps,kmax,true,false);
    [~,~,k(3),~,t(3),U3] = DecompositionPrediction(N,A,b,C,eps,kmax,true);
    
    figure(1)
    subplot(1,3,1);
    plot(1:(k(1)-1),U1)
    legend('u1','u2','u3','u4','u5')
    xlabel('Nombre d itérations (k)')
    ylabel('Solutions u_i')
    title('Décomposition par les Prix')
    
    subplot(1,3,2);
    plot(1:(k(2)-1),U2,'--')
    legend('u1','u2','u3','u4','u5')
    xlabel('Nombre d itérations (k)')
    ylabel('Solutions u_i')
    title('Décomposition par les Quantités')
    
    subplot(1,3,3);
    plot(1:(k(3)-1),U3,'.')
    legend('u1','u2','u3','u4','u5')
    xlabel('Nombre d itérations (k)')
    ylabel('Solutions u_i')
    title('Décomposition par Prédiction')
    
    
    % Comparaison nb d'itérations/temps d'execution pour N variant de 5 à 50
    N_vec = 5:5:50;
    k_vec = zeros(3,length(N_vec));
    t_vec = zeros(3,length(N_vec));
    i=1;
    for n = N_vec
        [A_temp,b_temp,C_temp] = CreateInstance(n);
        [~,~,k_vec(1,i),~,t_vec(1,i),~] = DecompositionPrix(n,A_temp,b_temp,C_temp,rho,eps,kmax,true);
        [~,~,k_vec(2,i),~,t_vec(2,i),~,~] = DecompositionQuantites(n,A_temp,b_temp,C_temp,rho,eps,kmax,true,false);
        [~,~,k_vec(3,i),~,t_vec(3,i),~] = DecompositionPrediction(n,A_temp,b_temp,C_temp,eps,kmax,true);
        i = i+1;
    end
    
    figure(2)
    subplot(1,2,1)
    semilogy(N_vec,k_vec)
    subplot(1,2,2)
    semilogy(N_vec,t_vec)

    
    
end

