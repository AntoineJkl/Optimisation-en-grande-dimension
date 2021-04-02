function [ U,t,k,complex ] = ComparMethode(rho,eps,kmax,plotNpetit,plotNbIt,plotTmpsEx,plotComplex)
    %%% Pas du tout fini
    
    if plotNpetit
        % Comparaison evolution des solutions pour chaque méthode pour N=5
        N=5;
        [A,b,C] = CreateInstance(N);
        t = zeros(3,1);
        k = zeros(3,1);
        complex = [];

        [~,~,k(1),~,t(1),U1] = DecompositionPrix(N,A,b,C,rho,eps,kmax,true);
        [~,~,k(2),~,t(2),U2] = DecompositionQuantites(N,A,b,C,rho,eps,kmax,true);
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
    end
    
    if plotNbIt || plotTmpsEx || plotComplex
        % Comparaison nb d'itérations/temps d'execution pour N variant de 5 à 50
        N_vec = 5:5:50;
        k_vec = zeros(3,length(N_vec));
        t_vec = zeros(3,length(N_vec));
        erreur=zeros(3,length(N_vec));
        i=1;
        
        options = optimoptions('fmincon','Display','off');
        
        for n = N_vec
            [A_temp,b_temp,C_temp] = CreateInstance(n);
            [u1,~,k_vec(1,i),~,t_vec(1,i),~] = DecompositionPrix(n,A_temp,b_temp,C_temp,rho,eps,kmax,true);
            [u2,~,k_vec(2,i),~,t_vec(2,i),~] = DecompositionQuantites(n,A_temp,b_temp,C_temp,rho,eps,kmax,true);
            [u3,~,k_vec(3,i),~,t_vec(3,i),~] = DecompositionPrediction(n,A_temp,b_temp,C_temp,eps,kmax,true);
            
            if plotComplex
                f_obj = @(u) 1/2*u'*A_temp*u - b_temp'*u;
                sol_exa = fmincon(f_obj,zeros(n,1),C_temp,zeros(n,1),[],[],[],[],[],options);
                erreur(1,i) = norm(sol_exa-u1);
                erreur(2,i) = norm(sol_exa-u2);
                erreur(3,i) = norm(sol_exa-u3);
            end
            i = i+1;
        end

        figure(2)
        if plotNbIt
            subplot(1,2,1)
            plot(N_vec,k_vec)
            title('Nombre d itérations')
            legend('Prix','Quantités','Prédiction')
            xlabel('Taille N')
            ylabel('Nombre d itérations (k)')
        end
        if plotTmpsEx
            subplot(1,2,2)
            plot(N_vec,t_vec)
            title('Temps d execution')
            legend('Prix','Quantités','Prédiction')
            xlabel('Taille N')
            ylabel('Temps d execution (en s)')
        end
        if plotComplex
            figure(3)
            plot(N_vec,erreur)
            title('Erreur')
            legend('Prix','Quantités','Prédiction')
            xlabel('Taille N')
            ylabel('Erreur')
        end
    end

    
    
end

