function [] = AffichageMultiplicateur(Multiplicateur,Multiplicateur_exact,titre)
    
    %Vecteur du nombre d'iterations:
    Iterations = 1:length(Multiplicateur);
    
    %Positionnement de la Figure:
    fig = figure();
    set(fig,'Position',[200 100 800 400]);
    
    %Affichage:
    plot(Iterations,Multiplicateur,'r','LineWidth',2);
    hold on;
    plot(Iterations,Multiplicateur_exact*ones(size(Iterations)),'LineStyle',':','Color','k');
    hold off;
    legend('Multiplicateur algorithme','Multiplicateur exact','Location','eastOutside');
    ylabel('Valeur du multiplicateur');
    xlabel('Iteration k');
    title(titre);

end

