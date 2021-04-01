function [] = AffichageSolution(P_tab,P_exact)

    Iterations = 1:length(P_tab(1,:));
    plots_tab = [];
    
    fig = figure(1);
    set(fig,'Position',[200 100 800 400]);
    
    %Centrale à charbon:
    plots_tab(1) = plot(Iterations,P_tab(1,:),'LineStyle','-','Color',[224 63 29]/255);
    hold on;
    plot(Iterations,P_exact(1)*ones(size(Iterations)),'LineStyle',':','Color','k');
    hold on;
    
    %Eolienne:
    for i=2:6
        plots_tab(i) = plot(Iterations,P_tab(i,:),'LineStyle','-','Color',[14 255 10]/255);
        hold on;
        plot(Iterations,P_exact(i)*ones(size(Iterations)),'LineStyle',':','Color','k');
        hold on;
    end
    
    %Barrage:
    plots_tab(7) = plot(Iterations,P_tab(7,:),'LineStyle','-','Color',[109 3 210]/255);
    hold on;
    plot(Iterations,P_exact(7)*ones(size(Iterations)),'LineStyle',':','Color','k');
    hold on;
    
    %Panneaux photovoltaiques:
    for i=8:9
        plots_tab(i) = plot(Iterations,P_tab(i,:),'LineStyle','-','Color',[10 177 0]/255);
        hold on;
        plot(Iterations,P_exact(i)*ones(size(Iterations)),'LineStyle',':','Color','k');
        hold on;
    end
    
    %DataCenter:
    plots_tab(10) = plot(Iterations,P_tab(10,:),'LineStyle','-','Color',[4 121 255]/255);
    hold on;
    plot(Iterations,P_exact(10)*ones(size(Iterations)),'LineStyle',':','Color','k');
    hold on;
    
    %Logements:
    plots_tab(11) = plot(Iterations,P_tab(11,:),'LineStyle','-','Color',[255 4 250]/255);
    hold on;
    plot(Iterations,P_exact(11)*ones(size(Iterations)),'LineStyle',':','Color','k');
    hold on;
    
    %Usine:
    plots_tab(12) = plot(Iterations,P_tab(12,:),'LineStyle','-','Color',[255 254 4]/255);
    hold on;
    plot(Iterations,P_exact(12)*ones(size(Iterations)),'LineStyle',':','Color','k');
    hold on;
    
    %Tramways:
    for i = 13:14
        plots_tab(i) = plot(Iterations,P_tab(i,:),'LineStyle','-','Color',[255 0 0]/255);
        hold on;
        plot(Iterations,P_exact(12)*ones(size(Iterations)),'LineStyle',':','Color','k');
        hold on;
    end
    
    %Hopital:
    plots_tab(15) = plot(Iterations,P_tab(15,:),'LineStyle','-','Color',[255 4 203]/255);
    hold on;
    exact = plot(Iterations,P_exact(15)*ones(size(Iterations)),'LineStyle',':','Color','k');
    hold off;
    
    %Affichage:
    legend([plots_tab([1,2,7,8,10,11,12,13,15]),exact],{'Centrale à charbon','Eolienne','Barrage','Panneau photovoltaique','DataCenter','Logement','Usine','Tramway','Hopital','Solution exacte'},'Location','eastOutside');
    xlim([1 length(Iterations)]);
    ylim([-14 14]);
    ylabel('Puissance fournie (MW)');
    xlabel('Iteration k');
    title({'Convergence des solutions','pour chaque agent'});
end

