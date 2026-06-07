% uses values obtained from simulation/experiments and plots the
% coefficient of variation. model_n stores the model being simulated.

function []= cv_plots(td, lb, ld, li, cd, gr, model_n)

    [cvlb, errlb]=calc_cv(lb);
    [cvld, errld]=calc_cv(ld);
    [cvtd, errtd]=calc_cv(td);
    [cvli,errli]=calc_cv(li);
    [cvcd, errcd]=calc_cv(cd);
    [cvgr, errgr]=calc_cv(gr);

    y=[1,2,3,4,5,6]; % each y corresponds to CV values
    x=[cvgr cvcd cvli cvtd cvld cvlb];
    errbr=[errgr; errcd; errli; errtd; errld; errlb];
    if model_n == 1 % this is CH model
        d=scatter(x, y, 140, 'filled', '^');
        d.MarkerFaceColor = 'g';
        d.MarkerEdgeColor = 'k';
    elseif model_n == 2 % parallel adder
        d=scatter(x, y, 140, 'filled', 's');
        d.MarkerFaceColor = 'r';
        d.MarkerEdgeColor = 'k';
    elseif model_n == 9 % concurrent (bir add init tim)
        d=scatter(x, y, 140, 'filled');
        d.MarkerFaceColor = 'm';
        d.MarkerEdgeColor = 'k';
    else
        d=scatter(x, y, 140, 'filled', 'd');
        d.MarkerFaceColor = [0 0.4470 0.7410];
        d.MarkerEdgeColor = 'k';
    end
    hold on
    a = get(gca,'YTickLabel');
    set(gca,'YTickLabel',a,'fontsize',15)
    xlabel('coefficient of variation', 'FontSize',15);
    yticks([1 2 3 4 5 6])
    yticklabels({'\lambda', 'T_{id}', 'L_i' , 'T_d', 'L_d', 'L_b'})
    if model_n==0
        eb = errorbar(x,y,errbr, 'horizontal', 'LineStyle', 'none', 'LineWidth', 3, 'Color', [0 0.4470 0.7410]);
        set(get(get(eb,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    end
    box on

    ylim([0.5 6.5])
end
