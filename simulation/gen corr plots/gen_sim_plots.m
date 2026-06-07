% uses values obtained from simulation/experiments and plots the
% correlation coefficients. model_n stores the model being simulated.

function []= gen_sim_plots(td, lb, ld, li, cd, gr, moth_dau, model_n)

if model_n~= 0
    mutter = moth_dau(3:length(moth_dau));
    tochter= 3:length(moth_dau);
    moth_dau = [mutter'; tochter];
end
[Rlbd, errlbd]=calc_err(lb, ld);
[Rtdb, errtdb]=calc_err(lb, td);
[Rcdli,errcdli]=calc_err(li,cd);
[Rcdtd,errcdtd]=calc_err(cd, td);
[Ridi, erridi]=calc_err(li, ld-li);
[Rbii, errbii]=calc_err(li, lb);
li_m_i=[]; li_d_i=[]; ld_m_i=[]; lb_m_i =[]; lb_d_i =[];
for i = 1:length(moth_dau(1,:))
    if lb(moth_dau(1,i))~=li(moth_dau(1,i)) && lb(moth_dau(2,i))~=li(moth_dau(2,i))
        li_m_i = [li_m_i; li(moth_dau(1,i))];
        li_d_i = [li_d_i; li(moth_dau(2,i))*2];
        lb_m_i = [lb_m_i; lb(moth_dau(1,i))];
        lb_d_i = [lb_d_i; lb(moth_dau(2,i))];
        ld_m_i = [ld_m_i; ld(moth_dau(1,i))];
    end
end
if model_n == 0
    [Rbiii,errbiii]=calc_err(li_m_i(li_m_i<5),li_d_i(li_m_i<5));
else
    [Rbiii,errbiii]=calc_err(li_m_i,li_d_i);
end
[Rbd_ig,errbd_ig]=cond_xy(lb,ld,[li,gr]);
[Rid_b ,errid_b]=cond_xy(li,ld,lb);
[Rgd_b,errgd_b]=cond_xy(gr,ld,lb);

    y=[1,2,3,4,5,6,7, 8, 9, 10 ]; % each y corresponds to corrl coeff
    x=[Rlbd Rtdb Rcdli Rcdtd Ridi Rbii Rbiii Rbd_ig Rid_b Rgd_b]; 
    errbr=[errlbd; errtdb; errcdli; errcdtd; erridi; errbii; errbiii; errbd_ig; errid_b; errgd_b];
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
    xlabel('correlation coefficients', 'FontSize',15);
    yticks([1 2 3 4 5 6 7 8 9 10]);
    yticklabels({'(L_b,L_d)','(L_b,T_d)','(L_i,T_{id})','(T_{id},T_d)','(L_i,\DeltaL_{id})','(L_i,L_{b})', '(L_i,L_{i+1})', '(L_b,L_d|(L_i,\lambda_{lin})', '(L_i,L_{d}|L_b)', '(\lambda_{lin},L_d|L_b)' })
    if model_n==0
        eb = errorbar(x,y,errbr(:,2), errbr(:,1), 'horizontal', 'LineStyle', 'none', 'LineWidth', 3, 'Color', [0 0.4470 0.7410]);
        set(get(get(eb,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    end
    box on
    ylim([0.5 10.5])
end
