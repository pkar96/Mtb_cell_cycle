% Analysing cell cycle data.

clear

[num,txt,raw] = xlsread('christin_cellcycletiming_all_new.xlsx');
% Various columns
mov= txt(2:end, 1);
cell_no=num(:,1);
pcell= num(:,2);
dcell_acc= num(:,3); 
dcell_alt= num(:,4);
ty= num(:,5);
lb= num(:, 6)*109.54/1024;
li= num(:, 7)*109.54/1024;
lt= num(:, 8)*109.54/1024;
ld= num(:, 9)*109.54/1024;
B= num(:, 10);
C= num(:, 11);
D= num(:, 12);
E= num(:, 13);
td= num(:,14);
ve= num(:,15)*109.54/1024;
b_or_e = num(:,16); % 0 if it means b = 0, 1 if it means e period. To be used to correct li for cells with lb=li and val =1
cd= C+D+E;
acc = find(ty==1);
alt = find(ty==0);
oth = find(ty==2);
gr= num(:,17);

% [length(cd), length(acc), length(alt), length(oth)]

% get daughter cell characteristics
moth_dau =[]; % 1st row stores the mother and second stores the daughter index
for i=1:length(lb)
    [row_acc, column, dataValues] = find(cell_no==dcell_acc(i) & strcmp(mov, mov{i}));
    if ~isempty(row_acc)
        moth_dau=[moth_dau [i; row_acc; 1]];
    end
    [row_alt, column, dataValues] = find(cell_no==dcell_alt(i) & strcmp(mov, mov{i}));
    if ~isempty(row_alt)
        moth_dau=[moth_dau [i; row_alt; 0]]; % stores mother index, daughter index and type of daughter 
    end
end

%%
% 
% % % find daughter cells if their B=0, check if the mother cell has e~=0
beq0 = find(B==0);
mutter= moth_dau(1,:);
tochter= moth_dau(2,:);
cdp=cd;
lip=li;
indm=[];
for ind1=1:length(beq0)
    if b_or_e(beq0(ind1))==1
        indd = find(tochter==beq0(ind1));
        if ~isempty(indd)
            indm =[indm; beq0(ind1) mutter(indd) E(mutter(indd)) ve(mutter(indd))];
            cdp(beq0(ind1))=cdp(beq0(ind1))+E(mutter(indd));
            lip(beq0(ind1))=ve(mutter(indd))/2;
        end
    elseif isnan(b_or_e(beq0(ind1))) || b_or_e(beq0(ind1))>1
        "problem"
    end
end
cd=cdp;
li=lip;


%%
% find cd<td for moth and dau cells
li_m_i=[]; li_d_i=[]; ld_m_i=[]; ld_d_i=[]; lt_m_i=[]; lb_m_i =[]; lb_d_i =[]; gr_m_i =[]; gr_d_i =[]; 
for i = 1:length(moth_dau(1,:))
    if lb(moth_dau(1,i))~=li(moth_dau(1,i)) && lb(moth_dau(2,i))~=li(moth_dau(2,i))
%     if lb(moth_dau(2,i))~=li(moth_dau(2,i)) %&& lt(moth_dau(1,i))~=ld(moth_dau(1,i))
        li_m_i = [li_m_i; li(moth_dau(1,i))];
        li_d_i = [li_d_i; li(moth_dau(2,i))*2];
        lb_m_i = [lb_m_i; lb(moth_dau(1,i))];
        lb_d_i = [lb_d_i; lb(moth_dau(2,i))];
        lt_m_i = [lt_m_i; lt(moth_dau(1,i))];
        ld_m_i = [ld_m_i; ld(moth_dau(1,i))];
        ld_d_i = [ld_d_i; ld(moth_dau(2,i))];
        gr_m_i = [gr_m_i; gr(moth_dau(1,i))];
        gr_d_i = [gr_d_i; gr(moth_dau(2,i))];
    end
end

% P=polyfit(ld_m_i, li_d_i*2, 1);
% (1-P(1))/(1-P(1)*0.08)

%%
% Plot binned data y vs x conditioned on condn
x = lb; y = li; condn= [];

if ~isempty(condn)
    condn=[condn ones(size(x))];
    
    by = regress(y,condn);
    bx = regress(x,condn);
    yfit=0; xfit=0;
    for icnt =1 : length(bx)
        yfit = by(icnt)*condn(:,icnt)+yfit;
        xfit = bx(icnt)*condn(:,icnt)+xfit;
    end
    
    ry= y-yfit;
    rx= x-xfit;
    
else
    ry= y;
    rx= x;
end

figure
d=scatter(rx,ry,60,'filled','MarkerFaceAlpha',0.6);
[bin,da,yfit, P]=binning(ry,rx);
hold on
scatter(bin,da, 100,'filled', 'MarkerEdgeColor', 'k');
plot(bin, P(1)*bin+P(2),'k', 'LineWidth',3)
yL=get(gca,'YLim');
xL=get(gca,'XLim');
[cor p]=corrcoef(rx,ry);
% p(1,2)
eq = sprintf('\\itr\\rm = %.2f, p=%.2d', cor(1,2),p(1,2));
% eq = sprintf('L_d = %.2f L_b + %.2f', P(1), P(2));
text((xL(1)+xL(2))/3,yL(2),eq,'HorizontalAlignment','left', 'VerticalAlignment','top','BackgroundColor',[1 1 1], 'FontSize',24);
ylabel('L_i (\mum)');
xlabel('L_b (\mum)');
box on
set(gca, 'FontSize', 28)
set(gcf, 'Position',[276,42,777,602])

%%

% find the percentage of cells with E period. How much % of these cells
% had daughter with B = 0. Cells with B = 0-> How much % of cells had an E
% period and how many of them had an E period greater than 1.

ind_e_n0 = find(E~=0); per_E = length(ind_e_n0)/length(B);
ind_dau_e_n0 = [];
for i =1:length(ind_e_n0)
    ind_dau_e_n0 = [ind_dau_e_n0 find(ind_e_n0(i)==moth_dau(1,:))];
end
ind_dau_e_n0 = ind_dau_e_n0';
dau_e_n0 = moth_dau(2,ind_dau_e_n0)'; % daughter cells of cells with E neq 0
% fraction of cells where mother has non-zero E and Li not equal to Lb
li_eq_lb = length(find(li(dau_e_n0)-lb(dau_e_n0)==0))/length(dau_e_n0) % not a good statistic

% find cells where B = 0 and mother cells has E 0/non zero 
beq0 = find(B==0);
ind_moth_beq0 =[];
for i=1:length(beq0)
    ind_moth_beq0 = [ind_moth_beq0; find(tochter==beq0(i))];
end
moth_beq0 = moth_dau(1,ind_moth_beq0)'; % daughter cells of cells with E neq 0
% fraction of cells where mother has non-zero E when B = 0
e_nz_b_0 = length(find(E(moth_beq0)~=0))/length(moth_beq0)

% fraction of cells where B not equal to 0
b_neq_0 = length(find(B(dau_e_n0)~=0))/length(dau_e_n0) 

