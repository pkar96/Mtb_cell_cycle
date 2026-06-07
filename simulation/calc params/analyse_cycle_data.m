% analysing cell cycle data.

function [td, lb, ld, li, cd, gr, moth_dau] = analyse_cycle_data()

[num,txt,raw] = xlsread('../../christin_cellcycletiming_all_new.xlsx');
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

% find daughter cells if their B=0, check if the mother cell has e~=0
beq0 = find(B==0);
mutter= moth_dau(1,:);
tochter= moth_dau(2,:);
cdp=cd;
lip=li;
% indm=[];
for ind1=1:length(beq0)
    if b_or_e(beq0(ind1))==1
        indd = find(tochter==beq0(ind1));
        if ~isempty(indd)
%             indm =[indm; beq0(ind1) mutter(indd) E(mutter(indd)) ve(mutter(indd))];
            cdp(beq0(ind1))=cdp(beq0(ind1))+E(mutter(indd));
            lip(beq0(ind1))=ve(mutter(indd))/2;
        end
    elseif isnan(b_or_e(beq0(ind1))) || b_or_e(beq0(ind1))>1
        "problem"
    end
end
cd=cdp;
li=lip;

