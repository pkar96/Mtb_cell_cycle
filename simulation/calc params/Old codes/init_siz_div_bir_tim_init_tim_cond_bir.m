% Concurrent process model to be simulated for multiple
% generations over a population. Initiation sizer and concurrent process -
% timer from birth and timer from initiation. More stringent conditions on
% initiation length i.e., it can happen after one of the concurrent 
% processes (from birth) is complete.

function [dist] = init_siz_div_bir_tim_init_tim_cond_bir(params, observed_data, num, num_med)

cells=[];

num_in=num_med;

% Input parameters
v_inp=2;
tau=num(num_in,1);% generation time in minutes
gr=num(num_in,2); % growth rate in min^-1
ngen=12; % number of generation
cd=params(1); % C+D period
bii= params(2); %num(num_in,23); % change in volume from initiation to initiation
tbd=params(3); % change in volume from birth to division
corbgr = num(num_in,14); mlb = num(num_in,18); stdlb = num(num_in,19)*mlb; % correlation between birth and growth rate; mean length at birth; std of length at birth

%-----Noise parameters
cvl= num(num_in,3)*gr; % std dev in growth rate
cvt= params(4)*num(num_in,23); %num(num_in,24)*bii; % std dev of vol additive initiation noise
sigcd= params(5)*num(num_in,4); % 0.35*num(num_in,4); % std dev of time additive division noise
cvtbd= params(6)*num(num_in,1); % 0.2*num(num_in,1); % std dev of vol additive division noise
cvr= num(num_in, 36); % std dev in r
cve=num(num_in,6); % std dev of experimental noise

cntr= 0; % counts number if times replication is limiting
tcnt = 0;
%-----For keeping track of events

for no_in=1:2
x= Cell(v_inp);
x.rate= gr;

x.oris = [1 1];
x.tNextDiv_in(1)= cd+randn()*sigcd;
while x.tNextDiv_in(1)<0
    x.tNextDiv_in(1)=cd+randn()*sigcd;
end

x.tNextDiv_dv = tbd+ randn()*cvtbd;
while x.tNextDiv_dv<0
    x.tNextDiv_dv = tbd+ randn()*cvtbd;
end

ii_in = bii+ randn()*cvt;
while ii_in<x.tNextDiv_dv*x.rate+x.v
    ii_in=bii+ randn()*cvt;
end
ii_in_1 = bii+ randn()*cvt;
while ii_in_1<x.tNextDiv_dv*x.rate+x.v
    ii_in_1=bii+ randn()*cvt;
end
x.vNextInit = [x.oris(1)*(ii_in) x.oris(2)*(ii_in_1)];

x.divr = 0.5 + randn()*cvr;

cells=[cells x];
end

%-----Outputs
vbs = []; % records volume at birth for the cell
vds = []; % records vol at division. Hence vds(1:end-1)= 2*vbs(2:end)
vis = []; % volume at initiation
vios= []; % volume at initiation per origin
tds = []; % generation time
cds = []; % initiation to division time
rts = []; % stores growth rate 
oriss=[]; % stores the number of origins just after division event
divrs=[]; % division ratio
moth_infos=[];% mother cell of the current cell

%-----Advance in time for %gens generations
gens = ngen*tau; % total time in mins
tStep = 0.2; % units of hr
simTime = 0;

ctit = 0;
flag = 0;

while simTime < gens
    %-----Step
    simTime = simTime + tStep;

    for ct=1:length(cells)
    x=cells(ct);
    x.t = x.t + tStep;
    x.v = x.v+tStep*x.rate; %grow linearly
    x.tNextDiv_in=x.tNextDiv_in-tStep;
    x.tNextDiv_dv=x.tNextDiv_dv-tStep;
    end    
    
    %-----Perform events

    cellsNew=[];
    for ct=1:length(cells)
    x=cells(ct);

    %-----Initiate
    %Upon accumulating enough volume
    ind = find(x.v - x.vNextInit>=0);
    if  ~isempty(ind)
        %-----Record event
        if length(x.vNextInit)==1
            x.vOfInits= [x.vOfInits x.v];
            x.tOfInits= [x.tOfInits x.t];
            x.oOfInits= [x.oOfInits 1];
            x.oris = [1 1];
            ii_in = bii+ randn()*cvt;
            while ii_in<x.vNextInit || ii_in<x.tNextDiv_dv*x.rate+x.v
                ii_in=bii+ randn()*cvt;
                ctit=ctit+1;
                if ctit>50
                    flag = 1;
                    break
                end
            end
            if flag == 1
                break
            else
                ctit = 0;
            end
            ii_in_1 = bii+ randn()*cvt;
            while ii_in_1<x.vNextInit || ii_in_1<x.tNextDiv_dv*x.rate+x.v
                ii_in_1=bii+ randn()*cvt;
                ctit=ctit+1;
                if ctit>50
                    flag = 1;
                    break
                end
            end
            if flag == 1
                break
            else
                ctit = 0;
            end
            x.vNextInit = [x.oris(1)*(ii_in) x.oris(2)*(ii_in_1)];
            xcd=cd+randn()*sigcd;
            while xcd<0
                xcd=cd+randn()*sigcd;
            end
            x.tNextDiv_in=[x.tNextDiv_in xcd];

        elseif length(x.vNextInit)==2
            x.vOfInits= [x.vOfInits x.v*ones(1,length(ind))];
            x.tOfInits= [x.tOfInits x.t*ones(1,length(ind))];
            x.oOfInits= [x.oOfInits 2*ones(1,length(ind))];
            if length(ind)==1
                x.oris = [1,2,2];
            elseif length(ind)==2
                x.oris = [2,2,2,2];
            else
                "Too many hits"
            end
            x.vNextInit(ind)=[];
            for i=1:length(ind)
                if length(ind)==2
                    ii_in = bii+ randn()*cvt;
                    while ii_in<bii/2 || ii_in<(x.v+x.tNextDiv_dv*x.rate) || ii_in<(x.v+x.tNextDiv_in(1)*x.rate)
                        ii_in=bii+ randn()*cvt;
                        ctit=ctit+1;
                        if ctit>50
                            flag = 1;
                            break
                        end
                    end
                    if flag == 1
                        break
                    else
                        ctit = 0;
                    end
                    ii_in_1 = bii+ randn()*cvt;
                    while ii_in_1<bii/2 || ii_in_1<(x.v+x.tNextDiv_dv*x.rate) || ii_in_1<(x.v+x.tNextDiv_in(1)*x.rate)
                        ii_in_1=bii+ randn()*cvt;
                        ctit=ctit+1;
                        if ctit>50
                            flag = 1;
                            break
                        end
                    end
                    if flag == 1
                        break
                    else
                        ctit = 0;
                    end
                    x.vNextInit = [x.vNextInit 2*(ii_in) 2*(ii_in_1)];
                elseif length(ind)==1
                    ii_in = bii+ randn()*cvt;
                    while ii_in*2<x.vNextInit(1) || ii_in<(x.v+x.tNextDiv_dv*x.rate) || ii_in<(x.v+x.tNextDiv_in(1)*x.rate)
                        ii_in=bii+ randn()*cvt;
                        ctit=ctit+1;
                        if ctit>50
                            flag = 1;
                            break
                        end
                    end
                    if flag == 1
                        break
                    else
                        ctit = 0;
                    end
                    ii_in_1 = bii+ randn()*cvt;
                    while ii_in_1*2<x.vNextInit(1) || ii_in_1<(x.v+x.tNextDiv_dv*x.rate) || ii_in_1<(x.v+x.tNextDiv_in(1)*x.rate)
                        ii_in_1=bii+ randn()*cvt;
                        ctit=ctit+1;
                        if ctit>50
                            flag = 1;
                            break
                        end
                    end
                    if flag == 1
                        break
                    else
                        ctit = 0;
                    end
                    x.vNextInit = [x.vNextInit 2*(ii_in) 2*(ii_in_1)];
                end
                xcd=cd+randn()*sigcd;
                while xcd<x.tNextDiv_in(1)
                    xcd=cd+randn()*sigcd;
                end
                x.tNextDiv_in=[x.tNextDiv_in xcd];
            end

        elseif length(x.vNextInit)==3
            x.vOfInits= [x.vOfInits x.v];
            x.tOfInits= [x.tOfInits x.t];
            x.oOfInits= [x.oOfInits 2];
            x.oris = [2,2,2,2];
            if ind~=1
                "x.vNextInit is not 1"
            end
            x.vNextInit(ind)=[];
            ii_in = bii+ randn()*cvt;
            while ii_in<bii/2 || ii_in<(x.v+x.tNextDiv_dv*x.rate) || ii_in<(x.v+x.tNextDiv_in(1)*x.rate)
                ii_in=bii+ randn()*cvt;
                ctit=ctit+1;
                if ctit>50
                    flag = 1;
                    break
                end
            end
            if flag == 1
                break
            else
                ctit = 0;
            end
            ii_in_1 = bii+ randn()*cvt;
            while ii_in_1<bii/2 || ii_in_1<(x.v+x.tNextDiv_dv*x.rate) || ii_in_1<(x.v+x.tNextDiv_in(1)*x.rate)
                ii_in_1=bii+ randn()*cvt;
                ctit=ctit+1;
                if ctit>50
                    flag = 1;
                    break
                end
            end
            if flag == 1
                break
            else
                ctit = 0;
            end
            x.vNextInit = [x.vNextInit 2*(ii_in) 2*(ii_in_1)];
            xcd=cd+randn()*sigcd;
            while xcd<x.tNextDiv_in(1)
                xcd=cd+randn()*sigcd;
            end
            x.tNextDiv_in=[x.tNextDiv_in xcd];
        else
            "Too many x.vNextInit"
        end
        
    end

    if flag == 1
        break
    end

    %-----Divide
    %Upon accumulating enough volume
    if ~isempty(x.tNextDiv_in) && x.tNextDiv_in(1)<=0 && x.tNextDiv_dv <=0

        % find no of cases in which replication is limiting 
        if x.tNextDiv_in>=-tStep
            cntr=cntr+1;
        end
        tcnt=tcnt+1;

        %-----Record event and output
        vds=[vds; (x.v+(round(x.t-x.tLastDiv)-(x.t-x.tLastDiv))*x.rate+randn()*2*cve)];
        tds=[tds; round(x.t-x.tLastDiv)];
        vbs=[vbs; (x.vb+(round(x.tLastDiv)-x.tLastDiv)*x.rate+randn()*cve)]; %with size additive measurement error
        if (x.t - x.tOfInits(1))<(x.t-x.tLastDiv+1) && (x.t - x.tOfInits(1))>(x.t-x.tLastDiv-1)
            vis=[vis; vbs(end)];
            vios=[vios; vbs(end)];
            cds=[cds; tds(end)];
        else
            vis=[vis; (x.vOfInits(1)+randn()*cve)];
            vios=[vios; vis(end)/x.oOfInits(1)+(round(x.t - x.tOfInits(1))-(x.t-x.tOfInits(1)))*x.rate+randn()*cve];    
            cds=[cds; round(x.t - x.tOfInits(1))];
        end
        rts=[rts; x.rate];
        if isnan(x.moth_info)
            divrs =[divrs; NaN];
            moth_infos=[moth_infos; NaN];
        else
            divrs =[divrs; vbs(end)/vds(x.moth_info)];
            moth_infos=[moth_infos; x.moth_info];
        end
%         conttime=[conttime, [vds(end); t; 2; oris/2], [vis(end); tOfInits(1); 1; oOfInits(1)]];
        
        x.vOfInits(1)=[];
        x.tOfInits(1)=[];
        x.oOfInits(1)=[];
        
        %-----Update cell
        x.tLastDiv = x.t;
        x.vd = x.v;
        x.v = x.vd*x.divr;
        x.vb = x.v;
        x.rate = max(gr + corbgr*(x.vb-mlb)/stdlb*cvl+ sqrt(1-corbgr^2)*randn()*cvl,gr/10);
        x.tNextDiv_in(1)=[];

        x.tNextDiv_dv = tbd+ randn()*cvtbd;
        while x.tNextDiv_dv<0
            x.tNextDiv_dv = tbd+ randn()*cvtbd;
        end

        %----------Update vNextInit and vNextDiv (see Methods)
        x.vNextInit = x.vNextInit/2;
        x.moth_info= length(vbs);

        %----------other daughter
        y=copy(x);
        y.v = y.v/y.divr*(1-y.divr);
        y.vb = y.v;
        y.rate = max(gr + corbgr*(y.vb-mlb)/stdlb*cvl+ sqrt(1-corbgr^2)*randn()*cvl,gr/10);
        y.divr = 0.5+ randn()*cvr;
        y.tNextDiv_dv = tbd+ randn()*cvtbd;
        while y.tNextDiv_dv<0
            y.tNextDiv_dv = tbd+ randn()*cvtbd;
        end

        if length(x.vNextInit)==2
            y.vNextInit = [x.vNextInit(1)];
            x.vNextInit(1) = [];
            x.oris = [1]; y.oris=[1];
        elseif length(x.vNextInit)==3
            y.vNextInit = [x.vNextInit(2) x.vNextInit(3)]; y.oris=[1,1];
            x.vNextInit(2:3) = []; x.oris=[1]; 
            x.vOfInits(1)=[]; x.tOfInits(1)=[]; x.oOfInits(1)=[]; x.tNextDiv_in(1) = [];
        elseif length(x.vNextInit)==4
            y.vNextInit = x.vNextInit(3:4); y.oris = [1,1];
            x.vNextInit(3:4) = []; x.oris = [1,1];
            x.vOfInits(2)=[]; x.tOfInits(2)=[]; x.oOfInits(2)=[]; x.tNextDiv_in(2) = [];
            y.vOfInits(1)=[]; y.tOfInits(1)=[]; y.oOfInits(1)=[]; y.tNextDiv_in(1) = [];
        elseif isempty(x.vNextInit)
            "Next Init is empty"
        else
            "Next Init has 1 or more than 4"
        end
        x.divr = 0.5+ randn()*cvr;
        cellsNew =[cellsNew y];
    end
    
    end
    if flag == 1
        break
    end
    cells = [cells cellsNew];
end

if flag == 1
    dist = 50000;
else
    simulated_data = [tds, vbs, vds, vios, cds, rts];
    mean_estimated = mean(simulated_data);
    cov_estimated = cov(simulated_data);

    dist = -sum(log(mvnpdf(observed_data, mean_estimated, cov_estimated)))*2;
end

