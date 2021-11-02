function [Solution,x] = PSOrepv19(parallel, experiment, plot, IC)
%This function runs a PSO optimization algorithm to parameterize the 
%3-gene model.

%INPUTS:
%- parallel: 0 or 1, for whether to run code in parallel.
%- experiment: 1, 2, or 3 for whether to fit exp 1, 2, or their average.
%- plot: 0 or 1, for whether to plot best cost function vs. iteration
%- IC: 0 or 1, for whether to use user-defined IC for particle swarm.

%Check if MATLAB is earlier than R2016b
if verLessThan('matlab', '9.1.0.441655')
    options = optimoptions(@particleswarm);
    options.Display = 'final';
    options.SwarmSize = 100; 
    options.MinFractionNeighbors = 0.1; 
    options.StallIterLimit = 200; 
    options.TolFun = 1e-6; 
    if parallel == true
        poolObj = parpool('local', 'AttachedFiles', ...
            {'repv19.m', 'heaviSideTrue.m', ...
            'getLambda_v19.m', 'getFluxv19.m', 'getCostZeng.m', ...
            'getCostv4.m', 'getCostThuv2.m', 'fv19.m', 'fv19_repv3.m'});
        options.UseParallel = true;
    end
else
    if plot == 1
        options = optimoptions(@particleswarm, 'PlotFcn', 'pswplotbestf');
    else
        options = optimoptions(@particleswarm);
    end
    options.Display = 'iter'; 
    options.SwarmSize = 100; 
    options.MinNeighborsFraction = 0.1; 
    options.MaxStallIterations = 200; 
    options.FunctionTolerance = 1e-6; 
    options.InertiaRange = [0.3, 1]; 
    if parallel == true
        poolObj = parpool('local', 'AttachedFiles', ...
            {'repv19.m', 'heaviSideTrue.m', 'getFluxv19.m', 'getCostZeng.m', ...
            'getCostv4.m', 'getCostThuv2.m', 'fv19.m', 'fv19_repv3.m'});
        options.UseParallel = true;
    end
end

%%
%DATA

%TY Infection Data Avg (cI857, 30C, O-)------------------------------------
%format: t cI cro cII ste_cI ste_cro ste_cII std_cI std_cro std_cII
load('033020_TYData.mat'); %load avg RNA data

if experiment == 1
    ty_m1 = ty_inf1_m1;
    ty_m2 = ty_inf1_m2;
    ty_m3 = ty_inf1_m3;
    ty_m4 = ty_inf1_m4;
    ty_m5 = ty_inf1_m5;
elseif experiment == 2
    ty_m1 = ty_inf2_m1;
    ty_m2 = ty_inf2_m2;
    ty_m3 = ty_inf2_m3;
    ty_m4 = ty_inf2_m4;
    ty_m5 = ty_inf2_m5;
else
    error('Exp. not correctly selected!');
end

%Thu qPCR Data (cI857 30C, WT)---------------------------------------------
%format: PCR(1:3).[time, mean, error]
%1: <MOI> = 0.4, 2: <MOI> = 3, 3: <MOI> = 8
%t: 0     5    10    20    30    45    60    75    90   120
load('042121_thu_qPCR.mat');

%Zeng's iScience 2018 Data (cI857, 30C, MOI = 1)---------------------------
%xdata: t, ydata: RNA
%data_cII_croP, data_cII_cro, data_cII_cI, data_cII_wt, data_cII_P
%data_cI_croP, data_cI_wt
load('032020_zeng.mat');

%Set Weights---------------------------------------------------------------
%P-
i_taus = 1:10; %t=0.5:60 min
dataMax = [
    max(ty_m1(i_taus, 2)), max(ty_m2(i_taus, 2)), max(ty_m3(i_taus, 2)), ...
    max(ty_m4(i_taus, 2)), max(ty_m5(i_taus, 2)), ...
    max(ty_m1(i_taus, 3)), max(ty_m2(i_taus, 3)), max(ty_m3(i_taus, 3)), ...
    max(ty_m4(i_taus, 3)), max(ty_m5(i_taus, 3)), ...
    max(ty_m1(i_taus, 4)), max(ty_m2(i_taus, 4)), max(ty_m3(i_taus, 4)), ...
    max(ty_m4(i_taus, 4)), max(ty_m5(i_taus, 4)), ...
    ];
normWeight = ones(15, length(i_taus))./dataMax';

%WT
i_taus_wt = 1:7;
dataMax_wt = log(max(PCR(1).mean(i_taus_wt))); %transform to log
normWeight_wt  = ones(length(i_taus_wt), 1)./dataMax_wt; %12
dW_wt = sqrt(8); %1
normWeight_wt = dW_wt.*normWeight_wt;

%Norm. for Zeng data relative to TY - MOI = 1 cI, cro, cII
normTYZeng = ones(1, 3); %1

%%
%PARAMS

%Lower Bounds
lb = [
            0.01;       %cI PRM basal prod
            5;         %ratio of cI PRM active / basal  
            20;         %ratio cI PRE/PRM basal prod (15)
            0.05;       %CI translation rate
            1;       %cro prod rate
            1.5;       %Cro translation rate (0.1, 0.5)
            1;       %cII prod rate
            1.5;       %CII translation rate (0.1, 0.5)
            1/5;      %lambda replication rate 
            
            1/80;   %kdil
            log(2)/7;    %cI RNA degr. rate (0.1)
            1/120;      %CI degr. rate (0)
            log(2)/7;    %cro degr. rate
            1/60;      %Cro degr. rate
            log(2)/7;    %cII degr. rate
            1/7;    %CII degr. rate
            
            2;      %nPRM,CI+
            6;      %nPRM,CI- 
            2;      %nPRM,Cro 
            2;      %nPRE
            2;      %nCro,Cro 
            2;      %nCro,CI 
            2;      %nCII,Cro 
            2;      %nCII,CI
            2;      %nM,Cro 
            2;      %nM,CI
            0;      %nDeg,CII
            
            10;     %KPRM,CI+
            10;     %KPRM,CI-
            10;     %KPRM,Cro
            10;     %KPRE
            10;     %KCro,Cro
            10;     %KCro,CI
            10;     %KCII,Cro
            10;     %KCII,CI
            10;     %KM,Cro
            10;     %KM,CI
            10;     %KDeg,CII
            
            0;     %Dt
            0;     %tau
     ];
 

%Upper Bounds
ub = [
            0.1;    %cI PRM basal prod 
            10;     %ratio of cI PRM active/basal 
            30;     %ratio cI PRE/PRM prod rate
            1;      %CI translation rate 
            10;     %cro prod rate
            4.5;    %Cro translation rate 
            10;     %cII prod rate
            4.5;    %CII translation rate 
            1/3;    %lambda replication rate
       
            1/40;   %kdil
            1;      %cI RNA degr rate
            1/60;   %CI degr. rate
            1;      %cro RNA degr rate
            1/30;   %Cro degr. rate
            1;      %cII RNA degr rate
            1;      %CII degradation
            
            4;      %nPRM,CI+
            12;     %nPRM,CI-
            4;      %nPRM,Cro
            4;      %nPRE
            4;      %nCro,Cro
            6;      %nCro,CI
            4;      %nCII,Cro
            6;      %nCII,CI
            6;      %nM,Cro
            6;      %nM,CI
            10;     %nDeg,CII 
            
            100;     %KPRM,CI+ 
            200;     %KPRM,CI- 
            200;     %KPRM,Cro
            400;     %KPRE 
            400;     %KCro,Cro
            200;     %KCro,CI
            400;     %KCII,Cro
            200;     %KCII,CI
            1e3;     %KM,Cro
            250;     %KM,CI 
            1e3;     %KDeg,CII
            
            5;    %Dt 
            8;    %replication offset
     ];
 
%Save bounds 
d = clock;
fileDate = strcat(num2str(d(3)), num2str(d(2)), num2str(d(1)));
bounds = [lb ub]; 
save(strcat('PSO_bounds', num2str(fileDate), '.txt'), 'bounds', ...
    '-ascii');

%Run PSO code nrep times
nrep = 100;
n_param=length(lb);

%Set initial distribution
if IC == 1
    load('PSO_IC.mat', 'x0');
    options.InitialSwarmMatrix = x0;
end

%%
%Run PSO
nsol = 1;
Solution=ones(nrep,n_param+1)*1e10; 
for j=1:nrep
    [x,fval,exitflag]=particleswarm(@(params)repv19(params, ty_m1, ...
        ty_m2, ty_m3, ty_m4, ty_m5, normWeight, normWeight_wt, normTYZeng, PCR, ...
        data_cII_P, data_cII_wt, data_cII_croP, data_cII_cro, ...
        data_cII_cI, data_cI_croP, data_cI_wt),...
        n_param,lb,ub, options); 
    if fval<Solution(nsol,n_param+1)
        ['Solution ' num2str(nsol) ', good fval is: ' num2str(fval)]
        x_fval=cat(2,x, fval);
        Solution(nsol,:)=x_fval;
        Solution=sortrows(Solution,n_param+1);
        save(strcat('exp', num2str(experiment), '_fits.txt'), ...
            'Solution', '-ascii');
        nsol = nsol + 1; 
    end
end

%If a parallel pool object exists, terminate.
if exist('poolObj')
    delete(poolObj);
end



