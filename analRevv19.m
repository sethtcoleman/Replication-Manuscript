%This script takes the fitted parameter sets and performs any necessary
%additional analysis necessary to produce revision figures

%This model includes:
%- Cro, CII, CI (RNA & protein)
%- Lambda 

%%
%GET DATA

experiment = 2;

%TY Infection Data Avg (cI857, 30C, O-)------------------------------------
%format: t cI cro cII ste_cI ste_cro ste_cII std_cI std_cro std_cII
load('033020_TYData.mat'); %load avg RNA data

%Normalize data (divide by the max in each column)
ty_inf1_m1_norm = getNormData(ty_inf1_m1);
ty_inf1_m2_norm = getNormData(ty_inf1_m2);
ty_inf1_m3_norm = getNormData(ty_inf1_m3);
ty_inf1_m4_norm = getNormData(ty_inf1_m4);
ty_inf1_m5_norm = getNormData(ty_inf1_m5);

ty_inf2_m1_norm = getNormData(ty_inf2_m1);
ty_inf2_m2_norm = getNormData(ty_inf2_m2);
ty_inf2_m3_norm = getNormData(ty_inf2_m3);
ty_inf2_m4_norm = getNormData(ty_inf2_m4);
ty_inf2_m5_norm = getNormData(ty_inf2_m5);

if experiment == 1
    ty_m1 = ty_inf1_m1;
    ty_m2 = ty_inf1_m2;
    ty_m3 = ty_inf1_m3;
    ty_m4 = ty_inf1_m4;
    ty_m5 = ty_inf1_m5;
    ty_AvgEar_m1 = ty_inf1_AvgEar_m1;
    ty_AvgEar_m2 = ty_inf1_AvgEar_m2;
    ty_AvgEar_m3 = ty_inf1_AvgEar_m3;
    ty_AvgEar_m4 = ty_inf1_AvgEar_m4;
    ty_AvgEar_m5 = ty_inf1_AvgEar_m5;
elseif experiment == 2
    ty_m1 = ty_inf2_m1;
    ty_m2 = ty_inf2_m2;
    ty_m3 = ty_inf2_m3;
    ty_m4 = ty_inf2_m4;
    ty_m5 = ty_inf2_m5;
    ty_AvgEar_m1 = ty_inf2_AvgEar_m1;
    ty_AvgEar_m2 = ty_inf2_AvgEar_m2;
    ty_AvgEar_m3 = ty_inf2_AvgEar_m3;
    ty_AvgEar_m4 = ty_inf2_AvgEar_m4;
    ty_AvgEar_m5 = ty_inf2_AvgEar_m5;
elseif experiment == 3
    ty_m1 = ty_infAvg_m1;
    ty_m2 = ty_infAvg_m2;
    ty_m3 = ty_infAvg_m3;
    ty_m4 = ty_infAvg_m4;
    ty_m5 = ty_infAvg_m5;
    ty_AvgEar_m1 = ty_infAvg_AvgEar_m1;
    ty_AvgEar_m2 = ty_infAvg_AvgEar_m2;
    ty_AvgEar_m3 = ty_infAvg_AvgEar_m3;
    ty_AvgEar_m4 = ty_infAvg_AvgEar_m4;
    ty_AvgEar_m5 = ty_infAvg_AvgEar_m5;
else
    error('Exp. not correctly selected!');
end;

%Thu qPCR Data (cI857 30C, WT)---------------------------------------------
%format: PCR(1:3).[time, mean, error]
%1: <MOI> = 0.4, 2: <MOI> = 3, 3: <MOI> = 8
%t: 0     5    10    20    30    45    60    75    90   120
load('042121_thu_qPCR.mat');

%%
%GET SOLUTIONS

filenames = {
    'exp2_fits.txt';
    %'exp1_fits.txt';
};

ncol = 41;

bigNum = 1e5;
data = zeros(bigNum, ncol);
for i = 1:length(filenames)
    fileID = fopen(filenames{i}, 'r');
    formatSpec = '%f'; %floating point numbers
    sizeA = [ncol Inf];
    A = fscanf(fileID, formatSpec, sizeA);
    data( (1 + (i-1)*size(A, 2)):(i*size(A, 2)), :) = A.';
    fclose(fileID);
end;

%Delete zero rows, keep only solutions
failNum = 1e10;
data(~any(data,2), :) = [];
data = data(~all(data == failNum, 2), :);

%Sort solutions by last column
data = sortrows(data, ncol);

%Keep solutions only up to cutoff
cutOff = 1/3;
numKeep = round(cutOff*size(data, 1));
data = data(1:numKeep, :);

%%
%GET DATA

%TY Infection Data Avg (cI857, 30C, O-)------------------------------------
%format: t cI cro cII ste_cI ste_cro ste_cII std_cI std_cro std_cII
load('033020_TYData.mat'); %load avg RNA data

experiment = 2;

if experiment == 1
    ty_m1 = ty_inf1_m1;
    ty_m2 = ty_inf1_m2;
    ty_m3 = ty_inf1_m3;
    ty_m4 = ty_inf1_m4;
    ty_m5 = ty_inf1_m5;
    ty_AvgEar_m1 = ty_inf1_AvgEar_m1;
    ty_AvgEar_m2 = ty_inf1_AvgEar_m2;
    ty_AvgEar_m3 = ty_inf1_AvgEar_m3;
    ty_AvgEar_m4 = ty_inf1_AvgEar_m4;
    ty_AvgEar_m5 = ty_inf1_AvgEar_m5;
elseif experiment == 2
    ty_m1 = ty_inf2_m1;
    ty_m2 = ty_inf2_m2;
    ty_m3 = ty_inf2_m3;
    ty_m4 = ty_inf2_m4;
    ty_m5 = ty_inf2_m5;
    ty_AvgEar_m1 = ty_inf2_AvgEar_m1;
    ty_AvgEar_m2 = ty_inf2_AvgEar_m2;
    ty_AvgEar_m3 = ty_inf2_AvgEar_m3;
    ty_AvgEar_m4 = ty_inf2_AvgEar_m4;
    ty_AvgEar_m5 = ty_inf2_AvgEar_m5;
elseif experiment == 3
    ty_m1 = ty_infAvg_m1;
    ty_m2 = ty_infAvg_m2;
    ty_m3 = ty_infAvg_m3;
    ty_m4 = ty_infAvg_m4;
    ty_m5 = ty_infAvg_m5;
    ty_AvgEar_m1 = ty_infAvg_AvgEar_m1;
    ty_AvgEar_m2 = ty_infAvg_AvgEar_m2;
    ty_AvgEar_m3 = ty_infAvg_AvgEar_m3;
    ty_AvgEar_m4 = ty_infAvg_AvgEar_m4;
    ty_AvgEar_m5 = ty_infAvg_AvgEar_m5;
else
    error('Exp. not correctly selected!');
end;

%Thu qPCR Data (cI857 30C, WT)---------------------------------------------
%format: PCR(1:3).[time, mean, error]
%1: <MOI> = 0.4, 2: <MOI> = 3, 3: <MOI> = 8
%t: 0     5    10    20    30    45    60    75    90   120
%load('122120_thu_qPCR.mat');
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
normTYZeng = ones(1, 3);

%%
%GET TRAJECTORIES

%Create data structure for each solution - parameters, protein & RNA trajectories (MOI
%1:5), norm. RNA trajectories (MOI 1:5), cost functions (MOI 1:5)

sol_Ind = 1:size(data, 1);
checkSS = zeros(length(sol_Ind), 1);
ssCut = 6;
checkDec = zeros(length(sol_Ind), 1); %1 = success
checkLyt1 = zeros(length(sol_Ind), 1); %1 = fail
checkLyt2 = zeros(length(sol_Ind), 1); %1 = fail
checkPRE = zeros(length(sol_Ind), 1); %1 = true
checkSS = zeros(length(sol_Ind), 1);
checkLysogeny = zeros(length(sol_Ind), 1);

sol(sol_Ind) = struct('parameters', [], ...
    'prod', [], 'degr', [], 'n', [], 'K', [], 'tau', [], 'dt', [], ...
    'm1', [], 'm2', [], 'm3', [], 'm4', [], 'm5', [], ...
    'm1Num', [], 'm2Num', [], 'm3Num', [], 'm4Num', [], 'm5Num', [], ...
    'm1_wt', [], 'm2_wt', [], 'm3_wt', [], 'm4_wt', [], 'm5_wt', [], ...
    'm1Num_wt', [], 'm2Num_wt', [], 'm3Num_wt', [], 'm4Num_wt', [], ...
    'm5Num_wt', [], ...
    'tauOn_m1', [], 'tauOn_m2', [], 'tauOn_m3', [], 'tauOn_m4', [], ...
    'tauOn_m5', [], ...
    'i_tauOn_m1', [], 'i_tauOn_m2', [], 'i_tauOn_m3', [], 'i_tauOn_m4', [], ...
    'i_tauOn_m5', [], ...
    'tauOn_m1_wt', [], 'tauOn_m2_wt', [], 'tauOn_m3_wt', [], 'tauOn_m4_wt', [], ...
    'tauOn_m5_wt', [], ...
    'i_tauOn_m1_wt', [], 'i_tauOn_m2_wt', [], 'i_tauOn_m3_wt', [], 'i_tauOn_m4_wt', [], ...
    'i_tauOn_m5_wt', [], ...
    'KCI', [], 'KCro', [], 'KQ', [], ...
    'KCIR', [], 'KCroR', [], 'KQR', []);

%Constants
convFac = (1e9*1e15)/(6.022e23); %# of um^3 -> nM
V0 = 1; %um^3
tauFail = 1e3;

%Grid spacings
drlam = 0.1; %0.01
dkCII = 0.1;
dMOI = 0.2; %0.05
dMOIThr = 0.01; %0.05
dtD = 0.1; % 

ss = 0;

%ICs
MOI = 1:5;
y0_m1 = [0 0 0 0 0 0 1*convFac/V0];
y0_m2 = [0 0 0 0 0 0 2*convFac/V0];
y0_m3 = [0 0 0 0 0 0 3*convFac/V0];
y0_m4 = [0 0 0 0 0 0 4*convFac/V0];
y0_m5 = [0 0 0 0 0 0 5*convFac/V0];

%Solve ODEs
options = odeset('Nonnegative', [], 'RelTol', 1e-6, ...
    'AbsTol', 1e-6);

for i = sol_Ind
    %ASSIGN PARAMS---------------------------------------------------------
    %cI PRM basal prod, ratio of cI PRM active / basal, ratio cI PRE/PRM basal prod
    %CI translation rate, cro prod rate, Cro translation rate, cII prod rate
    %CII translation rate, replication rate
    prod(1:9) = data(i, 1:9);
    prod(3) = prod(1)*prod(3);
    degr = zeros(8, 1);
    %kdil, kcI, kCI, kcro, kCro, kcII, kCII, kM
    [degr(1), degr(2), degr(3), degr(4), degr(5), degr(6), degr(7), degr(8)] = ...
        deal(data(i, 10), data(i, 11), data(i, 12), data(i, 13), data(i, 14), ...
        data(i, 15), data(i, 16), 0); 
    %nPRM_CIu, nPRM_CId, nPRM_Cro, nPRE, nPR_Cro, nPR_CI, nDeg_CII
    n = [
        data(i, 17);      %nPRM,CI+
        data(i, 18);      %nPRM,CI-
        data(i, 19);      %nPRM,Cro
        data(i, 20);      %nPRE
        data(i, 21);      %nCro,Cro
        data(i, 22);      %nCro,CI
        data(i, 23);      %nCII,Cro
        data(i, 24);      %nCII,CI
        data(i, 25);      %nM,Cro
        data(i, 26);      %nM,CI
        data(i, 27);      %nDeg,CII
         ];
    %KPRM_CIu, KPRM_CId, KPRM_Cro, KPRE, KCro_Cro, KCro_CI, 
    %KCII_Cro, KCII_CI, KM_Cro, KM_CI, KDeg_CII
    K = zeros(11, 1);
    [K(1), K(2), K(3), K(4), K(5), K(6), K(7), K(8), K(9), K(10), K(11)] = ...
        deal(data(i, 28), data(i, 29), data(i, 30), data(i, 31), data(i, 32), ...
        data(i, 33), data(i, 34), data(i, 35), data(i, 36), data(i, 37), ...
        data(i, 38));
    dt = data(i, end-2); %offset
    tau = data(i, end-1);
    kdil = degr(1);
    sol(i).parameters = data(i, :);
    sol(i).prod = prod;
    sol(i).degr = degr;
    sol(i).n = n;
    sol(i).K = K;
    sol(i).dt = dt;
    sol(i).tau = tau;
    
    tspan = 0:0.1:65;
    tspanLong = 0:0.1:240;
    
    V = V0.*exp(tspan.*kdil)';
    
    %SOLVE ODES ===========================================================
    
    %P- -------------------------------------------------------------------
    %MOI=1
    [tm1_OP, ym1_OP] = ode15s(@fv19, tspan, y0_m1(1:6), options, n([1:8, end]), ...
        prod(1:end-1), degr, K([1:8, end]), 1, V0, convFac, ss);
    sol(i).m1 = [tm1_OP, ym1_OP];
    sol(i).m1Num = [tm1_OP, ym1_OP.*V./convFac];
    %MOI=2
    [tm2_OP, ym2_OP] = ode15s(@fv19, tspan, y0_m2(1:6), options, n([1:8, end]), ...
        prod(1:end-1), degr, K([1:8, end]), 2, V0, convFac, ss);
    sol(i).m2 = [tm2_OP, ym2_OP];
    sol(i).m2Num = [tm2_OP, ym2_OP.*V./convFac];
    %MOI=3
    [tm3_OP, ym3_OP] = ode15s(@fv19, tspan, y0_m3(1:6), options, n([1:8, end]), ...
        prod(1:end-1), degr, K([1:8, end]), 3, V0, convFac, ss);
    sol(i).m3 = [tm3_OP, ym3_OP];
    sol(i).m3Num = [tm3_OP, ym3_OP.*V./convFac];
    %MOI=4
    [tm4_OP, ym4_OP] = ode15s(@fv19, tspan, y0_m4(1:6), options, n([1:8, end]), ...
        prod(1:end-1), degr, K([1:8, end]), 4, V0, convFac, ss);
    sol(i).m4 = [tm4_OP, ym4_OP];
    sol(i).m4Num = [tm4_OP, ym4_OP.*V./convFac];
    %MOI=5
    [tm5_OP, ym5_OP] = ode15s(@fv19, tspan, y0_m5(1:6), options, n([1:8, end]), ...
        prod(1:end-1), degr, K([1:8, end]), 5, V0, convFac, ss);
    sol(i).m5 = [tm5_OP, ym5_OP];
    sol(i).m5Num = [tm5_OP, ym5_OP.*V./convFac];
    
    %production rates
    [fm1, gm1] = getFluxv19(tm1_OP, [ym1_OP, 1.*convFac./(V0.*exp(kdil.*tm1_OP))], ...
        n, prod, degr, K, tau, V0, convFac);
    [fm5, gm5] = getFluxv19(tm5_OP, [ym5_OP, 5.*convFac./(V0.*exp(kdil.*tm5_OP))], ...
        n, prod, degr, K, tau, V0, convFac);
    sol(i).fm1 = fm1;
    sol(i).gm1 = gm1;
    sol(i).fm5 = fm5;
    sol(i).gm5 = gm5;
    
    %weights
    wm1 = getWeightsv19(tm1_OP, ym1_OP, ...
        n, prod, degr, K, tau, V0, convFac);
    wm2 = getWeightsv19(tm2_OP, ym2_OP, ...
        n, prod, degr, K, tau, V0, convFac);
    wm3 = getWeightsv19(tm3_OP, ym3_OP, ...
        n, prod, degr, K, tau, V0, convFac);
    wm4 = getWeightsv19(tm4_OP, ym4_OP, ...
        n, prod, degr, K, tau, V0, convFac);
    wm5 = getWeightsv19(tm5_OP, ym5_OP, ...
        n, prod, degr, K, tau, V0, convFac);
    sol(i).wm1 = wm1;
    sol(i).wm2 = wm2;
    sol(i).wm3 = wm3;
    sol(i).wm4 = wm4;
    sol(i).wm5 = wm5;
    
    [r, c] = find(wm1(:, 6) >= 0.1);
    if ~isempty(r)
        sol(i).tauOn_m1 = [tm1_OP(r(1)), tm1_OP(r(end))];
        sol(i).i_tauOn_m1 = [r(1), r(end), 1];
    end;
    
    [r, c] = find(wm2(:, 6) >= 0.1);
    if ~isempty(r)
        sol(i).tauOn_m2 = [tm2_OP(r(1)), tm2_OP(r(end))];
        sol(i).i_tauOn_m2 = [r(1), r(end), 1];
    end;
    
    [r, c] = find(wm3(:, 6) >= 0.1);
    if ~isempty(r)
        sol(i).tauOn_m3 = [tm3_OP(r(1)), tm3_OP(r(end))];
        sol(i).i_tauOn_m3 = [r(1), r(end), 1];
    end;
    
    [r, c] = find(wm4(:, 6) >= 0.1);
    if ~isempty(r)
        sol(i).tauOn_m4 = [tm4_OP(r(1)), tm4_OP(r(end))];
        sol(i).i_tauOn_m4 = [r(1), r(end), 1];
    end;
    
    [r, c] = find(wm5(:, 6) >= 0.1);
    if ~isempty(r)
        sol(i).tauOn_m5 = [tm5_OP(r(1)), tm5_OP(r(end))];
        sol(i).i_tauOn_m5 = [r(1), r(end), 1];
    end;
    
    %WT -------------------------------------------------------------------
    %MOI=1
    [tm1_wt, ym1_wt] = ode15s(@fv19_repv3, tspan, y0_m1, options, n, prod, ...
        degr, K, tau, V0, convFac);
    sol(i).m1_wt = [tm1_wt, ym1_wt];
    sol(i).m1Num_wt = [tm1_wt, ym1_wt.*V./convFac];
    %MOI=2
    [tm2_wt, ym2_wt] = ode15s(@fv19_repv3, tspan, y0_m2, options, n, prod, ...
        degr, K, tau, V0, convFac);
    sol(i).m2_wt = [tm2_wt, ym2_wt];
    sol(i).m2Num_wt = [tm2_wt, ym2_wt.*V./convFac];
    %MOI=3
    [tm3_wt, ym3_wt] = ode15s(@fv19_repv3, tspan, y0_m3, options, n, prod, ...
        degr, K, tau, V0, convFac);
    sol(i).m3_wt = [tm3_wt, ym3_wt];
    sol(i).m3Num_wt = [tm3_wt, ym3_wt.*V./convFac];
    %MOI=4
    [tm4_wt, ym4_wt] = ode15s(@fv19_repv3, tspan, y0_m4, options, n, prod, ...
        degr, K, tau, V0, convFac);
    sol(i).m4_wt = [tm4_wt, ym4_wt];
    sol(i).m4Num_wt = [tm4_wt, ym4_wt.*V./convFac];
    %MOI=5
    [tm5_wt, ym5_wt] = ode15s(@fv19_repv3, tspan, y0_m5, options, n, prod, ...
        degr, K, tau, V0, convFac);
    sol(i).m5_wt = [tm5_wt, ym5_wt];
    sol(i).m5Num_wt = [tm5_wt, ym5_wt.*V./convFac];
    
    %production rates
    [fm1_wt, gm1_wt] = getFluxv19(tm1_wt, ym1_wt, ...
        n, prod, degr, K, tau, V0, convFac);
    [fm5_wt, gm5_wt] = getFluxv19(tm5_wt, ym5_wt, ...
        n, prod, degr, K, tau, V0, convFac);
    sol(i).fm1_wt = fm1_wt;
    sol(i).gm1_wt = gm1_wt;
    sol(i).fm5_wt = fm5_wt;
    sol(i).gm5_wt = gm5_wt;
    
    %weights
    wm1_wt = getWeightsv19(tm1_wt, ym1_wt, ...
        n, prod, degr, K, tau, V0, convFac);
    wm2_wt = getWeightsv19(tm2_wt, ym2_wt, ...
        n, prod, degr, K, tau, V0, convFac);
    wm3_wt = getWeightsv19(tm3_wt, ym3_wt, ...
        n, prod, degr, K, tau, V0, convFac);
    wm4_wt = getWeightsv19(tm4_wt, ym4_wt, ...
        n, prod, degr, K, tau, V0, convFac);
    wm5_wt = getWeightsv19(tm5_wt, ym5_wt, ...
        n, prod, degr, K, tau, V0, convFac);
    sol(i).wm1_wt = wm1_wt;
    sol(i).wm2_wt = wm2_wt;
    sol(i).wm3_wt = wm3_wt;
    sol(i).wm4_wt = wm4_wt;
    sol(i).wm5_wt = wm5_wt;
    
    [r, c] = find(wm1_wt(:, 6) >= 0.1);
    if ~isempty(r)
        sol(i).tauOn_m1_wt = [tm1_wt(r(1)), tm1_wt(r(end))];
        sol(i).i_tauOn_m1_wt = [r(1), r(end), 1];
    end;
    
    [r, c] = find(wm2_wt(:, 6) >= 0.1);
    if ~isempty(r)
        sol(i).tauOn_m2_wt = [tm2_wt(r(1)), tm2_wt(r(end))];
        sol(i).i_tauOn_m2_wt = [r(1), r(end), 1];
    end;
    
    [r, c] = find(wm3_wt(:, 6) >= 0.1);
    if ~isempty(r)
        sol(i).tauOn_m3_wt = [tm3_wt(r(1)), tm3_wt(r(end))];
        sol(i).i_tauOn_m3_wt = [r(1), r(end), 1];
    end;
    
    [r, c] = find(wm4_wt(:, 6) >= 0.1);
    if ~isempty(r)
        sol(i).tauOn_m4_wt = [tm4_wt(r(1)), tm4_wt(r(end))];
        sol(i).i_tauOn_m4_wt = [r(1), r(end), 1];
    end;
    
    [r, c] = find(wm5_wt(:, 6) >= 0.1);
    if ~isempty(r)
        sol(i).tauOn_m5_wt = [tm5_wt(r(1)), tm5_wt(r(end))];
        sol(i).i_tauOn_m5_wt = [r(1), r(end), 1];
    end;
    
    %Lineage ==============================================================
    %P- -------------------------------------------------------------------
    %MOI=1
    [tm1_OP_l1, ym1_OP_l1] = ode15s(@fv19, tspan, y0_m1(1:6), options, n([1:8, end]), ...
        prod(1:end-1), degr, K([1:8, end]), 1, V0, convFac, ss);
    %tspan3 = tm1_OP_l1(end) + tspan;
    [tm1_OP_l2, ym1_OP_l2] = ode15s(@fv19, tspan, ...
        ym1_OP_l1(end, 1:6), options, n([1:8, end]), ...
        prod(1:end-1), degr, K([1:8, end]), 1, V0, convFac, ss);
    tm1_OP_l = [tm1_OP_l1(1:end-1); tm1_OP_l2 + tspan(end);];
    ym1_OP_l = [ym1_OP_l1(1:end-1, :); ym1_OP_l2;];
    sol(i).m1_l = [tm1_OP_l, ym1_OP_l];
    
    %Cro- =================================================================
    prodCro = prod;
    prodCro(5) = 0; 
    %MOI=1
    [tm1_cro, ym1_cro] = ode15s(@fv19_repv3, tspan, y0_m1, options, n, prodCro, ...
        degr, K, tau, V0, convFac);
    sol(i).m1_cro = [tm1_cro, ym1_cro];
    sol(i).m1Num_cro = [tm1_cro, ym1_cro.*V./convFac];
    
    %Thresholds============================================================
    [KCI, KCro, KCIR, KCroR, KCI_o, KCro_o, KCIR_o, KCroR_o] = ...
        getThresholdsv19(sol(i), MOI, dMOIThr, tspan, convFac, V0);
    
    sol(i).KCro_o = KCro_o;
    sol(i).KCI_o = KCI_o;
    sol(i).KCroR_o = KCroR_o;
    sol(i).KCIR_o = KCIR_o;
    sol(i).KCIR = KCIR;
    sol(i).KCroR = KCroR;
    sol(i).KCI = KCI;
    sol(i).KCro = KCro;
    
    %SIMULATE cII tx = cro tx, cII mRNA degr = 0.5 cro degr================
    %P- -------------------------------------------------------------------
    prodCIITX = prod;
    prodCIITX(7) = prodCIITX(5); %cII = cro tx
    degrCIITX = degr;
    degrCIITX(6) = 0.5*degrCIITX(4); %cII = 1/2 cro degr
    %MOI=1
    [tm1_OP_cIItx, ym1_OP_cIItx] = ode15s(@fv19, tspan, y0_m1(1:6), options, n([1:8, end]), ...
        prodCIITX(1:end-1), degrCIITX, K([1:8, end]), 1, V0, convFac, ss);
    sol(i).m1_cIItx = [tm1_OP_cIItx, ym1_OP_cIItx];
    sol(i).m1Num_cIItx = [tm1_OP_cIItx, ym1_OP_cIItx.*V./convFac];
    %MOI=2
    [tm2_OP_cIItx, ym2_OP_cIItx] = ode15s(@fv19, tspan, y0_m2(1:6), options, n([1:8, end]), ...
        prodCIITX(1:end-1), degrCIITX, K([1:8, end]), 2, V0, convFac, ss);
    sol(i).m2_cIItx = [tm2_OP_cIItx, ym2_OP_cIItx];
    sol(i).m2Num_cIItx = [tm2_OP_cIItx, ym2_OP_cIItx.*V./convFac];
    %MOI=3
    [tm3_OP_cIItx, ym3_OP_cIItx] = ode15s(@fv19, tspan, y0_m3(1:6), options, n([1:8, end]), ...
        prodCIITX(1:end-1), degrCIITX, K([1:8, end]), 3, V0, convFac, ss);
    sol(i).m3_cIItx = [tm3_OP_cIItx, ym3_OP_cIItx];
    sol(i).m3Num_cIItx = [tm3_OP_cIItx, ym3_OP_cIItx.*V./convFac];
    %MOI=4
    [tm4_OP_cIItx, ym4_OP_cIItx] = ode15s(@fv19, tspan, y0_m4(1:6), options, n([1:8, end]), ...
        prodCIITX(1:end-1), degrCIITX, K([1:8, end]), 4, V0, convFac, ss);
    sol(i).m4_cIItx = [tm4_OP_cIItx, ym4_OP_cIItx];
    sol(i).m4Num_cIItx = [tm4_OP_cIItx, ym4_OP_cIItx.*V./convFac];
    %MOI=5
    [tm5_OP_cIItx, ym5_OP_cIItx] = ode15s(@fv19, tspan, y0_m5(1:6), options, n([1:8, end]), ...
        prodCIITX(1:end-1), degrCIITX, K([1:8, end]), 5, V0, convFac, ss);
    sol(i).m5_cIItx = [tm5_OP_cIItx, ym5_OP_cIItx];
    sol(i).m5Num_cIItx = [tm5_OP_cIItx, ym5_OP_cIItx.*V./convFac];
    
    %cII tx 2: cII tx = 2 cro tx, cII degr = cro degr----------------------
    prodCIITX = prod;
    prodCIITX(7) = 2*prodCIITX(5); %cII = 2*cro tx
    degrCIITX = degr;
    degrCIITX(6) = degrCIITX(4); %cII = cro degr
    %MOI=1
    [tm1_OP_cIItx2, ym1_OP_cIItx2] = ode15s(@fv19, tspan, y0_m1(1:6), options, n([1:8, end]), ...
        prodCIITX(1:end-1), degrCIITX, K([1:8, end]), 1, V0, convFac, ss);
    sol(i).m1_cIItx2 = [tm1_OP_cIItx2, ym1_OP_cIItx2];
    sol(i).m1Num_cIItx2 = [tm1_OP_cIItx2, ym1_OP_cIItx2.*V./convFac];
    %MOI=2
    [tm2_OP_cIItx2, ym2_OP_cIItx2] = ode15s(@fv19, tspan, y0_m2(1:6), options, n([1:8, end]), ...
        prodCIITX(1:end-1), degrCIITX, K([1:8, end]), 2, V0, convFac, ss);
    sol(i).m2_cIItx2 = [tm2_OP_cIItx2, ym2_OP_cIItx2];
    sol(i).m2Num_cIItx2 = [tm2_OP_cIItx2, ym2_OP_cIItx2.*V./convFac];
    %MOI=3
    [tm3_OP_cIItx2, ym3_OP_cIItx2] = ode15s(@fv19, tspan, y0_m3(1:6), options, n([1:8, end]), ...
        prodCIITX(1:end-1), degrCIITX, K([1:8, end]), 3, V0, convFac, ss);
    sol(i).m3_cIItx2 = [tm3_OP_cIItx2, ym3_OP_cIItx2];
    sol(i).m3Num_cIItx2 = [tm3_OP_cIItx2, ym3_OP_cIItx2.*V./convFac];
    %MOI=4
    [tm4_OP_cIItx2, ym4_OP_cIItx2] = ode15s(@fv19, tspan, y0_m4(1:6), options, n([1:8, end]), ...
        prodCIITX(1:end-1), degrCIITX, K([1:8, end]), 4, V0, convFac, ss);
    sol(i).m4_cIItx2 = [tm4_OP_cIItx2, ym4_OP_cIItx2];
    sol(i).m4Num_cIItx2 = [tm4_OP_cIItx2, ym4_OP_cIItx2.*V./convFac];
    %MOI=5
    [tm5_OP_cIItx2, ym5_OP_cIItx2] = ode15s(@fv19, tspan, y0_m5(1:6), options, n([1:8, end]), ...
        prodCIITX(1:end-1), degrCIITX, K([1:8, end]), 5, V0, convFac, ss);
    sol(i).m5_cIItx2 = [tm5_OP_cIItx2, ym5_OP_cIItx2];
    sol(i).m5Num_cIItx2 = [tm5_OP_cIItx2, ym5_OP_cIItx2.*V./convFac];
    
    %SIMULATE 4HRS (MOI = 1, P-, Shao)=====================================
    %P-, MOI=1, MOI dilutions
    Vlong = V0.*exp(kdil.*tspanLong)';
    [tm1_OP_long, ym1_OP_long] = ode15s(@fv19, tspanLong, y0_m1(1:6), options, ...
        n([1:8, end]), prod(1:end-1), degr, K([1:8, end]), 1, V0, convFac, ss);
    sol(i).m1_long = [tm1_OP_long, ym1_OP_long];
    sol(i).m1Num_long = [tm1_OP_long, ym1_OP_long.*Vlong./convFac];
    %MOI = 2
    [tm2_OP_long, ym2_OP_long] = ode15s(@fv19, tspanLong, y0_m2(1:6), options, ...
        n([1:8, end]), prod(1:end-1), degr, K([1:8, end]), 2, V0, convFac, ss);
    sol(i).m2_long = [tm2_OP_long, ym2_OP_long];
    sol(i).m2Num_long = [tm2_OP_long, ym2_OP_long.*Vlong./convFac];
    %MOI = 3
    [tm3_OP_long, ym3_OP_long] = ode15s(@fv19, tspanLong, y0_m3(1:6), options, ...
        n([1:8, end]), prod(1:end-1), degr, K([1:8, end]), 3, V0, convFac, ss);
    sol(i).m3_long = [tm3_OP_long, ym3_OP_long];
    sol(i).m3Num_long = [tm3_OP_long, ym3_OP_long.*Vlong./convFac];
    %MOI = 4
    [tm4_OP_long, ym4_OP_long] = ode15s(@fv19, tspanLong, y0_m4(1:6), options, ...
        n([1:8, end]), prod(1:end-1), degr, K([1:8, end]), 4, V0, convFac, ss);
    sol(i).m4_long = [tm4_OP_long, ym4_OP_long];
    sol(i).m4Num_long = [tm4_OP_long, ym4_OP_long.*Vlong./convFac];
    %MOI = 5
    [tm5_OP_long, ym5_OP_long] = ode15s(@fv19, tspanLong, y0_m5(1:6), options, ...
        n([1:8, end]), prod(1:end-1), degr, K([1:8, end]), 5, V0, convFac, ss);
    sol(i).m5_long = [tm5_OP_long, ym5_OP_long];
    sol(i).m5Num_long = [tm5_OP_long, ym5_OP_long.*Vlong./convFac];
    
    %P-, MOI = 1, no MOI dilution
    [tm1_OP_longSS, ym1_OP_longSS] = ode15s(@fv19, tspanLong, y0_m1(1:6), options, ...
        n([1:8, end]), prod(1:end-1), degr, K([1:8, end]), 1, V0, convFac, 1);
    sol(i).m1_longSS = [tm1_OP_longSS, ym1_OP_longSS];
    sol(i).m1Num_longSS = [tm1_OP_longSS, ym1_OP_longSS.*Vlong./convFac];
    
    %PHASE PLANE TRAJ, ASYMPTOTIC==========================================
    tspanSS = 0:0.1:6000;
    %no MOI dilution-------------------------------------------------------
    %0
    [tm1_OP_longSS1, ym1_OP_longSS1] = ode15s(@fv19, tspanSS, ...
        [0, 0, 0, 0, 0, 0], options, ...
        n([1:8, end]), prod(1:end-1), degr, K([1:8, end]), 1, V0, convFac, 1);
    sol(i).m1_longSS1 = [tm1_OP_longSS1, ym1_OP_longSS1];
    %[0 0 0 1e3 0 0]
    [tm1_OP_longSS2, ym1_OP_longSS2] = ode15s(@fv19, tspanSS, ...
        [0, 0, 0, 1e3, 0, 0], options, ...
        n([1:8, end]), prod(1:end-1), degr, K([1:8, end]), 1, V0, convFac, 1);
    sol(i).m1_longSS2 = [tm1_OP_longSS2, ym1_OP_longSS2];
    %[0 0 0 1e3 1e3 0]
    [tm1_OP_longSS3, ym1_OP_longSS3] = ode15s(@fv19, tspanSS, ...
        [0, 0, 0, 1e3, 1e3, 0], options, ...
        n([1:8, end]), prod(1:end-1), degr, K([1:8, end]), 1, V0, convFac, 1);
    sol(i).m1_longSS3 = [tm1_OP_longSS3, ym1_OP_longSS3];
    %[0 0 0 0 1e3 0]
    [tm1_OP_longSS4, ym1_OP_longSS4] = ode15s(@fv19, tspanSS, ...
        [0, 0, 0, 0, 1e3, 0], options, ...
        n([1:8, end]), prod(1:end-1), degr, K([1:8, end]), 1, V0, convFac, 1);
    sol(i).m1_longSS4 = [tm1_OP_longSS4, ym1_OP_longSS4];
    %[0 0 0 1e3, 333, 0]
    [tm1_OP_longSS5, ym1_OP_longSS5] = ode15s(@fv19, tspanSS, ...
        [0, 0, 0, 1e3, 333, 0], options, ...
        n([1:8, end]), prod(1:end-1), degr, K([1:8, end]), 1, V0, convFac, 1);
    sol(i).m1_longSS5 = [tm1_OP_longSS5, ym1_OP_longSS5];
    %[0 0 0 333 1e3 0]
    [tm1_OP_longSS6, ym1_OP_longSS6] = ode15s(@fv19, tspanSS, ...
        [0, 0, 0, 333, 1e3, 0], options, ...
        n([1:8, end]), prod(1:end-1), degr, K([1:8, end]), 1, V0, convFac, 1);
    sol(i).m1_longSS6 = [tm1_OP_longSS6, ym1_OP_longSS6];
    
    %High CI Degradation---------------------------------------------------
    %0
    degrCI = degr; degrCI(3) = 2*degrCI(3);
    [tm1_OP_longSS1_d, ym1_OP_longSS1_d] = ode15s(@fv19, tspanSS, ...
        [0, 0, 0, 0, 0, 0], options, ...
        n([1:8, end]), prod(1:end-1), degrCI, K([1:8, end]), 1, V0, convFac, 1);
    sol(i).m1_longSS1_d = [tm1_OP_longSS1_d, ym1_OP_longSS1_d];
    %[0 0 0 1e3 0 0]
    [tm1_OP_longSS2_d, ym1_OP_longSS2_d] = ode15s(@fv19, tspanSS, ...
        [0, 0, 0, 1e3, 0, 0], options, ...
        n([1:8, end]), prod(1:end-1), degrCI, K([1:9, end]), 1, V0, convFac, 1);
    sol(i).m1_longSS2_d = [tm1_OP_longSS2_d, ym1_OP_longSS2_d];
    %[0 0 0 1e3 1e3 0]
    [tm1_OP_longSS3_d, ym1_OP_longSS3_d] = ode15s(@fv19, tspanSS, ...
        [0, 0, 0, 1e3, 1e3, 0], options, ...
        n([1:8, end]), prod(1:end-1), degrCI, K([1:8, end]), 1, V0, convFac, 1);
    sol(i).m1_longSS3_d = [tm1_OP_longSS3_d, ym1_OP_longSS3_d];
    %[0 0 0 0 1e3 0]
    [tm1_OP_longSS4_d, ym1_OP_longSS4_d] = ode15s(@fv19, tspanSS, ...
        [0, 0, 0, 0, 1e3, 0], options, ...
        n([1:8, end]), prod(1:end-1), degrCI, K([1:8, end]), 1, V0, convFac, 1);
    sol(i).m1_longSS4_d = [tm1_OP_longSS4_d, ym1_OP_longSS4_d];
    %[0 0 0 1e3, 333, 0]
    [tm1_OP_longSS5_d, ym1_OP_longSS5_d] = ode15s(@fv19, tspanSS, ...
        [0, 0, 0, 1e3, 333, 0], options, ...
        n([1:8, end]), prod(1:end-1), degrCI, K([1:8, end]), 1, V0, convFac, 1);
    sol(i).m1_longSS5_d = [tm1_OP_longSS5_d, ym1_OP_longSS5_d];
    %[0 0 0 333 1e3 0]
    [tm1_OP_longSS6_d, ym1_OP_longSS6_d] = ode15s(@fv19, tspanSS, ...
        [0, 0, 0, 333, 1e3, 0], options, ...
        n([1:8, end]), prod(1:end-1), degrCI, K([1:8, end]), 1, V0, convFac, 1);
    sol(i).m1_longSS6_d = [tm1_OP_longSS6_d, ym1_OP_longSS6_d];
    
    %MOI dilution----------------------------------------------------------
    %0
    [tm1_OP_long1, ym1_OP_long1] = ode15s(@fv19, tspanSS, ...
        [0, 0, 0, 0, 0, 0], options, ...
        n([1:8, end]), prod(1:end-1), degr, K([1:8, end]), 1, V0, convFac, 0);
    sol(i).m1_long1 = [tm1_OP_long1, ym1_OP_long1];
    %[0 0 0 1e3 0 0]
    [tm1_OP_long2, ym1_OP_long2] = ode15s(@fv19, tspanSS, ...
        [0, 0, 0, 1e3, 0, 0], options, ...
        n([1:8, end]), prod(1:end-1), degr, K([1:8, end]), 1, V0, convFac, 0);
    sol(i).m1_long2 = [tm1_OP_long2, ym1_OP_long2];
    %[0 0 0 1e3 1e3 0]
    [tm1_OP_long3, ym1_OP_long3] = ode15s(@fv19, tspanSS, ...
        [0, 0, 0, 1e3, 1e3, 0], options, ...
        n([1:8, end]), prod(1:end-1), degr, K([1:8, end]), 1, V0, convFac, 0);
    sol(i).m1_long3 = [tm1_OP_long3, ym1_OP_long3];
    %[0 0 0 0 1e3 0]
    [tm1_OP_long4, ym1_OP_long4] = ode15s(@fv19, tspanSS, ...
        [0, 0, 0, 0, 1e3, 0], options, ...
        n([1:8, end]), prod(1:end-1), degr, K([1:8, end]), 1, V0, convFac, 0);
    sol(i).m1_long4 = [tm1_OP_long4, ym1_OP_long4];
    %[0 0 0 1e3, 333, 0]
    [tm1_OP_long5, ym1_OP_long5] = ode15s(@fv19, tspanSS, ...
        [0, 0, 0, 1e3, 333, 0], options, ...
        n([1:8, end]), prod(1:end-1), degr, K([1:8, end]), 1, V0, convFac, 0);
    sol(i).m1_long5 = [tm1_OP_long5, ym1_OP_long5];
    %[0 0 0 333 1e3 0]
    [tm1_OP_long6, ym1_OP_long6] = ode15s(@fv19, tspanSS, ...
        [0, 0, 0, 333, 1e3, 0], options, ...
        n([1:8, end]), prod(1:end-1), degr, K([1:8, end]), 1, V0, convFac, 0);
    sol(i).m1_long6 = [tm1_OP_long6, ym1_OP_long6];
    
    %DETERMINE MOI = 1 BISTABILITY=========================================
    optionsSS = optimoptions('fsolve', 'Display', 'iter', 'OptimalityTolerance', ...
        1e-8);
    %P-, MOI = 1, no MOI dilution---lysogenic
    tspanSS = 0:0.1:6000;
    [tm1_OP_longSS1, ym1_OP_longSS1] = ode15s(@fv19, tspanSS, [0,0,0,1e3,0,0], options, ...
        n([1:8, end]), prod(1:end-1), degr, K([1:8, end]), 1, V0, convFac, 1);
    [x1, fval1, exitflag1, output1, jacobian1] = fsolve(@fv19_ss, ym1_OP_longSS1(end, :), ...
        optionsSS, n([1:8, end]), prod(1:end-1), degr, K([1:8, end]), 1, V0, convFac);
    sol(i).lysSS = x1;
    
    %P-, MOI = 1, no MOI dilution---lysogenic
    [tm1_OP_longSS2, ym1_OP_longSS2] = ode15s(@fv19, tspanSS, [0,0,0,0,1e3,0], options, ...
        n([1:8, end]), prod(1:end-1), degr, K([1:8, end]), 1, V0, convFac, 1);
    [x2, fval2, exitflag2, output2, jacobian2] = fsolve(@fv19_ss, ym1_OP_longSS2(end, :), ...
        optionsSS, n([1:8, end]), prod(1:end-1), degr, K([1:8, end]), 1, V0, convFac);
    sol(i).lytSS = x2;
    
    %Check---origin
    [tm1_OP_longSS3, ym1_OP_longSS3] = ode15s(@fv19, tspanSS, [0,0,0,0,0,0], options, ...
        n([1:8, end]), prod(1:end-1), degr, K([1:8, end]), 1, V0, convFac, 1);
    [x3, fval3, exitflag3, output3, jacobian3] = fsolve(@fv19_ss, ym1_OP_longSS3(end, :), ...
        optionsSS, n([1:8, end]), prod(1:end-1), degr, K([1:8, end]), 1, V0, convFac);
    
    if sum((x1-x2).^2) > ssCut %Distinct lytic and lysogenic states
       checkSS(i) = 1; %Distinct 
    end
    
    if x1(4) > x1(5) || x3(4) > x3(5)
       checkLysogeny(i) = 1; %Lysogenic state 
    end
    
    %CALCULATE BIFURCATION, BEST SOL ONLY==================================
    optionsSS2 = optimoptions('fsolve', 'Display', 'none', 'OptimalityTolerance', ...
        1e-8);
    x = 0:0.01:1;
    x_lyt = zeros(length(x), 6);
    x_lys = zeros(length(x), 6);
    if i == 1
       for j = 1:length(x)
           %Normal---------------------------------------------------------
           %lys
           [tm1_OP_SS1, ym1_OP_SS1] = ode15s(@fv19, tspanSS, [0,0,0,1e3,0,0], options, ...
               n([1:8, end]), prod(1:end-1), degr, K([1:8, end]), x(j), V0, convFac, 1);
           [x1, fval1, exitflag1, output1, jacobian1] = fsolve(@fv19_ss, ym1_OP_SS1(end, :), ...
               optionsSS2, n([1:8, end]), prod(1:end-1), degr, K([1:8, end]), x(j), V0, convFac);
           x_lys(j, :) = x1;
           
           %lyt
           [tm1_OP_SS2, ym1_OP_SS2] = ode15s(@fv19, tspanSS, [0,0,0,0,1e3,0], options, ...
               n([1:8, end]), prod(1:end-1), degr, K([1:8, end]), x(j), V0, convFac, 1);
           [x2, fval2, exitflag2, output2, jacobian2] = fsolve(@fv19_ss, ym1_OP_SS2(end, :), ...
               optionsSS2, n([1:8, end]), prod(1:end-1), degr, K([1:8, end]), x(j), V0, convFac);
           x_lyt(j, :) = x2;
       end
       sol(i).MOI_bf = x;
       sol(i).lytSS_bf = x_lyt;
       sol(i).lysSS_bf = x_lys;
    end
    
    %COMPUTE BIFURCATION DIAGRAM===========================================
    boolBif = 0; %start from last SS
    b = 1:0.01:10;
    [b, xSS] = getBifurc(sol(i), prod, degr, K, n, V0, [0 0 0 1e3 0 0], ...
        b, convFac, boolBif);
    sol(i).b_lh = b;
    sol(i).xSS_lh = xSS;
    
    [~, xSS] = getBifurc(sol(i), prod, degr, K, n, V0, [0 0 0 1e3 0 0], ...
        flip(b), convFac, boolBif);
    sol(i).b_hl = flip(b);
    sol(i).xSS_hl = xSS;
    
    %SIMULATE PRM TRANSLATION INEFFICIENCY=================================
    aPRMtr = 0.05;
    prodPRMtr = [prod(1:3), aPRMtr*prod(4), prod(4), prod(5:end)];
    %P- -------------------------------------------------------------------
    %MOI=1
    [tm1_OP_PRMtr, ym1_OP_PRMtr] = ode15s(@fv19_PRMtr, tspan, [y0_m1(1:6), 0], ...
        options, n([1:8, end]), prodPRMtr(1:end-1), degr, K([1:8, end]), 1, V0, convFac, ss);
    sol(i).m1_PRMtr = [tm1_OP_PRMtr, ym1_OP_PRMtr];
    sol(i).m1Num_PRMtr = [tm1_OP_PRMtr, ym1_OP_PRMtr.*V./convFac];
    %MOI=2
    [tm2_OP_PRMtr, ym2_OP_PRMtr] = ode15s(@fv19_PRMtr, tspan, [y0_m2(1:6), 0], ...
        options, n([1:8, end]), prodPRMtr(1:end-1), degr, K([1:8, end]), 2, V0, convFac, ss);
    sol(i).m2_PRMtr = [tm2_OP_PRMtr, ym2_OP_PRMtr];
    sol(i).m2Num_PRMtr = [tm2_OP_PRMtr, ym2_OP_PRMtr.*V./convFac];
    %MOI=3
    [tm3_OP_PRMtr, ym3_OP_PRMtr] = ode15s(@fv19_PRMtr, tspan, [y0_m3(1:6), 0], ...
        options, n([1:8, end]), prodPRMtr(1:end-1), degr, K([1:8, end]), 3, V0, convFac, ss);
    sol(i).m3_PRMtr = [tm3_OP_PRMtr, ym3_OP_PRMtr];
    sol(i).m3Num_PRMtr = [tm3_OP_PRMtr, ym3_OP_PRMtr.*V./convFac];
    %MOI=4
    [tm4_OP_PRMtr, ym4_OP_PRMtr] = ode15s(@fv19_PRMtr, tspan, [y0_m4(1:6), 0], ...
        options, n([1:8, end]), prodPRMtr(1:end-1), degr, K([1:8, end]), 4, V0, convFac, ss);
    sol(i).m4_PRMtr = [tm4_OP_PRMtr, ym4_OP_PRMtr];
    sol(i).m4Num_PRMtr = [tm4_OP_PRMtr, ym4_OP_PRMtr.*V./convFac];
    %MOI=5
    [tm5_OP_PRMtr, ym5_OP_PRMtr] = ode15s(@fv19_PRMtr, tspan, [y0_m5(1:6), 0], ...
        options, n([1:8, end]), prodPRMtr(1:end-1), degr, K([1:8, end]), 5, V0, convFac, ss);
    sol(i).m5_PRMtr = [tm5_OP_PRMtr, ym5_OP_PRMtr];
    sol(i).m5Num_PRMtr = [tm5_OP_PRMtr, ym5_OP_PRMtr.*V./convFac];
    
    %WT -------------------------------------------------------------------
    %MOI=1
    [tm1_wt_PRMtr, ym1_wt_PRMtr] = ode15s(@fv19_repv3_PRMtr, tspan, [0, y0_m1], ...
        options, n, prodPRMtr, degr, K, tau, V0, convFac);
    sol(i).m1_wt_PRMtr = [tm1_wt_PRMtr, ym1_wt_PRMtr];
    sol(i).m1Num_wt_PRMtr = [tm1_wt_PRMtr, ym1_wt_PRMtr.*V./convFac];
    %MOI=2
    [tm2_wt_PRMtr, ym2_wt_PRMtr] = ode15s(@fv19_repv3_PRMtr, tspan, [0, y0_m2], ...
        options, n, prodPRMtr, degr, K, tau, V0, convFac);
    sol(i).m2_wt_PRMtr = [tm2_wt_PRMtr, ym2_wt_PRMtr];
    sol(i).m2Num_wt_PRMtr = [tm2_wt_PRMtr, ym2_wt_PRMtr.*V./convFac];
    %MOI=3
    [tm3_wt_PRMtr, ym3_wt_PRMtr] = ode15s(@fv19_repv3_PRMtr, tspan, [0, y0_m3], ...
        options, n, prodPRMtr, degr, K, tau, V0, convFac);
    sol(i).m3_wt_PRMtr = [tm3_wt_PRMtr, ym3_wt_PRMtr];
    sol(i).m3Num_wt_PRMtr = [tm3_wt_PRMtr, ym3_wt_PRMtr.*V./convFac];
    %MOI=4
    [tm4_wt_PRMtr, ym4_wt_PRMtr] = ode15s(@fv19_repv3_PRMtr, tspan, [0, y0_m4], ...
        options, n, prodPRMtr, degr, K, tau, V0, convFac);
    sol(i).m4_wt_PRMtr = [tm4_wt_PRMtr, ym4_wt_PRMtr];
    sol(i).m4Num_wt_PRMtr = [tm4_wt_PRMtr, ym4_wt_PRMtr.*V./convFac];
    %MOI=5
    [tm5_wt_PRMtr, ym5_wt_PRMtr] = ode15s(@fv19_repv3_PRMtr, tspan, [0, y0_m5], ...
        options, n, prodPRMtr, degr, K, tau, V0, convFac);
    sol(i).m5_wt_PRMtr = [tm5_wt_PRMtr, ym5_wt_PRMtr];
    sol(i).m5Num_wt_PRMtr = [tm5_wt_PRMtr, ym5_wt_PRMtr.*V./convFac];
    
    %REPLICATION CONTROL===================================================
    %WT, no CI control
    KCIrep = K;
    KCIrep(10) = Inf;
    %MOI=1
    [tm1_wt_CIrep, ym1_wt_CIrep] = ode15s(@fv19_repv3, tspan, y0_m1, options, n, prod, ...
        degr, KCIrep, tau, V0, convFac);
    sol(i).m1_wt_CIrep = [tm1_wt_CIrep, ym1_wt_CIrep];
    sol(i).m1Num_wt_CIrep = [tm1_wt_CIrep, ym1_wt_CIrep.*V./convFac];
    %MOI=2
    [tm2_wt_CIrep, ym2_wt_CIrep] = ode15s(@fv19_repv3, tspan, y0_m2, options, n, prod, ...
        degr, KCIrep, tau, V0, convFac);
    sol(i).m2_wt_CIrep = [tm2_wt_CIrep, ym2_wt_CIrep];
    sol(i).m2Num_wt_CIrep = [tm2_wt_CIrep, ym2_wt_CIrep.*V./convFac];
    %MOI=3
    [tm3_wt_CIrep, ym3_wt_CIrep] = ode15s(@fv19_repv3, tspan, y0_m3, options, n, prod, ...
        degr, KCIrep, tau, V0, convFac);
    sol(i).m3_wt_CIrep = [tm3_wt_CIrep, ym3_wt_CIrep];
    sol(i).m3Num_wt_CIrep = [tm3_wt_CIrep, ym3_wt_CIrep.*V./convFac];
    %MOI=4
    [tm4_wt_CIrep, ym4_wt_CIrep] = ode15s(@fv19_repv3, tspan, y0_m4, options, n, prod, ...
        degr, KCIrep, tau, V0, convFac);
    sol(i).m4_wt_CIrep = [tm4_wt_CIrep, ym4_wt_CIrep];
    sol(i).m4Num_wt_CIrep = [tm4_wt_CIrep, ym4_wt_CIrep.*V./convFac];
    %MOI=5
    [tm5_wt_CIrep, ym5_wt_CIrep] = ode15s(@fv19_repv3, tspan, y0_m5, options, n, prod, ...
        degr, KCIrep, tau, V0, convFac);
    sol(i).m5_wt_CIrep = [tm5_wt_CIrep, ym5_wt_CIrep];
    sol(i).m5Num_wt_CIrep = [tm5_wt_CIrep, ym5_wt_CIrep.*V./convFac];
    
    %WT, no Cro control
    KCrorep = K;
    KCrorep(9) = Inf;
    %MOI=1
    [tm1_wt_Crorep, ym1_wt_Crorep] = ode15s(@fv19_repv3, tspan, y0_m1, options, n, prod, ...
        degr, KCrorep, tau, V0, convFac);
    sol(i).m1_wt_Crorep = [tm1_wt_Crorep, ym1_wt_Crorep];
    sol(i).m1Num_wt_Crorep = [tm1_wt_Crorep, ym1_wt_Crorep.*V./convFac];
    %MOI=2
    [tm2_wt_Crorep, ym2_wt_Crorep] = ode15s(@fv19_repv3, tspan, y0_m2, options, n, prod, ...
        degr, KCrorep, tau, V0, convFac);
    sol(i).m2_wt_Crorep = [tm2_wt_Crorep, ym2_wt_Crorep];
    sol(i).m2Num_wt_Crorep = [tm2_wt_Crorep, ym2_wt_Crorep.*V./convFac];
    %MOI=3
    [tm3_wt_Crorep, ym3_wt_Crorep] = ode15s(@fv19_repv3, tspan, y0_m3, options, n, prod, ...
        degr, KCrorep, tau, V0, convFac);
    sol(i).m3_wt_Crorep = [tm3_wt_Crorep, ym3_wt_Crorep];
    sol(i).m3Num_wt_Crorep = [tm3_wt_Crorep, ym3_wt_Crorep.*V./convFac];
    %MOI=4
    [tm4_wt_Crorep, ym4_wt_Crorep] = ode15s(@fv19_repv3, tspan, y0_m4, options, n, prod, ...
        degr, KCrorep, tau, V0, convFac);
    sol(i).m4_wt_Crorep = [tm4_wt_Crorep, ym4_wt_Crorep];
    sol(i).m4Num_wt_Crorep = [tm4_wt_Crorep, ym4_wt_Crorep.*V./convFac];
    %MOI=5
    [tm5_wt_Crorep, ym5_wt_Crorep] = ode15s(@fv19_repv3, tspan, y0_m5, options, n, prod, ...
        degr, KCrorep, tau, V0, convFac);
    sol(i).m5_wt_Crorep = [tm5_wt_Crorep, ym5_wt_Crorep];
    sol(i).m5Num_wt_Crorep = [tm5_wt_Crorep, ym5_wt_Crorep.*V./convFac];
    
    %WT, Neither Cro or CI control
    KCroCIrep = K;
    KCroCIrep([9, 10]) = Inf;
    %MOI=1
    [tm1_wt_CroCIrep, ym1_wt_CroCIrep] = ode15s(@fv19_repv3, tspan, y0_m1, options, n, prod, ...
        degr, KCroCIrep, tau, V0, convFac);
    sol(i).m1_wt_CroCIrep = [tm1_wt_CroCIrep, ym1_wt_CroCIrep];
    sol(i).m1Num_wt_CroCIrep = [tm1_wt_CroCIrep, ym1_wt_CroCIrep.*V./convFac];
    %MOI=2
    [tm2_wt_CroCIrep, ym2_wt_CroCIrep] = ode15s(@fv19_repv3, tspan, y0_m2, options, n, prod, ...
        degr, KCroCIrep, tau, V0, convFac);
    sol(i).m2_wt_CroCIrep = [tm2_wt_CroCIrep, ym2_wt_CroCIrep];
    sol(i).m2Num_wt_CroCIrep = [tm2_wt_CroCIrep, ym2_wt_CroCIrep.*V./convFac];
    %MOI=3
    [tm3_wt_CroCIrep, ym3_wt_CroCIrep] = ode15s(@fv19_repv3, tspan, y0_m3, options, n, prod, ...
        degr, KCroCIrep, tau, V0, convFac);
    sol(i).m3_wt_CroCIrep = [tm3_wt_CroCIrep, ym3_wt_CroCIrep];
    sol(i).m3Num_wt_CroCIrep = [tm3_wt_CroCIrep, ym3_wt_CroCIrep.*V./convFac];
    %MOI=4
    [tm4_wt_CroCIrep, ym4_wt_CroCIrep] = ode15s(@fv19_repv3, tspan, y0_m4, options, n, prod, ...
        degr, KCroCIrep, tau, V0, convFac);
    sol(i).m4_wt_CroCIrep = [tm4_wt_CroCIrep, ym4_wt_CroCIrep];
    sol(i).m4Num_wt_CroCIrep = [tm4_wt_CroCIrep, ym4_wt_CroCIrep.*V./convFac];
    %MOI=5
    [tm5_wt_CroCIrep, ym5_wt_CroCIrep] = ode15s(@fv19_repv3, tspan, y0_m5, options, n, prod, ...
        degr, KCroCIrep, tau, V0, convFac);
    sol(i).m5_wt_CroCIrep = [tm5_wt_CroCIrep, ym5_wt_CroCIrep];
    sol(i).m5Num_wt_CroCIrep = [tm5_wt_CroCIrep, ym5_wt_CroCIrep.*V./convFac];
    
    %NO CRO CONTROL OF CII=================================================
    %WT, no Cro control
    KCroCII = K;
    KCroCII(7) = Inf;
    %MOI=1
    [tm1_wt_CroCII, ym1_wt_CroCII] = ode15s(@fv19_repv3, tspan, y0_m1, options, n, prod, ...
        degr, KCroCII, tau, V0, convFac);
    sol(i).m1_wt_CroCII = [tm1_wt_CroCII, ym1_wt_CroCII];
    sol(i).m1Num_wt_CroCII = [tm1_wt_CroCII, ym1_wt_CroCII.*V./convFac];
    %MOI=2
    [tm2_wt_CroCII, ym2_wt_CroCII] = ode15s(@fv19_repv3, tspan, y0_m2, options, n, prod, ...
        degr, KCroCII, tau, V0, convFac);
    sol(i).m2_wt_CroCII = [tm2_wt_CroCII, ym2_wt_CroCII];
    sol(i).m2Num_wt_CroCII = [tm2_wt_CroCII, ym2_wt_CroCII.*V./convFac];
    %MOI=3
    [tm3_wt_CroCII, ym3_wt_CroCII] = ode15s(@fv19_repv3, tspan, y0_m3, options, n, prod, ...
        degr, KCroCII, tau, V0, convFac);
    sol(i).m3_wt_CroCII = [tm3_wt_CroCII, ym3_wt_CroCII];
    sol(i).m3Num_wt_CroCII = [tm3_wt_CroCII, ym3_wt_CroCII.*V./convFac];
    %MOI=4
    [tm4_wt_CroCII, ym4_wt_CroCII] = ode15s(@fv19_repv3, tspan, y0_m4, options, n, prod, ...
        degr, KCroCII, tau, V0, convFac);
    sol(i).m4_wt_CroCII = [tm4_wt_CroCII, ym4_wt_CroCII];
    sol(i).m4Num_wt_CroCII = [tm4_wt_CroCII, ym4_wt_CroCII.*V./convFac];
    %MOI=5
    [tm5_wt_CroCII, ym5_wt_CroCII] = ode15s(@fv19_repv3, tspan, y0_m5, options, n, prod, ...
        degr, KCroCII, tau, V0, convFac);
    sol(i).m5_wt_CroCII = [tm5_wt_CroCII, ym5_wt_CroCII];
    sol(i).m5Num_wt_CroCII = [tm5_wt_CroCII, ym5_wt_CroCII.*V./convFac];
    
    %weights
    wm1_wt_CroCII = getWeightsv19(tm1_wt_CroCII, ym1_wt_CroCII, ...
        n, prod, degr, KCroCII, tau, V0, convFac);
    wm2_wt_CroCII = getWeightsv19(tm2_wt_CroCII, ym2_wt_CroCII, ...
        n, prod, degr, KCroCII, tau, V0, convFac);
    wm3_wt_CroCII = getWeightsv19(tm3_wt_CroCII, ym3_wt_CroCII, ...
        n, prod, degr, KCroCII, tau, V0, convFac);
    wm4_wt_CroCII = getWeightsv19(tm4_wt_CroCII, ym4_wt_CroCII, ...
        n, prod, degr, KCroCII, tau, V0, convFac);
    wm5_wt_CroCII = getWeightsv19(tm5_wt_CroCII, ym5_wt_CroCII, ...
        n, prod, degr, KCroCII, tau, V0, convFac);
    sol(i).wm1_wt_CroCII = wm1_wt_CroCII;
    sol(i).wm2_wt_CroCII = wm2_wt_CroCII;
    sol(i).wm3_wt_CroCII = wm3_wt_CroCII;
    sol(i).wm4_wt_CroCII = wm4_wt_CroCII;
    sol(i).wm5_wt_CroCII = wm5_wt_CroCII;
    
    [r, c] = find(wm1_wt_CroCII(:, 6) >= 0.1);
    if ~isempty(r)
        sol(i).tauOn_m1_wt_CroCII = [tm1_wt_CroCII(r(1)), tm1_wt_CroCII(r(end))];
        sol(i).i_tauOn_m1_wt_CroCII = [r(1), r(end), 1];
    end;
    
    [r, c] = find(wm2_wt_CroCII(:, 6) >= 0.1);
    if ~isempty(r)
        sol(i).tauOn_m2_wt_CroCII = [tm2_wt_CroCII(r(1)), tm2_wt_CroCII(r(end))];
        sol(i).i_tauOn_m2_wt_CroCII = [r(1), r(end), 1];
    end;
    
    [r, c] = find(wm3_wt_CroCII(:, 6) >= 0.1);
    if ~isempty(r)
        sol(i).tauOn_m3_wt_CroCII = [tm3_wt_CroCII(r(1)), tm3_wt_CroCII(r(end))];
        sol(i).i_tauOn_m3_wt_CroCII = [r(1), r(end), 1];
    end;
    
    [r, c] = find(wm4_wt_CroCII(:, 6) >= 0.1);
    if ~isempty(r)
        sol(i).tauOn_m4_wt_CroCII = [tm4_wt_CroCII(r(1)), tm4_wt_CroCII(r(end))];
        sol(i).i_tauOn_m4_wt_CroCII = [r(1), r(end), 1];
    end;
    
    [r, c] = find(wm5_wt_CroCII(:, 6) >= 0.1);
    if ~isempty(r)
        sol(i).tauOn_m5_wt_CroCII = [tm5_wt_CroCII(r(1)), tm5_wt_CroCII(r(end))];
        sol(i).i_tauOn_m5_wt_CroCII = [r(1), r(end), 1];
    end;
    
    %NO CI CONTROL OF CII=================================================
    %WT, no Cro control
    KCICII = K;
    KCICII(8) = Inf;
    %MOI=1
    [tm1_wt_CICII, ym1_wt_CICII] = ode15s(@fv19_repv3, tspan, y0_m1, options, n, prod, ...
        degr, KCICII, tau, V0, convFac);
    sol(i).m1_wt_CICII = [tm1_wt_CICII, ym1_wt_CICII];
    sol(i).m1Num_wt_CICII = [tm1_wt_CICII, ym1_wt_CICII.*V./convFac];
    %MOI=2
    [tm2_wt_CICII, ym2_wt_CICII] = ode15s(@fv19_repv3, tspan, y0_m2, options, n, prod, ...
        degr, KCICII, tau, V0, convFac);
    sol(i).m2_wt_CICII = [tm2_wt_CICII, ym2_wt_CICII];
    sol(i).m2Num_wt_CICII = [tm2_wt_CICII, ym2_wt_CICII.*V./convFac];
    %MOI=3
    [tm3_wt_CICII, ym3_wt_CICII] = ode15s(@fv19_repv3, tspan, y0_m3, options, n, prod, ...
        degr, KCICII, tau, V0, convFac);
    sol(i).m3_wt_CICII = [tm3_wt_CICII, ym3_wt_CICII];
    sol(i).m3Num_wt_CICII = [tm3_wt_CICII, ym3_wt_CICII.*V./convFac];
    %MOI=4
    [tm4_wt_CICII, ym4_wt_CICII] = ode15s(@fv19_repv3, tspan, y0_m4, options, n, prod, ...
        degr, KCICII, tau, V0, convFac);
    sol(i).m4_wt_CICII = [tm4_wt_CICII, ym4_wt_CICII];
    sol(i).m4Num_wt_CICII = [tm4_wt_CICII, ym4_wt_CICII.*V./convFac];
    %MOI=5
    [tm5_wt_CICII, ym5_wt_CICII] = ode15s(@fv19_repv3, tspan, y0_m5, options, n, prod, ...
        degr, KCICII, tau, V0, convFac);
    sol(i).m5_wt_CICII = [tm5_wt_CICII, ym5_wt_CICII];
    sol(i).m5Num_wt_CICII = [tm5_wt_CICII, ym5_wt_CICII.*V./convFac];
    
    %weights
    wm1_wt_CICII = getWeightsv19(tm1_wt_CICII, ym1_wt_CICII, ...
        n, prod, degr, KCICII, tau, V0, convFac);
    wm2_wt_CICII = getWeightsv19(tm2_wt_CICII, ym2_wt_CICII, ...
        n, prod, degr, KCICII, tau, V0, convFac);
    wm3_wt_CICII = getWeightsv19(tm3_wt_CICII, ym3_wt_CICII, ...
        n, prod, degr, KCICII, tau, V0, convFac);
    wm4_wt_CICII = getWeightsv19(tm4_wt_CICII, ym4_wt_CICII, ...
        n, prod, degr, KCICII, tau, V0, convFac);
    wm5_wt_CICII = getWeightsv19(tm5_wt_CICII, ym5_wt_CICII, ...
        n, prod, degr, KCICII, tau, V0, convFac);
    sol(i).wm1_wt_CICII = wm1_wt_CICII;
    sol(i).wm2_wt_CICII = wm2_wt_CICII;
    sol(i).wm3_wt_CICII = wm3_wt_CICII;
    sol(i).wm4_wt_CICII = wm4_wt_CICII;
    sol(i).wm5_wt_CICII = wm5_wt_CICII;
    
    [r, c] = find(wm1_wt_CICII(:, 6) >= 0.1);
    if ~isempty(r)
        sol(i).tauOn_m1_wt_CICII = [tm1_wt_CICII(r(1)), tm1_wt_CICII(r(end))];
        sol(i).i_tauOn_m1_wt_CICII = [r(1), r(end), 1];
    end;
    
    [r, c] = find(wm2_wt_CICII(:, 6) >= 0.1);
    if ~isempty(r)
        sol(i).tauOn_m2_wt_CICII = [tm2_wt_CICII(r(1)), tm2_wt_CICII(r(end))];
        sol(i).i_tauOn_m2_wt_CICII = [r(1), r(end), 1];
    end;
    
    [r, c] = find(wm3_wt_CICII(:, 6) >= 0.1);
    if ~isempty(r)
        sol(i).tauOn_m3_wt_CICII = [tm3_wt_CICII(r(1)), tm3_wt_CICII(r(end))];
        sol(i).i_tauOn_m3_wt_CICII = [r(1), r(end), 1];
    end;
    
    [r, c] = find(wm4_wt_CICII(:, 6) >= 0.1);
    if ~isempty(r)
        sol(i).tauOn_m4_wt_CICII = [tm4_wt_CICII(r(1)), tm4_wt_CICII(r(end))];
        sol(i).i_tauOn_m4_wt_CICII = [r(1), r(end), 1];
    end;
    
    [r, c] = find(wm5_wt_CICII(:, 6) >= 0.1);
    if ~isempty(r)
        sol(i).tauOn_m5_wt_CICII = [tm5_wt_CICII(r(1)), tm5_wt_CICII(r(end))];
        sol(i).i_tauOn_m5_wt_CICII = [r(1), r(end), 1];
    end;
    
    %NO CI or Cro CONTROL OF CII=================================================
    %WT, no Cro control
    KCroCICII = K;
    KCroCICII(7:8) = Inf;
    %MOI=1
    [tm1_wt_CroCICII, ym1_wt_CroCICII] = ode15s(@fv19_repv3, tspan, y0_m1, options, n, prod, ...
        degr, KCroCICII, tau, V0, convFac);
    sol(i).m1_wt_CroCICII = [tm1_wt_CroCICII, ym1_wt_CroCICII];
    sol(i).m1Num_wt_CroCICII = [tm1_wt_CroCICII, ym1_wt_CroCICII.*V./convFac];
    %MOI=2
    [tm2_wt_CroCICII, ym2_wt_CroCICII] = ode15s(@fv19_repv3, tspan, y0_m2, options, n, prod, ...
        degr, KCroCICII, tau, V0, convFac);
    sol(i).m2_wt_CroCICII = [tm2_wt_CroCICII, ym2_wt_CroCICII];
    sol(i).m2Num_wt_CroCICII = [tm2_wt_CroCICII, ym2_wt_CroCICII.*V./convFac];
    %MOI=3
    [tm3_wt_CroCICII, ym3_wt_CroCICII] = ode15s(@fv19_repv3, tspan, y0_m3, options, n, prod, ...
        degr, KCroCICII, tau, V0, convFac);
    sol(i).m3_wt_CroCICII = [tm3_wt_CroCICII, ym3_wt_CroCICII];
    sol(i).m3Num_wt_CroCICII = [tm3_wt_CroCICII, ym3_wt_CroCICII.*V./convFac];
    %MOI=4
    [tm4_wt_CroCICII, ym4_wt_CroCICII] = ode15s(@fv19_repv3, tspan, y0_m4, options, n, prod, ...
        degr, KCroCICII, tau, V0, convFac);
    sol(i).m4_wt_CroCICII = [tm4_wt_CroCICII, ym4_wt_CroCICII];
    sol(i).m4Num_wt_CroCICII = [tm4_wt_CroCICII, ym4_wt_CroCICII.*V./convFac];
    %MOI=5
    [tm5_wt_CroCICII, ym5_wt_CroCICII] = ode15s(@fv19_repv3, tspan, y0_m5, options, n, prod, ...
        degr, KCroCICII, tau, V0, convFac);
    sol(i).m5_wt_CroCICII = [tm5_wt_CroCICII, ym5_wt_CroCICII];
    sol(i).m5Num_wt_CroCICII = [tm5_wt_CroCICII, ym5_wt_CroCICII.*V./convFac];
    
    %weights
    wm1_wt_CroCICII = getWeightsv19(tm1_wt_CroCICII, ym1_wt_CroCICII, ...
        n, prod, degr, KCroCICII, tau, V0, convFac);
    wm2_wt_CroCICII = getWeightsv19(tm2_wt_CroCICII, ym2_wt_CroCICII, ...
        n, prod, degr, KCroCICII, tau, V0, convFac);
    wm3_wt_CroCICII = getWeightsv19(tm3_wt_CroCICII, ym3_wt_CroCICII, ...
        n, prod, degr, KCroCICII, tau, V0, convFac);
    wm4_wt_CroCICII = getWeightsv19(tm4_wt_CroCICII, ym4_wt_CroCICII, ...
        n, prod, degr, KCroCICII, tau, V0, convFac);
    wm5_wt_CroCICII = getWeightsv19(tm5_wt_CroCICII, ym5_wt_CroCICII, ...
        n, prod, degr, KCroCICII, tau, V0, convFac);
    sol(i).wm1_wt_CroCICII = wm1_wt_CroCICII;
    sol(i).wm2_wt_CroCICII = wm2_wt_CroCICII;
    sol(i).wm3_wt_CroCICII = wm3_wt_CroCICII;
    sol(i).wm4_wt_CroCICII = wm4_wt_CroCICII;
    sol(i).wm5_wt_CroCICII = wm5_wt_CroCICII;
    
    [r, c] = find(wm1_wt_CroCICII(:, 6) >= 0.1);
    if ~isempty(r)
        sol(i).tauOn_m1_wt_CroCICII = [tm1_wt_CroCICII(r(1)), tm1_wt_CroCICII(r(end))];
        sol(i).i_tauOn_m1_wt_CroCICII = [r(1), r(end), 1];
    end;
    
    [r, c] = find(wm2_wt_CroCICII(:, 6) >= 0.1);
    if ~isempty(r)
        sol(i).tauOn_m2_wt_CroCICII = [tm2_wt_CroCICII(r(1)), tm2_wt_CroCICII(r(end))];
        sol(i).i_tauOn_m2_wt_CroCICII = [r(1), r(end), 1];
    end;
    
    [r, c] = find(wm3_wt_CroCICII(:, 6) >= 0.1);
    if ~isempty(r)
        sol(i).tauOn_m3_wt_CroCICII = [tm3_wt_CroCICII(r(1)), tm3_wt_CroCICII(r(end))];
        sol(i).i_tauOn_m3_wt_CroCICII = [r(1), r(end), 1];
    end;
    
    [r, c] = find(wm4_wt_CroCICII(:, 6) >= 0.1);
    if ~isempty(r)
        sol(i).tauOn_m4_wt_CroCICII = [tm4_wt_CroCICII(r(1)), tm4_wt_CroCICII(r(end))];
        sol(i).i_tauOn_m4_wt_CroCICII = [r(1), r(end), 1];
    end;
    
    [r, c] = find(wm5_wt_CroCICII(:, 6) >= 0.1);
    if ~isempty(r)
        sol(i).tauOn_m5_wt_CroCICII = [tm5_wt_CroCICII(r(1)), tm5_wt_CroCICII(r(end))];
        sol(i).i_tauOn_m5_wt_CroCICII = [r(1), r(end), 1];
    end;
    
    %NO CRO CONTROL OF CII, NO CI CONTROL REPLICATION======================
    %WT, no Cro control
    KCroCII_CIM = K;
    KCroCII_CIM(7) = Inf;
    KCroCII_CIM(10) = Inf;
    %MOI=1
    [tm1_wt_CroCII_CIM, ym1_wt_CroCII_CIM] = ode15s(@fv19_repv3, tspan, y0_m1, ...
        options, n, prod, degr, KCroCII_CIM, tau, V0, convFac);
    sol(i).m1_wt_CroCII_CIM = [tm1_wt_CroCII_CIM, ym1_wt_CroCII_CIM];
    sol(i).m1Num_wt_CroCII_CIM = [tm1_wt_CroCII_CIM, ym1_wt_CroCII_CIM.*V./convFac];
    
    %NO CRO CONTROL OF PR (+Rep.)==========================================
    %WT
    KCroPR = K;
    KCroPR([5, 7, 9]) = Inf;
    %MOI=1
    [tm1_wt_CroPR, ym1_wt_CroPR] = ode15s(@fv19_repv3, tspan, y0_m1, options, n, prod, ...
        degr, KCroPR, tau, V0, convFac);
    sol(i).m1_wt_CroPR = [tm1_wt_CroPR, ym1_wt_CroPR];
    sol(i).m1Num_wt_CroPR = [tm1_wt_CroPR, ym1_wt_CroPR.*V./convFac];
    %MOI=2
    [tm2_wt_CroPR, ym2_wt_CroPR] = ode15s(@fv19_repv3, tspan, y0_m2, options, n, prod, ...
        degr, KCroPR, tau, V0, convFac);
    sol(i).m2_wt_CroPR = [tm2_wt_CroPR, ym2_wt_CroPR];
    sol(i).m2Num_wt_CroPR = [tm2_wt_CroPR, ym2_wt_CroPR.*V./convFac];
    %MOI=3
    [tm3_wt_CroPR, ym3_wt_CroPR] = ode15s(@fv19_repv3, tspan, y0_m3, options, n, prod, ...
        degr, KCroPR, tau, V0, convFac);
    sol(i).m3_wt_CroPR = [tm3_wt_CroPR, ym3_wt_CroPR];
    sol(i).m3Num_wt_CroPR = [tm3_wt_CroPR, ym3_wt_CroPR.*V./convFac];
    %MOI=4
    [tm4_wt_CroPR, ym4_wt_CroPR] = ode15s(@fv19_repv3, tspan, y0_m4, options, n, prod, ...
        degr, KCroPR, tau, V0, convFac);
    sol(i).m4_wt_CroPR = [tm4_wt_CroPR, ym4_wt_CroPR];
    sol(i).m4Num_wt_CroPR = [tm4_wt_CroPR, ym4_wt_CroPR.*V./convFac];
    %MOI=5
    [tm5_wt_CroPR, ym5_wt_CroPR] = ode15s(@fv19_repv3, tspan, y0_m5, options, n, prod, ...
        degr, KCroPR, tau, V0, convFac);
    sol(i).m5_wt_CroPR = [tm5_wt_CroPR, ym5_wt_CroPR];
    sol(i).m5Num_wt_CroPR = [tm5_wt_CroPR, ym5_wt_CroPR.*V./convFac];
    
    %weights
    wm1_wt_CroPR = getWeightsv19(tm1_wt_CroPR, ym1_wt_CroPR, ...
        n, prod, degr, KCroPR, tau, V0, convFac);
    wm2_wt_CroPR = getWeightsv19(tm2_wt_CroPR, ym2_wt_CroPR, ...
        n, prod, degr, KCroPR, tau, V0, convFac);
    wm3_wt_CroPR = getWeightsv19(tm3_wt_CroPR, ym3_wt_CroPR, ...
        n, prod, degr, KCroPR, tau, V0, convFac);
    wm4_wt_CroPR = getWeightsv19(tm4_wt_CroPR, ym4_wt_CroPR, ...
        n, prod, degr, KCroPR, tau, V0, convFac);
    wm5_wt_CroPR = getWeightsv19(tm5_wt_CroPR, ym5_wt_CroPR, ...
        n, prod, degr, KCroPR, tau, V0, convFac);
    sol(i).wm1_wt_CroPR = wm1_wt_CroPR;
    sol(i).wm2_wt_CroPR = wm2_wt_CroPR;
    sol(i).wm3_wt_CroPR = wm3_wt_CroPR;
    sol(i).wm4_wt_CroPR = wm4_wt_CroPR;
    sol(i).wm5_wt_CroPR = wm5_wt_CroPR;
    
    [r, c] = find(wm1_wt_CroPR(:, 6) >= 0.1);
    if ~isempty(r)
        sol(i).tauOn_m1_wt_CroPR = [tm1_wt_CroPR(r(1)), tm1_wt_CroPR(r(end))];
        sol(i).i_tauOn_m1_wt_CroPR = [r(1), r(end), 1];
    end;
    
    [r, c] = find(wm2_wt_CroPR(:, 6) >= 0.1);
    if ~isempty(r)
        sol(i).tauOn_m2_wt_CroPR = [tm2_wt_CroPR(r(1)), tm2_wt_CroPR(r(end))];
        sol(i).i_tauOn_m2_wt_CroPR = [r(1), r(end), 1];
    end;
    
    [r, c] = find(wm3_wt_CroPR(:, 6) >= 0.1);
    if ~isempty(r)
        sol(i).tauOn_m3_wt_CroPR = [tm3_wt_CroPR(r(1)), tm3_wt_CroPR(r(end))];
        sol(i).i_tauOn_m3_wt_CroPR = [r(1), r(end), 1];
    end;
    
    [r, c] = find(wm4_wt_CroPR(:, 6) >= 0.1);
    if ~isempty(r)
        sol(i).tauOn_m4_wt_CroPR = [tm4_wt_CroPR(r(1)), tm4_wt_CroPR(r(end))];
        sol(i).i_tauOn_m4_wt_CroPR = [r(1), r(end), 1];
    end;
    
    [r, c] = find(wm5_wt_CroPR(:, 6) >= 0.1);
    if ~isempty(r)
        sol(i).tauOn_m5_wt_CroPR = [tm5_wt_CroPR(r(1)), tm5_wt_CroPR(r(end))];
        sol(i).i_tauOn_m5_wt_CroPR = [r(1), r(end), 1];
    end;
    
    %NO CI REPRESSION OF PR (+ REPL.) =====================================
    %WT
    KCIPR = K;
    KCIPR([6, 8, 10]) = Inf;
    %MOI=1
    [tm1_wt_CIPR, ym1_wt_CIPR] = ode15s(@fv19_repv3, tspan, y0_m1, options, n, prod, ...
        degr, KCIPR, tau, V0, convFac);
    sol(i).m1_wt_CIPR = [tm1_wt_CIPR, ym1_wt_CIPR];
    sol(i).m1Num_wt_CIPR = [tm1_wt_CIPR, ym1_wt_CIPR.*V./convFac];
    %MOI=2
    [tm2_wt_CIPR, ym2_wt_CIPR] = ode15s(@fv19_repv3, tspan, y0_m2, options, n, prod, ...
        degr, KCIPR, tau, V0, convFac);
    sol(i).m2_wt_CIPR = [tm2_wt_CIPR, ym2_wt_CIPR];
    sol(i).m2Num_wt_CIPR = [tm2_wt_CIPR, ym2_wt_CIPR.*V./convFac];
    %MOI=3
    [tm3_wt_CIPR, ym3_wt_CIPR] = ode15s(@fv19_repv3, tspan, y0_m3, options, n, prod, ...
        degr, KCIPR, tau, V0, convFac);
    sol(i).m3_wt_CIPR = [tm3_wt_CIPR, ym3_wt_CIPR];
    sol(i).m3Num_wt_CIPR = [tm3_wt_CIPR, ym3_wt_CIPR.*V./convFac];
    %MOI=4
    [tm4_wt_CIPR, ym4_wt_CIPR] = ode15s(@fv19_repv3, tspan, y0_m4, options, n, prod, ...
        degr, KCIPR, tau, V0, convFac);
    sol(i).m4_wt_CIPR = [tm4_wt_CIPR, ym4_wt_CIPR];
    sol(i).m4Num_wt_CIPR = [tm4_wt_CIPR, ym4_wt_CIPR.*V./convFac];
    %MOI=5
    [tm5_wt_CIPR, ym5_wt_CIPR] = ode15s(@fv19_repv3, tspan, y0_m5, options, n, prod, ...
        degr, KCIPR, tau, V0, convFac);
    sol(i).m5_wt_CIPR = [tm5_wt_CIPR, ym5_wt_CIPR];
    sol(i).m5Num_wt_CIPR = [tm5_wt_CIPR, ym5_wt_CIPR.*V./convFac];
    
    %Neither CI or Cro repression of PR====================================
    %WT
    KCroCIPR = K;
    KCroCIPR([5, 6, 7, 8, 9, 10]) = Inf;
    %MOI=1
    [tm1_wt_CroCIPR, ym1_wt_CroCIPR] = ode15s(@fv19_repv3, tspan, y0_m1, options, n, prod, ...
        degr, KCroCIPR, tau, V0, convFac);
    sol(i).m1_wt_CroCIPR = [tm1_wt_CroCIPR, ym1_wt_CroCIPR];
    sol(i).m1Num_wt_CroCIPR = [tm1_wt_CroCIPR, ym1_wt_CroCIPR.*V./convFac];
    %MOI=2
    [tm2_wt_CroCIPR, ym2_wt_CroCIPR] = ode15s(@fv19_repv3, tspan, y0_m2, options, n, prod, ...
        degr, KCroCIPR, tau, V0, convFac);
    sol(i).m2_wt_CroCIPR = [tm2_wt_CroCIPR, ym2_wt_CroCIPR];
    sol(i).m2Num_wt_CroCIPR = [tm2_wt_CroCIPR, ym2_wt_CroCIPR.*V./convFac];
    %MOI=3
    [tm3_wt_CroCIPR, ym3_wt_CroCIPR] = ode15s(@fv19_repv3, tspan, y0_m3, options, n, prod, ...
        degr, KCroCIPR, tau, V0, convFac);
    sol(i).m3_wt_CroCIPR = [tm3_wt_CroCIPR, ym3_wt_CroCIPR];
    sol(i).m3Num_wt_CroCIPR = [tm3_wt_CroCIPR, ym3_wt_CroCIPR.*V./convFac];
    %MOI=4
    [tm4_wt_CroCIPR, ym4_wt_CroCIPR] = ode15s(@fv19_repv3, tspan, y0_m4, options, n, prod, ...
        degr, KCroCIPR, tau, V0, convFac);
    sol(i).m4_wt_CroCIPR = [tm4_wt_CroCIPR, ym4_wt_CroCIPR];
    sol(i).m4Num_wt_CroCIPR = [tm4_wt_CroCIPR, ym4_wt_CroCIPR.*V./convFac];
    %MOI=5
    [tm5_wt_CroCIPR, ym5_wt_CroCIPR] = ode15s(@fv19_repv3, tspan, y0_m5, options, n, prod, ...
        degr, KCroCIPR, tau, V0, convFac);
    sol(i).m5_wt_CroCIPR = [tm5_wt_CroCIPR, ym5_wt_CroCIPR];
    sol(i).m5Num_wt_CroCIPR = [tm5_wt_CroCIPR, ym5_wt_CroCIPR.*V./convFac];
    
    %nPRE =================================================================
    %Set to 2
    nPRE = n;
    nPRE(4) = 2;
    %P- -------------------------------------------------------------------
    %MOI=1
    [tm1_OP_nPRE, ym1_OP_nPRE] = ode15s(@fv19, tspan, y0_m1(1:6), options, nPRE([1:8, end]), ...
        prod(1:end-1), degr, K([1:8, end]), 1, V0, convFac, ss);
    sol(i).m1_nPRE = [tm1_OP_nPRE, ym1_OP_nPRE];
    sol(i).m1Num_nPRE = [tm1_OP_nPRE, ym1_OP_nPRE.*V./convFac];
    %MOI=2
    [tm2_OP_nPRE, ym2_OP_nPRE] = ode15s(@fv19, tspan, y0_m2(1:6), options, nPRE([1:8, end]), ...
        prod(1:end-1), degr, K([1:8, end]), 2, V0, convFac, ss);
    sol(i).m2_nPRE = [tm2_OP_nPRE, ym2_OP_nPRE];
    sol(i).m2Num_nPRE = [tm2_OP_nPRE, ym2_OP_nPRE.*V./convFac];
    %MOI=3
    [tm3_OP_nPRE, ym3_OP_nPRE] = ode15s(@fv19, tspan, y0_m3(1:6), options, nPRE([1:8, end]), ...
        prod(1:end-1), degr, K([1:8, end]), 3, V0, convFac, ss);
    sol(i).m3_nPRE = [tm3_OP_nPRE, ym3_OP_nPRE];
    sol(i).m3Num_nPRE = [tm3_OP_nPRE, ym3_OP_nPRE.*V./convFac];
    %MOI=4
    [tm4_OP_nPRE, ym4_OP_nPRE] = ode15s(@fv19, tspan, y0_m4(1:6), options, nPRE([1:8, end]), ...
        prod(1:end-1), degr, K([1:8, end]), 4, V0, convFac, ss);
    sol(i).m4_nPRE = [tm4_OP_nPRE, ym4_OP_nPRE];
    sol(i).m4Num_nPRE = [tm4_OP_nPRE, ym4_OP_nPRE.*V./convFac];
    %MOI=1
    [tm5_OP_nPRE, ym5_OP_nPRE] = ode15s(@fv19, tspan, y0_m5(1:6), options, nPRE([1:8, end]), ...
        prod(1:end-1), degr, K([1:8, end]), 5, V0, convFac, ss);
    sol(i).m5_nPRE = [tm5_OP_nPRE, ym5_OP_nPRE];
    sol(i).m5Num_nPRE = [tm5_OP_nPRE, ym5_OP_nPRE.*V./convFac];
    
    %weights
    wm1_nPRE = getWeightsv19(tm1_OP, ym1_OP, ...
        nPRE, prod, degr, K, tau, V0, convFac);
    wm2_nPRE = getWeightsv19(tm2_OP, ym2_OP, ...
        nPRE, prod, degr, K, tau, V0, convFac);
    wm3_nPRE = getWeightsv19(tm3_OP, ym3_OP, ...
        nPRE, prod, degr, K, tau, V0, convFac);
    wm4_nPRE = getWeightsv19(tm4_OP, ym4_OP, ...
        nPRE, prod, degr, K, tau, V0, convFac);
    wm5_nPRE = getWeightsv19(tm5_OP, ym5_OP, ...
        nPRE, prod, degr, K, tau, V0, convFac);
    sol(i).wm1_nPRE = wm1_nPRE;
    sol(i).wm2_nPRE = wm2_nPRE;
    sol(i).wm3_nPRE = wm3_nPRE;
    sol(i).wm4_nPRE = wm4_nPRE;
    sol(i).wm5_nPRE = wm5_nPRE;
    
    [r, c] = find(wm1_nPRE(:, 6) >= 0.1);
    if ~isempty(r)
        sol(i).tauOn_m1_nPRE = [tm1_OP_nPRE(r(1)), tm1_OP_nPRE(r(end))];
        sol(i).i_tauOn_m1_nPRE = [r(1), r(end), 1];
    end;
    
    [r, c] = find(wm2_nPRE(:, 6) >= 0.1);
    if ~isempty(r)
        sol(i).tauOn_m2_nPRE = [tm2_OP_nPRE(r(1)), tm2_OP_nPRE(r(end))];
        sol(i).i_tauOn_m2_nPRE = [r(1), r(end), 1];
    end;
    
    [r, c] = find(wm3_nPRE(:, 6) >= 0.1);
    if ~isempty(r)
        sol(i).tauOn_m3_nPRE = [tm3_OP_nPRE(r(1)), tm3_OP_nPRE(r(end))];
        sol(i).i_tauOn_m3_nPRE = [r(1), r(end), 1];
    end;
    
    [r, c] = find(wm4_nPRE(:, 6) >= 0.1);
    if ~isempty(r)
        sol(i).tauOn_m4_nPRE = [tm4_OP_nPRE(r(1)), tm4_OP_nPRE(r(end))];
        sol(i).i_tauOn_m4_nPRE = [r(1), r(end), 1];
    end;
    
    [r, c] = find(wm5_nPRE(:, 6) >= 0.1);
    if ~isempty(r)
        sol(i).tauOn_m5_nPRE = [tm5_OP_nPRE(r(1)), tm5_OP_nPRE(r(end))];
        sol(i).i_tauOn_m5_nPRE = [r(1), r(end), 1];
    end;
    
    %WT -------------------------------------------------------------------
    %MOI=1
    [tm1_wt_nPRE, ym1_wt_nPRE] = ode15s(@fv19_repv3, tspan, y0_m1, options, nPRE, prod, ...
        degr, K, tau, V0, convFac);
    sol(i).m1_wt_nPRE = [tm1_wt_nPRE, ym1_wt_nPRE];
    sol(i).m1Num_wt_nPRE = [tm1_wt_nPRE, ym1_wt_nPRE.*V./convFac];
    %MOI=2
    [tm2_wt_nPRE, ym2_wt_nPRE] = ode15s(@fv19_repv3, tspan, y0_m2, options, nPRE, prod, ...
        degr, K, tau, V0, convFac);
    sol(i).m2_wt_nPRE = [tm2_wt_nPRE, ym2_wt_nPRE];
    sol(i).m2Num_wt_nPRE = [tm2_wt_nPRE, ym2_wt_nPRE.*V./convFac];
    %MOI=3
    [tm3_wt_nPRE, ym3_wt_nPRE] = ode15s(@fv19_repv3, tspan, y0_m3, options, nPRE, prod, ...
        degr, K, tau, V0, convFac);
    sol(i).m3_wt_nPRE = [tm3_wt_nPRE, ym3_wt_nPRE];
    sol(i).m3Num_wt_nPRE = [tm3_wt_nPRE, ym3_wt_nPRE.*V./convFac];
    %MOI=4
    [tm4_wt_nPRE, ym4_wt_nPRE] = ode15s(@fv19_repv3, tspan, y0_m4, options, nPRE, prod, ...
        degr, K, tau, V0, convFac);
    sol(i).m4_wt_nPRE = [tm4_wt_nPRE, ym4_wt_nPRE];
    sol(i).m4Num_wt_nPRE = [tm4_wt_nPRE, ym4_wt_nPRE.*V./convFac];
    %MOI=5
    [tm5_wt_nPRE, ym5_wt_nPRE] = ode15s(@fv19_repv3, tspan, y0_m5, options, nPRE, prod, ...
        degr, K, tau, V0, convFac);
    sol(i).m5_wt_nPRE = [tm5_wt_nPRE, ym5_wt_nPRE];
    sol(i).m5Num_wt_nPRE = [tm5_wt_nPRE, ym5_wt_nPRE.*V./convFac];
    
    %weights
    wm1_wt_nPRE = getWeightsv19(tm1_wt_nPRE, ym1_wt_nPRE, ...
        nPRE, prod, degr, K, tau, V0, convFac);
    wm2_wt_nPRE = getWeightsv19(tm2_wt_nPRE, ym2_wt_nPRE, ...
        nPRE, prod, degr, K, tau, V0, convFac);
    wm3_wt_nPRE = getWeightsv19(tm3_wt_nPRE, ym3_wt_nPRE, ...
        nPRE, prod, degr, K, tau, V0, convFac);
    wm4_wt_nPRE = getWeightsv19(tm4_wt_nPRE, ym4_wt_nPRE, ...
        nPRE, prod, degr, K, tau, V0, convFac);
    wm5_wt_nPRE = getWeightsv19(tm5_wt_nPRE, ym5_wt_nPRE, ...
        nPRE, prod, degr, K, tau, V0, convFac);
    sol(i).wm1_wt_nPRE = wm1_wt_nPRE;
    sol(i).wm2_wt_nPRE = wm2_wt_nPRE;
    sol(i).wm3_wt_nPRE = wm3_wt_nPRE;
    sol(i).wm4_wt_nPRE = wm4_wt_nPRE;
    sol(i).wm5_wt_nPRE = wm5_wt_nPRE;
    
    [r, c] = find(wm1_wt_nPRE(:, 6) >= 0.1);
    if ~isempty(r)
        sol(i).tauOn_m1_wt_nPRE = [tm1_wt_nPRE(r(1)), tm1_wt_nPRE(r(end))];
        sol(i).i_tauOn_m1_wt_nPRE = [r(1), r(end), 1];
    end;
    
    [r, c] = find(wm2_wt_nPRE(:, 6) >= 0.1);
    if ~isempty(r)
        sol(i).tauOn_m2_wt_nPRE = [tm2_wt_nPRE(r(1)), tm2_wt_nPRE(r(end))];
        sol(i).i_tauOn_m2_wt_nPRE = [r(1), r(end), 1];
    end;
    
    [r, c] = find(wm3_wt_nPRE(:, 6) >= 0.1);
    if ~isempty(r)
        sol(i).tauOn_m3_wt_nPRE = [tm3_wt_nPRE(r(1)), tm3_wt_nPRE(r(end))];
        sol(i).i_tauOn_m3_wt_nPRE = [r(1), r(end), 1];
    end;
    
    [r, c] = find(wm4_wt_nPRE(:, 6) >= 0.1);
    if ~isempty(r)
        sol(i).tauOn_m4_wt_nPRE = [tm4_wt_nPRE(r(1)), tm4_wt_nPRE(r(end))];
        sol(i).i_tauOn_m4_wt_nPRE = [r(1), r(end), 1];
    end;
    
    [r, c] = find(wm5_wt_nPRE(:, 6) >= 0.1);
    if ~isempty(r)
        sol(i).tauOn_m5_wt_nPRE = [tm5_wt_nPRE(r(1)), tm5_wt_nPRE(r(end))];
        sol(i).i_tauOn_m5_wt_nPRE = [r(1), r(end), 1];
    end;
    
    %EXTEND CII LIFETIME===================================================
    degrCIIext = degr;
    degrCIIext(end-1) = 0.5*degrCIIext(end-1);
    %WT--------------------------------------------------------------------
    %MOI=1
    [tm1_wt_CIIext, ym1_wt_CIIext] = ode15s(@fv19_repv3, tspan, y0_m1, ...
        options, n, prod, degrCIIext, K, tau, V0, convFac);
    sol(i).m1_wt_CIIext = [tm1_wt_CIIext, ym1_wt_CIIext];
    sol(i).m1Num_wt_CIIext = [tm1_wt_CIIext, ym1_wt_CIIext.*V./convFac];
    %MOI=2
    [tm2_wt_CIIext, ym2_wt_CIIext] = ode15s(@fv19_repv3, tspan, y0_m2, ...
        options, n, prod, degrCIIext, K, tau, V0, convFac);
    sol(i).m2_wt_CIIext = [tm2_wt_CIIext, ym2_wt_CIIext];
    sol(i).m2Num_wt_CIIext = [tm2_wt_CIIext, ym2_wt_CIIext.*V./convFac];
    %MOI=3
    [tm3_wt_CIIext, ym3_wt_CIIext] = ode15s(@fv19_repv3, tspan, y0_m3, ...
        options, n, prod, degrCIIext, K, tau, V0, convFac);
    sol(i).m3_wt_CIIext = [tm3_wt_CIIext, ym3_wt_CIIext];
    sol(i).m3Num_wt_CIIext = [tm3_wt_CIIext, ym3_wt_CIIext.*V./convFac];
    %MOI=4
    [tm4_wt_CIIext, ym4_wt_CIIext] = ode15s(@fv19_repv3, tspan, y0_m4, ...
        options, n, prod, degrCIIext, K, tau, V0, convFac);
    sol(i).m4_wt_CIIext = [tm4_wt_CIIext, ym4_wt_CIIext];
    sol(i).m4Num_wt_CIIext = [tm4_wt_CIIext, ym4_wt_CIIext.*V./convFac];
    %MOI=5
    [tm5_wt_CIIext, ym5_wt_CIIext] = ode15s(@fv19_repv3, tspan, y0_m5, ...
        options, n, prod, degrCIIext, K, tau, V0, convFac);
    sol(i).m5_wt_CIIext = [tm5_wt_CIIext, ym5_wt_CIIext];
    sol(i).m5Num_wt_CIIext = [tm5_wt_CIIext, ym5_wt_CIIext.*V./convFac];
    
    %weights
    wm1_wt_CIIext = getWeightsv19(tm1_wt_CIIext, ym1_wt_CIIext, ...
        n, prod, degrCIIext, K, tau, V0, convFac);
    wm2_wt_CIIext = getWeightsv19(tm2_wt_CIIext, ym2_wt_CIIext, ...
        n, prod, degrCIIext, K, tau, V0, convFac);
    wm3_wt_CIIext = getWeightsv19(tm3_wt_CIIext, ym3_wt_CIIext, ...
        n, prod, degrCIIext, K, tau, V0, convFac);
    wm4_wt_CIIext = getWeightsv19(tm4_wt_CIIext, ym4_wt_CIIext, ...
        n, prod, degrCIIext, K, tau, V0, convFac);
    wm5_wt_CIIext = getWeightsv19(tm5_wt_CIIext, ym5_wt_CIIext, ...
        n, prod, degrCIIext, K, tau, V0, convFac);
    sol(i).wm1_wt_CIIext = wm1_wt_CIIext;
    sol(i).wm2_wt_CIIext = wm2_wt_CIIext;
    sol(i).wm3_wt_CIIext = wm3_wt_CIIext;
    sol(i).wm4_wt_CIIext = wm4_wt_CIIext;
    sol(i).wm5_wt_CIIext = wm5_wt_CIIext;
    
    [r, c] = find(wm1_wt_CIIext(:, 6) >= 0.1);
    if ~isempty(r)
        sol(i).tauOn_m1_wt_CIIext = [tm1_wt_CIIext(r(1)), tm1_wt_CIIext(r(end))];
        sol(i).i_tauOn_m1_wt_CIIext = [r(1), r(end), 1];
    end;
    
    [r, c] = find(wm2_wt_CIIext(:, 6) >= 0.1);
    if ~isempty(r)
        sol(i).tauOn_m2_wt_CIIext = [tm2_wt_CIIext(r(1)), tm2_wt_CIIext(r(end))];
        sol(i).i_tauOn_m2_wt_CIIext = [r(1), r(end), 1];
    end;
    
    [r, c] = find(wm3_wt_CIIext(:, 6) >= 0.1);
    if ~isempty(r)
        sol(i).tauOn_m3_wt_CIIext = [tm3_wt_CIIext(r(1)), tm3_wt_CIIext(r(end))];
        sol(i).i_tauOn_m3_wt_CIIext = [r(1), r(end), 1];
    end;
    
    [r, c] = find(wm4_wt_CIIext(:, 6) >= 0.1);
    if ~isempty(r)
        sol(i).tauOn_m4_wt_CIIext = [tm4_wt_CIIext(r(1)), tm4_wt_CIIext(r(end))];
        sol(i).i_tauOn_m4_wt_CIIext = [r(1), r(end), 1];
    end;
    
    [r, c] = find(wm5_wt_CIIext(:, 6) >= 0.1);
    if ~isempty(r)
        sol(i).tauOn_m5_wt_CIIext = [tm5_wt_CIIext(r(1)), tm5_wt_CIIext(r(end))];
        sol(i).i_tauOn_m5_wt_CIIext = [r(1), r(end), 1];
    end;
    
    %nDegCII=5 ============================================================
    nDegCII = n;
    nDegCII(end) = 5;
    %WT -------------------------------------------------------------------
    %MOI=1
    [tm1_wt_nDegCII, ym1_wt_nDegCII] = ode15s(@fv19_repv3, tspan, y0_m1, options, ...
        nDegCII, prod, degr, K, tau, V0, convFac);
    sol(i).m1_wt_nDegCII = [tm1_wt_nDegCII, ym1_wt_nDegCII];
    sol(i).m1Num_wt_nDegCII = [tm1_wt_nDegCII, ym1_wt_nDegCII.*V./convFac];
    %MOI=2
    [tm2_wt_nDegCII, ym2_wt_nDegCII] = ode15s(@fv19_repv3, tspan, y0_m2, options, ...
        nDegCII, prod, degr, K, tau, V0, convFac);
    sol(i).m2_wt_nDegCII = [tm2_wt_nDegCII, ym2_wt_nDegCII];
    sol(i).m2Num_wt_nDegCII = [tm2_wt_nDegCII, ym2_wt_nDegCII.*V./convFac];
    %MOI=3
    [tm3_wt_nDegCII, ym3_wt_nDegCII] = ode15s(@fv19_repv3, tspan, y0_m3, options, ...
        nDegCII, prod, degr, K, tau, V0, convFac);
    sol(i).m3_wt_nDegCII = [tm3_wt_nDegCII, ym3_wt_nDegCII];
    sol(i).m3Num_wt_nDegCII = [tm3_wt_nDegCII, ym3_wt_nDegCII.*V./convFac];
    %MOI=4
    [tm4_wt_nDegCII, ym4_wt_nDegCII] = ode15s(@fv19_repv3, tspan, y0_m4, options, ...
        nDegCII, prod, degr, K, tau, V0, convFac);
    sol(i).m4_wt_nDegCII = [tm4_wt_nDegCII, ym4_wt_nDegCII];
    sol(i).m4Num_wt_nDegCII = [tm4_wt_nDegCII, ym4_wt_nDegCII.*V./convFac];
    %MOI=5
    [tm5_wt_nDegCII, ym5_wt_nDegCII] = ode15s(@fv19_repv3, tspan, y0_m5, options, ...
        nDegCII, prod, degr, K, tau, V0, convFac);
    sol(i).m5_wt_nDegCII = [tm5_wt_nDegCII, ym5_wt_nDegCII];
    sol(i).m5Num_wt_nDegCII = [tm5_wt_nDegCII, ym5_wt_nDegCII.*V./convFac];
    
    %weights
    wm1_wt_nDegCII = getWeightsv19(tm1_wt, ym1_wt, ...
        nDegCII, prod, degr, K, tau, V0, convFac);
    wm2_wt_nDegCII = getWeightsv19(tm2_wt, ym2_wt, ...
        nDegCII, prod, degr, K, tau, V0, convFac);
    wm3_wt_nDegCII = getWeightsv19(tm3_wt, ym3_wt, ...
        nDegCII, prod, degr, K, tau, V0, convFac);
    wm4_wt_nDegCII = getWeightsv19(tm4_wt, ym4_wt, ...
        nDegCII, prod, degr, K, tau, V0, convFac);
    wm5_wt_nDegCII = getWeightsv19(tm5_wt, ym5_wt, ...
        nDegCII, prod, degr, K, tau, V0, convFac);
    sol(i).wm1_wt_nDegCII = wm1_wt_nDegCII;
    sol(i).wm2_wt_nDegCII = wm2_wt_nDegCII;
    sol(i).wm3_wt_nDegCII = wm3_wt_nDegCII;
    sol(i).wm4_wt_nDegCII = wm4_wt_nDegCII;
    sol(i).wm5_wt_nDegCII = wm5_wt_nDegCII;
    
    [r, c] = find(wm1_wt_nDegCII(:, 6) >= 0.1);
    if ~isempty(r)
        sol(i).tauOn_m1_wt_nDegCII = [tm1_wt_nDegCII(r(1)), tm1_wt_nDegCII(r(end))];
        sol(i).i_tauOn_m1_wt_nDegCII = [r(1), r(end), 1];
    end;
    
    [r, c] = find(wm2_wt_nDegCII(:, 6) >= 0.1);
    if ~isempty(r)
        sol(i).tauOn_m2_wt_nDegCII = [tm2_wt_nDegCII(r(1)), tm2_wt_nDegCII(r(end))];
        sol(i).i_tauOn_m2_wt_nDegCII = [r(1), r(end), 1];
    end;
    
    [r, c] = find(wm3_wt_nDegCII(:, 6) >= 0.1);
    if ~isempty(r)
        sol(i).tauOn_m3_wt_nDegCII = [tm3_wt_nDegCII(r(1)), tm3_wt_nDegCII(r(end))];
        sol(i).i_tauOn_m3_wt_nDegCII = [r(1), r(end), 1];
    end;
    
    [r, c] = find(wm4_wt_nDegCII(:, 6) >= 0.1);
    if ~isempty(r)
        sol(i).tauOn_m4_wt_nDegCII = [tm4_wt_nDegCII(r(1)), tm4_wt_nDegCII(r(end))];
        sol(i).i_tauOn_m4_wt_nDegCII = [r(1), r(end), 1];
    end;
    
    [r, c] = find(wm5_wt_nDegCII(:, 6) >= 0.1);
    if ~isempty(r)
        sol(i).tauOn_m5_wt_nDegCII = [tm5_wt_nDegCII(r(1)), tm5_wt_nDegCII(r(end))];
        sol(i).i_tauOn_m5_wt_nDegCII = [r(1), r(end), 1];
    end;
    
    %CII SWITCH-LIE (nCII = 5, KCII = 2KCII) ==============================
    KDegCII = K;
    KDegCII(end) = 1.5*KDegCII(end);
    %P- -------------------------------------------------------------------
    %MOI=1
    [tm1_OP_CIIDeg, ym1_OP_CIIDeg] = ode15s(@fv19, tspan, y0_m1(1:6), options, ...
        nDegCII([1:8, end]), prod(1:end-1), degr, KDegCII([1:8, end]), 1, V0, convFac, ss);
    sol(i).m1_CIIDeg = [tm1_OP_CIIDeg, ym1_OP_CIIDeg];
    sol(i).m1Num_CIIDeg = [tm1_OP_CIIDeg, ym1_OP_CIIDeg.*V./convFac];
    %MOI=2
    [tm2_OP_CIIDeg, ym2_OP_CIIDeg] = ode15s(@fv19, tspan, y0_m2(1:6), options, ...
        nDegCII([1:8, end]), prod(1:end-1), degr, KDegCII([1:8, end]), 2, V0, convFac, ss);
    sol(i).m2_CIIDeg = [tm2_OP_CIIDeg, ym2_OP_CIIDeg];
    sol(i).m2Num_CIIDeg = [tm2_OP_CIIDeg, ym2_OP_CIIDeg.*V./convFac];
    %MOI=3
    [tm3_OP_CIIDeg, ym3_OP_CIIDeg] = ode15s(@fv19, tspan, y0_m3(1:6), options, ...
        nDegCII([1:8, end]), prod(1:end-1), degr, KDegCII([1:8, end]), 3, V0, convFac, ss);
    sol(i).m3_CIIDeg = [tm3_OP_CIIDeg, ym3_OP_CIIDeg];
    sol(i).m3Num_CIIDeg = [tm3_OP_CIIDeg, ym3_OP_CIIDeg.*V./convFac];
    %MOI=4
    [tm4_OP_CIIDeg, ym4_OP_CIIDeg] = ode15s(@fv19, tspan, y0_m4(1:6), options, ...
        nDegCII([1:8, end]), prod(1:end-1), degr, KDegCII([1:8, end]), 4, V0, convFac, ss);
    sol(i).m4_CIIDeg = [tm4_OP_CIIDeg, ym4_OP_CIIDeg];
    sol(i).m4Num_CIIDeg = [tm4_OP_CIIDeg, ym4_OP_CIIDeg.*V./convFac];
    %MOI=5
    [tm5_OP_CIIDeg, ym5_OP_CIIDeg] = ode15s(@fv19, tspan, y0_m5(1:6), options, ...
        nDegCII([1:8, end]), prod(1:end-1), degr, KDegCII([1:8, end]), 5, V0, convFac, ss);
    sol(i).m5_CIIDeg = [tm5_OP_CIIDeg, ym5_OP_CIIDeg];
    sol(i).m5Num_CIIDeg = [tm5_OP_CIIDeg, ym5_OP_CIIDeg.*V./convFac];
    
    %weights
    wm1_CIIDeg = getWeightsv19(tm1_OP_CIIDeg, ym1_OP_CIIDeg, ...
        nDegCII, prod, degr, KDegCII, tau, V0, convFac);
    wm2_CIIDeg = getWeightsv19(tm2_OP_CIIDeg, ym2_OP_CIIDeg, ...
        nDegCII, prod, degr, KDegCII, tau, V0, convFac);
    wm3_CIIDeg = getWeightsv19(tm3_OP_CIIDeg, ym3_OP_CIIDeg, ...
        nDegCII, prod, degr, KDegCII, tau, V0, convFac);
    wm4_CIIDeg = getWeightsv19(tm4_OP_CIIDeg, ym4_OP_CIIDeg, ...
        nDegCII, prod, degr, KDegCII, tau, V0, convFac);
    wm5_CIIDeg = getWeightsv19(tm5_OP_CIIDeg, ym5_OP_CIIDeg, ...
        nDegCII, prod, degr, KDegCII, tau, V0, convFac);
    sol(i).wm1_CIIDeg = wm1_CIIDeg;
    sol(i).wm2_CIIDeg = wm2_CIIDeg;
    sol(i).wm3_CIIDeg = wm3_CIIDeg;
    sol(i).wm4_CIIDeg = wm4_CIIDeg;
    sol(i).wm5_CIIDeg = wm5_CIIDeg;
    
    [r, c] = find(wm1_CIIDeg(:, 6) >= 0.1);
    if ~isempty(r)
        sol(i).tauOn_m1_CIIDeg = [tm1_OP_CIIDeg(r(1)), tm1_OP_CIIDeg(r(end))];
        sol(i).i_tauOn_m1_CIIDeg = [r(1), r(end), 1];
    end;
    
    [r, c] = find(wm2_CIIDeg(:, 6) >= 0.1);
    if ~isempty(r)
        sol(i).tauOn_m2_CIIDeg = [tm2_OP_CIIDeg(r(1)), tm2_OP_CIIDeg(r(end))];
        sol(i).i_tauOn_m2_CIIDeg = [r(1), r(end), 1];
    end;
    
    [r, c] = find(wm3_CIIDeg(:, 6) >= 0.1);
    if ~isempty(r)
        sol(i).tauOn_m3_CIIDeg = [tm3_OP_CIIDeg(r(1)), tm3_OP_CIIDeg(r(end))];
        sol(i).i_tauOn_m3_CIIDeg = [r(1), r(end), 1];
    end;
    
    [r, c] = find(wm4_CIIDeg(:, 6) >= 0.1);
    if ~isempty(r)
        sol(i).tauOn_m4_CIIDeg = [tm4_OP_CIIDeg(r(1)), tm4_OP_CIIDeg(r(end))];
        sol(i).i_tauOn_m4_CIIDeg = [r(1), r(end), 1];
    end;
    
    [r, c] = find(wm5_CIIDeg(:, 6) >= 0.1);
    if ~isempty(r)
        sol(i).tauOn_m5_CIIDeg = [tm5_OP_CIIDeg(r(1)), tm5_OP_CIIDeg(r(end))];
        sol(i).i_tauOn_m5_CIIDeg = [r(1), r(end), 1];
    end;
    
    %WT -------------------------------------------------------------------
    %MOI=1
    [tm1_wt_CIIDeg, ym1_wt_CIIDeg] = ode15s(@fv19_repv3, tspan, y0_m1, options, ...
        nDegCII, prod, degr, KDegCII, tau, V0, convFac);
    sol(i).m1_wt_CIIDeg = [tm1_wt_CIIDeg, ym1_wt_CIIDeg];
    sol(i).m1Num_wt_CIIDeg = [tm1_wt_CIIDeg, ym1_wt_CIIDeg.*V./convFac];
    %MOI=2
    [tm2_wt_CIIDeg, ym2_wt_CIIDeg] = ode15s(@fv19_repv3, tspan, y0_m2, options, ...
        nDegCII, prod, degr, KDegCII, tau, V0, convFac);
    sol(i).m2_wt_CIIDeg = [tm2_wt_CIIDeg, ym2_wt_CIIDeg];
    sol(i).m2Num_wt_CIIDeg = [tm2_wt_CIIDeg, ym2_wt_CIIDeg.*V./convFac];
    %MOI=3
    [tm3_wt_CIIDeg, ym3_wt_CIIDeg] = ode15s(@fv19_repv3, tspan, y0_m3, options, ...
        nDegCII, prod, degr, KDegCII, tau, V0, convFac);
    sol(i).m3_wt_CIIDeg = [tm3_wt_CIIDeg, ym3_wt_CIIDeg];
    sol(i).m3Num_wt_CIIDeg = [tm3_wt_CIIDeg, ym3_wt_CIIDeg.*V./convFac];
    %MOI=4
    [tm4_wt_CIIDeg, ym4_wt_CIIDeg] = ode15s(@fv19_repv3, tspan, y0_m4, options, ...
        nDegCII, prod, degr, KDegCII, tau, V0, convFac);
    sol(i).m4_wt_CIIDeg = [tm4_wt_CIIDeg, ym4_wt_CIIDeg];
    sol(i).m4Num_wt_CIIDeg = [tm4_wt_CIIDeg, ym4_wt_CIIDeg.*V./convFac];
    %MOI=5
    [tm5_wt_CIIDeg, ym5_wt_CIIDeg] = ode15s(@fv19_repv3, tspan, y0_m5, options, ...
        nDegCII, prod, degr, KDegCII, tau, V0, convFac);
    sol(i).m5_wt_CIIDeg = [tm5_wt_CIIDeg, ym5_wt_CIIDeg];
    sol(i).m5Num_wt_CIIDeg = [tm5_wt_CIIDeg, ym5_wt_CIIDeg.*V./convFac];
    
    %weights
    wm1_wt_CIIDeg = getWeightsv19(tm1_wt_CIIDeg, ym1_wt_CIIDeg, ...
        nDegCII, prod, degr, KDegCII, tau, V0, convFac);
    wm2_wt_CIIDeg = getWeightsv19(tm2_wt_CIIDeg, ym2_wt_CIIDeg, ...
        nDegCII, prod, degr, KDegCII, tau, V0, convFac);
    wm3_wt_CIIDeg = getWeightsv19(tm3_wt_CIIDeg, ym3_wt_CIIDeg, ...
        nDegCII, prod, degr, KDegCII, tau, V0, convFac);
    wm4_wt_CIIDeg = getWeightsv19(tm4_wt_CIIDeg, ym4_wt_CIIDeg, ...
        nDegCII, prod, degr, KDegCII, tau, V0, convFac);
    wm5_wt_CIIDeg = getWeightsv19(tm5_wt_CIIDeg, ym5_wt_CIIDeg, ...
        nDegCII, prod, degr, KDegCII, tau, V0, convFac);
    sol(i).wm1_wt_CIIDeg = wm1_wt_CIIDeg;
    sol(i).wm2_wt_CIIDeg = wm2_wt_CIIDeg;
    sol(i).wm3_wt_CIIDeg = wm3_wt_CIIDeg;
    sol(i).wm4_wt_CIIDeg = wm4_wt_CIIDeg;
    sol(i).wm5_wt_CIIDeg = wm5_wt_CIIDeg;
    
    [r, c] = find(wm1_wt_CIIDeg(:, 6) >= 0.1);
    if ~isempty(r)
        sol(i).tauOn_m1_wt_CIIDeg = [tm1_wt_CIIDeg(r(1)), tm1_wt_CIIDeg(r(end))];
        sol(i).i_tauOn_m1_wt_CIIDeg = [r(1), r(end), 1];
    end;
    
    [r, c] = find(wm2_wt_CIIDeg(:, 6) >= 0.1);
    if ~isempty(r)
        sol(i).tauOn_m2_wt_CIIDeg = [tm2_wt_CIIDeg(r(1)), tm2_wt_CIIDeg(r(end))];
        sol(i).i_tauOn_m2_wt_CIIDeg = [r(1), r(end), 1];
    end;
    
    [r, c] = find(wm3_wt_CIIDeg(:, 6) >= 0.1);
    if ~isempty(r)
        sol(i).tauOn_m3_wt_CIIDeg = [tm3_wt_CIIDeg(r(1)), tm3_wt_CIIDeg(r(end))];
        sol(i).i_tauOn_m3_wt_CIIDeg = [r(1), r(end), 1];
    end;
    
    [r, c] = find(wm4_wt_CIIDeg(:, 6) >= 0.1);
    if ~isempty(r)
        sol(i).tauOn_m4_wt_CIIDeg = [tm4_wt_CIIDeg(r(1)), tm4_wt_CIIDeg(r(end))];
        sol(i).i_tauOn_m4_wt_CIIDeg = [r(1), r(end), 1];
    end;
    
    [r, c] = find(wm5_wt_CIIDeg(:, 6) >= 0.1);
    if ~isempty(r)
        sol(i).tauOn_m5_wt_CIIDeg = [tm5_wt_CIIDeg(r(1)), tm5_wt_CIIDeg(r(end))];
        sol(i).i_tauOn_m5_wt_CIIDeg = [r(1), r(end), 1];
    end;
    
    %NO CII DEGR. PROT ====================================================
    KCIII = K;
    KCIII(end) = Inf;
    %WT -------------------------------------------------------------------
    %MOI=1
    [tm1_wt_CIII, ym1_wt_CIII] = ode15s(@fv19_repv3, tspan, y0_m1, options, n, prod, ...
        degr, KCIII, tau, V0, convFac);
    sol(i).m1_wt_CIII = [tm1_wt_CIII, ym1_wt_CIII];
    sol(i).m1Num_wt_CIII = [tm1_wt_CIII, ym1_wt_CIII.*V./convFac];
    %MOI=2
    [tm2_wt_CIII, ym2_wt_CIII] = ode15s(@fv19_repv3, tspan, y0_m2, options, n, prod, ...
        degr, KCIII, tau, V0, convFac);
    sol(i).m2_wt_CIII = [tm2_wt_CIII, ym2_wt_CIII];
    sol(i).m2Num_wt_CIII = [tm2_wt_CIII, ym2_wt_CIII.*V./convFac];
    %MOI=3
    [tm3_wt_CIII, ym3_wt_CIII] = ode15s(@fv19_repv3, tspan, y0_m3, options, n, prod, ...
        degr, KCIII, tau, V0, convFac);
    sol(i).m3_wt_CIII = [tm3_wt_CIII, ym3_wt_CIII];
    sol(i).m3Num_wt_CIII = [tm3_wt_CIII, ym3_wt_CIII.*V./convFac];
    %MOI=4
    [tm4_wt_CIII, ym4_wt_CIII] = ode15s(@fv19_repv3, tspan, y0_m4, options, n, prod, ...
        degr, KCIII, tau, V0, convFac);
    sol(i).m4_wt_CIII = [tm4_wt_CIII, ym4_wt_CIII];
    sol(i).m4Num_wt_CIII = [tm4_wt_CIII, ym4_wt_CIII.*V./convFac];
    %MOI=5
    [tm5_wt_CIII, ym5_wt_CIII] = ode15s(@fv19_repv3, tspan, y0_m5, options, n, prod, ...
        degr, KCIII, tau, V0, convFac);
    sol(i).m5_wt_CIII = [tm5_wt_CIII, ym5_wt_CIII];
    sol(i).m5Num_wt_CIII = [tm5_wt_CIII, ym5_wt_CIII.*V./convFac];
    
    %weights
    wm1_wt_CIII = getWeightsv19(tm1_wt_CIII, ym1_wt_CIII, ...
        n, prod, degr, KCIII, tau, V0, convFac);
    wm2_wt_CIII = getWeightsv19(tm2_wt_CIII, ym2_wt_CIII, ...
        n, prod, degr, KCIII, tau, V0, convFac);
    wm3_wt_CIII = getWeightsv19(tm3_wt_CIII, ym3_wt_CIII, ...
        n, prod, degr, KCIII, tau, V0, convFac);
    wm4_wt_CIII = getWeightsv19(tm4_wt_CIII, ym4_wt_CIII, ...
        n, prod, degr, KCIII, tau, V0, convFac);
    wm5_wt_CIII = getWeightsv19(tm5_wt_CIII, ym5_wt_CIII, ...
        n, prod, degr, KCIII, tau, V0, convFac);
    sol(i).wm1_wt_CIII = wm1_wt_CIII;
    sol(i).wm2_wt_CIII = wm2_wt_CIII;
    sol(i).wm3_wt_CIII = wm3_wt_CIII;
    sol(i).wm4_wt_CIII = wm4_wt_CIII;
    sol(i).wm5_wt_CIII = wm5_wt_CIII;
    
    [r, c] = find(wm1_wt_CIII(:, 6) >= 0.1);
    if ~isempty(r)
        sol(i).tauOn_m1_wt_CIII = [tm1_wt_CIII(r(1)), tm1_wt_CIII(r(end))];
        sol(i).i_tauOn_m1_wt_CIII = [r(1), r(end), 1];
    end;
    
    [r, c] = find(wm2_wt_CIII(:, 6) >= 0.1);
    if ~isempty(r)
        sol(i).tauOn_m2_wt_CIII = [tm2_wt_CIII(r(1)), tm2_wt_CIII(r(end))];
        sol(i).i_tauOn_m2_wt_CIII = [r(1), r(end), 1];
    end;
    
    [r, c] = find(wm3_wt_CIII(:, 6) >= 0.1);
    if ~isempty(r)
        sol(i).tauOn_m3_wt_CIII = [tm3_wt_CIII(r(1)), tm3_wt_CIII(r(end))];
        sol(i).i_tauOn_m3_wt_CIII = [r(1), r(end), 1];
    end;
    
    [r, c] = find(wm4_wt_CIII(:, 6) >= 0.1);
    if ~isempty(r)
        sol(i).tauOn_m4_wt_CIII = [tm4_wt_CIII(r(1)), tm4_wt_CIII(r(end))];
        sol(i).i_tauOn_m4_wt_CIII = [r(1), r(end), 1];
    end;
    
    [r, c] = find(wm5_wt_CIII(:, 6) >= 0.1);
    if ~isempty(r)
        sol(i).tauOn_m5_wt_CIII = [tm5_wt_CIII(r(1)), tm5_wt_CIII(r(end))];
        sol(i).i_tauOn_m5_wt_CIII = [r(1), r(end), 1];
    end;
    
end

%%
%FUNCTIONS

%--------------------------------------------------------------------------
function [b, xSS] = getBifurc(sol, prod, degr, K, n, V0, x0, b, convFac, boolBif)
%Computes steadystates (xSS) as a function of b (bifurcation parameter).
%Either uses last point of xSS as starting point of future sims, or uses x0

options = odeset('Nonnegative', [], 'RelTol', 1e-6, ...
    'AbsTol', 1e-6);

optionsSS = optimoptions('fsolve', 'Display', 'iter', 'OptimalityTolerance', ...
    1e-8);
tspanSS = 0:0.1:6000;
degrBif = degr;
xSS = zeros(length(b), 6);

for i = 1:length(b)
    degrBif(3) = b(i)*degr(3);
    if i == 1 || boolBif == 1
        [tm1, ym1] = ode15s(@fv19, tspanSS, x0, options, ...
            n([1:8, end]), prod(1:end-1), degrBif, K([1:8, end]), 1, V0, convFac, 1);
    else
        [tm1, ym1] = ode15s(@fv19, tspanSS, xSS(i-1, :), options, ...
            n([1:8, end]), prod(1:end-1), degrBif, K([1:8, end]), 1, V0, convFac, 1);
    end
    [x, fval, exitflag, output, jacobian] = fsolve(@fv19_ss, ym1(end, :), ...
        optionsSS, n([1:8, end]), prod(1:end-1), degrBif, K([1:8, end]), 1, V0, convFac);
    xSS(i, :) = x;
end


end
