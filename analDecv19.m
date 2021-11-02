%This script takes solutions to replication v19 and tests the decision and
%the role of replication in the decision. 

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
%data = sortrows(data, ncol);

%Keep solutions only up to cutoff
cutOff = 1/3;
numKeep = round(cutOff*size(data, 1));
data = data(1:numKeep, :);

%%
%GET TRAJECTORIES

%Create data structure for each solution - parameters, protein & RNA trajectories (MOI
%1:5), norm. RNA trajectories (MOI 1:5), cost functions (MOI 1:5)

sol_Ind = 1:size(data, 1);
decCheck = zeros(length(sol_Ind), 1); %1 = success
lytCheck1 = zeros(length(sol_Ind), 1); %1 = fail
lytCheck2 = zeros(length(sol_Ind), 1); %1 = fail
PRECheck = zeros(length(sol_Ind), 1); %1 = true

sol(sol_Ind) = struct('parameters', [], ...
    'prod', [], 'degr', [], 'n', [], 'K', [], 'tau', [], 'dt', [], ...
    'm1', [], 'm2', [], 'm3', [], 'm4', [], 'm5', [], 'm10', [], 'm20', [], ...
    'fm1', [], 'gm1', [], 'fm5', [], 'gm5', [], 'fm10', [], 'gm10', [], ...
    'm1Num', [], 'm2Num', [], 'm3Num', [], 'm4Num', [], 'm5Num', [], ...
    'm10Num', [], 'm20Num', [], ...
    'wm1', [], 'wm2', [], 'wm3', [], 'wm4', [], 'wm5', [], ...
    'm1_wt', [], 'm2_wt', [], 'm3_wt', [], 'm4_wt', [], 'm5_wt', [], ...
    'm10_wt', [], 'fm1_wt', [], 'gm1_wt', [], 'fm5_wt', [], 'gm5_wt', [], ...
    'fm10_wt', [], 'gm10_wt', [], ...
    'wm1_wt', [], 'wm2_wt', [], 'wm3_wt', [], 'wm4_wt', [], 'wm5_wt', [], ...
    'm1Num_wt', [], 'm2Num_wt', [], 'm3Num_wt', [], 'm4Num_wt', [], 'm5Num_wt', [], ...
    'm10Num_wt', [], 'm20Num_wt', [], ...  
    'm1_Q', [], 'm2_Q', [], 'm3_Q', [], 'm4_Q', [], 'm5_Q', [], 'm10_Q', [], ...
    'm1_Q_wt', [], 'm2_Q_wt', [], 'm3_Q_wt', [], 'm4_Q_wt', [], 'm5_Q_wt', [], ...
    'm10_Q_wt', [], ...
    'm2_d1_wt', [], 'm2_d2_wt', [], 'm2_d3_wt', [], ...
    'm2_d4_wt', [], 'm2_d5_wt', [], 'm2_d6_wt', [], 'm2_d7_wt', [], 'm1_a1_wt', [], ...
    'm2Num_d1_wt', [], 'm2Num_d2_wt', [], 'm2Num_d3_wt', [], ...
    'm2Num_d4_wt', [], 'm2Num_d5_wt', [], 'm2Num_d6_wt', [], 'm2Num_d7_wt', [], 'm1Num_a1_wt', [], ...
    'fm2_d1_wt', [], 'fm2_d2_wt', [], 'fm2_d3_wt', [], ...
    'fm2_d4_wt', [], 'fm2_d5_wt', [], 'fm2_d6_wt', [], 'fm2_d7_wt', [], 'fm1_a1_wt', [], ...
    'wm2_d1_wt', [], 'wm2_d2_wt', [], 'wm2_d3_wt', [], ...
    'wm2_d4_wt', [], 'wm2_d5_wt', [], 'wm2_d6_wt', [], 'wm2_d7_wt', [], ...
    'tauOn_m2_d1_wt', [], 'tauOn_m2_d2_wt', [], 'tauOn_m2_d3_wt', [], ...
    'tauOn_m2_d4_wt', [], 'tauOn_m2_d5_wt', [], 'tauOn_m2_d6_wt', [], 'tauOn_m2_d7_wt', [], ...
    'tauOn_m1_a1_wt', [], ...
    'i_tauOn_m2_d1_wt', [], 'i_tauOn_m2_d2_wt', [], 'i_tauOn_m2_d3_wt', [], ...
    'i_tauOn_m2_d4_wt', [], 'i_tauOn_m2_d5_wt', [], 'i_tauOn_m2_d6_wt', [], 'i_tauOn_m2_d7_wt', [], ...
    'i_tauOn_m1_a1_wt', [], ...
    'm1_croPRM_wt', [], 'm2_croPRM_wt', [], 'm3_croPRM_wt', [], ...
    'm4_croPRM_wt', [], 'm5_croPRM_wt', [], 'm10_croPRM_wt', [], ...
    'm1Num_croPRM_wt', [], 'm2Num_croPRM_wt', [], 'm3Num_croPRM_wt', [], ...
    'm4Num_croPRM_wt', [], 'm5Num_croPRM_wt', [], 'm10Num_croPRM_wt', [], ...
    'm1_croPR_wt', [], 'm2_croPR_wt', [], 'm3_croPR_wt', [], ...
    'm4_croPR_wt', [], 'm5_croPR_wt', [], 'm10_croPR_wt', [], ...
    'm1Num_croPR_wt', [], 'm2Num_croPR_wt', [], 'm3Num_croPR_wt', [], ...
    'm4Num_croPR_wt', [], 'm5Num_croPR_wt', [], 'm10Num_croPR_wt', [], ...
    'm1_cIIDeg_wt', [], 'm2_cIIDeg_wt', [], 'm3_cIIDeg_wt', [], ...
    'm4_cIIDeg_wt', [], 'm5_cIIDeg_wt', [], 'm10_cIIDeg_wt', [], ...
    'm1Num_cIIDeg_wt', [], 'm2Num_cIIDeg_wt', [], 'm3Num_cIIDeg_wt', [], ...
    'm4Num_cIIDeg_wt', [], 'm5Num_cIIDeg_wt', [], 'm10Num_cIIDeg_wt', [], ...
    'm1_cILambda_wt', [], 'm2_cILambda_wt', [], 'm3_cILambda_wt', [], ...
    'm4_cILambda_wt', [], 'm5_cILambda_wt', [], 'm10_cILambda_wt', [], ...
    'm1Num_cILambda_wt', [], 'm2Num_cILambda_wt', [], 'm3Num_cILambda_wt', [], ...
    'm4Num_cILambda_wt', [], 'm5Num_cILambda_wt', [], 'm10Num_cILambda_wt', [], ...
    'm1_pert1', [], 'm2_pert1', [], 'm3_pert1', [], 'm4_pert1', [], 'm5_pert1', [], ...
    'm10_pert1', [], ...
    'm1Num_pert1', [], 'm2Num_pert1', [], 'm3Num_pert1', [], 'm4Num_pert1', [], ...
    'm5Num_pert1', [], 'm10Num_pert1', [], ...
    'm1_pert1_wt', [], 'm2_pert1_wt', [], 'm3_pert1_wt', [], 'm4_pert1_wt', [], 'm5_pert1_wt', [], ...
    'm10_pert1_wt', [], ...
    'm1Num_pert1_wt', [], 'm2Num_pert1_wt', [], 'm3Num_pert1_wt', [], ...
    'm4Num_pert1_wt', [], 'm5Num_pert1_wt', [], 'm10Num_pert1_wt', [], ...
    'm1_pert2', [], 'm2_pert2', [], 'm3_pert2', [], 'm4_pert2', [], 'm5_pert2', [], ...
    'm10_pert2', [], ...
    'm1Num_pert2', [], 'm2Num_pert2', [], 'm3Num_pert2', [], 'm4Num_pert2', [], ...
    'm5Num_pert2', [], 'm10Num_pert2', [], ...
    'm1_pert2_wt', [], 'm2_pert2_wt', [], 'm3_pert2_wt', [], 'm4_pert2_wt', [], 'm5_pert2_wt', [], ...
    'm10_pert2_wt', [], ...
    'm1Num_pert2_wt', [], 'm2Num_pert2_wt', [], 'm3Num_pert2_wt', [], ...
    'm4Num_pert2_wt', [], 'm5Num_pert2_wt', [], 'm10Num_pert2_wt', [], ...
    'm1_pert3_wt', [], 'm2_pert3_wt', [], 'm3_pert3_wt', [], 'm4_pert3_wt', [], 'm5_pert3_wt', [], ...
    'm10_pert3_wt', [], ...
    'm1Num_pert3_wt', [], 'm2Num_pert3_wt', [], 'm3Num_pert3_wt', [], ...
    'm4Num_pert3_wt', [], 'm5Num_pert3_wt', [], 'm10Num_pert3_wt', [], ...
    'm1_tauRep0_wt', [], 'm2_tauRep0_wt', [], 'm3_tauRep0_wt', [], ...
    'm4_tauRep0_wt', [], 'm5_tauRep0_wt', [], ...
    'm1Num_tauRep0_wt', [], 'm2Num_tauRep0_wt', [], 'm3Num_tauRep0_wt', [], ...
    'm4Num_tauRep0_wt', [], 'm5Num_tauRep0_wt', [], ...
    'm1_tauRep2_wt', [], 'm2_tauRep2_wt', [], 'm3_tauRep2_wt', [], ...
    'm4_tauRep2_wt', [], 'm5_tauRep2_wt', [], ...
    'm1Num_tauRep2_wt', [], 'm2Num_tauRep2_wt', [], 'm3Num_tauRep2_wt', [], ...
    'm4Num_tauRep2_wt', [], 'm5Num_tauRep2_wt', [], ...
    'x_endTauRep_wt', [], 'x_maxTauRep_wt', [], ...
    'wm1_pert3_wt', [], 'wm2_pert3_wt', [], 'wm3_pert3_wt', [], ...
    'wm4_pert3_wt', [], 'wm5_pert3_wt', [], ...
    'fm1_pert3_wt', [], 'gm1_pert3_wt', [], ...
    'wm1_cIIDeg_wt', [], 'wm2_cIIDeg_wt', [], 'wm3_cIIDeg_wt', [], ...
    'wm4_cIIDeg_wt', [], 'wm5_cIIDeg_wt', [], ...
    'wm1_cILambda_wt', [], 'wm2_cILambda_wt', [], 'wm3_cILambda_wt', [], ...
    'wm4_cILambda_wt', [], 'wm5_cILambda_wt', [], ...
    'r_lam', [], 'x_end_lo', [], 'x_end_hi', [], 'x_max_lo', [], 'x_max_hi', [], ...
    'x_int_lo', [], 'x_int_hi', [], ...
    'MOI', [], 'x_int', [], 'x_int_wt', [], 'x_end', [], 'x_end_wt', [], ...
    'x_max', [], 'x_max_wt', [], 'x_max_cI', [], 'x_max_OP_cro', [], ...
    'x_max_totEar2', [], 'x_max_totEar2_wt', [], ...
    'x_max_totLate', [], 'x_max_totLate_wt', [], ...
    'x_max_pert1', [], 'x_max_pert1_wt', [], ...
    'x_avg', [], 'x_avg_wt', [], 'x_avg_totEar2', [], 'x_avg_totEar2_wt', [], ...
    'x_avg_totLate', [], 'x_avg_totLate_wt', [], ...
    'fm_tot', [], 'fm_tot_wt', [], 'fm_totEar1_wt', [], 'fm_totEar2_wt', [], ...
    'fm_totEar2_per_wt', [], 'fm_totEar2', [], 'fm_totLate', [], ...
    'fm_totLate_wt', [], 'fm_tot_per', [], 'fm_tot_per_wt', [], ...
    'wm_tot_wt', [], 'wm_totEar1_wt', [], 'wm_totEar2_wt', [], ...
    'wm_totEar2', [], ...
    'tauOn', [], 'tauOn_wt', [], 'i_tauOn', [], 'i_tauOn_wt', [], ...
    'tauOn_Ear2_wt', [], ...
    'tauOn_Ear2', [], 'tauOn_cIIDeg_Ear2', [], 'tauOn_lambda_wt', [], ...
    'tauOn_cro', [], 'tauOn_cro_wt', [], 'tauOn_cII', [], ...
    'tauOn_cII_wt', [], 'tauOn_PRM', [], 'tauOn_PRM_wt', [], ...
    'wm_fin_wt', [], 'wm_finEar1_wt', [], 'wm_finEar2_wt', [], ...
    'x_avgNum_wt', [], ...
    'KCro_o', [], 'KCI_o', [], 'KQ', [], ...
    'KCroR_o', [], 'KCIR_o', [], 'KQR', [], ...
    'KCro', [], 'KCI', [], ...
    'KCroR', [], 'KCIR', [], ...
    'rlam_PD', [], 'MOI_PD', [], 'PDMat', [], ...
    'tauDec', [], 'tauDec_wt', [], 'tauDec_PD', [], 'earDec_PD', [], ...
    'tauOn_PD', [], 'i_tauOn_PD', [], 'tauOn_PD_kCII', [], 'i_tauOn_PD_kCII', [], ...
    'tD_PD', [], 'dMOI_PD', [], 'PDMat_PD', [], 'tauDec_PD_tD', [], ...
    'kCII_PD', [], 'PDMat_kCII', [], 'tauDec_PD_kCII', [], ...
    'f_CI_1', [], 'f_CI_2', [], 'gof_CI_1', [], 'gof_CI_2', [], ...
    'f_CI_1_wt', [], 'f_CI_2_wt', [], 'gof_CI_1_wt', [], 'gof_CI_2_wt', [], ...
    'f_Cro_1', [], 'f_Cro_2', [], 'gof_Cro_1', [], 'gof_Cro_2', [], ...
    'f_Cro_1_wt', [], 'f_Cro_2_wt', [], 'gof_Cro_1_wt', [], 'gof_Cro_2_wt', [], ...
    'f_CII_1', [], 'f_CII_2', [], 'gof_CII_1', [], 'gof_CII_2', [], ...
    'f_CII_1_wt', [], 'f_CII_2_wt', [], 'gof_CII_1_wt', [], 'gof_CII_2_wt', [], ...
    'f_Q_1', [], 'f_Q_2', [], 'gof_Q_1', [], 'gof_Q_2', [], ...
    'f_Q_1_wt', [], 'f_Q_2_wt', [], 'gof_Q_1_wt', [], 'gof_Q_2_wt', [], ...
    'tauOn_m1', [], 'tauOn_m2', [], 'tauOn_m3', [], 'tauOn_m4', [], 'tauOn_m5', [], ...
    'i_tauOn_m1', [], 'i_tauOn_m2', [], 'i_tauOn_m3', [], 'i_tauOn_m4', [], 'i_tauOn_m5', [], ...
    'tauOn_m1_wt', [], 'tauOn_m2_wt', [], 'tauOn_m3_wt', [], 'tauOn_m4_wt', [], 'tauOn_m5_wt', [], ...
    'i_tauOn_m1_wt', [], 'i_tauOn_m2_wt', [], 'i_tauOn_m3_wt', [], ...
    'i_tauOn_m4_wt', [], 'i_tauOn_m5_wt', [], ...
    'm1_cro', [], 'm1Num_cro', [], 'm1_cI', [], 'm1Num_cI', [], ...
    'm1_OP_cro', [], 'm1Num_OP_cro', []);

%Constants
convFac = (1e9*1e15)/(6.022e23); %# of um^3 -> nM
V0 = 1; %um^3
tauFail = 1e3;

%Grid spacings
drlam = 0.1; %0.001
dkCII = 0.2; %0.005
dMOI = 0.1; %0.01
dMOIThr = 0.01; %0.05
dtD = 0.2; %0.001

ss = 0;

%ICs
MOI = 1:5;
y0_m1 = [0 0 0 0 0 0 1*convFac/V0];
y0_m2 = [0 0 0 0 0 0 2*convFac/V0];
y0_m3 = [0 0 0 0 0 0 3*convFac/V0];
y0_m4 = [0 0 0 0 0 0 4*convFac/V0];
y0_m5 = [0 0 0 0 0 0 5*convFac/V0];
y0_m10 = [0 0 0 0 0 0 10*convFac/V0];

%Solve ODEs
options = odeset('Nonnegative', [], 'RelTol', 1e-6, ...
    'AbsTol', 1e-6);

%tspan
tspan = 0:0.1:65;

for i = sol_Ind
    tic
    %ASSIGN PARAMS---------------------------------------------------------
    %cI PRM basal prod, ratio of cI PRM active / basal, ratio cI PRE/PRM basal prod
    %CI translation rate, cro prod rate, Cro translation rate, cII prod rate
    %CII translation rate, replication rate
    prod = data(i, 1:9);
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
    
    V = V0.*exp(tspan.*kdil)';
    
    %SOLVE ODES ===========================================================
    
    %P- & P+---------------------------------------------------------------
    [tm1_OP, ym1_OP, tm2_OP, ym2_OP, tm3_OP, ym3_OP, tm4_OP, ym4_OP, ...
    tm5_OP, ym5_OP, fm1_OP, fm5_OP, gm1_OP, gm5_OP, wm1_OP, wm2_OP, ...
    wm3_OP, wm4_OP, wm5_OP, tauOn_m1, i_tauOn_m1, tauOn_m2, i_tauOn_m2, ...
    tauOn_m3, i_tauOn_m3, tauOn_m4, i_tauOn_m4, tauOn_m5, i_tauOn_m5, ...
    tm1_wt, ym1_wt, tm2_wt, ym2_wt, tm3_wt, ym3_wt, tm4_wt, ym4_wt, ...
    tm5_wt, ym5_wt, fm1_wt, fm5_wt, gm1_wt, gm5_wt, wm1_wt, wm2_wt, ...
    wm3_wt, wm4_wt, wm5_wt, tauOn_m1_wt, i_tauOn_m1_wt, tauOn_m2_wt, ...
    i_tauOn_m2_wt, tauOn_m3_wt, i_tauOn_m3_wt, tauOn_m4_wt, i_tauOn_m4_wt, ...
    tauOn_m5_wt, i_tauOn_m5_wt] = getOutputv19(prod, degr, n, K, tau, options, ...
    tspan, convFac, V0);
    %P-
    sol(i).m1 = [tm1_OP, ym1_OP];
    sol(i).m1Num = [tm1_OP, ym1_OP.*V./convFac];
    sol(i).m2 = [tm2_OP, ym2_OP];
    sol(i).m2Num = [tm2_OP, ym2_OP.*V./convFac];
    sol(i).m3 = [tm3_OP, ym3_OP];
    sol(i).m3Num = [tm3_OP, ym3_OP.*V./convFac];
    sol(i).m4 = [tm4_OP, ym4_OP];
    sol(i).m4Num = [tm4_OP, ym4_OP.*V./convFac];
    sol(i).m5 = [tm5_OP, ym5_OP];
    sol(i).m5Num = [tm5_OP, ym5_OP.*V./convFac];
    sol(i).fm1_OP = fm1_OP;
    sol(i).fm5_OP = fm5_OP;
    sol(i).gm1_OP = gm1_OP;
    sol(i).gm5_OP = gm5_OP;
    sol(i).wm1_OP = wm1_OP;
    sol(i).wm2_OP = wm2_OP;
    sol(i).wm3_OP = wm3_OP;
    sol(i).wm4_OP = wm4_OP;
    sol(i).wm5_OP = wm5_OP;
    sol(i).tauOn_m1 = tauOn_m1;
    sol(i).i_tauOn_m1 = i_tauOn_m1;
    sol(i).tauOn_m2 = tauOn_m2;
    sol(i).i_tauOn_m2 = i_tauOn_m2;
    sol(i).tauOn_m3 = tauOn_m3;
    sol(i).i_tauOn_m3 = i_tauOn_m3;
    sol(i).tauOn_m4 = tauOn_m4;
    sol(i).i_tauOn_m4 = i_tauOn_m4;
    sol(i).tauOn_m5 = tauOn_m5;
    sol(i).i_tauOn_m5 = i_tauOn_m5;
    %P+
    sol(i).m1_wt = [tm1_wt, ym1_wt];
    sol(i).m1Num_wt = [tm1_wt, ym1_wt.*V./convFac];
    sol(i).m2_wt = [tm2_wt, ym2_wt];
    sol(i).m2Num_wt = [tm2_wt, ym2_wt.*V./convFac];
    sol(i).m3_wt = [tm3_wt, ym3_wt];
    sol(i).m3Num_wt = [tm3_wt, ym3_wt.*V./convFac];
    sol(i).m4_wt = [tm4_wt, ym4_wt];
    sol(i).m4Num_wt = [tm4_wt, ym4_wt.*V./convFac];
    sol(i).m5_wt = [tm5_wt, ym5_wt];
    sol(i).m5Num_wt = [tm5_wt, ym5_wt.*V./convFac];
    sol(i).fm1_wt = fm1_wt;
    sol(i).fm5_wt = fm5_wt;
    sol(i).gm1_wt = gm1_wt;
    sol(i).gm5_wt = gm5_wt;
    sol(i).wm1_wt = wm1_wt;
    sol(i).wm2_wt = wm2_wt;
    sol(i).wm3_wt = wm3_wt;
    sol(i).wm4_wt = wm4_wt;
    sol(i).wm5_wt = wm5_wt;
    sol(i).tauOn_m1_wt = tauOn_m1_wt;
    sol(i).i_tauOn_m1_wt = i_tauOn_m1_wt;
    sol(i).tauOn_m2_wt = tauOn_m2_wt;
    sol(i).i_tauOn_m2_wt = i_tauOn_m2_wt;
    sol(i).tauOn_m3_wt = tauOn_m3_wt;
    sol(i).i_tauOn_m3_wt = i_tauOn_m3_wt;
    sol(i).tauOn_m4_wt = tauOn_m4_wt;
    sol(i).i_tauOn_m4_wt = i_tauOn_m4_wt;
    sol(i).tauOn_m5_wt = tauOn_m5_wt;
    sol(i).i_tauOn_m5_wt = i_tauOn_m5_wt;

    %Delay (MOI = 2)-------------------------------------------------------
    %5 min
    [tm2_d1_1_wt, ym2_d1_1_wt] = ode15s(@fv19_repv3, 0:0.1:5, y0_m1, ...
        options, n, prod, degr, K, tau, V0, convFac);
    [tm2_d1_2_wt, ym2_d1_2_wt] = ode15s(@fv19_repv3, 5:0.1:tspan(end), ...
        [ym2_d1_1_wt(end, 1:end-1), ym2_d1_1_wt(end, end)+convFac/(V0*exp(kdil*5))], ...
        options, n, prod, degr, K, tau, V0, convFac);
    sol(i).m2_d1_wt = [[tm2_d1_1_wt(1:end-1); tm2_d1_2_wt(1:end);], ...
        [ym2_d1_1_wt(1:end-1, :); ym2_d1_2_wt;]];
    sol(i).m2Num_d1_wt = [[tm2_d1_1_wt(1:end-1); tm2_d1_2_wt(1:end);], ...
        [ym2_d1_1_wt(1:end-1, :).*(V0.*exp(kdil.*tm2_d1_1_wt(1:end-1)))./convFac; ...
        ym2_d1_2_wt.*(V0.*exp(kdil.*tm2_d1_2_wt))./convFac;]];
    
    [fm2_d1_wt, gm2_d1_wt] = getFluxv19([tm2_d1_1_wt(1:end-1); tm2_d1_2_wt(1:end);], ...
        [ym2_d1_1_wt(1:end-1, :); ym2_d1_2_wt;], ...
        n, prod, degr, K, tau, V0, convFac);
    sol(i).fm2_d1_wt = fm2_d1_wt;

    wm2_d1_wt = getWeightsv19([tm2_d1_1_wt(1:end-1); tm2_d1_2_wt(1:end);], ...
        [ym2_d1_1_wt(1:end-1, :); ym2_d1_2_wt;], ...
        n, prod, degr, K, tau, V0, convFac);
    sol(i).wm2_d1_wt = wm2_d1_wt;
    
    [r, c] = find(wm2_d1_wt(:, 6) >= 0.1);
    t = [tm2_d1_1_wt(1:end-1); tm2_d1_2_wt(1:end);];
    if ~isempty(r)
        sol(i).tauOn_m2_d1_wt = [t(r(1)), t(r(end))];
        sol(i).i_tauOn_m2_d1_wt = [r(1), r(end), 1];
    end
    
    %10 min
    [tm2_d2_1_wt, ym2_d2_1_wt] = ode15s(@fv19_repv3, 0:0.1:10, y0_m1, ...
        options, n, prod, degr, K, tau, V0, convFac);
    [tm2_d2_2_wt, ym2_d2_2_wt] = ode15s(@fv19_repv3, 10:0.1:tspan(end), ...
        [ym2_d2_1_wt(end, 1:end-1), ym2_d2_1_wt(end, end)+convFac/(V0*exp(kdil*10))], ...
        options, n, prod, degr, K, tau, V0, convFac);
    sol(i).m2_d2_wt = [[tm2_d2_1_wt(1:end-1); tm2_d2_2_wt(1:end);], ...
        [ym2_d2_1_wt(1:end-1, :); ym2_d2_2_wt;]];
    sol(i).m2Num_d2_wt = [[tm2_d2_1_wt(1:end-1); tm2_d2_2_wt(1:end);], ...
        [ym2_d2_1_wt(1:end-1, :).*(V0.*exp(kdil.*tm2_d2_1_wt(1:end-1)))./convFac; ...
        ym2_d2_2_wt.*(V0.*exp(kdil.*tm2_d2_2_wt))./convFac;]];
    
    [fm2_d2_wt, gm2_d2_wt] = getFluxv19([tm2_d2_1_wt(1:end-1); tm2_d2_2_wt(1:end);], ...
        [ym2_d2_1_wt(1:end-1, :); ym2_d2_2_wt;], ...
        n, prod, degr, K, tau, V0, convFac);
    sol(i).fm2_d2_wt = fm2_d2_wt;

    wm2_d2_wt = getWeightsv19([tm2_d2_1_wt(1:end-1); tm2_d2_2_wt(1:end);], ...
        [ym2_d2_1_wt(1:end-1, :); ym2_d2_2_wt;], ...
        n, prod, degr, K, tau, V0, convFac);
    sol(i).wm2_d2_wt = wm2_d2_wt;
    
    [r, c] = find(wm2_d2_wt(:, 6) >= 0.1);
    t = [tm2_d2_1_wt(1:end-1); tm2_d2_2_wt(1:end);];
    if ~isempty(r)
        sol(i).tauOn_m2_d2_wt = [t(r(1)), t(r(end))];
        sol(i).i_tauOn_m2_d2_wt = [r(1), r(end), 1];
    end
    
    %15 minutes
    [tm2_d3_1_wt, ym2_d3_1_wt] = ode15s(@fv19_repv3, 0:0.1:15, y0_m1, ...
        options, n, prod, degr, K, tau, V0, convFac);
    [tm2_d3_2_wt, ym2_d3_2_wt] = ode15s(@fv19_repv3, 15:0.1:tspan(end), ...
        [ym2_d3_1_wt(end, 1:end-1), ym2_d3_1_wt(end, end)+convFac/(V0*exp(kdil*15))], ...
        options, n, prod, degr, K, tau, V0, convFac);
    sol(i).m2_d3_wt = [[tm2_d3_1_wt(1:end-1); tm2_d3_2_wt(1:end);], ...
        [ym2_d3_1_wt(1:end-1, :); ym2_d3_2_wt;]];
    sol(i).m2Num_d3_wt = [[tm2_d3_1_wt(1:end-1); tm2_d3_2_wt(1:end);], ...
        [ym2_d3_1_wt(1:end-1, :).*(V0.*exp(kdil.*tm2_d3_1_wt(1:end-1)))./convFac; ...
        ym2_d3_2_wt.*(V0.*exp(kdil.*tm2_d3_2_wt))./convFac;]];
    
    [fm2_d3_wt, gm2_d3_wt] = getFluxv19([tm2_d3_1_wt(1:end-1); tm2_d3_2_wt(1:end);], ...
        [ym2_d3_1_wt(1:end-1, :); ym2_d3_2_wt;], ...
        n, prod, degr, K, tau, V0, convFac);
    sol(i).fm2_d3_wt = fm2_d3_wt;

    wm2_d3_wt = getWeightsv19([tm2_d3_1_wt(1:end-1); tm2_d3_2_wt(1:end);], ...
        [ym2_d3_1_wt(1:end-1, :); ym2_d3_2_wt;], ...
        n, prod, degr, K, tau, V0, convFac);
    sol(i).wm2_d3_wt = wm2_d3_wt;
    
    [r, c] = find(wm2_d3_wt(:, 6) >= 0.1);
    t = [tm2_d3_1_wt(1:end-1); tm2_d3_2_wt(1:end);];
    if ~isempty(r)
        sol(i).tauOn_m2_d3_wt = [t(r(1)), t(r(end))];
        sol(i).i_tauOn_m2_d3_wt = [r(1), r(end), 1];
    end
    
    %Full window
    [tm2_d4_1_wt, ym2_d4_1_wt] = ode15s(@fv19_repv3, 0:0.1:1*sol(i).tauOn_m1_wt(2), y0_m1, ...
        options, n, prod, degr, K, tau, V0, convFac);
    [tm2_d4_2_wt, ym2_d4_2_wt] = ode15s(@fv19_repv3, 1*sol(i).tauOn_m1_wt(2):0.1:tspan(end), ...
        [ym2_d4_1_wt(end, 1:end-1), ym2_d4_1_wt(end, end)+convFac/(V0*exp(kdil*1*sol(i).tauOn_m1_wt(2)))], ...
        options, n, prod, degr, K, tau, V0, convFac);
    sol(i).m2_d4_wt = [[tm2_d4_1_wt(1:end-1); tm2_d4_2_wt(1:end);], ...
        [ym2_d4_1_wt(1:end-1, :); ym2_d4_2_wt;]];
    sol(i).m2Num_d4_wt = [[tm2_d4_1_wt(1:end-1); tm2_d4_2_wt(1:end);], ...
        [ym2_d4_1_wt(1:end-1, :).*(V0.*exp(kdil.*tm2_d4_1_wt(1:end-1)))./convFac; ...
        ym2_d4_2_wt.*(V0.*exp(kdil.*tm2_d4_2_wt))./convFac;]];
    
    [fm2_d4_wt, gm2_d4_wt] = getFluxv19([tm2_d4_1_wt(1:end-1); tm2_d4_2_wt(1:end);], ...
        [ym2_d4_1_wt(1:end-1, :); ym2_d4_2_wt;], ...
        n, prod, degr, K, tau, V0, convFac);
    sol(i).fm2_d4_wt = fm2_d4_wt;

    wm2_d4_wt = getWeightsv19([tm2_d4_1_wt(1:end-1); tm2_d4_2_wt(1:end);], ...
        [ym2_d4_1_wt(1:end-1, :); ym2_d4_2_wt;], ...
        n, prod, degr, K, tau, V0, convFac);
    sol(i).wm2_d4_wt = wm2_d4_wt;
    
    [r, c] = find(wm2_d4_wt(:, 6) >= 0.1);
    t = [tm2_d4_1_wt(1:end-1); tm2_d4_2_wt(1:end);];
    if ~isempty(r)
        sol(i).tauOn_m2_d4_wt = [t(r(1)), t(r(end))];
        sol(i).i_tauOn_m2_d4_wt = [r(1), r(end), 1];
    end
    
    %1.5*Full window
    [tm2_d5_1_wt, ym2_d5_1_wt] = ode15s(@fv19_repv3, 0:0.1:1.5*sol(i).tauOn_m1_wt(2), y0_m1, ...
        options, n, prod, degr, K, tau, V0, convFac);
    [tm2_d5_2_wt, ym2_d5_2_wt] = ode15s(@fv19_repv3, 1.5*sol(i).tauOn_m1_wt(2):0.1:tspan(end), ...
        [ym2_d5_1_wt(end, 1:end-1), ym2_d5_1_wt(end, end)+convFac/(V0*exp(kdil*1.5*sol(i).tauOn_m1_wt(2)))], ...
        options, n, prod, degr, K, tau, V0, convFac);
    sol(i).m2_d5_wt = [[tm2_d5_1_wt(1:end-1); tm2_d5_2_wt(1:end);], ...
        [ym2_d5_1_wt(1:end-1, :); ym2_d5_2_wt;]];
    sol(i).m2Num_d5_wt = [[tm2_d5_1_wt(1:end-1); tm2_d5_2_wt(1:end);], ...
        [ym2_d5_1_wt(1:end-1, :).*(V0.*exp(kdil.*tm2_d5_1_wt(1:end-1)))./convFac; ...
        ym2_d5_2_wt.*(V0.*exp(kdil.*tm2_d5_2_wt))./convFac;]];
    
    [fm2_d5_wt, gm2_d5_wt] = getFluxv19([tm2_d5_1_wt(1:end-1); tm2_d5_2_wt(1:end);], ...
        [ym2_d5_1_wt(1:end-1, :); ym2_d5_2_wt;], ...
        n, prod, degr, K, tau, V0, convFac);
    sol(i).fm2_d5_wt = fm2_d5_wt;

    wm2_d5_wt = getWeightsv19([tm2_d5_1_wt(1:end-1); tm2_d5_2_wt(1:end);], ...
        [ym2_d5_1_wt(1:end-1, :); ym2_d5_2_wt;], ...
        n, prod, degr, K, tau, V0, convFac);
    sol(i).wm2_d5_wt = wm2_d5_wt;
    
    [r, c] = find(wm2_d5_wt(:, 6) >= 0.1);
    t = [tm2_d5_1_wt(1:end-1); tm2_d5_2_wt(1:end);];
    if ~isempty(r)
        sol(i).tauOn_m2_d5_wt = [t(r(1)), t(r(end))];
        sol(i).i_tauOn_m2_d5_wt = [r(1), r(end), 1];
    end
    
    %Quarter of Window
    [tm2_d6_1_wt, ym2_d6_1_wt] = ode15s(@fv19_repv3, 0:0.1:(sol(i).tauOn_m1_wt(2)/4), y0_m1, ...
        options, n, prod, degr, K, tau, V0, convFac);
    [tm2_d6_2_wt, ym2_d6_2_wt] = ode15s(@fv19_repv3, (sol(i).tauOn_m1_wt(2)/4):0.1:tspan(end), ...
        [ym2_d6_1_wt(end, 1:end-1), ym2_d6_1_wt(end, end)+convFac/(V0*exp(kdil*1*sol(i).tauOn_m1_wt(2)/4))], ...
        options, n, prod, degr, K, tau, V0, convFac);
    sol(i).m2_d6_wt = [[tm2_d6_1_wt(1:end-1); tm2_d6_2_wt(1:end);], ...
        [ym2_d6_1_wt(1:end-1, :); ym2_d6_2_wt;]];
    sol(i).m2Num_d6_wt = [[tm2_d6_1_wt(1:end-1); tm2_d6_2_wt(1:end);], ...
        [ym2_d6_1_wt(1:end-1, :).*(V0.*exp(kdil.*tm2_d6_1_wt(1:end-1)))./convFac; ...
        ym2_d6_2_wt.*(V0.*exp(kdil.*tm2_d6_2_wt))./convFac;]];
    
    [fm2_d6_wt, gm2_d6_wt] = getFluxv19([tm2_d6_1_wt(1:end-1); tm2_d6_2_wt(1:end);], ...
        [ym2_d6_1_wt(1:end-1, :); ym2_d6_2_wt;], ...
        n, prod, degr, K, tau, V0, convFac);
    sol(i).fm2_d6_wt = fm2_d6_wt;

    wm2_d6_wt = getWeightsv19([tm2_d6_1_wt(1:end-1); tm2_d6_2_wt(1:end);], ...
        [ym2_d6_1_wt(1:end-1, :); ym2_d6_2_wt;], ...
        n, prod, degr, K, tau, V0, convFac);
    sol(i).wm2_d6_wt = wm2_d6_wt;
    
    [r, c] = find(wm2_d6_wt(:, 6) >= 0.1);
    t = [tm2_d6_1_wt(1:end-1); tm2_d6_2_wt(1:end);];
    if ~isempty(r)
        sol(i).tauOn_m2_d6_wt = [t(r(1)), t(r(end))];
        sol(i).i_tauOn_m2_d6_wt = [r(1), r(end), 1];
    end
    
    %Third of Window
    [tm2_d7_1_wt, ym2_d7_1_wt] = ode15s(@fv19_repv3, 0:0.1:(sol(i).tauOn_m1_wt(2)/3), y0_m1, ...
        options, n, prod, degr, K, tau, V0, convFac);
    [tm2_d7_2_wt, ym2_d7_2_wt] = ode15s(@fv19_repv3, (sol(i).tauOn_m1_wt(2)/3):0.1:tspan(end), ...
        [ym2_d7_1_wt(end, 1:end-1), ym2_d7_1_wt(end, end)+convFac/(V0*exp(kdil*1*sol(i).tauOn_m1_wt(2)/3))], ...
        options, n, prod, degr, K, tau, V0, convFac);
    sol(i).m2_d7_wt = [[tm2_d7_1_wt(1:end-1); tm2_d7_2_wt(1:end);], ...
        [ym2_d7_1_wt(1:end-1, :); ym2_d7_2_wt;]];
    sol(i).m2Num_d7_wt = [[tm2_d7_1_wt(1:end-1); tm2_d7_2_wt(1:end);], ...
        [ym2_d7_1_wt(1:end-1, :).*(V0.*exp(kdil.*tm2_d7_1_wt(1:end-1)))./convFac; ...
        ym2_d7_2_wt.*(V0.*exp(kdil.*tm2_d7_2_wt))./convFac;]];
    
    [fm2_d7_wt, gm2_d7_wt] = getFluxv19([tm2_d7_1_wt(1:end-1); tm2_d7_2_wt(1:end);], ...
        [ym2_d7_1_wt(1:end-1, :); ym2_d7_2_wt;], ...
        n, prod, degr, K, tau, V0, convFac);
    sol(i).fm2_d7_wt = fm2_d7_wt;

    wm2_d7_wt = getWeightsv19([tm2_d7_1_wt(1:end-1); tm2_d7_2_wt(1:end);], ...
        [ym2_d7_1_wt(1:end-1, :); ym2_d7_2_wt;], ...
        n, prod, degr, K, tau, V0, convFac);
    sol(i).wm2_d7_wt = wm2_d7_wt;
    
    [r, c] = find(wm2_d7_wt(:, 6) >= 0.1);
    t = [tm2_d7_1_wt(1:end-1); tm2_d7_2_wt(1:end);];
    if ~isempty(r)
        sol(i).tauOn_m2_d7_wt = [t(r(1)), t(r(end))];
        sol(i).i_tauOn_m2_d7_wt = [r(1), r(end), 1];
    end
    
    
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
    
    
    %CHECK DECISION========================================================
    Cro_check = [max(ym1_wt(:, 5)), max(ym2_wt(:, 5)), max(ym3_wt(:, 5)), ...
        max(ym4_wt(:, 5)), max(ym5_wt(:, 5))];
    CI_check = [max(ym1_wt(:, 4)), max(ym2_wt(:, 4)), max(ym3_wt(:, 4)), ...
        max(ym4_wt(:, 4)), max(ym5_wt(:, 4))];
    if any(Cro_check(2:end) < Cro_check(1))
        sol(i).dec = 1;
        decCheck(i) = 1;
    else
        sol(i).dec = 0;
    end
    
    %Check PRE=============================================================
    fm5_tot = trapz(tm5_wt, fm5_wt);
    fracPRE = fm5_tot(3)/(fm5_tot(1) + fm5_tot(2) + fm5_tot(3));
    
    if fracPRE >= 0.5
        PRECheck(i) = 1;
    end
    
    %Test MOI==============================================================
    MOI = 1:dMOIThr:5; %10
    sol(i).MOI = MOI;
    prod_pert1 = prod;
    prod_pert2 = prod;
    prod_pert1(2) = 0;
    prod_pert2(3) = 0;
    K_pert1 = K;
    K_pert1(3) = 1e10;
    K_pert2 = K;
    K_pert2([5, 7, 9]) = 1e10;
    K_pert3 = K;
    K_pert3(end) = 1e10;
    for j = 1:length(MOI)
        %Simulate----------------------------------------------------------
        [t_wt, y_wt] = ode15s(@fv19_repv3, tspan, [0 0 0 0 0 0 MOI(j)*convFac/V0], ...
            options, n, prod, degr, K, tau, V0, convFac);
        [t, y] = ode15s(@fv19, tspan, [0 0 0 0 0 0], options, n([1:8, end]), ...
            prod(1:end-1), degr, K([1:8, end]), MOI(j), V0, convFac);
        
        wm = getWeightsv19(t, [y, MOI(j).*convFac./V], ...
            n, prod, degr, K, tau, V0, convFac);
        wm_wt = getWeightsv19(t_wt, y_wt, ...
            n, prod, degr, K, tau, V0, convFac);

        %PRE on time-------------------------------------------------------
        [r, c] = find(wm(:, 6) >= 0.1);
        if ~isempty(r)
            tauOn = [t(r(1)), t(r(end))];
            i_tauOn = [r(1), r(end), MOI(j)];
        else
            tauOn = [0, 0];
            i_tauOn = [1, 1, MOI(j)];
        end
        
        [r, c] = find(wm_wt(:, 6) >= 0.1);
        if ~isempty(r)
            tauOn_wt = [t_wt(r(1)), t_wt(r(end))];
            i_tauOn_wt = [r(1), r(end), MOI(j)];
        else
            tauOn_wt = [0, 0];
            i_tauOn_wt = [1, 1, MOI(j)];
        end
        
        sol(i).PREInt_wt = trapz(t_wt, wm_wt(:, 6)); 

        %Production rates--------------------------------------------------
        [fm, gm] = getFluxv19(t, [y, MOI(j).*convFac./V], ...
            n, prod, degr, K, tau, V0, convFac);
        [fm_wt, gm_wt] = getFluxv19(t_wt, y_wt, ...
            n, prod, degr, K, tau, V0, convFac);
        fm_tot = trapz(t, fm);
        fm_tot_wt = trapz(t_wt, fm_wt);
        fm_totEar2 = trapz(t(1:i_tauOn(end-1)), fm(1:i_tauOn(end-1), :));
        fm_tot_per = trapz(t, fm./(MOI(j).*convFac./(V0.*exp(degr(1).*tspan')) ) );
        fm_totEar2_wt = trapz(t_wt(1:i_tauOn_wt(end-1)), fm_wt(1:i_tauOn_wt(end-1), :));
        if i_tauOn(end-1) ~= length(t)
            fm_totLate = trapz(t(i_tauOn(end-1):end), fm(i_tauOn(end-1):end, :));
        else
            fm_totLate = fm(i_tauOn(end-1):end, :);
        end
        if i_tauOn_wt(end-1) ~= length(t_wt)
            fm_totLate_wt = trapz(t_wt(i_tauOn_wt(end-1):end), fm_wt(i_tauOn_wt(end-1):end, :));
        else
            fm_totLate_wt = fm_wt(i_tauOn_wt(end-1):end, :);
        end
        fm_totEar2_per_wt = trapz(t_wt(1:i_tauOn_wt(end-1)), fm_wt(1:i_tauOn_wt(end-1), :)./...
            (y_wt(1:i_tauOn_wt(end-1), end)).*convFac./(V0.*exp(degr(1).*tspan(1:i_tauOn_wt(end-1))')));
        fm_tot_per_wt = trapz(t_wt, fm_wt./(y_wt(:, end).*convFac./(V0.*exp(degr(1).*tspan'))));
        
        %Weights-----------------------------------------------------------
        wm_tot_wt = trapz(t_wt, wm_wt);
        
        wm_totEar2 = trapz(t(1:i_tauOn(end-1)), wm(1:i_tauOn(end-1), :));
        wm_totEar2_wt = trapz(t_wt(1:i_tauOn_wt(end-1)), wm_wt(1:i_tauOn_wt(end-1), :));
        
        wm_fin_wt = wm_wt(end, :);
        wm_finEar2_wt = wm_wt(i_tauOn_wt(end-1), :);
        
        %TauOn for cro, cII, PRM-------------------------------------------
        %cro
        [r, c] = find(wm(:, 7) < 0.1);
        if ~isempty(r)
            tauOn_cro = t(r(1));
        else
            tauOn_cro = 0;
        end
        
        [r, c] = find(wm_wt(:, 7) < 0.1);
        if ~isempty(r)
            tauOn_cro_wt = t_wt(r(1));
        else
            tauOn_cro_wt = 0;
        end
        
        %cII
        [r, c] = find(wm(:, 10) < 0.1);
        if ~isempty(r)
            tauOn_cII = t(r(1));
        else
            tauOn_cII = 0;
        end
        
        [r, c] = find(wm_wt(:, 10) < 0.1);
        if ~isempty(r)
            tauOn_cII_wt = t_wt(r(1));
        else
            tauOn_cII_wt = 0;
        end
        
        %PRM
        [r, c] = find(wm(:, 3) + wm(:, 4) >= 0.1);
        if ~isempty(r)
            tauOn_PRM = t(r(1));
        else
            tauOn_PRM = 0;
        end
        
        [r, c] = find(wm_wt(:, 3) + wm_wt(:, 4) >= 0.1);
        if ~isempty(r)
            tauOn_PRM_wt = t_wt(r(1));
        else
            tauOn_PRM_wt = 0;
        end
        
        [r, c] = find(wm_wt(:, 16) + wm_wt(:, 17) >= 0.1);
        if ~isempty(r)
            tauOn_lambda_wt = t_wt(r(1));
        else
            tauOn_lambda_wt = 0;
        end
               
        %Decision Times----------------------------------------------------
        [r1, ~] = find(y(:, 4) >= sol(i).KCI);
        [r2, ~] = find(y(:, 5) >= sol(i).KCro);
        if ~isempty(r1) && ~isempty(r2)
            tauDec = [t(r1(1)), t(r2(1))];
        elseif isempty(r1) && ~isempty(r2)
            tauDec = [tauFail, t(r2(1))];
        elseif ~isempty(r1) && isempty(r2)
            tauDec = [t(r1(1)), tauFail];
        else
            tauDec = [tauFail, tauFail];
        end
        
        [r1_wt, ~] = find(y_wt(:, 4) >= sol(i).KCI);
        [r2_wt, ~] = find(y_wt(:, 5) >= sol(i).KCro);
        if ~isempty(r1_wt) && ~isempty(r2_wt)
            tauDec_wt = [t_wt(r1_wt(1)), t_wt(r2_wt(1))];
        elseif isempty(r1_wt) && ~isempty(r2_wt)
            tauDec_wt = [tauFail, t_wt(r2_wt(1))];
        elseif ~isempty(r1_wt) && isempty(r2_wt)
            tauDec_wt = [t_wt(r1_wt(1)), tauFail];
        else
            tauDec_wt = [tauFail, tauFail];
        end
        
        %Assignment--------------------------------------------------------
        x_end = y(end, :);
        x_end_wt = y_wt(end, :);
        x_max = max(y);
        x_max_wt = max(y_wt);
        x_max_totEar2 = max(y(1:i_tauOn(end-1), :), [], 1);
        x_max_totEar2_wt = max(y_wt(1:i_tauOn_wt(end-1), :), [], 1);
        x_max_totLate = max(y(i_tauOn(end-1):end, :), [], 1);
        x_max_totLate_wt = max(y_wt(i_tauOn_wt(end-1):end, :), [], 1);
        sol(i).x_end = cat(1, sol(i).x_end, x_end);
        sol(i).x_end_wt = cat(1, sol(i).x_end_wt, x_end_wt);
        sol(i).x_max = cat(1, sol(i).x_max, x_max);
        sol(i).x_max_wt = cat(1, sol(i).x_max_wt, x_max_wt);
        sol(i).x_max_totEar2 = cat(1, sol(i).x_max_totEar2, x_max_totEar2);
        sol(i).x_max_totEar2_wt = cat(1, sol(i).x_max_totEar2_wt, x_max_totEar2_wt);
        sol(i).x_max_totLate = cat(1, sol(i).x_max_totLate, x_max_totLate);
        sol(i).x_max_totLate_wt = cat(1, sol(i).x_max_totLate_wt, x_max_totLate_wt);
        sol(i).x_int = cat(1, sol(i).x_int, trapz(t, fm));
        sol(i).x_int_wt = cat(1, sol(i).x_int_wt, trapz(t, fm_wt));
        sol(i).x_avg = cat(1, sol(i).x_avg, trapz(t, y)./t(end));
        sol(i).x_avg_wt = cat(1, sol(i).x_avg_wt, trapz(t_wt, y_wt)./t_wt(end));
        if ~isempty(i_tauOn)
            sol(i).x_avg_totEar2 = cat(1, sol(i).x_avg_totEar2, ...
                trapz(t(1:i_tauOn(2)), y(1:i_tauOn(2), :))./t(i_tauOn(2)));
            sol(i).x_avg_totEar2_wt = cat(1, sol(i).x_avg_totEar2_wt, ...
                trapz(t_wt(1:i_tauOn_wt(2)), y_wt(1:i_tauOn_wt(2), :))./t_wt(i_tauOn_wt(2)));
        end
        if i_tauOn(2) ~= length(t) 
            sol(i).x_avg_totLate = cat(1, sol(i).x_avg_totLate, ...
                trapz(t(i_tauOn(2):end), y(i_tauOn(2):end, :))./...
                (t(end) - t(i_tauOn(2))));
        elseif ~isempty(i_tauOn)
            sol(i).x_avg_totLate = cat(1, sol(i).x_avg_totLate, ...
                y(i_tauOn(2):end, :));
        end
        if i_tauOn_wt(2) ~= length(t_wt)
            sol(i).x_avg_totLate_wt = cat(1, sol(i).x_avg_totLate_wt, ...
                trapz(t_wt(i_tauOn_wt(2):end), y_wt(i_tauOn_wt(2):end, :))./...
                (t_wt(end) - t_wt(i_tauOn_wt(2))));
        else
            sol(i).x_avg_totLate_wt = cat(1, sol(i).x_avg_totLate_wt, ...
                y_wt(i_tauOn_wt(2):end, :));
        end
        sol(i).fm_tot = cat(1, sol(i).fm_tot, fm_tot);
        sol(i).fm_tot_wt = cat(1, sol(i).fm_tot_wt, fm_tot_wt);
        %sol(i).fm_totEar1_wt = cat(1, sol(i).fm_totEar1_wt, fm_totEar1_wt);
        sol(i).fm_totEar2_wt = cat(1, sol(i).fm_totEar2_wt, fm_totEar2_wt);
        sol(i).fm_totEar2_per_wt = cat(1, sol(i).fm_totEar2_per_wt, ...
            fm_totEar2_per_wt);
        sol(i).fm_totEar2 = cat(1, sol(i).fm_totEar2, fm_totEar2);
        sol(i).fm_totLate = cat(1, sol(i).fm_totLate, fm_totLate);
        sol(i).fm_totLate_wt = cat(1, sol(i).fm_totLate_wt, fm_totLate_wt);
        sol(i).fm_tot_per = cat(1, sol(i).fm_tot_per, fm_tot_per);
        sol(i).fm_tot_per_wt = cat(1, sol(i).fm_tot_per_wt, fm_tot_per_wt);
        
        sol(i).tauDec = cat(1, sol(i).tauDec, tauDec);
        sol(i).tauDec_wt = cat(1, sol(i).tauDec_wt, tauDec_wt);
        
        sol(i).wm_tot_wt = cat(1, sol(i).wm_tot_wt, wm_tot_wt);
        sol(i).wm_totEar2 = cat(1, sol(i).wm_totEar2, wm_totEar2);
        sol(i).wm_totEar2_wt = cat(1, sol(i).wm_totEar2_wt, wm_totEar2_wt);
        sol(i).wm_fin_wt = cat(1, sol(i).wm_fin_wt, wm_fin_wt);
        sol(i).wm_finEar2_wt = cat(1, sol(i).wm_finEar2_wt, wm_finEar2_wt);
        
        sol(i).tauOn = cat(1, sol(i).tauOn, tauOn);
        sol(i).tauOn_wt = cat(1, sol(i).tauOn_wt, tauOn_wt);
        sol(i).i_tauOn = cat(1, sol(i).i_tauOn, i_tauOn);
        sol(i).i_tauOn_wt = cat(1, sol(i).i_tauOn_wt, i_tauOn_wt);
        
        sol(i).tauOn_cro = cat(1, sol(i).tauOn_cro, tauOn_cro);
        sol(i).tauOn_cro_wt = cat(1, sol(i).tauOn_cro_wt, tauOn_cro_wt);
        sol(i).tauOn_cII = cat(1, sol(i).tauOn_cII, tauOn_cII);
        sol(i).tauOn_cII_wt = cat(1, sol(i).tauOn_cII_wt, tauOn_cII_wt);
        sol(i).tauOn_PRM = cat(1, sol(i).tauOn_PRM, tauOn_PRM);
        sol(i).tauOn_PRM_wt = cat(1, sol(i).tauOn_PRM_wt, tauOn_PRM_wt);
        sol(i).tauOn_lambda_wt = cat(1, sol(i).tauOn_lambda_wt, tauOn_lambda_wt);
    end
    
    %Create Phase diagram - repl. rate vs. MOI-----------------------------
    rlam_PD = 0:drlam:2; %0.01
    MOI_PD = 1:dMOI:7; %0.05
    prodPD = prod;
    tauOn_PD = zeros(length(rlam_PD), 2);
    i_tauOn_PD = zeros(length(rlam_PD), 3);
    for j = 1:length(rlam_PD)
        prodPD(end) = rlam_PD(j)*prod(end);
        for k = 1:length(MOI_PD)           
            [t_PD, y_PD] = ode15s(@fv19_repv3, tspan, [0 0 0 0 0 0 MOI_PD(k)*convFac/V0], ...
                options, n, prodPD, degr, K, tau, V0, convFac);
            [fPD, gPD] = getFluxv19(t_PD, y_PD(:, [1:6, end]), ...
                n, prodPD, degr, K, tau, V0, convFac);
            [r1_PD, ~] = find(y_PD(:, 4) >= sol(i).KCI);
            [r2_PD, ~] = find(y_PD(:, 5) >= sol(i).KCro);
            %MOI = 1 PRE on time
            if MOI_PD(k) == 1
                wm_PD = getWeightsv19(t_PD, y_PD, n, prod, degr, K, tau, V0, convFac);
                [r, c] = find(wm_PD(:, 6) >= 0.1);
                if ~isempty(r)
                    tauOn_PD(j, :) = [t_PD(r(1)), t_PD(r(end))];
                    i_tauOn_PD(j, :) = [r(1), r(end), MOI(k)];
                else
                    tauOn_PD(j, :) = [0, 0];
                    i_tauOn_PD(j, :) = [1, 1, MOI(k)];
                end
            end
            %Threshold crossing criteria
            tauDec_PD = [tauFail, tauFail];
            if ~isempty(r1_PD)
                tauDec_PD(1) = t_PD(r1_PD(1));
            end
            if ~isempty(r2_PD)
                tauDec_PD(2) = t_PD(r2_PD(1));
            end
            sol(i).PDMat(j, k, :) = [max(y_PD(:, 4))/sol(i).KCI, ...
                max(y_PD(:, 5))/sol(i).KCro];
            sol(i).tauDec_PD(j, k, :) = tauDec_PD;
        end
    end
    sol(i).rlam_PD = rlam_PD;
    sol(i).MOI_PD = MOI_PD;
    sol(i).tauOn_PD = tauOn_PD;
    sol(i).i_tauOn_PD = i_tauOn_PD;
    
    %Create Phase diagram - CII degr. rate vs. MOI-------------------------    
    kCII_PD = 0:dkCII:2; %0.01
    MOI_PD = 1:dMOI:7; %0.05
    tauOn_PD_kCII = zeros(length(kCII_PD), 2);
    i_tauOn_PD_kCII = zeros(length(kCII_PD), 3);
    degrPD = degr;
    for j = 1:length(kCII_PD)
        degrPD(end-1) = kCII_PD(j)*degr(end-1);
        for k = 1:length(MOI_PD)           
            [t_PD, y_PD] = ode15s(@fv19_repv3, tspan, [0 0 0 0 0 0 MOI_PD(k)*convFac/V0], ...
                options, n, prod, degrPD, K, tau, V0, convFac);
            [fPD, gPD] = getFluxv19(t_PD, y_PD(:, [1:6, end]), ...
                n, prod, degrPD, K, tau, V0, convFac);
            [r1_PD, ~] = find(y_PD(:, 4) >= sol(i).KCI);
            [r2_PD, ~] = find(y_PD(:, 5) >= sol(i).KCro);
            %MOI = 1 PRE on time
            if MOI_PD(k) == 1
                wm_PD = getWeightsv19(t_PD, y_PD, n, prod, degrPD, K, tau, V0, convFac);
                [r, c] = find(wm_PD(:, 6) >= 0.1);
                if ~isempty(r)
                    tauOn_PD_kCII(j, :) = [t_PD(r(1)), t_PD(r(end))];
                    i_tauOn_PD_kCII(j, :) = [r(1), r(end), MOI(k)];
                else
                    tauOn_PD_kCII(j, :) = [0, 0];
                    i_tauOn_PD_kCII(j, :) = [1, 1, MOI(k)];
                end
            end
            %Threshold crossing criteria
            tauDec_PD = [tauFail, tauFail];
            if ~isempty(r1_PD)
                tauDec_PD(1) = t_PD(r1_PD(1));
            end
            if ~isempty(r2_PD)
                tauDec_PD(2) = t_PD(r2_PD(1));
            end
            sol(i).PDMat_kCII(j, k, :) = [max(y_PD(:, 4))/sol(i).KCI, ...
                max(y_PD(:, 5))/sol(i).KCro];
            sol(i).tauDec_PD_kCII(j, k, :) = tauDec_PD;
        end
    end
    sol(i).kCII_PD = kCII_PD;
    sol(i).tauOn_PD_kCII = tauOn_PD_kCII;
    sol(i).i_tauOn_PD_kCII = i_tauOn_PD_kCII;
    
    %Create Phase diagram - delay vs. dMOI---------------------------------
    tD_PD = 0:dtD:1; %0.005
    dMOI_PD = 1:dMOI:6; %0.05
    PDMat = zeros(length(tD_PD), length(dMOI_PD), 2); %CI/KCI, Cro/KCro
    for j = 1:length(tD_PD)
        for k = 1:length(dMOI_PD)           
            if tD_PD(j) > 0
                [t_1_PD, y_1_PD] = ode15s(@fv19_repv3, 0:0.01:(tD_PD(j)*sol(i).tauOn_m1_wt(2)), ...
                    y0_m1, options, n, prod, degr, K, tau, V0, convFac);
                [t_2_PD, y_2_PD] = ode15s(@fv19_repv3, (tD_PD(j)*sol(i).tauOn_m1_wt(2)):0.01:tspan(end), ...
                    [y_1_PD(end, 1:end-1), y_1_PD(end, end) + ...
                    dMOI_PD(k)*convFac/(V0*exp(kdil*tD_PD(j)*sol(i).tauOn_m1_wt(2)))], ...
                    options, n, prod, degr, K, tau, ...
                    (V0*exp(kdil*tD_PD(j)*sol(i).tauOn_m1_wt(2))), convFac);
                t_PD = [t_1_PD(1:end-1); t_2_PD(1:end);];
                y_PD = [y_1_PD(1:end-1, :); y_2_PD;];
            else
                [t_PD, y_PD] = ode15s(@fv19_repv3, 0:0.01:tspan(end), ...
                    (1+dMOI_PD(k))*y0_m1, options, n, prod, degr, ...
                    K, tau, V0, convFac);
            end
            [fPD, gPD] = getFluxv19(t_PD, y_PD(:, [1:6, end]), ...
                n, prodPD, degr, K, tau, V0, convFac);
            [r1_PD, ~] = find(y_PD(:, 4) >= sol(i).KCI);
            [r2_PD, ~] = find(y_PD(:, 5) >= sol(i).KCro);
            %Threshold crossing criteria
            tauDec_PD = [tauFail, tauFail, tauFail];
            if ~isempty(r1_PD)
                tauDec_PD(1) = t_PD(r1_PD(1));
            end
            if ~isempty(r2_PD)
                tauDec_PD(2) = t_PD(r2_PD(1));
            end
            sol(i).PDMat_tD(j, k, :) = [max(y_PD(:, 4))/sol(i).KCI, ...
                max(y_PD(:, 5))/sol(i).KCro];
            sol(i).tauDec_PD_tD(j, k, :) = tauDec_PD;
        end
    end
    sol(i).tD_PD = tD_PD;
    sol(i).dMOI_PD = dMOI_PD;
    
    %PRM- =================================================================
    prod_pert1 = prod;
    prod_pert1(2) = 0;
    
    %WT -------------------------------------------------------------------
    %MOI=1
    [tm1_pert1_wt, ym1_pert1_wt] = ode15s(@fv19_repv3, tspan, y0_m1, ...
        options, n, prod_pert1, degr, K, tau, V0, convFac);
    sol(i).m1_pert1_wt = [tm1_pert1_wt, ym1_pert1_wt];
    sol(i).m1Num_pert1_wt = [tm1_pert1_wt, ym1_pert1_wt.*V./convFac];
    %MOI=2
    [tm2_pert1_wt, ym2_pert1_wt] = ode15s(@fv19_repv3, tspan, y0_m2, ...
        options, n, prod_pert1, degr, K, tau, V0, convFac);
    sol(i).m2_pert1_wt = [tm2_pert1_wt, ym2_pert1_wt];
    sol(i).m2Num_pert1_wt = [tm2_pert1_wt, ym2_pert1_wt.*V./convFac];
    %MOI=3
    [tm3_pert1_wt, ym3_pert1_wt] = ode15s(@fv19_repv3, tspan, y0_m3, ...
        options, n, prod_pert1, degr, K, tau, V0, convFac);
    sol(i).m3_pert1_wt = [tm3_pert1_wt, ym3_pert1_wt];
    sol(i).m3Num_pert1_wt = [tm3_pert1_wt, ym3_pert1_wt.*V./convFac];
    %MOI=4
    [tm4_pert1_wt, ym4_pert1_wt] = ode15s(@fv19_repv3, tspan, y0_m4, ...
        options, n, prod_pert1, degr, K, tau, V0, convFac);
    sol(i).m4_pert1_wt = [tm4_pert1_wt, ym4_pert1_wt];
    sol(i).m4Num_pert1_wt = [tm4_pert1_wt, ym4_pert1_wt.*V./convFac];
    %MOI=5
    [tm5_pert1_wt, ym5_pert1_wt] = ode15s(@fv19_repv3, tspan, y0_m5, ...
        options, n, prod_pert1, degr, K, tau, V0, convFac);
    sol(i).m5_pert1_wt = [tm5_pert1_wt, ym5_pert1_wt];
    sol(i).m5Num_pert1_wt = [tm5_pert1_wt, ym5_pert1_wt.*V./convFac];
    %MOI=10
    [tm10_pert1_wt, ym10_pert1_wt] = ode15s(@fv19_repv3, tspan, ...
        y0_m10, options, n, prod_pert1, degr, K, tau, V0, convFac);
    sol(i).m10_pert1_wt = [tm10_pert1_wt, ym10_pert1_wt];
    sol(i).m10Num_pert1_wt = [tm10_pert1_wt, ym10_pert1_wt.*V./convFac];
    
    %PRE- =================================================================
    prod_pert2 = prod;
    prod_pert2(3) = 0;
    
    %WT -------------------------------------------------------------------
    %MOI=1
    [tm1_pert2_wt, ym1_pert2_wt] = ode15s(@fv19_repv3, tspan, y0_m1, ...
        options, n, prod_pert2, degr, K, tau, V0, convFac);
    sol(i).m1_pert2_wt = [tm1_pert2_wt, ym1_pert2_wt];
    sol(i).m1Num_pert2_wt = [tm1_pert2_wt, ym1_pert2_wt.*V./convFac];
    %MOI=2
    [tm2_pert2_wt, ym2_pert2_wt] = ode15s(@fv19_repv3, tspan, y0_m2, ...
        options, n, prod_pert2, degr, K, tau, V0, convFac);
    sol(i).m2_pert2_wt = [tm2_pert2_wt, ym2_pert2_wt];
    sol(i).m2Num_pert2_wt = [tm2_pert2_wt, ym2_pert2_wt.*V./convFac];
    %MOI=3
    [tm3_pert2_wt, ym3_pert2_wt] = ode15s(@fv19_repv3, tspan, y0_m3, ...
        options, n, prod_pert2, degr, K, tau, V0, convFac);
    sol(i).m3_pert2_wt = [tm3_pert2_wt, ym3_pert2_wt];
    sol(i).m3Num_pert2_wt = [tm3_pert2_wt, ym3_pert2_wt.*V./convFac];
    %MOI=4
    [tm4_pert2_wt, ym4_pert2_wt] = ode15s(@fv19_repv3, tspan, y0_m4, ...
        options, n, prod_pert2, degr, K, tau, V0, convFac);
    sol(i).m4_pert2_wt = [tm4_pert2_wt, ym4_pert2_wt];
    sol(i).m4Num_pert2_wt = [tm4_pert2_wt, ym4_pert2_wt.*V./convFac];
    %MOI=5
    [tm5_pert2_wt, ym5_pert2_wt] = ode15s(@fv19_repv3, tspan, y0_m5, ...
        options, n, prod_pert2, degr, K, tau, V0, convFac);
    sol(i).m5_pert2_wt = [tm5_pert2_wt, ym5_pert2_wt];
    sol(i).m5Num_pert2_wt = [tm5_pert2_wt, ym5_pert2_wt.*V./convFac];
    %MOI=10
    [tm10_pert2_wt, ym10_pert2_wt] = ode15s(@fv19_repv3, tspan, ...
        y0_m10, options, n, prod_pert2, degr, K, tau, V0, convFac);
    sol(i).m10_pert2_wt = [tm10_pert2_wt, ym10_pert2_wt];
    sol(i).m10Num_pert2_wt = [tm10_pert2_wt, ym10_pert2_wt.*V./convFac];
    
    %CroPRM- ==============================================================
    K_pert1 = K;
    K_pert1(3) = Inf;
    %WT -------------------------------------------------------------------
    %MOI=1
    [tm1_croPRM_wt, ym1_croPRM_wt] = ode15s(@fv19_repv3, tspan, y0_m1, ...
        options, n, prod, degr, K_pert1, tau, V0, convFac);
    sol(i).m1_croPRM_wt = [tm1_croPRM_wt, ym1_croPRM_wt];
    sol(i).m1Num_croPRM_wt = [tm1_croPRM_wt, ym1_croPRM_wt.*V./convFac];
    %MOI=2
    [tm2_croPRM_wt, ym2_croPRM_wt] = ode15s(@fv19_repv3, tspan, y0_m2, ...
        options, n, prod, degr, K_pert1, tau, V0, convFac);
    sol(i).m2_croPRM_wt = [tm2_croPRM_wt, ym2_croPRM_wt];
    sol(i).m2Num_croPRM_wt = [tm2_croPRM_wt, ym2_croPRM_wt.*V./convFac];
    %MOI=3
    [tm3_croPRM_wt, ym3_croPRM_wt] = ode15s(@fv19_repv3, tspan, y0_m3, ...
        options, n, prod, degr, K_pert1, tau, V0, convFac);
    sol(i).m3_croPRM_wt = [tm3_croPRM_wt, ym3_croPRM_wt];
    sol(i).m3Num_croPRM_wt = [tm3_croPRM_wt, ym3_croPRM_wt.*V./convFac];
    %MOI=4
    [tm4_croPRM_wt, ym4_croPRM_wt] = ode15s(@fv19_repv3, tspan, y0_m4, ...
        options, n, prod, degr, K_pert1, tau, V0, convFac);
    sol(i).m4_croPRM_wt = [tm4_croPRM_wt, ym4_croPRM_wt];
    sol(i).m4Num_croPRM_wt = [tm4_croPRM_wt, ym4_croPRM_wt.*V./convFac];
    %MOI=5
    [tm5_croPRM_wt, ym5_croPRM_wt] = ode15s(@fv19_repv3, tspan, y0_m5, ...
        options, n, prod, degr, K_pert1, tau, V0, convFac);
    sol(i).m5_croPRM_wt = [tm5_croPRM_wt, ym5_croPRM_wt];
    sol(i).m5Num_croPRM_wt = [tm5_croPRM_wt, ym5_croPRM_wt.*V./convFac];
    %MOI=10
    [tm10_croPRM_wt, ym10_croPRM_wt] = ode15s(@fv19_repv3, tspan, ...
        y0_m10, options, n, prod, degr, K_pert1, tau, V0, convFac);
    sol(i).m10_croPRM_wt = [tm10_croPRM_wt, ym10_croPRM_wt];
    sol(i).m10Num_croPRM_wt = [tm10_croPRM_wt, ym10_croPRM_wt.*V./convFac];
    
    %Cro-PR ===============================================================
    K_pert2 = K;
    K_pert2([5, 7, 9]) = Inf;
    %WT -------------------------------------------------------------------
    %MOI=1
    [tm1_croPR_wt, ym1_croPR_wt] = ode15s(@fv19_repv3, tspan, y0_m1, ...
        options, n, prod, degr, K_pert2, tau, V0, convFac);
    sol(i).m1_croPR_wt = [tm1_croPR_wt, ym1_croPR_wt];
    sol(i).m1Num_croPR_wt = [tm1_croPR_wt, ym1_croPR_wt.*V./convFac];
    %MOI=2
    [tm2_croPR_wt, ym2_croPR_wt] = ode15s(@fv19_repv3, tspan, y0_m2, ...
        options, n, prod, degr, K_pert2, tau, V0, convFac);
    sol(i).m2_croPR_wt = [tm2_croPR_wt, ym2_croPR_wt];
    sol(i).m2Num_croPR_wt = [tm2_croPR_wt, ym2_croPR_wt.*V./convFac];
    %MOI=3
    [tm3_croPR_wt, ym3_croPR_wt] = ode15s(@fv19_repv3, tspan, y0_m3, ...
        options, n, prod, degr, K_pert2, tau, V0, convFac);
    sol(i).m3_croPR_wt = [tm3_croPR_wt, ym3_croPR_wt];
    sol(i).m3Num_croPR_wt = [tm3_croPR_wt, ym3_croPR_wt.*V./convFac];
    %MOI=4
    [tm4_croPR_wt, ym4_croPR_wt] = ode15s(@fv19_repv3, tspan, y0_m4, ...
        options, n, prod, degr, K_pert2, tau, V0, convFac);
    sol(i).m4_croPR_wt = [tm4_croPR_wt, ym4_croPR_wt];
    sol(i).m4Num_croPR_wt = [tm4_croPR_wt, ym4_croPR_wt.*V./convFac];
    %MOI=5
    [tm5_croPR_wt, ym5_croPR_wt] = ode15s(@fv19_repv3, tspan, y0_m5, ...
        options, n, prod, degr, K_pert2, tau, V0, convFac);
    sol(i).m5_croPR_wt = [tm5_croPR_wt, ym5_croPR_wt];
    sol(i).m5Num_croPR_wt = [tm5_croPR_wt, ym5_croPR_wt.*V./convFac];
    %MOI=10
    [tm10_croPR_wt, ym10_croPR_wt] = ode15s(@fv19_repv3, tspan, ...
        y0_m10, options, n, prod, degr, K_pert2, tau, V0, convFac);
    sol(i).m10_croPR_wt = [tm10_croPR_wt, ym10_croPR_wt];
    sol(i).m10Num_croPR_wt = [tm10_croPR_wt, ym10_croPR_wt.*V./convFac];
    
    %CII-Deg ===============================================================
    K_pert3 = K;
    K_pert3(end-1) = Inf;
    %WT -------------------------------------------------------------------
    %MOI=1
    [tm1_cIIDeg_wt, ym1_cIIDeg_wt] = ode15s(@fv19_repv3, tspan, y0_m1, ...
        options, n, prod, degr, K_pert3, tau, V0, convFac);
    sol(i).m1_cIIDeg_wt = [tm1_cIIDeg_wt, ym1_cIIDeg_wt];
    sol(i).m1Num_cIIDeg_wt = [tm1_cIIDeg_wt, ym1_cIIDeg_wt.*V./convFac];
    %MOI=2
    [tm2_cIIDeg_wt, ym2_cIIDeg_wt] = ode15s(@fv19_repv3, tspan, y0_m2, ...
        options, n, prod, degr, K_pert3, tau, V0, convFac);
    sol(i).m2_cIIDeg_wt = [tm2_cIIDeg_wt, ym2_cIIDeg_wt];
    sol(i).m2Num_cIIDeg_wt = [tm2_cIIDeg_wt, ym2_cIIDeg_wt.*V./convFac];
    %MOI=3
    [tm3_cIIDeg_wt, ym3_cIIDeg_wt] = ode15s(@fv19_repv3, tspan, y0_m3, ...
        options, n, prod, degr, K_pert3, tau, V0, convFac);
    sol(i).m3_cIIDeg_wt = [tm3_cIIDeg_wt, ym3_cIIDeg_wt];
    sol(i).m3Num_cIIDeg_wt = [tm3_cIIDeg_wt, ym3_cIIDeg_wt.*V./convFac];
    %MOI=4
    [tm4_cIIDeg_wt, ym4_cIIDeg_wt] = ode15s(@fv19_repv3, tspan, y0_m4, ...
        options, n, prod, degr, K_pert3, tau, V0, convFac);
    sol(i).m4_cIIDeg_wt = [tm4_cIIDeg_wt, ym4_cIIDeg_wt];
    sol(i).m4Num_cIIDeg_wt = [tm4_cIIDeg_wt, ym4_cIIDeg_wt.*V./convFac];
    %MOI=5
    [tm5_cIIDeg_wt, ym5_cIIDeg_wt] = ode15s(@fv19_repv3, tspan, y0_m5, ....
        options, n, prod, degr, K_pert3, tau, V0, convFac);
    sol(i).m5_cIIDeg_wt = [tm5_cIIDeg_wt, ym5_cIIDeg_wt];
    sol(i).m5Num_cIIDeg_wt = [tm5_cIIDeg_wt, ym5_cIIDeg_wt.*V./convFac];
    %MOI=10
    [tm10_cIIDeg_wt, ym10_cIIDeg_wt] = ode15s(@fv19_repv3, tspan, ...
        y0_m10, options, n, prod, degr, K_pert3, tau, V0, convFac);
    sol(i).m10_cIIDeg_wt = [tm10_cIIDeg_wt, ym10_cIIDeg_wt];
    sol(i).m10Num_cIIDeg_wt = [tm10_cIIDeg_wt, ym10_cIIDeg_wt.*V./convFac];
    
    %weights
    wm1_cIIDeg_wt = getWeightsv19(tm1_cIIDeg_wt, ym1_cIIDeg_wt, ...
        n, prod, degr, K_pert3, tau, V0, convFac);
    wm2_cIIDeg_wt = getWeightsv19(tm2_cIIDeg_wt, ym2_cIIDeg_wt, ...
        n, prod, degr, K_pert3, tau, V0, convFac);
    wm3_cIIDeg_wt = getWeightsv19(tm3_cIIDeg_wt, ym3_cIIDeg_wt, ...
        n, prod, degr, K_pert3, tau, V0, convFac);
    wm4_cIIDeg_wt = getWeightsv19(tm4_cIIDeg_wt, ym4_cIIDeg_wt, ...
        n, prod, degr, K_pert3, tau, V0, convFac);
    wm5_cIIDeg_wt = getWeightsv19(tm5_cIIDeg_wt, ym5_cIIDeg_wt, ...
        n, prod, degr, K_pert3, tau, V0, convFac);
    sol(i).wm1_cIIDeg_wt = wm1_cIIDeg_wt;
    sol(i).wm2_cIIDeg_wt = wm2_cIIDeg_wt;
    sol(i).wm3_cIIDeg_wt = wm3_cIIDeg_wt;
    sol(i).wm4_cIIDeg_wt = wm4_cIIDeg_wt;
    sol(i).wm5_cIIDeg_wt = wm5_cIIDeg_wt;
    
    %CI Repl===============================================================
    K_pert4 = K;
    K_pert4(10) = Inf;
    %WT -------------------------------------------------------------------
    %MOI=1
    [tm1_cILambda_wt, ym1_cILambda_wt] = ode15s(@fv19_repv3, tspan, y0_m1, ...
        options, n, prod, degr, K_pert4, tau, V0, convFac);
    sol(i).m1_cILambda_wt = [tm1_cILambda_wt, ym1_cILambda_wt];
    sol(i).m1Num_cILambda_wt = [tm1_cILambda_wt, ym1_cILambda_wt.*V./convFac];
    %MOI=2
    [tm2_cILambda_wt, ym2_cILambda_wt] = ode15s(@fv19_repv3, tspan, y0_m2, ...
        options, n, prod, degr, K_pert4, tau, V0, convFac);
    sol(i).m2_cILambda_wt = [tm2_cILambda_wt, ym2_cILambda_wt];
    sol(i).m2Num_cILambda_wt = [tm2_cILambda_wt, ym2_cILambda_wt.*V./convFac];
    %MOI=3
    [tm3_cILambda_wt, ym3_cILambda_wt] = ode15s(@fv19_repv3, tspan, y0_m3, ...
        options, n, prod, degr, K_pert4, tau, V0, convFac);
    sol(i).m3_cILambda_wt = [tm3_cILambda_wt, ym3_cILambda_wt];
    sol(i).m3Num_cILambda_wt = [tm3_cILambda_wt, ym3_cILambda_wt.*V./convFac];
    %MOI=4
    [tm4_cILambda_wt, ym4_cILambda_wt] = ode15s(@fv19_repv3, tspan, y0_m4, ...
        options, n, prod, degr, K_pert4, tau, V0, convFac);
    sol(i).m4_cILambda_wt = [tm4_cILambda_wt, ym4_cILambda_wt];
    sol(i).m4Num_cILambda_wt = [tm4_cILambda_wt, ym4_cILambda_wt.*V./convFac];
    %MOI=5
    [tm5_cILambda_wt, ym5_cILambda_wt] = ode15s(@fv19_repv3, tspan, y0_m5, ...
        options, n, prod, degr, K_pert4, tau, V0, convFac);
    sol(i).m5_cILambda_wt = [tm5_cILambda_wt, ym5_cILambda_wt];
    sol(i).m5Num_cILambda_wt = [tm5_cILambda_wt, ym5_cILambda_wt.*V./convFac];
    
    %weights
    wm1_cILambda_wt = getWeightsv19(tm1_cILambda_wt, ym1_cILambda_wt, ...
        n, prod, degr, K_pert4, tau, V0, convFac);
    wm2_cILambda_wt = getWeightsv19(tm2_cILambda_wt, ym2_cILambda_wt, ...
        n, prod, degr, K_pert4, tau, V0, convFac);
    wm3_cILambda_wt = getWeightsv19(tm3_cILambda_wt, ym3_cILambda_wt, ...
        n, prod, degr, K_pert4, tau, V0, convFac);
    wm4_cILambda_wt = getWeightsv19(tm4_cILambda_wt, ym4_cILambda_wt, ...
        n, prod, degr, K_pert4, tau, V0, convFac);
    wm5_cILambda_wt = getWeightsv19(tm5_cILambda_wt, ym5_cILambda_wt, ...
        n, prod, degr, K_pert4, tau, V0, convFac);
    sol(i).wm1_cILambda_wt = wm1_cILambda_wt;
    sol(i).wm2_cILambda_wt = wm2_cILambda_wt;
    sol(i).wm3_cILambda_wt = wm3_cILambda_wt;
    sol(i).wm4_cILambda_wt = wm4_cILambda_wt;
    sol(i).wm5_cILambda_wt = wm5_cILambda_wt;

    %No repl. delay =======================================================
    %Set tau to 0
    %WT -------------------------------------------------------------------
    %MOI=1
    [tm1_tauRep0_wt, ym1_tauRep0_wt] = ode15s(@fv19_repv3, tspan, y0_m1, ...
        options, n, prod, ...
        degr, K, 0, V0, convFac);
    sol(i).m1_tauRep0_wt = [tm1_tauRep0_wt, ym1_tauRep0_wt];
    sol(i).m1Num_tauRep0_wt = [tm1_tauRep0_wt, ym1_tauRep0_wt.*V./convFac];
    %MOI=2
    [tm2_tauRep0_wt, ym2_tauRep0_wt] = ode15s(@fv19_repv3, tspan, y0_m2, ...
        options, n, prod, ...
        degr, K, 0, V0, convFac);
    sol(i).m2_tauRep0_wt = [tm2_tauRep0_wt, ym2_tauRep0_wt];
    sol(i).m2Num_tauRep0_wt = [tm2_tauRep0_wt, ym2_tauRep0_wt.*V./convFac];
    %MOI=3
    [tm3_tauRep0_wt, ym3_tauRep0_wt] = ode15s(@fv19_repv3, tspan, y0_m3, ...
        options, n, prod, ...
        degr, K, 0, V0, convFac);
    sol(i).m3_tauRep0_wt = [tm3_tauRep0_wt, ym3_tauRep0_wt];
    sol(i).m3Num_tauRep0_wt = [tm3_tauRep0_wt, ym3_tauRep0_wt.*V./convFac];
    %MOI=4
    [tm4_tauRep0_wt, ym4_tauRep0_wt] = ode15s(@fv19_repv3, tspan, y0_m4, ...
        options, n, prod, ...
        degr, K, 0, V0, convFac);
    sol(i).m4_tauRep0_wt = [tm4_tauRep0_wt, ym4_tauRep0_wt];
    sol(i).m4Num_tauRep0_wt = [tm4_tauRep0_wt, ym4_tauRep0_wt.*V./convFac];
    %MOI=5
    [tm5_tauRep0_wt, ym5_tauRep0_wt] = ode15s(@fv19_repv3, tspan, y0_m5, ...
        options, n, prod, ...
        degr, K, 0, V0, convFac);
    sol(i).m5_tauRep0_wt = [tm5_tauRep0_wt, ym5_tauRep0_wt];
    sol(i).m5Num_tauRep0_wt = [tm5_tauRep0_wt, ym5_tauRep0_wt.*V./convFac];

    %2x repl. delay =======================================================
    %Set tau to 2x
    %WT -------------------------------------------------------------------
    %MOI=1
    [tm1_tauRep2_wt, ym1_tauRep2_wt] = ode15s(@fv19_repv3, tspan, y0_m1, ...
        options, n, prod, ...
        degr, K, 2*tau, V0, convFac);
    sol(i).m1_tauRep2_wt = [tm1_tauRep2_wt, ym1_tauRep2_wt];
    sol(i).m1Num_tauRep2_wt = [tm1_tauRep2_wt, ym1_tauRep2_wt.*V./convFac];
    %MOI=2
    [tm2_tauRep2_wt, ym2_tauRep2_wt] = ode15s(@fv19_repv3, tspan, y0_m2, ...
        options, n, prod, ...
        degr, K, 2*tau, V0, convFac);
    sol(i).m2_tauRep2_wt = [tm2_tauRep2_wt, ym2_tauRep2_wt];
    sol(i).m2Num_tauRep2_wt = [tm2_tauRep2_wt, ym2_tauRep2_wt.*V./convFac];
    %MOI=3
    [tm3_tauRep2_wt, ym3_tauRep2_wt] = ode15s(@fv19_repv3, tspan, y0_m3, ...
        options, n, prod, ...
        degr, K, 2*tau, V0, convFac);
    sol(i).m3_tauRep2_wt = [tm3_tauRep2_wt, ym3_tauRep2_wt];
    sol(i).m3Num_tauRep2_wt = [tm3_tauRep2_wt, ym3_tauRep2_wt.*V./convFac];
    %MOI=4
    [tm4_tauRep2_wt, ym4_tauRep2_wt] = ode15s(@fv19_repv3, tspan, y0_m4, ...
        options, n, prod, ...
        degr, K, 2*tau, V0, convFac);
    sol(i).m4_tauRep2_wt = [tm4_tauRep2_wt, ym4_tauRep2_wt];
    sol(i).m4Num_tauRep2_wt = [tm4_tauRep2_wt, ym4_tauRep2_wt.*V./convFac];
    %MOI=5
    [tm5_tauRep2_wt, ym5_tauRep2_wt] = ode15s(@fv19_repv3, tspan, y0_m5, ...
        options, n, prod, ...
        degr, K, 2*tau, V0, convFac);
    sol(i).m5_tauRep2_wt = [tm5_tauRep2_wt, ym5_tauRep2_wt];
    sol(i).m5Num_tauRep2_wt = [tm5_tauRep2_wt, ym5_tauRep2_wt.*V./convFac];
    
    %Cro- -----------------------------------------------------------------
    prodCro = prod;
    prodCro(5) = 0;
    %MOI=1
    [tm1_cro, ym1_cro] = ode15s(@fv19_repv3, tspan, y0_m1, ...
        options, n, prodCro, degr, K, tau, V0, convFac);
    sol(i).m1_cro = [tm1_cro, ym1_cro];
    sol(i).m1Num_cro = [tm1_cro, ym1_cro.*V./convFac];
    
    %Cro-P- ---------------------------------------------------------------
    %MOI = 1
    [tm1_OP_cro, ym1_OP_cro] = ode15s(@fv19, tspan, y0_m1(1:6), options, n([1:8, end]), ...
        prodCro(1:end-1), degr, K([1:8, end]), 1, V0, convFac);
    sol(i).m1_OP_cro = [tm1_OP_cro, ym1_OP_cro];
    sol(i).m1Num_OP_cro = [tm1_OP_cro, ym1_OP_cro.*V./convFac];
    
    %CI- ------------------------------------------------------------------
    prodCI = prod;
    prodCI(1) = 0;
    %MOI=1
    [tm1_cI, ym1_cI] = ode15s(@fv19_repv3, tspan, y0_m1, ...
        options, n, prodCI, degr, K, tau, V0, convFac);
    sol(i).m1_cI = [tm1_cI, ym1_cI];
    sol(i).m1Num_cI = [tm1_cI, ym1_cI.*V./convFac];

    %Scan replication delay (MOI = 1) =====================================
    tauRep = 0:0.1:2;
    for j = 1:length(tauRep)
        [t_tauRep_wt, y_tauRep_wt] = ode15s(@fv19_repv3, tspan, y0_m1, ...
            options, n, prod, ...
            degr, K, tau*tauRep(j), V0, convFac);
        x_endTauRep_wt = y_tauRep_wt(end, :);
        x_maxTauRep_wt = max(y_tauRep_wt);
        sol(i).x_endTauRep_wt = cat(1, sol(i).x_endTauRep_wt, x_endTauRep_wt);
        sol(i).x_maxTauRep_wt = cat(1, sol(i).x_maxTauRep_wt, x_maxTauRep_wt);
    end
    toc
    disp(string(strcat('Solution', {' '}, num2str(i), {' '}, 'finished!')));
    
end


%%
%CLOSE FILES

save(strcat('sol_exp', num2str(experiment), '.mat'), 'sol');
fclose('all');