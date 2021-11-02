%This script plots the new figures for our revisions to the PNAS manuscript

%%
%LOAD DATA

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
end

%Thu qPCR Data (cI857 30C, WT)---------------------------------------------
%format: PCR(1:3).[time, mean, error]
%1: <MOI> = 0.4, 2: <MOI> = 3, 3: <MOI> = 8
%t: 0     5    10    20    30    45    60    75    90   120
%load('121720_thu_qPCR.mat');
load('042121_thu_qPCR.mat');

%Zeng's iScience 2018 Data (cI857, 30C, MOI = 1)---------------------------
%xdata: t, ydata: RNA
%data_cII_croP, data_cII_cro, data_cII_cI, data_cII_wt, data_cII_P
%data_cI_croP, data_cI_wt
load('032020_zeng.mat');

cm = [
    168 222 233;
    101 185 222;
    75 123 190;
    88 78 160;
    7 24 50;
]./255;

cLyt = [108, 190, 69]./255;
cLys = [237, 28, 36]./255;
cShade = [229, 228, 228]./255;
cCII = [185, 83, 159]./255;

mSize = 5; %marker size
dXMax = 62; %max length, x-axis


%%
%LOAD SOLUTION

%Load revision solution
load('solRev_exp2.mat', 'sol');
sol2 = sol;
clear sol

%Load main text solution
load('sol_exp2.mat');

%CI, Cro normalizations
CINorm = sol(1).K(6); %KCro_CI
CroNorm = sol(1).K(3); %KPRM_Cro

%%
%ENSEMBLE TABLE

x_t = zeros(25, 1);
str_t = {'Lysis-to-lysogeny transition requires replication', ...
    'The critical MOI for lysogeny is lower for P+ infection', ...
    'Infection at MOI = 1 during P- results in failed infection', ...
    'PRE is required for a transition to lysogeny, while PRM is not'};

tableLength = length(sol2);

for i = 1:tableLength
    %OUTCOMES:
    %0 - failed infection
    %1 - lysogeny
    %2 - lysis
    %3 - mixed outcome
    
   %Main Table (Predictions)-----------------------------------------------
   %1. Lysis-to-lysogeny transition requires replication
   a = getOutCome(sol(i), {'m1', 'm2', 'm3', 'm4', 'm5'});
   [m_1a, i_1a] = find(a == 2);
   [m_2a, i_2a] = find(a == 1);
   if isempty(m_1a) || isempty(m_2a) || ~any(m_2a > m_1a)
       x_t(1) = x_t(1) + 1;
   end
   %2. The critical MOI for lysogeny is lower for P+ infection
   a = getOutCome(sol(i), {'m1', 'm2', 'm3', 'm4', 'm5'});
   b = getOutCome(sol(i), {'m1_wt', 'm2_wt', 'm3_wt', 'm4_wt', 'm5_wt'});
   [m_2a, i_2a] = find(a == 1);
   [m_2b, i_2b] = find(b == 1);
   if ~isempty(m_2b) && ~isempty(m_2a) && m_2a(1) > m_2b(1)
       x_t(2) = x_t(2) + 1;
   elseif ~isempty(m_2b) && isempty(m_2a)
       x_t(2) = x_t(2) + 1;
   end
   %3. Infection at MOI = 1 during P- results in failed infection
   if max(sol(i).m1(:, 5))/sol(i).KCI < 1 && max(sol(i).m1(:, 6))/sol(i).KCro
       x_t(3) = x_t(3) + 1;
   end
   %4. PRM transcription is not required for a lysis-to-lysogeny
   %transition, while PRE is
   a = getOutCome(sol(i), {'m1_pert1_wt', 'm2_pert1_wt', 'm3_pert1_wt', ...
       'm4_pert1_wt', 'm5_pert1_wt'}); %PRM-
   b = getOutCome(sol(i), {'m1_pert2_wt', 'm2_pert2_wt', 'm3_pert2_wt', ...
       'm4_pert2_wt', 'm5_pert2_wt'}); %PRE-
   [m_1a, i_1a] = find(a == 2); %lysis
   [m_1b, i_b1] = find(b == 2); %lysis
   if ~any(b == 1) && any(a == 1) %as long as no lysogeny in PRE-, but lysog in PRM-
        x_t(4) = x_t(4) + 1; 
   end
   %5. PRE activity occurs within the first 35 minutes during P+ infection
   a = [sol(i).tauOn_m1_wt(2), sol(i).tauOn_m2_wt(2), sol(i).tauOn_m3_wt(2), ...
       sol(i).tauOn_m4_wt(2), sol(i).tauOn_m5_wt(2)];
   if ~any(a > 35)
        x_t(5) = x_t(5) + 1;
   end
   %6. CI accumulation during the CII activity window scales with viral copy number
   a = [max(sol(i).m1_wt(1:sol(i).i_tauOn_m1_wt(2), 5)), ...
       max(sol(i).m2_wt(1:sol(i).i_tauOn_m2_wt(2), 5)), ...
       max(sol(i).m3_wt(1:sol(i).i_tauOn_m3_wt(2), 5)), ...
       max(sol(i).m4_wt(1:sol(i).i_tauOn_m4_wt(2), 5)), ...
       max(sol(i).m5_wt(1:sol(i).i_tauOn_m5_wt(2), 5))]/...
       max(sol(i).m1_wt(1:sol(i).i_tauOn_m1_wt(2), 5));
   b = [1, 2, 3, 4, 5];
   if sum(a(2:end) > b(2:end)) == 4
      x_t(6) = x_t(6) + 1; 
   end
   %7. CII-activated PRE transcription is not ultrasensitive
   a = diff(sol(i).PREInt_wt)/sol(i).PREInt_wt(1:end-1);
   b = diff(sol(i).MOI)/sol(i).MOI(1:end-1);
   if any(a./b >= 2)
       x_t(7) = x_t(7) + 1;
   end
   %10. Changes in viral copy number outside the CII activity window do not
   %flip the decision to lysogeny
   if ~(max(sol(i).m2_d5_wt(:, 5))/sol(i).KCI > 1 && max(sol(i).m2_d5_wt(:, 6))/sol(i).KCro < 1)
      x_t(10) = x_t(10) + 1; 
   else
      keyboard 
   end
   %11. High enough replication rate at MOI = 1 yields lysogeny
   a = sol(i).PDMat; %rlam, MOI, (CI, Cro, CII)
   [m_PD, i_PD] = find(a(:, 1, 1) >= 1);
   if any(a(m_PD, 1, 2) < 1)
      x_t(11) = x_t(11) + 1; 
   end
   %12. Extending CII lifetime yields lysogeny at MOI = 1
   a = sol(i).PDMat_kCII;
   [m_PD, i_PD] = find(a(:, 1, 1) >= 1);
   if any(a(m_PD, 1, 2) < 1)
      x_t(12) = x_t(12) + 1; 
   else
      %keyboard 
   end
   %13. Mutations inhibiting Cro repression of cII result in lysogeny at
   %MOI = 1
   if max(sol2(i).m1_wt_CroCII(:, 5))/sol2(i).KCI >= 1 && ...
           max(sol2(i).m1_wt_CroCII(:, 6))/sol2(i).KCro < 1
       x_t(13) = x_t(13) + 1;
   end
   
   %Reviewer comment table items-------------------------------------------
   %No lysis-to-lysogeny transition during P+ infection (Removal of CI repr. of PR)
   a = getOutCome(sol2(i), {'m1_wt_CIPR', 'm2_wt_CIPR', 'm3_wt_CIPR', ...
       'm4_wt_CIPR', 'm5_wt_CIPR'}); 
   [m_CIPR, i_CIPR] = find(a == 2);
   if ~isempty(i_CIPR) && m_CIPR(end) < 5 && a(m_CIPR(end) + 1) == 1
      x_t(14) = x_t(14) + 1; 
   end
   %No lysis-to-lysogeny transition during P+ infection (Removal of Cro
   %repr. of PR
   a = getOutCome(sol2(i), {'m1_wt_CroPR', 'm2_wt_CroPR', 'm3_wt_CroPR', ...
       'm4_wt_CroPR', 'm5_wt_CroPR'}); 
   [m_CroPR, i_CroPR] = find(a == 2);
   if ~isempty(i_CroPR) && m_CroPR(end) < 5 && a(m_CroPR(end) + 1) == 1
      x_t(15) = x_t(15) + 1; 
   end
   %No lysis-to-lysogeny transtion during P+ infection (removal of Cro and
   %CI repr. of PR)
   a = getOutCome(sol2(i), {'m1_wt_CroCIPR', 'm2_wt_CroCIPR', 'm3_wt_CroCIPR', ...
       'm4_wt_CroCIPR', 'm5_wt_CroCIPR'}); 
   [m_CroCIPR, i_CroCIPR] = find(a == 2);
   if ~isempty(i_CroCIPR) && m_CroCIPR(end) < 5 && a(m_CroCIPR(end) + 1) == 1
      x_t(16) = x_t(16) + 1; 
   end
   %No lysis-to-lysogeny transition during P+ infection (removal of CI
   %repression of replication
   a = getOutCome(sol2(i), {'m1_wt_CIrep', 'm2_wt_CIrep', 'm3_wt_CIrep', ...
       'm4_wt_CIrep', 'm5_wt_CIrep'}); 
   [m_CIrep, i_CIrep] = find(a == 2);
   if ~isempty(i_CIrep) && m_CIrep(end) < 5 && a(m_CIrep(end) + 1) == 1
      x_t(17) = x_t(17) + 1; 
   end
   %No lysis-to-lysogeny transition during P+ infection (removal of Cro
   %repression of replication)
   a = getOutCome(sol2(i), {'m1_wt_Crorep', 'm2_wt_Crorep', 'm3_wt_Crorep', ...
       'm4_wt_Crorep', 'm5_wt_Crorep'}); 
   [m_Crorep, i_Crorep] = find(a == 2);
   if ~isempty(i_Crorep) && m_Crorep(end) < 5 && a(m_Crorep(end) + 1) == 1
      x_t(18) = x_t(18) + 1; 
   end
   %No lysis-to-lysogeny transition during P+ infection (removal of CI and
   %Cro repression of replication)
   a = getOutCome(sol2(i), {'m1_wt_CroCIrep', 'm2_wt_CroCIrep', 'm3_wt_CroCIrep', ...
       'm4_wt_CroCIrep', 'm5_wt_CroCIrep'}); 
   [m_CroCIrep, i_CroCIrep] = find(a == 2);
   if ~isempty(i_CroCIrep) && m_CroCIrep(end) < 5 && a(m_CroCIrep(end) + 1) == 1
      x_t(19) = x_t(19) + 1; 
   end
   %No lysis-to-lysogeny transition during P+ infection (removal of Cro
   %repression of cII)
    a = getOutCome(sol2(i), {'m1_wt_CroCII', 'm2_wt_CroCII', 'm3_wt_CroCII', ...
       'm4_wt_CroCII', 'm5_wt_CroCII'}); 
   [m_CroCII, i_CroCII] = find(a == 2);
   if ~isempty(i_CroCII) && m_CroCII(end) < 5 && a(m_CroCII(end) + 1) == 1
      x_t(20) = x_t(20) + 1; 
   end
   %No lysis-to-lysogeny transition during P+ infection (removal of CI
   %autoactivation)
   a = getOutCome(sol(i), {'m1_pert1_wt', 'm2_pert1_wt', 'm3_pert1_wt', ...
       'm4_pert1_wt', 'm5_pert1_wt'}); 
   [m_pert1_wt, i_pert1_wt] = find(a == 2);
   if ~isempty(i_pert1_wt) && m_pert1_wt(end) < 5 && a(m_pert1_wt(end) + 1) == 1
      x_t(21) = x_t(21) + 1; 
   end
   %No lysis-to-lysogeny transition during P+ infection (removal of CII
   %activation of PRE)
   a = getOutCome(sol(i), {'m1_pert2_wt', 'm2_pert2_wt', 'm3_pert2_wt', ...
       'm4_pert2_wt', 'm5_pert2_wt'}); 
   [m_pert2_wt, i_pert2_wt] = find(a == 2);
   if ~isempty(i_pert2_wt) && m_pert2_wt(end) < 5 && a(m_pert2_wt(end) + 1) == 1
      x_t(22) = x_t(22) + 1; 
   end
   %Lysis-to-lysogeny transition during P+ infection if CII degradation
   %rate is doubled
   a = sol(i).PDMat_kCII;
   [~, i_PD] = find(sol(i).kCII_PD == 2);
   [~, i_PD2] = find(sol(i).MOI_PD <= 5);
   if any(a(i_PD, i_PD2, 2) < 1 & a(i_PD, i_PD2, 1) > 1)
      x_t(23) = x_t(23) + 1; 
   else
      %keyboard 
   end
   %Lysis-to-lysogeny transition during P+ infection if CII degradation
   %rate is halved
   a = sol(i).PDMat_kCII;
   [~, i_PD] = find(sol(i).kCII_PD == 0.5);
   [~, i_PD2] = find(a(i_PD, :, 2) > 1 & sol(i).MOI_PD <= 5);
   if ~isempty(i_PD2) && any(a(i_PD, i_PD(end):end, 1) > 1)
       x_t(24) = x_t(24) + 1;
      %keyboard 
   end
   %No lysis-to-lysogeny transition when nPRE = 2
   a = getOutCome(sol2(i), {'m1_wt_nPRE', 'm2_wt_nPRE', 'm3_wt_nPRE', ...
       'm4_wt_nPRE', 'm5_wt_nPRE'}); 
   [m_wt_nPRE, i_wt_nPRE] = find(a == 2);
   if ~isempty(i_wt_nPRE) && m_wt_nPRE(end) < 5 && a(m_wt_nPRE(end) + 1) == 1
      x_t(25) = x_t(25) + 1; 
   end
end
    
%%
%PRM TRANSLATION

fig3a = figure('Name', 'PRM, 20-fold lower translational rate');
dim = [.2 .5 .3 .3];
str_lett = {'A', 'B', 'C', 'D', 'E', 'F'};
numLet = 4;
for i = 1:numLet
    a_lett = annotation('textbox',dim,'String',str_lett{i},'FitBoxToText','on', ...
        'Color', 'k', 'EdgeColor', 'none', ...
        'BackgroundColor', 'none', 'Fontsize', 16);
end
%Lower PRM translational rate----------------------------------------------
subplot(2, 2, 1); %P-
    p1 = plot(sol2(1).m1_PRMtr(:, 6)./CINorm, ...
        sol2(1).m1_PRMtr(:, 7)./CroNorm, '-', 'Color', cm(1, :)); hold on;
    p2 = plot(sol2(1).m2_PRMtr(:, 6)./CINorm, ...
        sol2(1).m2_PRMtr(:, 7)./CroNorm, '-', 'Color', cm(2, :));
    p3 = plot(sol2(1).m3_PRMtr(:, 6)./CINorm, ...
        sol2(1).m3_PRMtr(:, 7)./CroNorm, '-', 'Color', cm(3, :));
    p4 = plot(sol2(1).m4_PRMtr(:, 6)./CINorm, ...
        sol2(1).m4_PRMtr(:, 7)./CroNorm, '-', 'Color', cm(4, :));
    p5 = plot(sol2(1).m5_PRMtr(:, 6)./CINorm, ...
        sol2(1).m5_PRMtr(:, 7)./CroNorm, '-', 'Color', cm(5, :));   
    %Plot thresholds
    xl = xlim;
    yl = ylim;
    cro_thresh = plot(xl, [sol2(1).KCro/CroNorm; sol2(1).KCro/CroNorm;], '--', 'Color', cLyt, ...
        'LineWidth', 1);
    cI_thresh = plot([sol2(1).KCI/CINorm; sol2(1).KCI/CINorm;], yl, '--', 'Color', cLys, ...
        'LineWidth', 1);
    dim = [.2 .5 .3 .3];
    str1 = 'Lytic threshold';
    a1 = annotation('textbox',dim,'String',str1,'FitBoxToText','on', ...
        'Color', cLyt, 'EdgeColor', 'none', ...
        'BackgroundColor', 'none');
    str2 = 'Lysogenic threshold';
    a2 = annotation('textbox',dim,'String',str2,'FitBoxToText','on', ...
        'Color', cLys, 'EdgeColor', 'none', ...
        'BackgroundColor', 'none');
    xlabel('CI concentration (normalized)');
    ylabel('Cro concentration (normalized)');
    legend([p1, p2, p3, p4, p5], 'MOI = 1 (P-)', 'MOI = 2 (P-)', 'MOI = 3 (P-)', ...
        'MOI = 4 (P-)', 'MOI = 5 (P-)', 'Color', 'none', 'EdgeColor', 'none');
    axis square
subplot(2, 2, 2); %P+
    p1 = plot(sol2(1).m1_wt_PRMtr(:, 6)./CINorm, ...
        sol2(1).m1_wt_PRMtr(:, 7)./CroNorm, '-', 'Color', cm(1, :)); hold on;
    p2 = plot(sol2(1).m2_wt_PRMtr(:, 6)./CINorm, ...
        sol2(1).m2_wt_PRMtr(:, 7)./CroNorm, '-', 'Color', cm(2, :));
    p3 = plot(sol2(1).m3_wt_PRMtr(:, 6)./CINorm, ...
        sol2(1).m3_wt_PRMtr(:, 7)./CroNorm, '-', 'Color', cm(3, :));
    p4 = plot(sol2(1).m4_wt_PRMtr(:, 6)./CINorm, ...
        sol2(1).m4_wt_PRMtr(:, 7)./CroNorm, '-', 'Color', cm(4, :));
    p5 = plot(sol2(1).m5_wt_PRMtr(:, 6)./CINorm, ...
        sol2(1).m5_wt_PRMtr(:, 7)./CroNorm, '-', 'Color', cm(5, :));   
    %Plot thresholds
    xl = xlim;
    yl = ylim;
    cro_thresh = plot(xl, [sol2(1).KCro/CroNorm; sol2(1).KCro/CroNorm;], '--', 'Color', cLyt, ...
        'LineWidth', 1);
    cI_thresh = plot([sol2(1).KCI/CINorm; sol2(1).KCI/CINorm;], yl, '--', 'Color', cLys, ...
        'LineWidth', 1);
    dim = [.2 .5 .3 .3];
    str1 = 'Lytic threshold';
    a1 = annotation('textbox',dim,'String',str1,'FitBoxToText','on', ...
        'Color', cLyt, 'EdgeColor', 'none', ...
        'BackgroundColor', 'none');
    str2 = 'Lysogenic threshold';
    a2 = annotation('textbox',dim,'String',str2,'FitBoxToText','on', ...
        'Color', cLys, 'EdgeColor', 'none', ...
        'BackgroundColor', 'none');
    xlabel('CI concentration (normalized)');
    ylabel('Cro concentration (normalized)');
    legend([p1, p2, p3, p4, p5], 'MOI = 1 (P+)', 'MOI = 2 (P+)', 'MOI = 3 (P+)', ...
        'MOI = 4 (P+)', 'MOI = 5 (P+)', 'Color', 'none', 'EdgeColor', 'none');
    axis square
%Normal PRM translation rate-----------------------------------------------
subplot(2, 2, 3); %P-   
    %Plot best solution
    p1 = plot(sol(1).m1(:, 5)./CINorm, sol(1).m1(:, 6)./CroNorm, ...
        '-', 'Color', cm(1, :), 'LineWidth', 1.5); hold on;
    p2 = plot(sol(1).m2(:, 5)./CINorm, sol(1).m2(:, 6)./CroNorm, ...
        '-', 'Color', cm(2, :), 'LineWidth', 1.5);
    p3 = plot(sol(1).m3(:, 5)./CINorm, sol(1).m3(:, 6)./CroNorm, ...
        '-', 'Color', cm(3, :), 'LineWidth', 1.5);
    p4 = plot(sol(1).m4(:, 5)./CINorm, sol(1).m4(:, 6)./CroNorm, ...
        '-', 'Color', cm(4, :), 'LineWidth', 1.5);
    p5 = plot(sol(1).m5(:, 5)./CINorm, sol(1).m5(:, 6)./CroNorm, ...
        '-', 'Color', cm(5, :), 'LineWidth', 1.5);
    %Plot thresholds
    xl = [0; 30;];
    yl = [0; 30;];
    cro_thresh = plot(xl, [sol(1).KCro/CroNorm; sol(1).KCro/CroNorm;], ...
        '--', 'Color', cLyt, 'LineWidth', 1);
    cI_thresh = plot([sol(1).KCI/CINorm; sol(1).KCI/CINorm;], yl, ...
        '--', 'Color', cLys, 'LineWidth', 1);
    dim = [.2 .5 .3 .3];
    str1 = 'Lytic threshold';
    a1 = annotation('textbox',dim,'String',str1,'FitBoxToText','on', ...
        'Color', cLyt, 'EdgeColor', 'none', ...
        'BackgroundColor', 'none');
    str2 = 'Lysogenic threshold';
    a2 = annotation('textbox',dim,'String',str2,'FitBoxToText','on', ...
        'Color', cLys, 'EdgeColor', 'none', ...
        'BackgroundColor', 'none');
    xlabel('CI concentration (normalized)');
    ylabel('Cro concentration (normalized)');
    legend([p1, p2, p3, p4, p5], 'MOI = 1 (P-)', 'MOI = 2 (P-)', 'MOI = 3 (P-)', ...
        'MOI = 4 (P-)', 'MOI = 5 (P-)', 'Color', 'none', 'EdgeColor', 'none');
    %set(legend, 'NumColumns', 2);
    axis square
subplot(2, 2, 4); %WT
    %Plot best solution
    p1 = plot(sol(1).m1_wt(:, 5)./CINorm, sol(1).m1_wt(:, 6)./CroNorm, ...
        '-', 'Color', cm(1, :), 'LineWidth', 1.5); hold on;
    p2 = plot(sol(1).m2_wt(:, 5)./CINorm, sol(1).m2_wt(:, 6)./CroNorm, ...
        '-', 'Color', cm(2, :), 'LineWidth', 1.5);
    p3 = plot(sol(1).m3_wt(:, 5)./CINorm, sol(1).m3_wt(:, 6)./CroNorm, ...
        '-', 'Color', cm(3, :), 'LineWidth', 1.5);
    p4 = plot(sol(1).m4_wt(:, 5)./CINorm, sol(1).m4_wt(:, 6)./CroNorm, ...
        '-', 'Color', cm(4, :), 'LineWidth', 1.5);
    p5 = plot(sol(1).m5_wt(:, 5)./CINorm, sol(1).m5_wt(:, 6)./CroNorm, ...
        '-', 'Color', cm(5, :), 'LineWidth', 1.5);
    %Plot thresholds
    xl = [0; 30;];
    yl = [0; 30;];
    cro_thresh = plot(xl, [sol(1).KCro/CroNorm; sol(1).KCro/CroNorm;], ...
        '--', 'Color', cLyt, 'LineWidth', 1);
    cI_thresh = plot([sol(1).KCI/CINorm; sol(1).KCI/CINorm;], yl, ...
        '--', 'Color', cLys, 'LineWidth', 1);
    dim = [.2 .5 .3 .3];
    str1 = 'Lytic threshold';
    a1 = annotation('textbox',dim,'String',str1,'FitBoxToText','on', ...
        'Color', cLyt, 'EdgeColor', 'none', ...
        'BackgroundColor', 'none');
    str2 = 'Lysogenic threshold';
    a2 = annotation('textbox',dim,'String',str2,'FitBoxToText','on', ...
        'Color', cLys, 'EdgeColor', 'none', ...
        'BackgroundColor', 'none');
    legend([p1, p2, p3, p4, p5], 'MOI = 1 (P+)', 'MOI = 2 (P+)', 'MOI = 3 (P+)', ...
        'MOI = 4 (P+)', 'MOI = 5 (P+)', 'Color', 'none', 'EdgeColor', 'none');
    set(legend, 'NumColumns', 2);
    xlabel('CI concentration (normalized)');
    ylabel('Cro concentration (normalized)');
    axis square

%%
%CRUCIAL MATHEMATICAL STRUCTURE
    
%CII REG ==================================================================
fig7c = figure('Name', 'CII Reg. phase planes');
dim = [.2 .5 .3 .3];
str_lett = {'A', 'B', 'C', 'D', 'E', 'F'};
numLet = 2;
for i = 1:numLet
    a_lett = annotation('textbox',dim,'String',str_lett{i},'FitBoxToText','on', ...
        'Color', 'k', 'EdgeColor', 'none', ...
        'BackgroundColor', 'none', 'Fontsize', 16);
end
subplot(2, 1, 1); %No Cro repression of cII
    p1 = plot(sol2(1).m1_wt_CroCII(:, 5)./CINorm, ...
        sol2(1).m1_wt_CroCII(:, 6)./CroNorm, '-', 'Color', cm(1, :)); hold on;
    p2 = plot(sol2(1).m2_wt_CroCII(:, 5)./CINorm, ...
        sol2(1).m2_wt_CroCII(:, 6)./CroNorm, '-', 'Color', cm(2, :));
    p3 = plot(sol2(1).m3_wt_CroCII(:, 5)./CINorm, ...
        sol2(1).m3_wt_CroCII(:, 6)./CroNorm, '-', 'Color', cm(3, :));
    p4 = plot(sol2(1).m4_wt_CroCII(:, 5)./CINorm, ...
        sol2(1).m4_wt_CroCII(:, 6)./CroNorm, '-', 'Color', cm(4, :));
    p5 = plot(sol2(1).m5_wt_CroCII(:, 5)./CINorm, ...
        sol2(1).m5_wt_CroCII(:, 6)./CroNorm, '-', 'Color', cm(5, :));
    %Plot thresholds
    yl = ylim;
    ylim([yl(1), 2*yl(2)]);
    xl = xlim;
    yl = ylim;
    cro_thresh = plot(xl, [sol2(1).KCro; sol2(1).KCro;]/CroNorm, '--', 'Color', cLyt, ...
        'LineWidth', 1);
    cI_thresh = plot([sol2(1).KCI; sol2(1).KCI;]/CINorm, yl, '--', 'Color', cLys, ...
        'LineWidth', 1);
    dim = [.2 .5 .3 .3];
    str1 = 'Lytic threshold';
    a1 = annotation('textbox',dim,'String',str1,'FitBoxToText','on', ...
        'Color', cLyt, 'EdgeColor', 'none', ...
        'BackgroundColor', 'none');
    str2 = 'Lysogenic threshold';
    a2 = annotation('textbox',dim,'String',str2,'FitBoxToText','on', ...
        'Color', cLys, 'EdgeColor', 'none', ...
        'BackgroundColor', 'none');
    xlabel('CI concentration (normalized)');
    ylabel('Cro concentration (normalized)');
    legend([p1, p2, p3, p4, p5], 'MOI = 1 (P+, K_{cII,Cro}=\infty)', ...
        'MOI = 2 (P+, K_{cII,Cro}=\infty)', ...
        'MOI = 3 (P+, K_{cII,Cro}=\infty)', ...
        'MOI = 4 (P+, K_{cII,Cro}=\infty)', ...
        'MOI = 5 (P+, K_{cII,Cro}=\infty)', 'Color', 'none', 'EdgeColor', 'none');
    axis square
subplot(2, 1, 2); %"", PRE
    p1 = plot(sol2(1).m1_wt_CroCII(:, 1), sol2(1).wm1_wt_CroCII(:, 6), ...
        '-', 'Color', cm(1, :)); hold on;
    p2 = plot(sol2(1).m2_wt_CroCII(:, 1), sol2(1).wm2_wt_CroCII(:, 6), ...
        '-', 'Color', cm(2, :));
    p3 = plot(sol2(1).m3_wt_CroCII(:, 1), sol2(1).wm3_wt_CroCII(:, 6), ...
        '-', 'Color', cm(3, :));
    p4 = plot(sol2(1).m4_wt_CroCII(:, 1), sol2(1).wm4_wt_CroCII(:, 6), ...
        '-', 'Color', cm(4, :));
    p5 = plot(sol2(1).m5_wt_CroCII(:, 1), sol2(1).wm5_wt_CroCII(:, 6), ...
        '-', 'Color', cm(5, :));
    xlabel('Time (min)');
    ylabel('PRE activity (per phage)');
    legend([p1, p2, p3, p4, p5], 'MOI = 1 (P+, K_{cII,Cro}=\infty)', ...
        'MOI = 2 (P+, K_{cII,Cro}=\infty)', ...
        'MOI = 3 (P+, K_{cII,Cro}=\infty)', ...
        'MOI = 4 (P+, K_{cII,Cro}=\infty)', ...
        'MOI = 5 (P+, K_{cII,Cro}=\infty)', 'Color', 'none', 'EdgeColor', 'none');
    pbaspect([2.5, 1, 1]);
    
    
%%
%nPRE DEPENDENCE

fig7b = figure('Name', 'nPRE = 2');
dim = [.2 .5 .3 .3];
str_lett = {'A', 'B', 'C', 'D', 'E', 'F'};
numLet = 6;
for i = 1:numLet
    a_lett = annotation('textbox',dim,'String',str_lett{i},'FitBoxToText','on', ...
        'Color', 'k', 'EdgeColor', 'none', ...
        'BackgroundColor', 'none', 'Fontsize', 16);
end 
subplot(3, 2, 1); %P+
    p1 = plot(sol2(1).m1_wt_nPRE(:, 1), sol2(1).wm1_wt_nPRE(:, 6), '-', 'Color', cm(1, :)); hold on;
    p2 = plot(sol2(1).m2_wt_nPRE(:, 1), sol2(1).wm2_wt_nPRE(:, 6), '-', 'Color', cm(2, :));
    p3 = plot(sol2(1).m3_wt_nPRE(:, 1), sol2(1).wm3_wt_nPRE(:, 6), '-', 'Color', cm(3, :));
    p4 = plot(sol2(1).m4_wt_nPRE(:, 1), sol2(1).wm4_wt_nPRE(:, 6), '-', 'Color', cm(4, :));
    p5 = plot(sol2(1).m5_wt_nPRE(:, 1), sol2(1).wm5_wt_nPRE(:, 6), '-', 'Color', cm(5, :));
    xlabel('Time (min)');
    ylabel('PRE activity (per phage)');
    legend([p1, p2, p3, p4, p5], 'MOI = 1 (P+, n_{PRE}=2)', 'MOI = 2 (P+, n_{PRE}=2)', ...
        'MOI = 3 (P+, n_{PRE}=2)', 'MOI = 4 (P+, n_{PRE}=2)', 'MOI = 5 (P+, n_{PRE}=2)');
    pbaspect([2.5, 1, 1]);
    xlim([0, 65]);
subplot(3, 2, 3); %P+
    p1 = plot(sol2(1).m1_wt(:, 1), sol2(1).wm1_wt(:, 6), '-', 'Color', cm(1, :)); hold on;
    p2 = plot(sol2(1).m2_wt(:, 1), sol2(1).wm2_wt(:, 6), '-', 'Color', cm(2, :));
    p3 = plot(sol2(1).m3_wt(:, 1), sol2(1).wm3_wt(:, 6), '-', 'Color', cm(3, :));
    p4 = plot(sol2(1).m4_wt(:, 1), sol2(1).wm4_wt(:, 6), '-', 'Color', cm(4, :));
    p5 = plot(sol2(1).m5_wt(:, 1), sol2(1).wm5_wt(:, 6), '-', 'Color', cm(5, :));
    xlabel('Time (min)');
    ylabel('PRE activity (per phage)');
    legend([p1, p2, p3, p4, p5], 'MOI = 1 (P+)', 'MOI = 2 (P+)', ...
        'MOI = 3 (P+)', 'MOI = 4 (P+)', 'MOI = 5 (P+)');
    pbaspect([2.5, 1, 1]);
    xlim([0, 65]);
subplot(3, 2, 5); %P+
    p1 = plot(1:5, [sol2(1).tauOn_m1_wt(2)-sol2(1).tauOn_m1_wt(1), ...
        sol2(1).tauOn_m2_wt(2)-sol2(1).tauOn_m2_wt(1), ...
        sol2(1).tauOn_m3_wt(2)-sol2(1).tauOn_m3_wt(1), ...
        sol2(1).tauOn_m4_wt(2)-sol2(1).tauOn_m4_wt(1), ...
        sol2(1).tauOn_m5_wt(2)-sol2(1).tauOn_m5_wt(1)], 'ksq', ...
        'MarkerFaceColor', 'k'); hold on;
    p2 = plot(1:5, [sol2(1).tauOn_m1_wt_nPRE(2)-sol2(1).tauOn_m1_wt_nPRE(1), ...
        sol2(1).tauOn_m2_wt_nPRE(2)-sol2(1).tauOn_m2_wt_nPRE(1), ...
        sol2(1).tauOn_m3_wt_nPRE(2)-sol2(1).tauOn_m3_wt_nPRE(1), ...
        sol2(1).tauOn_m4_wt_nPRE(2)-sol2(1).tauOn_m4_wt_nPRE(1), ...
        sol2(1).tauOn_m5_wt_nPRE(2)-sol2(1).tauOn_m5_wt_nPRE(1)], 'rsq', ...
        'MarkerFaceColor', 'r'); 
    xlabel('MOI');
    ylabel('Turn-on duration (min)');
    legend([p1, p2], 'P+, Fitted n_{PRE}', 'P+, n_{PRE}=2');
    pbaspect([2.5, 1, 1]);   
subplot(3, 2, 2); %P+, nPRE = 2
    p1 = plot(sol2(1).m1_wt_nPRE(:, 5)./CINorm, ...
        sol2(1).m1_wt_nPRE(:, 6)./CroNorm, '-', 'Color', cm(1, :)); hold on;
    p2 = plot(sol2(1).m2_wt_nPRE(:, 5)./CINorm, ...
        sol2(1).m2_wt_nPRE(:, 6)./CroNorm, '-', 'Color', cm(2, :));
    p3 = plot(sol2(1).m3_wt_nPRE(:, 5)./CINorm, ...
        sol2(1).m3_wt_nPRE(:, 6)./CroNorm, '-', 'Color', cm(3, :));
    p4 = plot(sol2(1).m4_wt_nPRE(:, 5)./CINorm, ...
        sol2(1).m4_wt_nPRE(:, 6)./CroNorm, '-', 'Color', cm(4, :));
    p5 = plot(sol2(1).m5_wt_nPRE(:, 5)./CINorm, ...
        sol2(1).m5_wt_nPRE(:, 6)./CroNorm, '-', 'Color', cm(5, :));
    %Plot thresholds
    xl = xlim;
    yl = ylim;
    cro_thresh = plot(xl, [sol2(1).KCro; sol2(1).KCro;]/CroNorm, '--', 'Color', cLyt, ...
            'LineWidth', 1);
    cI_thresh = plot([sol2(1).KCI; sol2(1).KCI;]/CINorm, yl, '--', 'Color', cLys, ...
            'LineWidth', 1);
    dim = [.2 .5 .3 .3];
    str1 = 'Lytic threshold';
    a1 = annotation('textbox',dim,'String',str1,'FitBoxToText','on', ...
        'Color', cLyt, 'EdgeColor', 'none', ...
        'BackgroundColor', 'none');
    str2 = 'Lysogenic threshold';
    a2 = annotation('textbox',dim,'String',str2,'FitBoxToText','on', ...
        'Color', cLys, 'EdgeColor', 'none', ...
        'BackgroundColor', 'none');
    xlabel('CI concentration (normalized)');
    ylabel('Cro concentration (normalized)');
    legend([p1, p2, p3, p4, p5], 'MOI = 1 (P+, n_{PRE}=2)', 'MOI = 2 (P+, n_{PRE}=2)', ...
        'MOI = 3 (P+, n_{PRE}=2)', 'MOI = 4 (P+, n_{PRE}=2)', 'MOI = 5 (P+, n_{PRE}=2)');
    axis square
subplot(3, 2, 4); %P+
    p1 = plot(sol2(1).m1_wt(:, 5)./CINorm, ...
        sol2(1).m1_wt(:, 6)./CroNorm, '-', 'Color', cm(1, :)); hold on;
    p2 = plot(sol2(1).m2_wt(:, 5)./CINorm, ...
        sol2(1).m2_wt(:, 6)./CroNorm, '-', 'Color', cm(2, :));
    p3 = plot(sol2(1).m3_wt(:, 5)./CINorm, ...
        sol2(1).m3_wt(:, 6)./CroNorm, '-', 'Color', cm(3, :));
    p4 = plot(sol2(1).m4_wt(:, 5)./CINorm, ...
        sol2(1).m4_wt(:, 6)./CroNorm, '-', 'Color', cm(4, :));
    p5 = plot(sol2(1).m5_wt(:, 5)./CINorm, ...
        sol2(1).m5_wt(:, 6)./CroNorm, '-', 'Color', cm(5, :));
    %Plot thresholds
    xl = xlim;
    yl = ylim;
    cro_thresh = plot(xl, [sol2(1).KCro; sol2(1).KCro;]/CroNorm, '--', 'Color', cLyt, ...
            'LineWidth', 1);
    cI_thresh = plot([sol2(1).KCI; sol2(1).KCI;]/CINorm, yl, '--', 'Color', cLys, ...
            'LineWidth', 1);
    dim = [.2 .5 .3 .3];
    str1 = 'Lytic threshold';
    a1 = annotation('textbox',dim,'String',str1,'FitBoxToText','on', ...
        'Color', cLyt, 'EdgeColor', 'none', ...
        'BackgroundColor', 'none');
    str2 = 'Lysogenic threshold';
    a2 = annotation('textbox',dim,'String',str2,'FitBoxToText','on', ...
        'Color', cLys, 'EdgeColor', 'none', ...
        'BackgroundColor', 'none');
    xlabel('CI concentration (normalized)');
    ylabel('Cro concentration (normalized)');
    legend([p1, p2, p3, p4, p5], 'MOI = 1 (P+)', 'MOI = 2 (P+)', ...
        'MOI = 3 (P+)', 'MOI = 4 (P+)', 'MOI = 5 (P+)');
    axis square
    
%%
%cII mRNA higher during P+ (no 2nd PRE pulse)

%cII mRNA (MOI=1)==========================================================
fig9a = figure('Name', 'cII mRNA, MOI = 1, P+, P-');
dim = [.2 .5 .3 .3];
str_lett = {'A', 'B', 'C', 'D', 'E', 'F'};
numLet = 3;
for i = 1:numLet
    a_lett = annotation('textbox',dim,'String',str_lett{i},'FitBoxToText','on', ...
        'Color', 'k', 'EdgeColor', 'none', ...
        'BackgroundColor', 'none', 'Fontsize', 16);
end
subplot(3, 1, 1);
    simNorm_cII = max(sol(i).m1Num(:, 4));
    zengNorm_cII = max(data_cII_P(:, 2));
    p1 = plot(sol(1).m1Num(:, 1), sol(1).m1Num(:, 4)./simNorm_cII, '-', 'Color', cm(1, :)); hold on;
    p2 = plot(sol(1).m1Num_wt(:, 1), sol(1).m1Num_wt(:, 4)./simNorm_cII, ...
            '--', 'Color', cm(1, :)); 
    s1 = plot(data_cII_wt(:, 1), data_cII_wt(:, 2)./zengNorm_cII, ...
        'o', 'Color', cm(1, :), 'MarkerFaceColor', cm(1, :), 'MarkerSize', mSize);
    xlabel('Time (min)');
    ylabel('cII mRNA per cell (normalized)');
    legend([p1, p2, s1], 'MOI = 1 (P-)', 'MOI = 1 (P+)', 'Data');
    pbaspect([2.5, 1, 1]);
    xlim([0, 65]);
subplot(3, 1, 2); %CII concentration
    p1 = plot(sol(1).m1(:, 1), sol(1).m1(:, 7)./sol(1).K(4), '-', 'Color', cm(1, :)); hold on;
    p2 = plot(sol(1).m1_wt(:, 1), sol(1).m1_wt(:, 7)./sol(1).K(4), ...
        '--', 'Color', cm(1, :));
    %Plot thresholds
    xl = xlim;
    cII_thresh = plot(xl, [1; 1;], '--', 'Color', cCII, ...
        'LineWidth', 1);
    dim = [.2 .5 .3 .3];
    str1 = 'K_{PRE}';
    a1 = annotation('textbox',dim,'String',str1,'FitBoxToText','on', ...
        'Color', cCII, 'EdgeColor', 'none', ...
        'BackgroundColor', 'none');
    xlabel('Time (min)');
    ylabel('CII concentration (normalized)');
    pbaspect([2.5, 1, 1]);
    xlim([0, 65]);
subplot(3, 1, 3); %PRE activity
    p1 = plot(sol(1).m1Num(:, 1), sol(1).wm1_OP(:, 6), '-', 'Color', cm(1, :)); hold on;
    p2 = plot(sol(1).m1Num_wt(:, 1), sol(1).wm1_wt(:, 6), '--', 'Color', cm(1, :));
    xlabel('Time (min)');
    ylabel('PRE activity (per phage)');
    legend([p1, p2], 'MOI = 1 (P-)', 'MOI = 1 (P+)');
    pbaspect([2.5, 1, 1]);
    xlim([0, 65]);
    
%%
%Cro- Viral replication ratio, CI repression of repl. (KcIICro = Inf)

tMark = 40;

fig10a = figure('Name', 'MOI=1, P+/cro- viral copy number');
p1 = plot(sol2(1).m1Num_wt(:, 1), sol2(1).m1Num_wt(:, end)./sol2(1).m1Num_cro(:, end), ...
    '-', 'Color', cm(1, :)); hold on;
yl = ylim;
p2 = plot([tMark, tMark], yl, 'r--');
xlabel('Time (min)');
ylabel('Viral copy number ratio');
legend('MOI = 1, cro+P+/cro-P+');
pbaspect([2.5, 1, 1]);
    
    
%%
%FUNCTIONS

%--------------------------------------------------------------------------
function outCome = getOutCome(sol, fieldNames)
%Returns the outcome:
%0 - failed infection
%1 - lysogeny
%2 - lysis
%3 - mixed outcome

outCome = zeros(length(fieldNames), 1);
for i = 1:length(fieldNames)
    if isfield(sol, fieldNames{i}) && ~isempty(sol.(fieldNames{i}))
       if max(sol.(fieldNames{i})(:, 5))/sol.KCI < 1 && ...
               max(sol.(fieldNames{i})(:, 6))/sol.KCro < 1
           outCome(i) = 0; %fail
       elseif max(sol.(fieldNames{i})(:, 5))/sol.KCI >= 1 && ...
               max(sol.(fieldNames{i})(:, 6))/sol.KCro >= 1
           outCome(i) = 3; %mixed
       elseif max(sol.(fieldNames{i})(:, 5))/sol.KCI >= 1 && ...
               max(sol.(fieldNames{i})(:, 6))/sol.KCro < 1
           outCome(i) = 1; %lysogeny
       elseif max(sol.(fieldNames{i})(:, 5))/sol.KCI < 1 && ...
               max(sol.(fieldNames{i})(:, 6))/sol.KCro >= 1
           outCome(i) = 2; %lysis
       end
    elseif isfield(sol, fieldNames{i}) && isempty(sol.(fieldNames{i}))
        outCome(i) = -1;
        warning('Field name is empty');
    elseif ~isfield(sol, fieldNames{i})
        outCome(i) = -2;
        warning('Field name does not exist');
    else
        outCome(i) = -3;
        warning('Unknown error');
    end
end

end
