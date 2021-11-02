%This script takes a solution data structure from one of the analDec
%scripts and creates figures for the replication paper.

%Load data structure
load('sol_exp2.mat', 'sol'); %sol_exp2_1_hiRes.mat

%Use best fit
solInd = 1;
sol = sol(solInd);

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
%load('121720_thu_qPCR.mat');
load('042121_thu_qPCR.mat');

%%
%PLOT SETTINGS

%Color map
cm = [
    168 222 233;
    101 185 222;
    75 123 190;
    88 78 160;
    7 24 50;
]./255;

%Other colors
cLyt = [108, 190, 69]./255;
cLys = [237, 28, 36]./255;
cShade = [229, 228, 228]./255;
cCII = [185, 83, 159]./255;

%Letters & Markers
dim = [.2 .5 .3 .3];
str_lett = {'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J'};
mSize = 5; %marker size
dXMax = 65; %max length, x-axis

%Normalizations
CINorm = sol(1).K(6); %KCro_CI
CroNorm = sol(1).K(3); %KPRM_Cro
CIINorm = sol(1).K(4); %KPRE

%%
%FIG 2: Data

t_ind = 1:7;

fig2 = figure('Name', 'Data plots');
numLet = 2;
for i = 1:numLet
    a_lett = annotation('textbox',dim,'String',str_lett{i},'FitBoxToText','on', ...
        'Color', 'k', 'EdgeColor', 'none', ...
        'BackgroundColor', 'none', 'Fontsize', 16);
end
%RNA + Viral copy number---------------------------------------------------
subplot(3, 2, 1); %cI
    %Plot best solution
    s1 = plot(sol(1).m1Num(:, 1), sol(1).m1Num(:, 2), '-', ...
        'Color', cm(1, :), 'LineWidth', 1.5); hold on;
    s2 = plot(sol(1).m2Num(:, 1), sol(1).m2Num(:, 2), '-', ...
        'Color', cm(2, :), 'LineWidth', 1.5);
    s3 = plot(sol(1).m3Num(:, 1), sol(1).m3Num(:, 2), '-', ...
        'Color', cm(3, :), 'LineWidth', 1.5);
    s4 = plot(sol(1).m4Num(:, 1), sol(1).m4Num(:, 2), '-', ...
        'Color', cm(4, :), 'LineWidth', 1.5);
    s5 = plot(sol(1).m5Num(:, 1), sol(1).m5Num(:, 2), '-', ...
        'Color', cm(5, :), 'LineWidth', 1.5);
    %Plot data
    p1 = errorbar(ty_AvgEar_m1(t_ind, 1), ty_AvgEar_m1(t_ind, 2), ty_AvgEar_m1(t_ind, 5), ...
        'o', 'Color', cm(1, :), 'MarkerFaceColor', cm(1, :), 'MarkerSize', mSize); 
    p2 = errorbar(ty_AvgEar_m2(t_ind, 1), ty_AvgEar_m2(t_ind, 2), ty_AvgEar_m2(t_ind, 5), ...
        'o', 'Color', cm(2, :), 'MarkerFaceColor', cm(2, :), 'MarkerSize', mSize);
    p3 = errorbar(ty_AvgEar_m3(t_ind, 1), ty_AvgEar_m3(t_ind, 2), ty_AvgEar_m3(t_ind, 5), ...
        'o', 'Color', cm(3, :), 'MarkerFaceColor', cm(3, :), 'MarkerSize', mSize);
    p4 = errorbar(ty_AvgEar_m4(t_ind, 1), ty_AvgEar_m4(t_ind, 2), ty_AvgEar_m4(t_ind, 5), ...
        'o', 'Color', cm(4, :), 'MarkerFaceColor', cm(4, :), 'MarkerSize', mSize);
    p5 = errorbar(ty_AvgEar_m5(t_ind, 1), ty_AvgEar_m5(t_ind, 2), ty_AvgEar_m5(t_ind, 5), ...
        'o', 'Color', cm(5, :), 'MarkerFaceColor', cm(5, :), 'MarkerSize', mSize);
    legend([s1, s2, s3, s4, s5], 'MOI = 1 (P-)', 'MOI = 2 (P-)', 'MOI = 3 (P-)', ...
        'MOI = 4 (P-)', 'MOI = 5 (P-)', 'Color', 'none', 'EdgeColor', 'none');
    %set(legend, 'NumColumns', 2);
    xlabel('Time (min)');
    ylabel('cI RNA (#)');
    xlim([0, dXMax]);
    pbaspect([2.5, 1, 1]);
subplot(3, 2, 3); %cro
    %Plot best solution
    s1 = plot(sol(1).m1Num(:, 1), sol(1).m1Num(:, 3), '-', ...
        'Color', cm(1, :), 'LineWidth', 1.5); hold on;
    s2 = plot(sol(1).m2Num(:, 1), sol(1).m2Num(:, 3), '-', ...
        'Color', cm(2, :), 'LineWidth', 1.5);
    s3 = plot(sol(1).m3Num(:, 1), sol(1).m3Num(:, 3), '-', ...
        'Color', cm(3, :), 'LineWidth', 1.5);
    s4 = plot(sol(1).m4Num(:, 1), sol(1).m4Num(:, 3), '-', ...
        'Color', cm(4, :), 'LineWidth', 1.5);
    s5 = plot(sol(1).m5Num(:, 1), sol(1).m5Num(:, 3), '-', ...
        'Color', cm(5, :), 'LineWidth', 1.5);
    %Plot data
    p1 = errorbar(ty_AvgEar_m1(t_ind, 1), ty_AvgEar_m1(t_ind, 3), ty_AvgEar_m1(t_ind, 6), ...
        'o', 'Color', cm(1, :), 'MarkerFaceColor', cm(1, :), 'MarkerSize', mSize); 
    p2 = errorbar(ty_AvgEar_m2(t_ind, 1), ty_AvgEar_m2(t_ind, 3), ty_AvgEar_m2(t_ind, 6), ...
        'o', 'Color', cm(2, :), 'MarkerFaceColor', cm(2, :), 'MarkerSize', mSize);
    p3 = errorbar(ty_AvgEar_m3(t_ind, 1), ty_AvgEar_m3(t_ind, 3), ty_AvgEar_m3(t_ind, 6), ...
        'o', 'Color', cm(3, :), 'MarkerFaceColor', cm(3, :), 'MarkerSize', mSize);
    p4 = errorbar(ty_AvgEar_m4(t_ind, 1), ty_AvgEar_m4(t_ind, 3), ty_AvgEar_m4(t_ind, 6), ...
        'o', 'Color', cm(4, :), 'MarkerFaceColor', cm(4, :), 'MarkerSize', mSize);
    p5 = errorbar(ty_AvgEar_m5(t_ind, 1), ty_AvgEar_m5(t_ind, 3), ty_AvgEar_m5(t_ind, 6), ...
        'o', 'Color', cm(5, :), 'MarkerFaceColor', cm(5, :), 'MarkerSize', mSize);
    xlabel('Time (min)');
    ylabel('cro RNA (#)');
    xlim([0, dXMax]);
    pbaspect([2.5, 1, 1]);
subplot(3, 2, 5); %cII
    %Plot best solution
    s1 = plot(sol(1).m1Num(:, 1), sol(1).m1Num(:, 4), '-', ...
        'Color', cm(1, :), 'LineWidth', 1.5); hold on;
    s2 = plot(sol(1).m2Num(:, 1), sol(1).m2Num(:, 4), '-', ...
        'Color', cm(2, :), 'LineWidth', 1.5);
    s3 = plot(sol(1).m3Num(:, 1), sol(1).m3Num(:, 4), '-', ...
        'Color', cm(3, :), 'LineWidth', 1.5);
    s4 = plot(sol(1).m4Num(:, 1), sol(1).m4Num(:, 4), '-', ...
        'Color', cm(4, :), 'LineWidth', 1.5);
    s5 = plot(sol(1).m5Num(:, 1), sol(1).m5Num(:, 4), '-', ...
        'Color', cm(5, :), 'LineWidth', 1.5);
    %Plot data
    p1 = errorbar(ty_AvgEar_m1(t_ind, 1), ty_AvgEar_m1(t_ind, 4), ty_AvgEar_m1(t_ind, 7), ...
        'o', 'Color', cm(1, :), 'MarkerFaceColor', cm(1, :), 'MarkerSize', mSize); 
    p2 = errorbar(ty_AvgEar_m2(t_ind, 1), ty_AvgEar_m2(t_ind, 4), ty_AvgEar_m2(t_ind, 7), ...
        'o', 'Color', cm(2, :), 'MarkerFaceColor', cm(2, :), 'MarkerSize', mSize);
    p3 = errorbar(ty_AvgEar_m3(t_ind, 1), ty_AvgEar_m3(t_ind, 4), ty_AvgEar_m3(t_ind, 7), ...
        'o', 'Color', cm(3, :), 'MarkerFaceColor', cm(3, :), 'MarkerSize', mSize);
    p4 = errorbar(ty_AvgEar_m4(t_ind, 1), ty_AvgEar_m4(t_ind, 4), ty_AvgEar_m4(t_ind, 7), ...
        'o', 'Color', cm(4, :), 'MarkerFaceColor', cm(4, :), 'MarkerSize', mSize);
    p5 = errorbar(ty_AvgEar_m5(t_ind, 1), ty_AvgEar_m5(t_ind, 4), ty_AvgEar_m5(t_ind, 7), ...
        'o', 'Color', cm(5, :), 'MarkerFaceColor', cm(5, :), 'MarkerSize', mSize);
    xlabel('Time (min)');
    ylabel('cII RNA (#)');
    xlim([0, dXMax]);
    pbaspect([2.5, 1, 1]);
subplot(3, 2, 2); %lambda
    %Plot best solution
    s1 = plot(sol(1).m1Num_wt(:, 1), sol(1).m1Num_wt(:, end), '--', ...
        'Color', cm(1, :), 'LineWidth', 1); hold on;
    xl = xlim;
    s2 = plot(xl, [1, 1], '-', 'Color', cm(1, :));
    %Plot data
    p1 = errorbar(PCR(1).time(1:7)+sol(1).dt, PCR(1).mean(1:7), PCR(1).error(1:7), ...
        'sq', 'Color', cm(1, :), 'MarkerFaceColor', cm(1, :), 'MarkerSize', mSize);
    p2 = errorbar(PCR(4).time+sol(1).dt, PCR(4).mean, PCR(4).error, ...
        'o', 'Color', cm(1, :), 'MarkerFaceColor', cm(1, :), 'MarkerSize', mSize);
    xlabel('Time (min)');
    ylabel('Viral copy (#)');
    legend('MOI = 1 (P+)', 'Color', 'none', 'EdgeColor', 'none');
    xlim([0, dXMax]);
    pbaspect([2.5, 1, 1]);
    
%%
%Fig 3: Phase plane traj.
   
%WITH THRESHOLDS===========================================================
fig3a = figure('Name', 'Phase plane traj.');
numLet = 8;
for i = 1:numLet
    a_lett = annotation('textbox',dim,'String',str_lett{i},'FitBoxToText','on', ...
        'Color', 'k', 'EdgeColor', 'none', ...
        'BackgroundColor', 'none', 'Fontsize', 16);
end
%Phase plane traj ---------------------------------------------------------
subplot(5, 2, 9); %P-   
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
subplot(5, 2, 10); %WT
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
%CI traj-------------------------------------------------------------------
subplot(5, 2, 1); %P-
    p1 = plot(sol(1).m1(:, 1), sol(1).m1(:, 5)/CINorm, '-', 'Color', cm(1, :)); hold on;
    p2 = plot(sol(1).m2(:, 1), sol(1).m2(:, 5)/CINorm, '-', 'Color', cm(2, :));
    p3 = plot(sol(1).m3(:, 1), sol(1).m3(:, 5)/CINorm, '-', 'Color', cm(3, :));
    p4 = plot(sol(1).m4(:, 1), sol(1).m4(:, 5)/CINorm, '-', 'Color', cm(4, :));
    p5 = plot(sol(1).m5(:, 1), sol(1).m5(:, 5)/CINorm, '-', 'Color', cm(5, :));
    %Plot decision thresholds
    xl = xlim;
    yl = ylim;
    cI_thresh = plot(xl, [sol(1).KCI/CINorm; sol(1).KCI/CINorm;], '--', 'Color', cLys, ...
        'LineWidth', 1);
    str1 = 'Lysogenic threshold';
    a1 = annotation('textbox',dim,'String',str1,'FitBoxToText','on', ...
        'Color', cLys, 'EdgeColor', 'none', ...
        'BackgroundColor', 'none');
    xlabel('Time (min)');
    ylabel('CI concentration (normalized)');
    legend([p1, p2, p3, p4, p5], 'MOI = 1 (P-)', 'MOI = 2 (P-)', 'MOI = 3 (P-)', ...
        'MOI = 4 (P-)', 'MOI = 5 (P-)', 'Color', 'none', 'EdgeColor', 'none');
    xlim([0, 65]);
    pbaspect([2, 1, 1]);
subplot(5, 2, 2); %P+
    p1 = plot(sol(1).m1_wt(:, 1), sol(1).m1_wt(:, 5)/CINorm, '-', 'Color', cm(1, :)); hold on;
    p2 = plot(sol(1).m2_wt(:, 1), sol(1).m2_wt(:, 5)/CINorm, '-', 'Color', cm(2, :));
    p3 = plot(sol(1).m3_wt(:, 1), sol(1).m3_wt(:, 5)/CINorm, '-', 'Color', cm(3, :));
    p4 = plot(sol(1).m4_wt(:, 1), sol(1).m4_wt(:, 5)/CINorm, '-', 'Color', cm(4, :));
    p5 = plot(sol(1).m5_wt(:, 1), sol(1).m5_wt(:, 5)/CINorm, '-', 'Color', cm(5, :));
    %Plot decision thresholds
    xl = xlim;
    yl = ylim;
    cI_thresh = plot(xl, [sol(1).KCI/CINorm; sol(1).KCI/CINorm;], '--', 'Color', cLys, ...
        'LineWidth', 1);
    dim = [.2 .5 .3 .3];
    str1 = 'Lysogenic threshold';
    a1 = annotation('textbox',dim,'String',str1,'FitBoxToText','on', ...
        'Color', cLys, 'EdgeColor', 'none', ...
        'BackgroundColor', 'none');
    xlabel('Time (min)');
    ylabel('CI concentration (normalized)');
    legend([p1, p2, p3, p4, p5], 'MOI = 1 (P+)', 'MOI = 2 (P+)', 'MOI = 3 (P+)', ...
        'MOI = 4 (P+)', 'MOI = 5 (P+)', 'Color', 'none', 'EdgeColor', 'none');
    xlim([0, 65]);
    pbaspect([2, 1, 1]);
%Cro traj.-----------------------------------------------------------------
subplot(5, 2, 3); %P-
    p1 = plot(sol(1).m1(:, 1), sol(1).m1(:, 6)/CroNorm, '-', 'Color', cm(1, :)); hold on;
    p2 = plot(sol(1).m2(:, 1), sol(1).m2(:, 6)/CroNorm, '-', 'Color', cm(2, :));
    p3 = plot(sol(1).m3(:, 1), sol(1).m3(:, 6)/CroNorm, '-', 'Color', cm(3, :));
    p4 = plot(sol(1).m4(:, 1), sol(1).m4(:, 6)/CroNorm, '-', 'Color', cm(4, :));
    p5 = plot(sol(1).m5(:, 1), sol(1).m5(:, 6)/CroNorm, '-', 'Color', cm(5, :));
    %Plot decision thresholds
    xl = xlim;
    yl = ylim;
    cI_thresh = plot(xl, [sol(1).KCro/CroNorm; sol(1).KCro/CroNorm;], '--', 'Color', cLyt, ...
        'LineWidth', 1);
    dim = [.2 .5 .3 .3];
    str1 = 'Lytic threshold';
    a1 = annotation('textbox',dim,'String',str1,'FitBoxToText','on', ...
        'Color', cLyt, 'EdgeColor', 'none', ...
        'BackgroundColor', 'none');
    xlabel('Time (min)');
    ylabel('Cro concentration (normalized)');
    xlim([0, 65]);
    pbaspect([2, 1, 1]);
subplot(5, 2, 4); %P+
    p1 = plot(sol(1).m1_wt(:, 1), sol(1).m1_wt(:, 6)/CroNorm, '-', 'Color', cm(1, :)); hold on;
    p2 = plot(sol(1).m2_wt(:, 1), sol(1).m2_wt(:, 6)/CroNorm, '-', 'Color', cm(2, :));
    p3 = plot(sol(1).m3_wt(:, 1), sol(1).m3_wt(:, 6)/CroNorm, '-', 'Color', cm(3, :));
    p4 = plot(sol(1).m4_wt(:, 1), sol(1).m4_wt(:, 6)/CroNorm, '-', 'Color', cm(4, :));
    p5 = plot(sol(1).m5_wt(:, 1), sol(1).m5_wt(:, 6)/CroNorm, '-', 'Color', cm(5, :));
     %Plot decision thresholds
    xl = xlim;
    yl = ylim;
    cI_thresh = plot(xl, [sol(1).KCro/CroNorm; sol(1).KCro/CroNorm;], '--', 'Color', cLyt, ...
        'LineWidth', 1);
    dim = [.2 .5 .3 .3];
    str1 = 'Lytic threshold';
    a1 = annotation('textbox',dim,'String',str1,'FitBoxToText','on', ...
        'Color', cLyt, 'EdgeColor', 'none', ...
        'BackgroundColor', 'none');
    xlabel('Time (min)');
    ylabel('Cro concentration (normalized)');
    xlim([0, 65]);
    pbaspect([2, 1, 1]);
%CII-----------------------------------------------------------------------
subplot(5, 2, 5); %P-
    p1 = plot(sol(1).m1(:, 1), sol(1).m1(:, 7)/CIINorm, '-', 'Color', cm(1, :)); hold on;
    p2 = plot(sol(1).m2(:, 1), sol(1).m2(:, 7)/CIINorm, '-', 'Color', cm(2, :));
    p3 = plot(sol(1).m3(:, 1), sol(1).m3(:, 7)/CIINorm, '-', 'Color', cm(3, :));
    p4 = plot(sol(1).m4(:, 1), sol(1).m4(:, 7)/CIINorm, '-', 'Color', cm(4, :));
    p5 = plot(sol(1).m5(:, 1), sol(1).m5(:, 7)/CIINorm, '-', 'Color', cm(5, :));
    xlabel('Time (min)');
    ylabel('CII concentration (normalized)');
    xlim([0, 65]);
    pbaspect([2, 1, 1]);
subplot(5, 2, 6); %P+
    p1 = plot(sol(1).m1_wt(:, 1), sol(1).m1_wt(:, 7)/CIINorm, '-', 'Color', cm(1, :)); hold on;
    p2 = plot(sol(1).m2_wt(:, 1), sol(1).m2_wt(:, 7)/CIINorm, '-', 'Color', cm(2, :));
    p3 = plot(sol(1).m3_wt(:, 1), sol(1).m3_wt(:, 7)/CIINorm, '-', 'Color', cm(3, :));
    p4 = plot(sol(1).m4_wt(:, 1), sol(1).m4_wt(:, 7)/CIINorm, '-', 'Color', cm(4, :));
    p5 = plot(sol(1).m5_wt(:, 1), sol(1).m5_wt(:, 7)/CIINorm, '-', 'Color', cm(5, :));
    xlabel('Time (min)');
    ylabel('CII concentration (normalized)');
    xlim([0, 65]);
    pbaspect([2, 1, 1]);
%Viral copy number---------------------------------------------------------
subplot(5, 2, 7); %P-
    p1 = plot([0, 65], [1, 1], '-', 'Color', cm(1, :)); hold on;
    p1 = plot([0, 65], [2, 2], '-', 'Color', cm(2, :));
    p1 = plot([0, 65], [3, 3], '-', 'Color', cm(3, :));
    p1 = plot([0, 65], [4, 4], '-', 'Color', cm(4, :));
    p1 = plot([0, 65], [5, 5], '-', 'Color', cm(5, :));
    xlabel('Time (min)');
    ylabel('Viral copy (#)');
    xlim([0, 65]);
    pbaspect([2, 1, 1]);
subplot(5, 2, 8); %P+
    p1 = plot(sol(1).m1Num_wt(:, 1), sol(1).m1Num_wt(:, end), '-', 'Color', cm(1, :)); hold on;
    p2 = plot(sol(1).m2Num_wt(:, 1), sol(1).m2Num_wt(:, end), '-', 'Color', cm(2, :));
    p3 = plot(sol(1).m3Num_wt(:, 1), sol(1).m3Num_wt(:, end), '-', 'Color', cm(3, :));
    p4 = plot(sol(1).m4Num_wt(:, 1), sol(1).m4Num_wt(:, end), '-', 'Color', cm(4, :));
    p5 = plot(sol(1).m5Num_wt(:, 1), sol(1).m5Num_wt(:, end), '-', 'Color', cm(5, :));
    xlabel('Time (min)');
    ylabel('Viral copy (#)');
    xlim([0, 65]);
    pbaspect([2, 1, 1]);
        
%%
%FIG 4: CI counts due to CII, turns off lytic functions at MOI >= 2 

fig4 = figure('Name', 'CI counts');
numLet = 5;
for i = 1:numLet
    a_lett = annotation('textbox',dim,'String',str_lett{i},'FitBoxToText','on', ...
        'Color', 'k', 'EdgeColor', 'none', ...
        'BackgroundColor', 'none', 'Fontsize', 16);
end
subplot(3, 2, 1); %PRE activation (per copy)
    %Plot tau_{PRE}
    yl = [0, 3];
    tauPRE1 = [sol(1).tauOn_m1_wt(1), sol(1).tauOn_m2_wt(1), ...
        sol(1).tauOn_m3_wt(1), sol(1).tauOn_m4_wt(1), sol(1).tauOn_m5_wt(1)];
    tauPRE2 = [sol(1).tauOn_m1_wt(2), sol(1).tauOn_m2_wt(2), ...
        sol(1).tauOn_m3_wt(2), sol(1).tauOn_m4_wt(2), sol(1).tauOn_m5_wt(2)];
    gtau = fill([mean(tauPRE1), mean(tauPRE1), mean(tauPRE2), mean(tauPRE2)], ...
        [0, yl(2), yl(2), 0], cShade, 'LineStyle', 'none'); hold on;
    dim = [.2 .5 .3 .3];
    str1 = 'CII activity window';
    a1 = annotation('textbox',dim,'String',str1,'FitBoxToText','on', ...
        'Color', cCII, 'EdgeColor', 'none', 'BackgroundColor', 'none');
    %Plot best solution
    s1 = plot(sol(1).m1_wt(:, 1), sol(1).wm1_wt(:, 6), '-', 'Color', cm(1, :), ...
        'LineWidth', 1.5); hold on;
    s2 = plot(sol(1).m2_wt(:, 1), sol(1).wm2_wt(:, 6), '-', 'Color', cm(2, :), ...
        'LineWidth', 1.5);
    s3 = plot(sol(1).m3_wt(:, 1), sol(1).wm3_wt(:, 6), '-', 'Color', cm(3, :), ...
        'LineWidth', 1.5);
    s4 = plot(sol(1).m4_wt(:, 1), sol(1).wm4_wt(:, 6), '-', 'Color', cm(4, :), ...
        'LineWidth', 1.5);
    s5 = plot(sol(1).m5_wt(:, 1), sol(1).wm5_wt(:, 6), '-', 'Color', cm(5, :), ...
        'LineWidth', 1.5);
    xlabel('Time (min)');
    ylabel('P_{RE} activity');
    legend([s1, s2, s3, s4, s5], ...
        'MOI = 1 (P+)', 'MOI = 2 (P+)', 'MOI = 3 (P+)', ...
        'MOI = 4 (P+)', 'MOI = 5 (P+)', 'Color', 'none', 'EdgeColor', 'none');
    ylim([0, 1]);
    pbaspect([2, 1, 1]);
subplot(3, 2, 3); %CI turns cro off on tau_PRE timescale
    %Plot tau_{PRE}
    yl = [0, 1];
    gtau = fill([mean(tauPRE1), mean(tauPRE1), mean(tauPRE2), mean(tauPRE2)], ...
        [0, yl(2), yl(2), 0], cShade, ...
        'LineStyle', 'none'); hold on;
    dim = [.2 .5 .3 .3];
    str1 = 'CII activity window';
    a1 = annotation('textbox',dim,'String',str1,'FitBoxToText','on', ...
        'Color', cCII, 'EdgeColor', 'none', ...
        'BackgroundColor', 'none');
    %Best solution
    p1 = plot(sol(1).m1_wt(:, 1), ...
        sol(1).wm1_wt(:, 9), ...
        '-', 'Color', cm(1, :), 'LineWidth', 1.5); hold on;
    p2 = plot(sol(1).m2_wt(:, 1), ...
        sol(1).wm2_wt(:, 9), ...
        '-', 'Color', cm(2, :), 'LineWidth', 1.5);
    p3 = plot(sol(1).m3_wt(:, 1), ...
        sol(1).wm3_wt(:, 9), ...
        '-', 'Color', cm(3, :), 'LineWidth', 1.5);
    p4 = plot(sol(1).m4_wt(:, 1), ...
        sol(1).wm4_wt(:, 9), ...
        '-', 'Color', cm(4, :), 'LineWidth', 1.5);
    p5 = plot(sol(1).m5_wt(:, 1), ...
        sol(1).wm5_wt(:, 9), ...
        '-', 'Color', cm(5, :), 'LineWidth', 1.5);
    ylim([0, 1]);
    xlabel('Time (min)');
    ylabel({'Repression of cro'; 'by CI'});
    pbaspect([2, 1, 1]);
subplot(3, 2, 5); %CI turns off replication
    %Plot tau_{PRE}
    yl = [0, 1];
    gtau = fill([mean(tauPRE1), mean(tauPRE1), mean(tauPRE2), mean(tauPRE2)], ...
        [0, yl(2), yl(2), 0], cShade, ...
        'LineStyle', 'none'); hold on;
    dim = [.2 .5 .3 .3];
    str1 = 'CII activity window';
    a1 = annotation('textbox',dim,'String',str1,'FitBoxToText','on', ...
        'Color', cCII, 'EdgeColor', 'none', ...
        'BackgroundColor', 'none');
    %Best solution
    p1 = plot(sol(1).m1_wt(:, 1), ...
        sol(1).wm1_wt(:, 17), ...
        '-', 'Color', cm(1, :), 'LineWidth', 1.5); 
    p2 = plot(sol(1).m2_wt(:, 1), ...
        sol(1).wm2_wt(:, 17), ...
        '-', 'Color', cm(2, :), 'LineWidth', 1.5);
    p3 = plot(sol(1).m3_wt(:, 1), ...
        sol(1).wm3_wt(:, 17), ...
        '-', 'Color', cm(3, :), 'LineWidth', 1.5);
    p4 = plot(sol(1).m4_wt(:, 1), ...
        sol(1).wm4_wt(:, 17), ...
        '-', 'Color', cm(4, :), 'LineWidth', 1.5);
    p5 = plot(sol(1).m5_wt(:, 1), ...
        sol(1).wm5_wt(:, 17), ...
        '-', 'Color', cm(5, :), 'LineWidth', 1.5);
    ylim([0, 1]);
    xlabel('Time (min)');
    ylabel({'Repression of replication'; 'by CI';});
    pbaspect([2, 1, 1]);    
%Cro repression------------------------------------------------------------
subplot(3, 2, 4); %cro repression of cro transcription
    %Plot tau_{PRE}
    yl = [0, 3];
    tauPRE1 = [sol(1).tauOn_m1_wt(1), sol(1).tauOn_m2_wt(1), ...
        sol(1).tauOn_m3_wt(1), sol(1).tauOn_m4_wt(1), sol(1).tauOn_m5_wt(1)];
    tauPRE2 = [sol(1).tauOn_m1_wt(2), sol(1).tauOn_m2_wt(2), ...
        sol(1).tauOn_m3_wt(2), sol(1).tauOn_m4_wt(2), sol(1).tauOn_m5_wt(2)];
    gtau = fill([mean(tauPRE1), mean(tauPRE1), mean(tauPRE2), mean(tauPRE2)], ...
        [0, yl(2), yl(2), 0], cShade, 'LineStyle', 'none'); hold on;
    dim = [.2 .5 .3 .3];
    str1 = 'CII activity window';
    a1 = annotation('textbox',dim,'String',str1,'FitBoxToText','on', ...
        'Color', cCII, 'EdgeColor', 'none', 'BackgroundColor', 'none');
    p1 = plot(sol(1).m1_wt(:, 1), sol(1).wm1_wt(:, 8), '-', 'Color', cm(1, :), ...
        'LineWidth', 1.5); hold on;
    p2 = plot(sol(1).m2_wt(:, 1), sol(1).wm2_wt(:, 8), '-', 'Color', cm(2, :), ...
        'LineWidth', 1.5);
    p3 = plot(sol(1).m3_wt(:, 1), sol(1).wm3_wt(:, 8), '-', 'Color', cm(3, :), ...
        'LineWidth', 1.5);
    p4 = plot(sol(1).m4_wt(:, 1), sol(1).wm4_wt(:, 8), '-', 'Color', cm(4, :), ...
        'LineWidth', 1.5);
    p5 = plot(sol(1).m5_wt(:, 1), sol(1).wm5_wt(:, 8), '-', 'Color', cm(5, :), ...
        'LineWidth', 1.5);
    xlabel('Time (min)');
    ylabel('Repression of cro by Cro');
    pbaspect([2, 1, 1]);  
subplot(3, 2, 6); %cro repression of viral replication
    %Plot tau_{PRE}
    yl = [0, 3];
    tauPRE1 = [sol(1).tauOn_m1_wt(1), sol(1).tauOn_m2_wt(1), ...
        sol(1).tauOn_m3_wt(1), sol(1).tauOn_m4_wt(1), sol(1).tauOn_m5_wt(1)];
    tauPRE2 = [sol(1).tauOn_m1_wt(2), sol(1).tauOn_m2_wt(2), ...
        sol(1).tauOn_m3_wt(2), sol(1).tauOn_m4_wt(2), sol(1).tauOn_m5_wt(2)];
    gtau = fill([mean(tauPRE1), mean(tauPRE1), mean(tauPRE2), mean(tauPRE2)], ...
        [0, yl(2), yl(2), 0], cShade, 'LineStyle', 'none'); hold on;
    dim = [.2 .5 .3 .3];
    str1 = 'CII activity window';
    a1 = annotation('textbox',dim,'String',str1,'FitBoxToText','on', ...
        'Color', cCII, 'EdgeColor', 'none', 'BackgroundColor', 'none');
    p1 = plot(sol(1).m1_wt(:, 1), sol(1).wm1_wt(:, 19), '-', 'Color', cm(1, :), ...
        'LineWidth', 1.5); hold on;
    p2 = plot(sol(1).m2_wt(:, 1), sol(1).wm2_wt(:, 19), '-', 'Color', cm(2, :), ...
        'LineWidth', 1.5);
    p3 = plot(sol(1).m3_wt(:, 1), sol(1).wm3_wt(:, 19), '-', 'Color', cm(3, :), ...
        'LineWidth', 1.5);
    p4 = plot(sol(1).m4_wt(:, 1), sol(1).wm4_wt(:, 19), '-', 'Color', cm(4, :), ...
        'LineWidth', 1.5);
    p5 = plot(sol(1).m5_wt(:, 1), sol(1).wm5_wt(:, 19), '-', 'Color', cm(5, :), ...
        'LineWidth', 1.5);
    xlabel('Time (min)');
    ylabel('Repression of replication by Cro');
    pbaspect([2, 1, 1]);

%%
%Fig 5: The system is insensitive to changes in viral copy number after tau_PRE

[i_delay1, ~] = find(sol(1).m2_d7_wt(:, 1) == sol(1).tauOn_m1_wt(2)/3);
[i_delay2, ~] = find(sol(1).m2_d4_wt(:, 1) == sol(1).tauOn_m1_wt(2));

fig5 = figure('Name', 'delays');
dim = [.2 .5 .3 .3];
str_lett = {'A', 'B', 'C', 'D', 'E', 'F'};
numLet = 3;
for i = 1:numLet
    a_lett = annotation('textbox',dim,'String',str_lett{i},'FitBoxToText','on', ...
        'Color', 'k', 'EdgeColor', 'none', ...
        'BackgroundColor', 'none', 'Fontsize', 16);
end
subplot(1, 3, 1); %CI-Cro - delay sensitivity
    %Best solution
    p1 = plot(sol(1).m1_wt(:, 5)./CINorm, sol(1).m1_wt(:, 6)./CroNorm, ...
        '-', 'Color', cm(1, :), 'LineWidth', 1.5); hold on;
    p2 = plot(sol(1).m2_wt(:, 5)./CINorm, sol(1).m2_wt(:, 6)./CroNorm, ...
        '-', 'Color', cm(4, :), 'LineWidth', 1.5);
    p3 = plot(sol(1).m2_d7_wt(:, 5)./CINorm, sol(1).m2_d7_wt(:, 6)./CroNorm, ...
        '--', 'Color', cm(4, :), 'LineWidth', 1.5); 
    p4 = plot(sol(1).m2_d4_wt(:, 5)./CINorm, sol(1).m2_d4_wt(:, 6)./CroNorm, ...
        ':', 'Color', cm(4, :), 'LineWidth', 1.5);
    %Plot thresholds
    xl = xlim;
    yl = ylim;
    cro_thresh = plot(xl, [sol(1).KCro/CroNorm; sol(1).KCro/CroNorm;], '--', 'Color', cLyt, ...
        'LineWidth', 1);
    cI_thresh = plot([sol(1).KCI/CINorm; sol(1).KCI/CINorm;], yl, '--', 'Color', cLys, ...
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
    ylabel('Cro concentratin (normalized)');
    legend([p3, p4, p1, p2], ...
        'MOI = 2, \tau_d = \tau_{PRE}/3 (P+)', ...
        'MOI = 2, \tau_d = \tau_{PRE} (P+)', ...
        'MOI = 1 (P+)', 'MOI = 2 (P+)', ...
        'Color', 'none', 'EdgeColor', 'none');
    axis square
subplot(1, 3, 2); %PRE per copy activation
    %Plot tau_{PRE}
    yl = [0, 5];
    tauPRE1 = [sol(1).tauOn_m1_wt(1), sol(1).tauOn_m2_wt(1), ...
        sol(1).tauOn_m3_wt(1), sol(1).tauOn_m4_wt(1), sol(1).tauOn_m5_wt(1)];
    tauPRE2 = [sol(1).tauOn_m1_wt(2), sol(1).tauOn_m2_wt(2), ...
        sol(1).tauOn_m3_wt(2), sol(1).tauOn_m4_wt(2), sol(1).tauOn_m5_wt(2)];
    gtau = fill([mean(tauPRE1), mean(tauPRE1), mean(tauPRE2), mean(tauPRE2)], ...
        [0, yl(2), yl(2), 0], cShade, ...
        'LineStyle', 'none'); hold on;
    dim = [.2 .5 .3 .3];
    str1 = 'CII activity window';
    a1 = annotation('textbox',dim,'String',str1,'FitBoxToText','on', ...
        'Color', cCII, 'EdgeColor', 'none', ...
        'BackgroundColor', 'none');
    %Plot best solution
    s1 = plot(sol(1).m1_wt(:, 1), sol(1).wm1_wt(:, 6), '-', 'Color', ...
        cm(1, :), 'LineWidth', 1.5); hold on;
    s2 = plot(sol(1).m2_wt(:, 1), sol(1).wm2_wt(:, 6), '-', 'Color', ...
        cm(4, :), 'LineWidth', 1.5); 
    s3 = plot(sol(1).m2_d7_wt(:, 1), ...
        [zeros(i_delay1-1, 1); ...
        sol(1).wm2_d7_wt(i_delay1:end, 6);], '--', ...
        'Color', cm(4, :), 'LineWidth', 1.5); hold on;
    s4 = plot(sol(1).m2_d4_wt(:, 1), [zeros(i_delay2-1, 1); ...
        sol(1).wm2_d4_wt(i_delay2:end, 6);], ':', ...
        'Color', cm(4, :), 'LineWidth', 1.5);
    xlabel('Time (min)');
    ylabel('P_{RE} activity');
    ylim([0, 1]);
    legend([s1, s2, s3, s4], 'MOI = 1', 'MOI = 2', 'MOI = 1 + 1', 'MOI = 1 + 2');
    pbaspect([2, 1, 1]);
subplot(1, 3, 3); %viral copy number - delay sensitivity
    %Plot tau_{PRE}
    yl = [0, 600];
    tauPRE1 = [sol(1).tauOn_m1_wt(1), sol(1).tauOn_m2_wt(1), ...
        sol(1).tauOn_m3_wt(1), sol(1).tauOn_m4_wt(1), sol(1).tauOn_m5_wt(1)];
    tauPRE2 = [sol(1).tauOn_m1_wt(2), sol(1).tauOn_m2_wt(2), ...
        sol(1).tauOn_m3_wt(2), sol(1).tauOn_m4_wt(2), sol(1).tauOn_m5_wt(2)];
    gtau = fill([mean(tauPRE1), mean(tauPRE1), mean(tauPRE2), mean(tauPRE2)], ...
        [eps, yl(2), yl(2), eps], cShade, ...
        'LineStyle', 'none'); hold on;
    dim = [.2 .5 .3 .3];
    str1 = 'CII activity window';
    a1 = annotation('textbox',dim,'String',str1,'FitBoxToText','on', ...
        'Color', cCII, 'EdgeColor', 'none', ...
        'BackgroundColor', 'none');
    %Best solution
    p1 = plot(sol(1).m1Num_wt(:, 1), sol(1).m1Num_wt(:, end), ...
        '-', 'Color', cm(1, :), 'LineWidth', 1.5); 
    p2 = plot(sol(1).m2Num_wt(:, 1), sol(1).m2Num_wt(:, end), ...
        '-', 'Color', cm(4, :), 'LineWidth', 1.5);
    p3 = plot(sol(1).m2Num_d7_wt(:, 1), sol(1).m2Num_d7_wt(:, end), ...
        '--', 'Color', cm(4, :), 'LineWidth', 1.5); 
    p4 = plot(sol(1).m2Num_d4_wt(:, 1), sol(1).m2Num_d4_wt(:, end), ...
        ':', 'Color', cm(4, :), 'LineWidth', 1.5);
    xlabel('Time (min)');
    ylabel('Viral copy (#)');
    pbaspect([2, 1, 1]);
    
%%
%Fig 6ABCE: P- fails to implement the decision

PD = getPD(sol, 1, 1);

fig6a = figure();
numLet = 5;
for i = 1:numLet
    a_lett = annotation('textbox',dim,'String',str_lett{i},'FitBoxToText','on', ...
        'Color', 'k', 'EdgeColor', 'none', ...
        'BackgroundColor', 'none', 'Fontsize', 16);
end
%lytic---------------------------------------------------------------------
subplot(3, 2, 1); %per copy cro expr.
    %Best solution
    p1 = plot(sol(1).m1(:, 1), ...
        sol(1).wm1_OP(:, 7), ...
        '-', 'Color', cm(1, :), 'LineWidth', 1.5); hold on;
    p2 = plot(sol(1).m1_wt(:, 1), ...
        sol(1).wm1_wt(:, 7), ...
        '--', 'Color', cm(1, :), 'LineWidth', 1.5);
    ylim([0, 1]);
    xlabel('Time (min)');
    ylabel('P_{R} activity');
    legend([p1, p2], 'MOI=1 (P-)', ...
        'MOI=1 (P+)', 'Color', 'none', 'EdgeColor', 'none');
    xlim([0, 65]);
    axis square
subplot(3, 2, 3); %Cro vs. t (MOI = 1)
    %Best solution
    p1 = plot(sol(1).m1(:, 1), ...
        sol(1).m1(:, 6)./CroNorm, ...
        '-', 'Color', cm(1, :), 'LineWidth', 1.5); hold on;
    p2 = plot(sol(1).m1_wt(:, 1), ...
        sol(1).m1_wt(:, 6)./CroNorm, ...
        '--', 'Color', cm(1, :), 'LineWidth', 1.5);
    %Plot thresholds
    xl = xlim;
    cro_thresh = plot(xl, [sol(1).KCro/CroNorm; sol(1).KCro/CroNorm;], '--', 'Color', cLyt, ...
        'LineWidth', 1);
    dim = [.2 .5 .3 .3];
    str = 'Lytic threshold';
    annotation('textbox',dim,'String',str,'FitBoxToText','on', ...
        'Color', cLyt, 'EdgeColor', 'none', ...
        'BackgroundColor', 'none');
    xlabel('Time (min)');
    ylabel('Cro concentration (normalized)');
    xlim([0, 65]);
    axis square
%cI------------------------------------------------------------------------
subplot(3, 2, 5); % [CI]/KCI vs t
    %Best solution
    p1 = plot(sol(1).m2(:, 1), sol(1).m2(:, 5)./CINorm, ...
        '-', 'Color', cm(2, :), 'LineWidth', 1.5); hold on; %P- MOI = 1
    p2 = plot(sol(1).m4(:, 1), sol(1).m4(:, 5)./CINorm, ...
        '-', 'Color', cm(4, :), 'LineWidth', 1.5);
    p3 = plot(sol(1).m2_wt(:, 1), sol(1).m2_wt(:, 5)./CINorm, ...
        '--', 'Color', cm(2, :), 'LineWidth', 1.5); %WT MOI = 2
    %Threshold
    xl = xlim;
    cI_thresh = plot(xl, [sol(1).KCI/CINorm; sol(1).KCI/CINorm;], '--', 'Color', cLys, ...
        'LineWidth', 1);
    dim = [.2 .5 .3 .3];
    str = 'Lysogenic threshold';
    annotation('textbox',dim,'String',str,'FitBoxToText','on', ...
        'Color', cLys, 'EdgeColor', 'none', ...
        'BackgroundColor', 'none');
    xlabel('Time (min)');
    ylabel('CI concentration (normalized)');
    legend([p1, p2, p3], 'MOI=2 (P-)', ...
        'MOI=4 (P-)', 'MOI=2 (P+)', ...
        'Color', 'none', 'EdgeColor', 'none');
    xlim([0, 65]);
    axis square 
%Phase diagram-------------------------------------------------------------
subplot(3, 2, 6); 
    h3 = imagesc(sol(1).MOI_PD, sol(1).rlam_PD, PD);
    colormap([cLys; cLyt; [1, 1, 1]; [0.5, 0.5, 0.5];]);
    caxis([1, 4]);
    str1 = 'Lysogeny';
    annotation('textbox',dim,'String',str1,'FitBoxToText','on', ...
        'EdgeColor', 'none', ...
        'BackgroundColor', 'none');
    str2 = 'Lysis';
    annotation('textbox',dim,'String',str2,'FitBoxToText','on', ...
        'EdgeColor', 'none', ...
        'BackgroundColor', 'none');
    str3 = 'Failed infection';
    annotation('textbox',dim,'String',str3,'FitBoxToText','on', ...
        'EdgeColor', 'none', ...
        'BackgroundColor', 'none');
    str4 = 'Mixed outcome';
    annotation('textbox',dim,'String',str4,'FitBoxToText','on', ...
        'EdgeColor', 'none', ...
        'BackgroundColor', 'none');
    xlabel('MOI');
    ylabel('Replication rate (normalized)');
    ax = gca;
    ax.YDir = 'normal';
    axis square;