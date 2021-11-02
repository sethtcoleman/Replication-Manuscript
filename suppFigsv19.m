%This script takes a sol2ution data structure from one of the analDec
%scripts and creates figures for the replication paper.

%Load Exp 1
load('sol_exp1.mat');
sol1 = sol;

%Load Exp 2
load('sol_exp2.mat'); %sol_exp2_1_hiRes.mat
sol2 = sol;

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

%Letters & markers
dim = [.2 .5 .3 .3];
str_lett = {'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J'};
mSize = 5; %marker size
dXMax = 62; %max length, x-axis

%CI, Cro normalizations
CINorm = sol2(1).K(6); %KCro_CI
CroNorm = sol2(1).K(3); %KPRM_Cro

%%
%FIG 1: Ensemble of P- fits (Exp 1 & 2)

t_ind = 1:7; %for early average

fig1a = figure('Name', 'P- Ensemble fits');
for i = 1:2
    subplot(3, 2, i); %cI--------------------------------------------------
        %Exp 1
        if i == 1
            for j = 1:length(sol1)
                p1 = plot(sol1(j).m1Num(:, 1), sol1(j).m1Num(:, 2), ...
                    '-', 'Color', cm(1, :)); hold on;
                p2 = plot(sol1(j).m2Num(:, 1), sol1(j).m2Num(:, 2), ...
                    '-', 'Color', cm(2, :));
                p3 = plot(sol1(j).m3Num(:, 1), sol1(j).m3Num(:, 2), ...
                    '-', 'Color', cm(3, :));
                p4 = plot(sol1(j).m4Num(:, 1), sol1(j).m4Num(:, 2), ...
                    '-', 'Color', cm(4, :));
                p5 = plot(sol1(j).m5Num(:, 1), sol1(j).m5Num(:, 2), ...
                    '-', 'Color', cm(5, :));
            end
            s1 = errorbar(ty_inf1_AvgEar_m1(t_ind, 1), ty_inf1_AvgEar_m1(t_ind, 2), ...
                ty_inf1_AvgEar_m1(t_ind, 5), 'o', 'Color', cm(1, :), ...
                'MarkerFaceColor', cm(1, :), 'MarkerSize', mSize);
            s2 = errorbar(ty_inf1_AvgEar_m2(t_ind, 1), ty_inf1_AvgEar_m2(t_ind, 2), ...
                ty_inf1_AvgEar_m2(t_ind, 5), 'o', 'Color', cm(2, :), ...
                'MarkerFaceColor', cm(2, :), 'MarkerSize', mSize);
            s3 = errorbar(ty_inf1_AvgEar_m3(t_ind, 1), ty_inf1_AvgEar_m3(t_ind, 2), ...
                ty_inf1_AvgEar_m3(t_ind, 5), 'o', 'Color', cm(3, :), ...
                'MarkerFaceColor', cm(3, :), 'MarkerSize', mSize);
            s4 = errorbar(ty_inf1_AvgEar_m4(t_ind, 1), ty_inf1_AvgEar_m4(t_ind, 2), ...
                ty_inf1_AvgEar_m4(t_ind, 5), 'o', 'Color', cm(4, :), ...
                'MarkerFaceColor', cm(4, :), 'MarkerSize', mSize);
            s5 = errorbar(ty_inf1_AvgEar_m5(t_ind, 1), ty_inf1_AvgEar_m5(t_ind, 2), ...
                ty_inf1_AvgEar_m5(t_ind, 5), 'o', 'Color', cm(5, :), ...
                'MarkerFaceColor', cm(5, :), 'MarkerSize', mSize);
            title('Exp. 1');
        %Exp 2
        elseif i == 2
            for j = 1:length(sol2)
                p1 = plot(sol2(j).m1Num(:, 1), sol2(j).m1Num(:, 2), ...
                    '-', 'Color', cm(1, :)); hold on;
                p2 = plot(sol2(j).m2Num(:, 1), sol2(j).m2Num(:, 2), ...
                    '-', 'Color', cm(2, :));
                p3 = plot(sol2(j).m3Num(:, 1), sol2(j).m3Num(:, 2), ...
                    '-', 'Color', cm(3, :));
                p4 = plot(sol2(j).m4Num(:, 1), sol2(j).m4Num(:, 2), ...
                    '-', 'Color', cm(4, :));
                p5 = plot(sol2(j).m5Num(:, 1), sol2(j).m5Num(:, 2), ...
                    '-', 'Color', cm(5, :));
            end
            s1 = errorbar(ty_inf2_AvgEar_m1(t_ind, 1), ty_inf2_AvgEar_m1(t_ind, 2), ...
                ty_inf2_AvgEar_m1(t_ind, 5), 'o', 'Color', cm(1, :), ...
                'MarkerFaceColor', cm(1, :), 'MarkerSize', mSize);
            s2 = errorbar(ty_inf2_AvgEar_m2(t_ind, 1), ty_inf2_AvgEar_m2(t_ind, 2), ...
                ty_inf2_AvgEar_m2(t_ind, 5), 'o', 'Color', cm(2, :), ...
                'MarkerFaceColor', cm(2, :), 'MarkerSize', mSize);
            s3 = errorbar(ty_inf2_AvgEar_m3(t_ind, 1), ty_inf2_AvgEar_m3(t_ind, 2), ...
                ty_inf2_AvgEar_m3(t_ind, 5), 'o', 'Color', cm(3, :), ...
                'MarkerFaceColor', cm(3, :), 'MarkerSize', mSize);
            s4 = errorbar(ty_inf2_AvgEar_m4(t_ind, 1), ty_inf2_AvgEar_m4(t_ind, 2), ...
                ty_inf2_AvgEar_m4(t_ind, 5), 'o', 'Color', cm(4, :), ...
                'MarkerFaceColor', cm(4, :), 'MarkerSize', mSize);
            s5 = errorbar(ty_inf2_AvgEar_m5(t_ind, 1), ty_inf2_AvgEar_m5(t_ind, 2), ...
                ty_inf2_AvgEar_m5(t_ind, 5), 'o', 'Color', cm(5, :), ...
                'MarkerFaceColor', cm(5, :), 'MarkerSize', mSize);
            title('Exp. 2');
        end
        legend([p1, p2, p3, p4, p5], 'MOI = 1 (P-)', 'MOI = 2 (P-)', 'MOI = 3 (P-)', ...
        'MOI = 4 (P-)', 'MOI = 5 (P-)', 'Color', 'none', 'EdgeColor', 'none');
        %set(legend, 'NumColumns', 2);
        xlabel('Time (min)');
        ylabel('cI RNA');
        pbaspect([2.5, 1, 1]);
    subplot(3, 2, i+2); %cro-------------------------------------------------
        %Exp 1
        if i == 1
            for j = 1:length(sol1)
                p1 = plot(sol1(j).m1Num(:, 1), sol1(j).m1Num(:, 3), ...
                    '-', 'Color', cm(1, :)); hold on;
                p2 = plot(sol1(j).m2Num(:, 1), sol1(j).m2Num(:, 3), ...
                    '-', 'Color', cm(2, :));
                p3 = plot(sol1(j).m3Num(:, 1), sol1(j).m3Num(:, 3), ...
                    '-', 'Color', cm(3, :));
                p4 = plot(sol1(j).m4Num(:, 1), sol1(j).m4Num(:, 3), ...
                    '-', 'Color', cm(4, :));
                p5 = plot(sol1(j).m5Num(:, 1), sol1(j).m5Num(:, 3), ...
                    '-', 'Color', cm(5, :));
            end
            s1 = errorbar(ty_inf1_AvgEar_m1(t_ind, 1), ty_inf1_AvgEar_m1(t_ind, 3), ...
                ty_inf1_AvgEar_m1(t_ind, 6), 'o', 'Color', cm(1, :), ...
                'MarkerFaceColor', cm(1, :), 'MarkerSize', mSize);
            s2 = errorbar(ty_inf1_AvgEar_m2(t_ind, 1), ty_inf1_AvgEar_m2(t_ind, 3), ...
                ty_inf1_AvgEar_m2(t_ind, 6), 'o', 'Color', cm(2, :), ...
                'MarkerFaceColor', cm(2, :), 'MarkerSize', mSize);
            s3 = errorbar(ty_inf1_AvgEar_m3(t_ind, 1), ty_inf1_AvgEar_m3(t_ind, 3), ...
                ty_inf1_AvgEar_m3(t_ind, 6), 'o', 'Color', cm(3, :), ...
                'MarkerFaceColor', cm(3, :), 'MarkerSize', mSize);
            s4 = errorbar(ty_inf1_AvgEar_m4(t_ind, 1), ty_inf1_AvgEar_m4(t_ind, 3), ...
                ty_inf1_AvgEar_m4(t_ind, 6), 'o', 'Color', cm(4, :), ...
                'MarkerFaceColor', cm(4, :), 'MarkerSize', mSize);
            s5 = errorbar(ty_inf1_AvgEar_m5(t_ind, 1), ty_inf1_AvgEar_m5(t_ind, 3), ...
                ty_inf1_AvgEar_m5(t_ind, 6), 'o', 'Color', cm(5, :), ...
                'MarkerFaceColor', cm(5, :), 'MarkerSize', mSize);
        %Exp 2
        elseif i == 2
            for j = 1:length(sol2)
                p1 = plot(sol2(j).m1Num(:, 1), sol2(j).m1Num(:, 3), ...
                    '-', 'Color', cm(1, :)); hold on;
                p2 = plot(sol2(j).m2Num(:, 1), sol2(j).m2Num(:, 3), ...
                    '-', 'Color', cm(2, :));
                p3 = plot(sol2(j).m3Num(:, 1), sol2(j).m3Num(:, 3), ...
                    '-', 'Color', cm(3, :));
                p4 = plot(sol2(j).m4Num(:, 1), sol2(j).m4Num(:, 3), ...
                    '-', 'Color', cm(4, :));
                p5 = plot(sol2(j).m5Num(:, 1), sol2(j).m5Num(:, 3), ...
                    '-', 'Color', cm(5, :));
            end
            s1 = errorbar(ty_inf2_AvgEar_m1(t_ind, 1), ty_inf2_AvgEar_m1(t_ind, 3), ...
                ty_inf2_AvgEar_m1(t_ind, 6), 'o', 'Color', cm(1, :), ...
                'MarkerFaceColor', cm(1, :), 'MarkerSize', mSize);
            s2 = errorbar(ty_inf2_AvgEar_m2(t_ind, 1), ty_inf2_AvgEar_m2(t_ind, 3), ...
                ty_inf2_AvgEar_m2(t_ind, 6), 'o', 'Color', cm(2, :), ...
                'MarkerFaceColor', cm(2, :), 'MarkerSize', mSize);
            s3 = errorbar(ty_inf2_AvgEar_m3(t_ind, 1), ty_inf2_AvgEar_m3(t_ind, 3), ...
                ty_inf2_AvgEar_m3(t_ind, 6), 'o', 'Color', cm(3, :), ...
                'MarkerFaceColor', cm(3, :), 'MarkerSize', mSize);
            s4 = errorbar(ty_inf2_AvgEar_m4(t_ind, 1), ty_inf2_AvgEar_m4(t_ind, 3), ...
                ty_inf2_AvgEar_m4(t_ind, 6), 'o', 'Color', cm(4, :), ...
                'MarkerFaceColor', cm(4, :), 'MarkerSize', mSize);
            s5 = errorbar(ty_inf2_AvgEar_m5(t_ind, 1), ty_inf2_AvgEar_m5(t_ind, 3), ...
                ty_inf2_AvgEar_m5(t_ind, 6), 'o', 'Color', cm(5, :), ...
                'MarkerFaceColor', cm(5, :), 'MarkerSize', mSize);
        end
        %set(legend, 'NumColumns', 2);
        xlabel('Time (min)');
        ylabel('cro RNA');
        pbaspect([2.5, 1, 1]);
  subplot(3, 2, i+4); %cII-------------------------------------------------
        %Exp 1
        if i == 1
            for j = 1:length(sol1)
                p1 = plot(sol1(j).m1Num(:, 1), sol1(j).m1Num(:, 4), ...
                    '-', 'Color', cm(1, :)); hold on;
                p2 = plot(sol1(j).m2Num(:, 1), sol1(j).m2Num(:, 4), ...
                    '-', 'Color', cm(2, :));
                p3 = plot(sol1(j).m3Num(:, 1), sol1(j).m3Num(:, 4), ...
                    '-', 'Color', cm(3, :));
                p4 = plot(sol1(j).m4Num(:, 1), sol1(j).m4Num(:, 4), ...
                    '-', 'Color', cm(4, :));
                p5 = plot(sol1(j).m5Num(:, 1), sol1(j).m5Num(:, 4), ...
                    '-', 'Color', cm(5, :));
            end
            s1 = errorbar(ty_inf1_AvgEar_m1(t_ind, 1), ty_inf1_AvgEar_m1(t_ind, 4), ...
                ty_inf1_AvgEar_m1(t_ind, 7), 'o', 'Color', cm(1, :), ...
                'MarkerFaceColor', cm(1, :), 'MarkerSize', mSize);
            s2 = errorbar(ty_inf1_AvgEar_m2(t_ind, 1), ty_inf1_AvgEar_m2(t_ind, 4), ...
                ty_inf1_AvgEar_m2(t_ind, 7), 'o', 'Color', cm(2, :), ...
                'MarkerFaceColor', cm(2, :), 'MarkerSize', mSize);
            s3 = errorbar(ty_inf1_AvgEar_m3(t_ind, 1), ty_inf1_AvgEar_m3(t_ind, 4), ...
                ty_inf1_AvgEar_m3(t_ind, 7), 'o', 'Color', cm(3, :), ...
                'MarkerFaceColor', cm(3, :), 'MarkerSize', mSize);
            s4 = errorbar(ty_inf1_AvgEar_m4(t_ind, 1), ty_inf1_AvgEar_m4(t_ind, 4), ...
                ty_inf1_AvgEar_m4(t_ind, 7), 'o', 'Color', cm(4, :), ...
                'MarkerFaceColor', cm(4, :), 'MarkerSize', mSize);
            s5 = errorbar(ty_inf1_AvgEar_m5(t_ind, 1), ty_inf1_AvgEar_m5(t_ind, 4), ...
                ty_inf1_AvgEar_m5(t_ind, 7), 'o', 'Color', cm(5, :), ...
                'MarkerFaceColor', cm(5, :), 'MarkerSize', mSize);
        %Exp 2
        elseif i == 2
            for j = 1:length(sol2)
                p1 = plot(sol2(j).m1Num(:, 1), sol2(j).m1Num(:, 4), ...
                    '-', 'Color', cm(1, :)); hold on;
                p2 = plot(sol2(j).m2Num(:, 1), sol2(j).m2Num(:, 4), ...
                    '-', 'Color', cm(2, :));
                p3 = plot(sol2(j).m3Num(:, 1), sol2(j).m3Num(:, 4), ...
                    '-', 'Color', cm(3, :));
                p4 = plot(sol2(j).m4Num(:, 1), sol2(j).m4Num(:, 4), ...
                    '-', 'Color', cm(4, :));
                p5 = plot(sol2(j).m5Num(:, 1), sol2(j).m5Num(:, 4), ...
                    '-', 'Color', cm(5, :));
            end
            s1 = errorbar(ty_inf2_AvgEar_m1(t_ind, 1), ty_inf2_AvgEar_m1(t_ind, 4), ...
                ty_inf2_AvgEar_m1(t_ind, 7), 'o', 'Color', cm(1, :), ...
                'MarkerFaceColor', cm(1, :), 'MarkerSize', mSize);
            s2 = errorbar(ty_inf2_AvgEar_m2(t_ind, 1), ty_inf2_AvgEar_m2(t_ind, 4), ...
                ty_inf2_AvgEar_m2(t_ind, 7), 'o', 'Color', cm(2, :), ...
                'MarkerFaceColor', cm(2, :), 'MarkerSize', mSize);
            s3 = errorbar(ty_inf2_AvgEar_m3(t_ind, 1), ty_inf2_AvgEar_m3(t_ind, 4), ...
                ty_inf2_AvgEar_m3(t_ind, 7), 'o', 'Color', cm(3, :), ...
                'MarkerFaceColor', cm(3, :), 'MarkerSize', mSize);
            s4 = errorbar(ty_inf2_AvgEar_m4(t_ind, 1), ty_inf2_AvgEar_m4(t_ind, 4), ...
                ty_inf2_AvgEar_m4(t_ind, 7), 'o', 'Color', cm(4, :), ...
                'MarkerFaceColor', cm(4, :), 'MarkerSize', mSize);
            s5 = errorbar(ty_inf2_AvgEar_m5(t_ind, 1), ty_inf2_AvgEar_m5(t_ind, 4), ...
                ty_inf2_AvgEar_m5(t_ind, 7), 'o', 'Color', cm(5, :), ...
                'MarkerFaceColor', cm(5, :), 'MarkerSize', mSize);
        end
        %set(legend, 'NumColumns', 2);
        xlabel('Time (min)');
        ylabel('cII RNA');
        pbaspect([2.5, 1, 1]);
end

set(gcf,'renderer','Painters');
    
%%
%Ensemble P+ 

fig2a = figure('Name', 'P+ ens. (Exp 2)');
%lambda--------------------------------------------------------------------
subplot(1, 1, 1); 
    for i = 1:length(sol2)
        p1 = plot(sol2(i).m1Num_wt(:, 1)-sol2(i).dt, sol2(i).m1Num_wt(:, end), ...
            '-', 'Color', cm(1, :)); hold on;
    end
    s1 = errorbar(PCR(1).time(1:7), PCR(1).mean(1:7), PCR(1).error(1:7), ...
        'o', 'Color', cm(1, :), 'MarkerFaceColor', cm(1, :), 'MarkerSize', mSize);
    xlabel('Time after triggering (min)');
    ylabel('Viral copy (#)');
    legend('MOI = 1 (P+)', 'Color', 'none', 'EdgeColor', 'none');
    pbaspect([2.5, 1, 1]);
    xlim([0, 65]);

%%
%Ensemble Zeng fits (Exp 1 & 2)

zengNorm_cII = max(data_cII_P(:, 2));
zengNorm_cI = max(data_cI_croP(:, 2));
    
fig2b = figure('Name', 'Zeng ens. (Exp 2)');
%cII-----------------------------------------------------------------------
subplot(3, 1, 1); 
    %P-
    for i = 1:length(sol2)
        simNorm_cII = max(sol2(i).m1Num(:, 4));
        p2 = plot(sol2(i).m1Num(:, 1), sol2(i).m1Num(:, 4)./simNorm_cII, ...
            'b-'); hold on;
    end
    %cro-
    for i = 1:length(sol2)
        simNorm_cII = max(sol2(i).m1Num(:, 4));
        p3 = plot(sol2(i).m1Num_cro(:, 1), sol2(i).m1Num_cro(:, 4)./simNorm_cII, ...
            'g-'); 
    end
    %cI-
    for i = 1:length(sol2)
        simNorm_cII = max(sol2(i).m1Num(:, 4));
        p5 = plot(sol2(i).m1Num_cI(:, 1), sol2(i).m1Num_cI(:, 4)./simNorm_cII, ...
            'm-');
    end
    %P-
    s2 = plot(data_cII_P(:, 1), data_cII_P(:, 2)./zengNorm_cII, ...
        'o', 'Color', 'b', 'MarkerFaceColor', 'b', 'MarkerSize', mSize);
    %Cro-
    s3 = plot(data_cII_cro(:, 1), data_cII_cro(:, 2)./zengNorm_cII, ...
        'o', 'Color', 'g', 'MarkerFaceColor', 'g', 'MarkerSize', mSize);
    %CI-
    s5 = plot(data_cII_cI(:, 1), data_cII_cI(:, 2)./zengNorm_cII, ...
        'o', 'Color', 'm', 'MarkerFaceColor', 'm', 'MarkerSize', mSize);
    xlim([0, 50]);
    xlabel('Time (min)');
    ylabel('cII mRNA (norm. by P-)');
    legend([p2, p3, p5], 'P-', 'Cro-', 'CI-', ...
        'EdgeColor', 'none', 'Color', 'none');
    pbaspect([2.5, 1, 1]);
%cII RNA 2-----------------------------------------------------------------
subplot(3, 1, 2);
    %P+
    for i = 1:length(sol2)
        simNorm_cII = max(sol2(i).m1Num(:, 4));
        p1 = plot(sol2(i).m1Num_wt(:, 1), sol2(i).m1Num_wt(:, 4)./simNorm_cII, ...
            '-', 'Color', cm(1, :)); hold on;
    end
    %cro-P-
    for i = 1:length(sol2)
        simNorm_cII = max(sol2(i).m1Num(:, 4));
        p4 = plot(sol2(i).m1Num_OP_cro(:, 1), sol2(i).m1Num_OP_cro(:, 4)./simNorm_cII, ...
            'r-');
    end
    %WT
    s1 = plot(data_cII_wt(:, 1), data_cII_wt(:, 2)./zengNorm_cII, ...
        'o', 'Color', cm(1, :), 'MarkerFaceColor', cm(1, :), 'MarkerSize', mSize);
    %Cro-P-
    s4 = plot(data_cII_croP(:, 1), data_cII_croP(:, 2)./zengNorm_cII, ...
        'o', 'Color', 'r', 'MarkerFaceColor', 'r', 'MarkerSize', mSize);
    xlim([0, 50]);
    xlabel('Time (min)');
    ylabel('cII mRNA (norm. by P-)');
    legend([p1, p4], 'P+', 'Cro-P-', ...
        'EdgeColor', 'none', 'Color', 'none');
    pbaspect([2.5, 1, 1]);
%cI RNA--------------------------------------------------------------------
subplot(3, 1, 3);
    %P+
    for i = 1:length(sol2)
        simNorm_cI = max(sol2(i).m1Num_OP_cro(:, 2));
        p1 = plot(sol2(i).m1Num_wt(:, 1), sol2(i).m1Num_wt(:, 2)./simNorm_cI, ...
            '-', 'Color', cm(1, :)); hold on;
    end
    %cro-P-
    for i = 1:length(sol2)
        simNorm_cI = max(sol2(i).m1Num_OP_cro(:, 2));
        p2 = plot(sol2(i).m1Num_OP_cro(:, 1), sol2(i).m1Num_OP_cro(:, 2)./simNorm_cI, ...
            'r-');
    end
    %WT
    s1 = plot(data_cI_wt(:, 1), data_cI_wt(:, 2)./zengNorm_cI, ...
        'o', 'Color', cm(1, :), 'MarkerFaceColor', cm(1, :), 'MarkerSize', mSize);
    %Cro-P-
    s2 = plot(data_cI_croP(:, 1), data_cI_croP(:, 2)./zengNorm_cI, ...
        'o', 'Color', 'r', 'MarkerFaceColor', 'r', 'MarkerSize', mSize);
    xlim([0, 50]);
    xlabel('Time (min)');
    ylabel('cI mRNA (norm. by cro-P-)');
    pbaspect([2.5, 1, 1]);
    
    
%%
%REMOVAL OF PRM

fig3a = figure();
subplot(1, 3, 1); %fraction of cI expr.
    p1 = plot(sol2(1).MOI, sol2(1).fm_tot_wt(:, 3)./...
        (sol2(1).fm_tot_wt(:, 1) + sol2(1).fm_tot_wt(:, 2) + sol2(1).fm_tot_wt(:, 3)), ...
        '-', 'Color', cCII, 'LineWidth', 2); hold on;
    p2 = plot(sol2(1).MOI, (sol2(1).fm_tot_wt(:, 1) + sol2(1).fm_tot_wt(:, 2))./...
        (sol2(1).fm_tot_wt(:, 1) + sol2(1).fm_tot_wt(:, 2) + sol2(1).fm_tot_wt(:, 3)), ...
        '-', 'Color', cLys, 'LineWidth', 2);
    title('P+, frac. cI expr.');
    xlim([sol2(i).MOI(1), sol2(i).MOI(end)]);
    xlabel('MOI');
    ylabel('fraction of cI expression');
    legend([p1, p2], 'P_{RE}', 'P_{RM}');
    axis square
%PRM- ---------------------------------------------------------------------
subplot(1, 3, 2);
    p1 = plot(sol2(1).m1_pert1_wt(:, 5)./CINorm, sol2(1).m1_pert1_wt(:, 6)./CroNorm, ...
        '-', 'Color', cm(1, :), 'LineWidth', 2); hold on;
    p2 = plot(sol2(1).m2_pert1_wt(:, 5)./CINorm, sol2(1).m2_pert1_wt(:, 6)./CroNorm, ...
        '-', 'Color', cm(2, :), 'LineWidth', 2);
    p3 = plot(sol2(1).m3_pert1_wt(:, 5)./CINorm, sol2(1).m3_pert1_wt(:, 6)./CroNorm, ...
        '-', 'Color', cm(3, :), 'LineWidth', 2);
    p4 = plot(sol2(1).m4_pert1_wt(:, 5)./CINorm, sol2(1).m4_pert1_wt(:, 6)./CroNorm, ...
        '-', 'Color', cm(4, :), 'LineWidth', 2);
    p5 = plot(sol2(1).m5_pert1_wt(:, 5)./CINorm, sol2(1).m5_pert1_wt(:, 6)./CroNorm, ...
        '-', 'Color', cm(5, :), 'LineWidth', 2);
    %Plot thresholds
    xl = xlim; 
    yl = ylim;
    cro_thresh = plot(xl, [sol2(1).KCro/CroNorm; sol2(1).KCro/CroNorm;], '--', 'Color', cLyt, ...
        'LineWidth', 1);
    cI_thresh = plot([sol2(1).KCI/CINorm; sol2(1).KCI/CINorm;], yl, '--', 'Color', cLys, ...
        'LineWidth', 1);
    xlabel('CI concentration (normalized)');
    ylabel('Cro concentration (normalized)');
    title('P+ (PRM-)');
    legend([p1, p2, p3, p4, p5], 'MOI = 1 (PRM-)', 'MOI = 2 (PRM-)', 'MOI = 3 (PRM-)', ...
        'MOI = 4 (PRM-)', 'MOI = 5 (PRM-)', 'Color', 'none', 'EdgeColor', 'none');
    axis square
%PRE- ---------------------------------------------------------------------
subplot(1, 3, 3); 
    p1 = plot(sol2(i).m1_pert2_wt(:, 5)./CINorm, sol2(i).m1_pert2_wt(:, 6)./CroNorm, ...
        '-', 'Color', cm(1, :), 'LineWidth', 2); hold on;
    p2 = plot(sol2(i).m2_pert2_wt(:, 5)./CINorm, sol2(i).m2_pert2_wt(:, 6)./CroNorm, ...
        '-', 'Color', cm(2, :), 'LineWidth', 2);
    p3 = plot(sol2(i).m3_pert2_wt(:, 5)./CINorm, sol2(i).m3_pert2_wt(:, 6)./CroNorm, ...
        '-', 'Color', cm(3, :), 'LineWidth', 2);
    p4 = plot(sol2(i).m4_pert2_wt(:, 5)./CINorm, sol2(i).m4_pert2_wt(:, 6)./CroNorm, ...
        '-', 'Color', cm(4, :), 'LineWidth', 2);
    p5 = plot(sol2(i).m5_pert2_wt(:, 5)./CINorm, sol2(i).m5_pert2_wt(:, 6)./CroNorm, ...
        '-', 'Color', cm(5, :), 'LineWidth', 2);
    %Plot thresholds
    xl = xlim; xl(1) = eps; xl(2) = 15;
    yl = ylim; yl(1) = eps;
    cro_thresh = plot(xl, [sol2(1).KCro/CroNorm; sol2(1).KCro/CroNorm;], '--', 'Color', cLyt, ...
        'LineWidth', 1);
    cI_thresh = plot([sol2(1).KCI/CINorm; sol2(1).KCI/CINorm;], yl, '--', 'Color', cLys, ...
        'LineWidth', 1);
    xlabel('CI concentration (normalized)');
    ylabel('Cro concentration (normalized)');
    title('P+ (PRE-)');
    axis square
    
%%
%PRE is MOI-independent (P+)

fig4 = figure('Name', 'PRE MOI-independence, P+');
subplot(2, 2, 1); %Maximum PRE activity (per phage)
    y = [max(sol2(1).wm1_wt(:, 6)), max(sol2(1).wm2_wt(:, 6)), ...
        max(sol2(1).wm3_wt(:, 6)), max(sol2(1).wm4_wt(:, 6)), max(sol2(1).wm5_wt(:, 6))];
    for i = 1:5
        p1 = bar(i, y(i)); hold on;
        set(p1, 'FaceColor', cm(i, :));
    end
    ylabel('Maximum PRE activity (per phage)');
    xlabel('MOI');
    pbaspect([2, 1, 1]);
subplot(2, 2, 2); %Activation time
    y = [sol2(1).tauOn_m1_wt(1), sol2(1).tauOn_m2_wt(1), ...
        sol2(1).tauOn_m3_wt(1), sol2(1).tauOn_m4_wt(1), sol2(1).tauOn_m5_wt(1)];
    for i = 1:5
        p1 = bar(i, y(i)); hold on;
        set(p1, 'FaceColor', cm(i, :));
    end
    ylabel('Activation time (min)');
    xlabel('MOI');
    pbaspect([2, 1, 1]);
subplot(2, 2, 3); %Deactivation time
    y = [sol2(1).tauOn_m1_wt(2), sol2(1).tauOn_m2_wt(2), ...
        sol2(1).tauOn_m3_wt(2), sol2(1).tauOn_m4_wt(2), sol2(1).tauOn_m5_wt(2)];
    for i = 1:5
        p1 = bar(i, y(i)); hold on;
        set(p1, 'FaceColor', cm(i, :));
    end
    ylabel('Deactivation time (min)');
    xlabel('MOI');
    pbaspect([2, 1, 1]);
subplot(2, 2, 4); %Deactivation time
    y = [sol2(1).tauOn_m1_wt(2)-sol2(1).tauOn_m1_wt(1), ...
        sol2(1).tauOn_m2_wt(2)-sol2(1).tauOn_m2_wt(1), ...
        sol2(1).tauOn_m3_wt(2)-sol2(1).tauOn_m3_wt(1), ...
        sol2(1).tauOn_m4_wt(2)-sol2(1).tauOn_m4_wt(1), ...
        sol2(1).tauOn_m5_wt(2)-sol2(1).tauOn_m5_wt(1)];
    for i = 1:5
        p1 = bar(i, y(i)); hold on;
        set(p1, 'FaceColor', cm(i, :));
    end
    ylabel('Turn-on duration (min)');
    xlabel('MOI');
    pbaspect([2, 1, 1]);
    
%%
%DELAY PD

dim = [.2 .5 .3 .3];

%Phase Diagram stuff
tauFail = 1e3;
tD_PD = sol2(1).tD_PD;
dMOI_PD = sol2(1).dMOI_PD;
PD_CI = zeros(length(tD_PD), length(dMOI_PD));
PD_Cro = zeros(length(tD_PD), length(dMOI_PD));
PD = zeros(length(tD_PD), length(dMOI_PD));
for j = 1:length(tD_PD)
    for k = 1:length(dMOI_PD)
        yPD = zeros(1, 4); %[CI]/KCI, [Cro]/KCro, tauCI, tauCro
        yPD(1, 1) = sol2(1).PDMat_tD(j, k, 1);
        yPD(1, 2) = sol2(1).PDMat_tD(j, k, 2);
        yPD(1, end-1) = sol2(1).tauDec_PD_tD(j, k, 1);
        yPD(1, end) = sol2(1).tauDec_PD_tD(j, k, 2);
        PD_CI(j, k) = yPD(:, 1); 
        PD_Cro(j, k) = yPD(:, 2); 
        %Find sol2utions that crossed thresholds
        [r_cI, c_cI] = find(yPD(:, 3) ~= tauFail);
        [r_cro, c_cro] = find(yPD(:, 4) ~= tauFail);
        %CI Cro
        if PD_CI(j, k) < 1 && PD_Cro(j, k) < 1
            PD(j, k) = 3; %fail
        elseif PD_CI(j, k) >= 1 && PD_Cro(j, k) < 1
            PD(j, k) = 1; %lys.
        elseif PD_CI(j, k) < 1 && PD_Cro(j, k) >= 1
            PD(j, k) = 2; %lyt.
        elseif PD_CI(j, k) >= 1 && PD_Cro(j, k) >= 1 
            PD(j, k) = 4; %mixed
        end
    end
end

%Mixed outcomes
figDelayPD = figure('Name', 'Delay PD');
subplot(1, 1, 1); %phase diagram
    h3 = imagesc(dMOI_PD, tD_PD, PD);
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
    xlabel('dMOI');
    ylabel('\tau_{d}/\tau_{PRE}');
    ax = gca;
    ax.YDir = 'normal';
    axis square;

%%
%CII Activity

[i_delay1, ~] = find(sol2(1).m2_d7_wt(:, 1) == sol2(1).tauOn_m1_wt(2)/3);
[i_delay2, ~] = find(sol2(1).m2_d4_wt(:, 1) == sol2(1).tauOn_m1_wt(2));
   
fig5b = figure('Name', 'CII activity/trans.');
    tauPRE1 = [sol2(1).tauOn_m1_wt(1), sol2(1).tauOn_m2_wt(1), ...
        sol2(1).tauOn_m3_wt(1), sol2(1).tauOn_m4_wt(1), sol2(1).tauOn_m5_wt(1)];
    tauPRE2 = [sol2(1).tauOn_m1_wt(2), sol2(1).tauOn_m2_wt(2), ...
        sol2(1).tauOn_m3_wt(2), sol2(1).tauOn_m4_wt(2), sol2(1).tauOn_m5_wt(2)];
gtau = fill([mean(tauPRE1), mean(tauPRE1), mean(tauPRE2), mean(tauPRE2)], ...
    [0, yl(2), yl(2), 0], cShade, ...
    'LineStyle', 'none'); hold on;
str1 = 'CII activity window';
a1 = annotation('textbox',dim,'String',str1,'FitBoxToText','on', ...
    'Color', cCII, 'EdgeColor', 'none', ...
    'BackgroundColor', 'none');
%Plot best solution
s1 = plot(sol2(1).m1_wt(:, 1), sol2(1).m1_wt(:, 7)./sol2(1).K(4), '-', 'Color', ...
    cm(1, :), 'LineWidth', 1.5); hold on;
s2 = plot(sol2(1).m2_wt(:, 1), sol2(1).m2_wt(:, 7)./sol2(1).K(4), '-', 'Color', ...
    cm(4, :), 'LineWidth', 1.5);
s3 = plot(sol2(1).m2_d7_wt(:, 1), sol2(1).m2_d7_wt(:, 7)./sol2(1).K(4), '--', ...
    'Color', cm(4, :), 'LineWidth', 1.5); hold on;
s4 = plot(sol2(1).m2_d4_wt(:, 1), sol2(1).m2_d4_wt(:, 7)./sol2(1).K(4), ':', ...
    'Color', cm(4, :), 'LineWidth', 1.5);
xl = xlim;
cII_thresh = plot(xl, [1; 1;], '--', 'Color', 'k', ...
    'LineWidth', 1);
legend([s1, s2, s3, s4, s5], ...
    'MOI = 1', 'MOI = 2', 'MOI = 1 + 1', 'MOI = 1 + 2');
xlabel('Time (min)');
ylabel('CII concentration (normalized)');
pbaspect([2.5, 1, 1]);
    
%%
%P- PRE ACTIVITY

figPMinPRE = figure('Name', 'P- PRE activity');
subplot(3, 2, 1); %PRE activation (per copy)
    %Plot tau_{PRE}
    yl = [0, 3];
    tauPRE1 = [sol2(1).tauOn_m1(1), sol2(1).tauOn_m2(1), ...
        sol2(1).tauOn_m3(1), sol2(1).tauOn_m4(1), sol2(1).tauOn_m5(1)];
    tauPRE2 = [sol2(1).tauOn_m1(2), sol2(1).tauOn_m2(2), ...
        sol2(1).tauOn_m3(2), sol2(1).tauOn_m4(2), sol2(1).tauOn_m5(2)];
    gtau = fill([mean(tauPRE1), mean(tauPRE1), mean(tauPRE2), mean(tauPRE2)], ...
        [0, yl(2), yl(2), 0], cShade, 'LineStyle', 'none'); hold on;
    dim = [.2 .5 .3 .3];
    str1 = 'CII activity window';
    a1 = annotation('textbox',dim,'String',str1,'FitBoxToText','on', ...
        'Color', cCII, 'EdgeColor', 'none', 'BackgroundColor', 'none');
    %Plot best sol2ution
    s1 = plot(sol2(1).m1(:, 1), sol2(1).wm1_OP(:, 6), '-', 'Color', cm(1, :), ...
        'LineWidth', 1.5); hold on;
    s2 = plot(sol2(1).m2(:, 1), sol2(1).wm2_OP(:, 6), '-', 'Color', cm(2, :), ...
        'LineWidth', 1.5);
    s3 = plot(sol2(1).m3(:, 1), sol2(1).wm3_OP(:, 6), '-', 'Color', cm(3, :), ...
        'LineWidth', 1.5);
    s4 = plot(sol2(1).m4(:, 1), sol2(1).wm4_OP(:, 6), '-', 'Color', cm(4, :), ...
        'LineWidth', 1.5);
    s5 = plot(sol2(1).m5(:, 1), sol2(1).wm5_OP(:, 6), '-', 'Color', cm(5, :), ...
        'LineWidth', 1.5);
    xlabel('Time (min)');
    ylabel('P_{RE} activity');
    legend([s1, s2, s3, s4, s5], ...
        'MOI = 1 (P-)', 'MOI = 2 (P-)', 'MOI = 3 (P-)', ...
        'MOI = 4 (P-)', 'MOI = 5 (P-)', 'Color', 'none', 'EdgeColor', 'none');
    ylim([0, 1]);
    pbaspect([2, 1, 1]);
subplot(3, 2, 2); %Maximum PRE activity (per phage)
    y = [max(sol2(1).wm1_OP(:, 6)), max(sol2(1).wm2_OP(:, 6)), ...
        max(sol2(1).wm3_OP(:, 6)), max(sol2(1).wm4_OP(:, 6)), ...
        max(sol2(1).wm5_OP(:, 6))];
    for i = 1:5
        p1 = bar(i, y(i)); hold on;
        set(p1, 'FaceColor', cm(i, :));
    end
    ylabel('Maximum PRE activity (per phage)');
    xlabel('MOI');
    pbaspect([2, 1, 1]);
subplot(3, 2, 3); %Activation time
    y = [sol2(1).tauOn_m1(1), sol2(1).tauOn_m2(1), ...
        sol2(1).tauOn_m3(1), sol2(1).tauOn_m4(1), sol2(1).tauOn_m5(1)];
    for i = 1:5
        p1 = bar(i, y(i)); hold on;
        set(p1, 'FaceColor', cm(i, :));
    end
    ylabel('Activation time (min)');
    xlabel('MOI');
    pbaspect([2, 1, 1]);
subplot(3, 2, 4); %Deactivation time
    y = [sol2(1).tauOn_m1(2), sol2(1).tauOn_m2(2), ...
        sol2(1).tauOn_m3(2), sol2(1).tauOn_m4(2), sol2(1).tauOn_m5(2)];
    for i = 1:5
        p1 = bar(i, y(i)); hold on;
        set(p1, 'FaceColor', cm(i, :));
    end
    ylabel('Deactivation time (min)');
    xlabel('MOI');
    pbaspect([2, 1, 1]);
subplot(3, 2, 5); %Deactivation time
    y = [sol2(1).tauOn_m1(2)-sol2(1).tauOn_m1(1), ...
        sol2(1).tauOn_m2(2)-sol2(1).tauOn_m2(1), ...
        sol2(1).tauOn_m3(2)-sol2(1).tauOn_m3(1), ...
        sol2(1).tauOn_m4(2)-sol2(1).tauOn_m4(1), ...
        sol2(1).tauOn_m5(2)-sol2(1).tauOn_m5(1)];
    for i = 1:5
        p1 = bar(i, y(i)); hold on;
        set(p1, 'FaceColor', cm(i, :));
    end
    ylabel('Turn-on duration (min)');
    xlabel('MOI');
    pbaspect([2, 1, 1]);
    
%%
%EQUIVALENCE(?) OF rlambda/kCII

dim = [.2 .5 .3 .3];

%Phase Diagram 2 - rlambda/kCII, scan kCII
tauFail = 1e3;
kCII_PD = sol2(1).kCII_PD;
MOI_PD = sol2(1).MOI_PD;
PD_CI = zeros(length(kCII_PD), length(MOI_PD));
PD_Cro = zeros(length(kCII_PD), length(MOI_PD));
PD_1 = zeros(length(kCII_PD), length(MOI_PD));
for j = 1:length(kCII_PD)
    for k = 1:length(MOI_PD)
        yPD = zeros(1, 4); %[CI]/KCI, [Cro]/KCro, tauCI, tauCro
        yPD(1, 1) = sol2(1).PDMat_kCII(j, k, 1);
        yPD(1, 2) = sol2(1).PDMat_kCII(j, k, 2);
        yPD(1, end-1) = sol2(1).tauDec_PD_kCII(j, k, 1);
        yPD(1, end) = sol2(1).tauDec_PD_kCII(j, k, 2);
        PD_CI(j, k) = yPD(:, 1); 
        PD_Cro(j, k) = yPD(:, 2); 
        %CI Cro
        if PD_CI(j, k) < 1 && PD_Cro(j, k) < 1
            PD_1(j, k) = 3; %fail
        elseif PD_CI(j, k) >= 1 && PD_Cro(j, k) < 1
            PD_1(j, k) = 1; %lys.
        elseif PD_CI(j, k) < 1 && PD_Cro(j, k) >= 1
            PD_1(j, k) = 2; %lyt.
        elseif PD_CI(j, k) >= 1 && PD_Cro(j, k) >= 1 
            PD_1(j, k) = 4; %mixed outcome
        end;
    end;
end;

%Plot
figrLamkCII = figure('Name', 'kCII PD');
subplot(1, 1, 1); %normal phase diagram
    h1 = imagesc(sol2(1).MOI_PD, kCII_PD, PD_1); hold on;
    colormap([cLys; cLyt; [1, 1, 1]; [0.5, 0.5, 0.5];]);
    caxis([1, 4]);
    c3 = colorbar('Ticks', [], 'TickLabels', {});
    xlabel('MOI');
    ylabel('CII degradation rate (normalized)');
    ax = gca;
    ax.YDir = 'normal';
    axis square;
    
%%
%PR ACTIVITY (P+, P-)

figPRActivity = figure('Name', 'PR Activity, P+, P-');
subplot(2, 1, 1); %per copy cro expr.
    %Best solution
    p1 = plot(sol(1).m1(:, 1), ...
        sol(1).wm1_OP(:, 7), ...
        '-', 'Color', cm(1, :), 'LineWidth', 1.5); hold on;
    p2 = plot(sol(1).m2(:, 1), ...
        sol(1).wm2_OP(:, 7), ...
        '-', 'Color', cm(2, :), 'LineWidth', 1.5);
    p3 = plot(sol(1).m3(:, 1), ...
        sol(1).wm3_OP(:, 7), ...
        '-', 'Color', cm(3, :), 'LineWidth', 1.5);
    p4 = plot(sol(1).m4(:, 1), ...
        sol(1).wm4_OP(:, 7), ...
        '-', 'Color', cm(4, :), 'LineWidth', 1.5);
    p5 = plot(sol(1).m5(:, 1), ...
        sol(1).wm5_OP(:, 7), ...
        '-', 'Color', cm(5, :), 'LineWidth', 1.5);
    ylim([0, 1]);
    xlabel('Time (min)');
    ylabel('P_{R} activity');
    legend([p1, p2, p3, p4, p5], 'MOI=1 (P-)', 'MOI=2 (P-)', ...
        'MOI=3 (P-)', 'MOI=4 (P-)', 'MOI=5 (P-)', ...
        'Color', 'none', 'EdgeColor', 'none');
    xlim([0, 65]);
    pbaspect([2.5, 1, 1]);
subplot(2, 1, 2); %per copy cro expr.
    %Best solution
    p1 = plot(sol(1).m1_wt(:, 1), ...
        sol(1).wm1_wt(:, 7), ...
        '-', 'Color', cm(1, :), 'LineWidth', 1.5); hold on;
    p2 = plot(sol(1).m2_wt(:, 1), ...
        sol(1).wm2_wt(:, 7), ...
        '-', 'Color', cm(2, :), 'LineWidth', 1.5);
    p3 = plot(sol(1).m3_wt(:, 1), ...
        sol(1).wm3_wt(:, 7), ...
        '-', 'Color', cm(3, :), 'LineWidth', 1.5);
    p4 = plot(sol(1).m4_wt(:, 1), ...
        sol(1).wm4_wt(:, 7), ...
        '-', 'Color', cm(4, :), 'LineWidth', 1.5);
    p5 = plot(sol(1).m5_wt(:, 1), ...
        sol(1).wm5_wt(:, 7), ...
        '-', 'Color', cm(5, :), 'LineWidth', 1.5);
    ylim([0, 1]);
    xlabel('Time (min)');
    ylabel('P_{R} activity');
    legend([p1, p2, p3, p4, p5], 'MOI=1 (P+)', 'MOI=2 (P+)', ...
        'MOI=3 (P+)', 'MOI=4 (P+)', 'MOI=5 (P+)', ...
        'Color', 'none', 'EdgeColor', 'none');
    xlim([0, 65]);
    axis square
    pbaspect([2.5, 1, 1]);  