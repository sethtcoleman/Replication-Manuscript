function cost = repv18_S(params, ty_m1, ...
    ty_m2, ty_m3, ty_m4, ty_m5, normWeight, normWeight_wt, normTYZeng, PCR, ...
    data_cII_P, data_cII_wt, data_cII_croP, data_cII_cro, ...
    data_cII_cI, data_cI_croP, data_cI_wt)
%Runs PSO on repv18

cost = 0;

%%
%Assign parameters

prod = params(1:9);
degr = zeros(8, 1);
[degr(1), degr(2), degr(3), degr(4), degr(5), degr(6), degr(7), degr(8)] = ...
    deal(params(10), params(11), params(12), params(13), params(14), params(15), ...
    params(16), 0); %no viral degr.
n = params(17:27);
K = zeros(11, 1);
[K(1), K(2), K(3), K(4), K(5), K(6), K(7), K(8), K(9), K(10), K(11)] = ...
    deal(params(28), params(29), params(30), params(31), params(32), ...
    params(33), params(34), params(35), params(36), params(37), params(38));
tau = params(end);
kdil = degr(1);

tspan = 0:0.1:60; 
tspan2 = 0:0.1:60; 
convFac = (1e9*1e15)/(6.022e23); %# of um^3 -> nM
V0 = 1; %um^3
V = V0.*exp(kdil.*tspan);
V2 = V0.*exp(kdil.*tspan2);

zengNorm_cII = max(data_cII_P(:, 2));
zengNorm_cI = max(data_cI_croP(:, 2));

%%
%SIMULATE NONREPL.

%P- =======================================================================

%ICs
y0_m1 = [0 0 0 0 0 0];
y0_m2 = [0 0 0 0 0 0];
y0_m3 = [0 0 0 0 0 0];
y0_m4 = [0 0 0 0 0 0];
y0_m5 = [0 0 0 0 0 0];

%Sparsity matrix
S = [
    1 0 0 1 1 1; %cI
    0 1 0 1 1 0; %cro
    0 0 1 1 1 0; %cII
    1 0 0 1 0 0; %CI
    0 1 0 0 1 0; %Cro
    0 0 1 0 0 1; %CII
];

%Solve ODEs
options = odeset('RelTol', 1e-4, 'AbsTol', 1e-6, ...
    'InitialStep', 1e-2, 'JPattern', S); %RelTol 1e-6

%MOI=1
try
    [tm1_OP, ym1_OP] = ode15s(@fv18_S, tspan, y0_m1, options, n([1:8, end]), ...
        prod(1:8), degr, K([1:8, end]), 1, V0, convFac);
%     if any(~isreal(ym1_OP))
%         cost = cost + 1e10;
%     end;
catch
    cost = cost + 1e10;
end;

%MOI=2
try
    [~, ym2_OP] = ode15s(@fv18_S, tspan, y0_m2, options, n([1:8, end]), ...
        prod(1:8), degr, K([1:8, end]), 2, V0, convFac);
%     if any(~isreal(ym2_OP))
%         cost = cost + 1e10;
%     end;
catch
    cost = cost + 1e10;
end;

%MOI=3
try
    [~, ym3_OP] = ode15s(@fv18_S, tspan, y0_m3, options, n([1:8, end]), ...
        prod(1:8), degr, K([1:8, end]), 3, V0, convFac);
%     if any(~isreal(ym3_OP))
%         cost = cost + 1e10;
%     end;
catch
    cost = cost + 1e10;
end;

%MOI=4
try
    [~, ym4_OP] = ode15s(@fv18_S, tspan, y0_m4, options, n([1:8, end]), ...
        prod(1:8), degr, K([1:8, end]), 4, V0, convFac);
%     if any(~isreal(ym4_OP))
%         cost = cost + 1e10;
%     end;
catch
    cost = cost + 1e10;
end;

%MOI=5
try
    [~, ym5_OP] = ode15s(@fv18_S, tspan, y0_m5, options, n([1:8, end]), ...
        prod(1:8), degr, K([1:8, end]), 5, V0, convFac);
%     if any(~isreal(ym5_OP))
%         cost = cost + 1e10;
%     end;
catch
    cost = cost + 1e10;
end;

%Cro-P- ===================================================================
prodCro = prod;
prodCro(5) = 0;
%MOI=1
try
    tic
    [tm1_OP_cro, ym1_OP_cro] = ode15s(@fv18_S, tspan, y0_m1, options, n([1:8, end]), ...
        prodCro(1:8), degr, K([1:8, end]), 1, V0, convFac);
%     if any(~isreal(ym1_OP_cro))
%         cost = cost + 1e10;
%     end;
catch
    cost = cost + 1e10;
end;

if cost >= 1e10
    return
end;

%%
%CALCULATE NONREPL. COST

simNorm_cII = max(ym1_OP(:, 3).*V'./convFac);
simNorm_cI = max(ym1_OP_cro(:, 1).*V'./convFac);

%P- =======================================================================
i_taus = 1:(size(ty_m1, 1) - 2); %t=0.5:60 min leave out 120, 180

[cost_OP, ~] = getCostv4(tm1_OP, ym1_OP, ym2_OP, ym3_OP, ym4_OP, ym5_OP, ...
    [ty_m1(:, 1), ty_m1(:, 2:end)], ...
    [ty_m2(:, 1), ty_m2(:, 2:end)], ...
    [ty_m3(:, 1), ty_m3(:, 2:end)], ...
    [ty_m4(:, 1), ty_m4(:, 2:end)], ...
    [ty_m5(:, 1), ty_m5(:, 2:end)], V, ...
    convFac, i_taus, normWeight);

%cII croP
i_taus_zeng = 1:size(data_cII_croP, 1);
normWeight_zeng = (normTYZeng(3)./(max(data_cII_croP(i_taus_zeng, 2))./zengNorm_cII)).*...
    ones(1, length(i_taus_zeng));
[cost_zeng_cII_croP, ~] = getCostZeng(tm1_OP_cro, ym1_OP_cro(:, 3)./simNorm_cII, ...
    [data_cII_croP(:, 1), data_cII_croP(:, 2)./zengNorm_cII], V, convFac, ...
    i_taus_zeng, normWeight_zeng);

%cI croP
i_taus_zeng = 1:size(data_cI_croP, 1);
normWeight_zeng = (normTYZeng(1)/(max(data_cI_croP(i_taus_zeng, 2))./zengNorm_cI)).*...
    ones(1, length(i_taus_zeng));
[cost_zeng_cI_croP, ~] = getCostZeng(tm1_OP_cro, ym1_OP_cro(:, 1)./simNorm_cI, ...
    [data_cI_croP(:, 1), data_cI_croP(:, 2)./zengNorm_cI], V, convFac, ...
    i_taus_zeng, normWeight_zeng);

cost = cost + cost_OP + cost_zeng_cI_croP + cost_zeng_cII_croP; 

if cost >= 1e10
    return
end;

%%
%SIMULATE REPLICATING

%WT =======================================================================
S_wt = [
    1 0 0 1 1 1 1; %cI
    0 1 0 1 1 0 1; %cro
    0 0 1 1 1 0 1; %cII
    1 0 0 1 0 0 0; %CI
    0 1 0 0 1 0 0; %Cro
    0 0 1 0 0 1 0; %CII
    0 0 0 1 1 0 1; %lambda
];

options2 = odeset('RelTol', 1e-4, 'AbsTol', 1e-6, ...
    'InitialStep', 1e-2, 'JPattern', S_wt); %RelTol 1e-6

%MOI=1
try
    [tm1_wt, ym1_wt] = ode15s(@fv18_repv3_S, tspan2, [y0_m1 1*convFac/V0], ...
        options2, n, prod, degr, K, tau, V0, convFac);
%     if any(~isreal(ym1_wt))
%         cost = cost + 1e10;
%     end;
catch
    cost = cost + 1e10;
end;

%MOI = 5
try
    [tm5_wt, ym5_wt] = ode15s(@fv18_repv3_S, tspan2, [y0_m5 5*convFac/V0], ...
        options2, n, prod, degr, K, tau, V0, convFac);
%     if any(~isreal(ym5_wt))
%         cost = cost + 1e10;
%     end;
catch
    cost = cost + 1e10;
end;

%Cro- =====================================================================
%MOI=1
try
    [tm1_cro, ym1_cro] = ode15s(@fv18_repv3_S, tspan2, [y0_m1 1*convFac/V0], ...
        options2, n, prodCro, degr, K, tau, V0, convFac);
%     if any(~isreal(ym1_cro))
%         cost = cost + 1e10;
%     end;
catch
    cost = cost + 1e10;
end;

%CI- ======================================================================
prodCI = prod;
prodCI(1) = 0;
%MOI=1
try
    [tm1_cI, ym1_cI] = ode15s(@fv18_repv3_S, tspan2, [y0_m1 1*convFac/V0], ...
        options2, n, prodCI, degr, K, tau, V0, convFac);
%     if any(~isreal(ym1_cI))
%         cost = cost + 1e10;
%     end;
catch
    cost = cost + 1e10;
end;

if cost >= 1e10
    return
end;

%%
%CALCULATE REPLICATING COSTS

%WT =======================================================================
i_taus_wt = 1:7;
PCR_fit = PCR;
PCR_fit(1).time = PCR_fit(1).time;
[cost_wt, ~] = getCostThu(tm1_wt, ym1_wt, PCR_fit, V2, convFac, ...
    i_taus_wt, normWeight_wt);

% if ~isreal(cost_wt)
%     cost = 1e10;
%     return
% end;

%cII WT
i_taus_zeng = 1:size(data_cII_wt, 1);
normWeight_zeng = (normTYZeng(3)./(max(data_cII_wt(i_taus_zeng, 2))./zengNorm_cII)).*...
    ones(1, length(i_taus_zeng));
[cost_zeng_cII_wt, ~] = getCostZeng(tm1_wt, ym1_wt(:, 3)./simNorm_cII, ...
    [data_cII_wt(:, 1), data_cII_wt(:, 2)./zengNorm_cII], V2, convFac, ...
    i_taus_zeng, normWeight_zeng);

%cII cro
i_taus_zeng = 1:size(data_cII_cro, 1);
normWeight_zeng = (normTYZeng(3)./(max(data_cII_cro(i_taus_zeng, 2))./zengNorm_cII)).*...
    ones(1, length(i_taus_zeng)); 
[cost_zeng_cII_cro, ~] = getCostZeng(tm1_cro, ym1_cro(:, 3)./simNorm_cII, ...
    [data_cII_cro(:, 1), data_cII_cro(:, 2)./zengNorm_cII], V2, convFac, ...
    i_taus_zeng, normWeight_zeng);

%cII cI
i_taus_zeng = 1:size(data_cII_cI, 1);
normWeight_zeng = (normTYZeng(3)./(max(data_cII_cI(i_taus_zeng, 2))./zengNorm_cII)).*...
    ones(1, length(i_taus_zeng));
[cost_zeng_cII_cI, ~] = getCostZeng(tm1_cI, ym1_cI(:, 3)./simNorm_cII, ...
    [data_cII_cI(:, 1), data_cII_cI(:, 2)./zengNorm_cII], V2, convFac, ...
    i_taus_zeng, normWeight_zeng);

%cI wt
i_taus_zeng = 1:size(data_cI_wt, 1);
normWeight_zeng = (normTYZeng(1)./(max(data_cI_wt(i_taus_zeng, 2))./zengNorm_cI)).*...
    ones(1, length(i_taus_zeng));
[cost_zeng_cI_wt, ~] = getCostZeng(tm1_wt, ym1_wt(:, 1)./simNorm_cI, ...
    [data_cI_wt(:, 1), data_cI_wt(:, 2)./zengNorm_cI], V2, convFac, ...
    i_taus_zeng, normWeight_zeng);

cost = cost + cost_wt + cost_zeng_cII_wt + cost_zeng_cII_cI + cost_zeng_cII_cro + ...
    cost_zeng_cI_wt;

if cost >= 1e10
    return
end;

%%
%REGULARIZED PENALTIES

t_ind1 = 401; %tspan2 < tspan
reg = zeros(26, 1);
regNorm = 0.1*ones(length(reg), 1); %1e-3

%Make sure for cro-, MOI = 1 viral copy number 
if ym1_wt(t_ind1, end)/ym1_cro(t_ind1, end) < 5
    reg(1) = regNorm(1)*( (ym1_wt(t_ind1, end)/ym1_cro(t_ind1, end) - 5)/...
        (ym1_wt(t_ind1, end)/ym1_cro(t_ind1, end)) )^2;
end;

%Make sure rCI < rCII
if prod(4) > prod(8)
    reg(2) = regNorm(2)*( (prod(4) - prod(8))/prod(8) )^2;
end;

%Make sure rPRM <= rPRE
if prod(1)*prod(2) > prod(3)
    reg(3) = regNorm(3)*( (prod(1)*prod(2) - prod(3))/prod(3) )^2;
end;

%Make sure KM,Cro >= KCro,Cro
if K(9) < K(5)
    reg(4) = regNorm(4)*( (K(9) - K(5))/K(9) )^2;
end;

%Make sure KM,CI >= KCro,CI
if K(10) < K(6)
    reg(5) = regNorm(5)*( (K(10) - K(6))/K(10) )^2;
end;

%Make sure rPRM,active < rcro
if prod(1)*prod(2) > prod(5)
    reg(6) = regNorm(6)*( (prod(1)*prod(2) - prod(5))/prod(5) )^2;
end;

%Make sure rPRE < rcro
if prod(3) > prod(5)
    reg(7) = regNorm(7)*( (prod(3) - prod(5))/prod(5) )^2;
end;

%Make sure rCI < rCro
if prod(4) > prod(6)
    reg(8) = regNorm(8)*( (prod(4) - prod(6))/prod(6) )^2;
end;

%Make sure KPRM+ >= KCROCI
if K(1) < K(6)
    reg(9) = regNorm(9)*( (K(1) - K(6))/K(1) )^2;
end;

%Make sure KPRM+ >= KCII,CI
if K(1) < K(8)
    reg(10) = regNorm(10)*( (K(1) - K(8))/K(1) )^2;
end;

%Make sure KPRM,Cro <= KCro,Cro
if K(3) > K(5)
    reg(11) = regNorm(11)*( (K(3) - K(5))/K(5) )^2;
end;

%Make sure KPRM,Cro <= KCII,Cro
if K(3) > K(7)
    reg(12) = regNorm(12)*( (K(3) - K(7))/K(7) )^2;
end;

%Make sure KPRM- > KPRM+
if K(1) > K(2)
    reg(13) = regNorm(13)*( (K(2) - K(1))/K(2) )^2;
end;

%Add point penalty, per Algorithms for Optimization
p = 0.1; %0.1
pointPen = p*sum(reg~=0);

% if ~isreal(sum(reg)) || any(reg) < 0
%     cost = penalty(3);
%     return
% end;

cost = cost + sum(reg) + pointPen;

% if cost >= 1e10
%     return
% end;

% if cost < 40
%    disp(strcat('cost:', num2str(cost_OP), ', cost_{wt}:', num2str(cost_wt), ...
%        ', reg:', num2str(sum(reg)), ', point:', num2str(pointPen))); 
% end

end
