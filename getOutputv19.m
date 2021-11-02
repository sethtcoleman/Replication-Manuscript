function [tm1_OP, ym1_OP, tm2_OP, ym2_OP, tm3_OP, ym3_OP, tm4_OP, ym4_OP, ...
    tm5_OP, ym5_OP, fm1_OP, fm5_OP, gm1_OP, gm5_OP, wm1_OP, wm2_OP, ...
    wm3_OP, wm4_OP, wm5_OP, tauOn_m1, i_tauOn_m1, tauOn_m2, i_tauOn_m2, ...
    tauOn_m3, i_tauOn_m3, tauOn_m4, i_tauOn_m4, tauOn_m5, i_tauOn_m5, ...
    tm1_wt, ym1_wt, tm2_wt, ym2_wt, tm3_wt, ym3_wt, tm4_wt, ym4_wt, ...
    tm5_wt, ym5_wt, fm1_wt, fm5_wt, gm1_wt, gm5_wt, wm1_wt, wm2_wt, ...
    wm3_wt, wm4_wt, wm5_wt, tauOn_m1_wt, i_tauOn_m1_wt, tauOn_m2_wt, ...
    i_tauOn_m2_wt, tauOn_m3_wt, i_tauOn_m3_wt, tauOn_m4_wt, i_tauOn_m4_wt, ...
    tauOn_m5_wt, i_tauOn_m5_wt] = getOutputv19(prod, degr, n, K, tau, ...
    options, tspan, convFac, V0)
%Returns model output, both with and without replication

kdil = degr(1);

%P- -------------------------------------------------------------------
%MOI=1
[tm1_OP, ym1_OP] = ode15s(@fv19, tspan, [0 0 0 0 0 0], options, n([1:8, end]), ...
    prod(1:end-1), degr, K([1:8, end]), 1, V0, convFac, 0);

%MOI=2
[tm2_OP, ym2_OP] = ode15s(@fv19, tspan, [0 0 0 0 0 0], options, n([1:8, end]), ...
    prod(1:end-1), degr, K([1:8, end]), 2, V0, convFac, 0);

%MOI=3
[tm3_OP, ym3_OP] = ode15s(@fv19, tspan, [0 0 0 0 0 0], options, n([1:8, end]), ...
    prod(1:end-1), degr, K([1:8, end]), 3, V0, convFac, 0);

%MOI=4
[tm4_OP, ym4_OP] = ode15s(@fv19, tspan, [0 0 0 0 0 0], options, n([1:8, end]), ...
    prod(1:end-1), degr, K([1:8, end]), 4, V0, convFac, 0);

%MOI=5
[tm5_OP, ym5_OP] = ode15s(@fv19, tspan, [0 0 0 0 0 0], options, n([1:8, end]), ...
    prod(1:end-1), degr, K([1:8, end]), 5, V0, convFac, 0);


%production rates
[fm1_OP, gm1_OP] = getFluxv19(tm1_OP, [ym1_OP, 1.*convFac./(V0.*exp(kdil.*tm1_OP))], ...
    n, prod, degr, K, tau, V0, convFac);
[fm5_OP, gm5_OP] = getFluxv19(tm5_OP, [ym5_OP, 5.*convFac./(V0.*exp(kdil.*tm5_OP))], ...
    n, prod, degr, K, tau, V0, convFac);

%weights
wm1_OP = getWeightsv19(tm1_OP, ym1_OP, ...
    n, prod, degr, K, tau, V0, convFac);
wm2_OP = getWeightsv19(tm2_OP, ym2_OP, ...
    n, prod, degr, K, tau, V0, convFac);
wm3_OP = getWeightsv19(tm3_OP, ym3_OP, ...
    n, prod, degr, K, tau, V0, convFac);
wm4_OP = getWeightsv19(tm4_OP, ym4_OP, ...
    n, prod, degr, K, tau, V0, convFac);
wm5_OP = getWeightsv19(tm5_OP, ym5_OP, ...
    n, prod, degr, K, tau, V0, convFac);

[r, c] = find(wm1_OP(:, 6) >= 0.1);
if ~isempty(r)
    tauOn_m1 = [tm1_OP(r(1)), tm1_OP(r(end))];
    i_tauOn_m1 = [r(1), r(end), 1];
end

[r, c] = find(wm2_OP(:, 6) >= 0.1);
if ~isempty(r)
    tauOn_m2 = [tm2_OP(r(1)), tm2_OP(r(end))];
    i_tauOn_m2 = [r(1), r(end), 1];
end

[r, c] = find(wm3_OP(:, 6) >= 0.1);
if ~isempty(r)
    tauOn_m3 = [tm3_OP(r(1)), tm3_OP(r(end))];
    i_tauOn_m3 = [r(1), r(end), 1];
end

[r, c] = find(wm4_OP(:, 6) >= 0.1);
if ~isempty(r)
    tauOn_m4 = [tm4_OP(r(1)), tm4_OP(r(end))];
    i_tauOn_m4 = [r(1), r(end), 1];
end

[r, c] = find(wm5_OP(:, 6) >= 0.1);
if ~isempty(r)
    tauOn_m5 = [tm5_OP(r(1)), tm5_OP(r(end))];
    i_tauOn_m5 = [r(1), r(end), 1];
end

%WT -------------------------------------------------------------------
%MOI=1
[tm1_wt, ym1_wt] = ode15s(@fv19_repv3, tspan, [0 0 0 0 0 0 1*convFac/V0], options, n, prod, ...
    degr, K, tau, V0, convFac);

%MOI=2
[tm2_wt, ym2_wt] = ode15s(@fv19_repv3, tspan, [0 0 0 0 0 0 2*convFac/V0], options, n, prod, ...
    degr, K, tau, V0, convFac);

%MOI=3
[tm3_wt, ym3_wt] = ode15s(@fv19_repv3, tspan, [0 0 0 0 0 0 3*convFac/V0], options, n, prod, ...
    degr, K, tau, V0, convFac);

%MOI=4
[tm4_wt, ym4_wt] = ode15s(@fv19_repv3, tspan, [0 0 0 0 0 0 4*convFac/V0], options, n, prod, ...
    degr, K, tau, V0, convFac);

%MOI=5
[tm5_wt, ym5_wt] = ode15s(@fv19_repv3, tspan, [0 0 0 0 0 0 5*convFac/V0], options, n, prod, ...
    degr, K, tau, V0, convFac);

%production rates
[fm1_wt, gm1_wt] = getFluxv19(tm1_wt, ym1_wt, ...
    n, prod, degr, K, tau, V0, convFac);
[fm5_wt, gm5_wt] = getFluxv19(tm5_wt, ym5_wt, ...
    n, prod, degr, K, tau, V0, convFac);

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

[r, c] = find(wm1_wt(:, 6) >= 0.1);
if ~isempty(r)
    tauOn_m1_wt = [tm1_wt(r(1)), tm1_wt(r(end))];
    i_tauOn_m1_wt = [r(1), r(end), 1];
end

[r, c] = find(wm2_wt(:, 6) >= 0.1);
if ~isempty(r)
    tauOn_m2_wt = [tm2_wt(r(1)), tm2_wt(r(end))];
    i_tauOn_m2_wt = [r(1), r(end), 1];
end

[r, c] = find(wm3_wt(:, 6) >= 0.1);
if ~isempty(r)
    tauOn_m3_wt = [tm3_wt(r(1)), tm3_wt(r(end))];
    i_tauOn_m3_wt = [r(1), r(end), 1];
end

[r, c] = find(wm4_wt(:, 6) >= 0.1);
if ~isempty(r)
    tauOn_m4_wt = [tm4_wt(r(1)), tm4_wt(r(end))];
    i_tauOn_m4_wt = [r(1), r(end), 1];
end

[r, c] = find(wm5_wt(:, 6) >= 0.1);
if ~isempty(r)
    tauOn_m5_wt = [tm5_wt(r(1)), tm5_wt(r(end))];
    i_tauOn_m5_wt = [r(1), r(end), 1];
end


end