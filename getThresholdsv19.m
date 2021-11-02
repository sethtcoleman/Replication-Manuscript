function [KCI, KCro, KCIR, KCroR, KCI_o, KCro_o, KCIR_o, KCroR_o] = ...
    getThresholdsv19(sol, MOI, dMOIThr, tspan, convFac, V0)
%Returns decision thresholds and their ranges, both with and without the
%robustness algorithm.

[m_Cro1, ind_Cro1] = max(sol.m5(:, 6)); %min Cro
[m_CI1, ind_CI1] = max(sol.m1_wt(:, 5)); %min CI
[m_Cro2, ind_Cro2] = max(sol.m1_wt(:, 6)); %max Cro
[m_CI2, ind_CI2] = max(sol.m2_wt(:, 5)); %max CI

KCro_o = (m_Cro2 - m_Cro1)/2 + m_Cro1; %2
KCI_o = (m_CI2 - m_CI1)/2 + m_CI1;
KCroR_o = [m_Cro1, m_Cro2];
KCIR_o = [m_CI1, m_CI2];

MOI = 1:dMOIThr:2; %10
maxKCI = KCIR_o(2);
maxKCro = KCroR_o(2);
options = odeset('Nonnegative', [], 'RelTol', 1e-6, ...
    'AbsTol', 1e-6);
for j = 1:length(MOI)
    %Simulate----------------------------------------------------------
    [t_wt, y_wt] = ode15s(@fv19_repv3, tspan, [0 0 0 0 0 0 MOI(j)*convFac/V0], ...
        options, sol.n, sol.prod, sol.degr, sol.K, sol.tau, V0, convFac);
    if max(y_wt(:, 4)) < KCIR_o(1) && max(y_wt(:, 5)) > KCroR_o(1) %has to be lysis
        if max(y_wt(:, 5))/maxKCro < 1
            maxKCro = max(y_wt(:, 5));
        end
    elseif max(y_wt(:, 5)) < KCroR_o(1) && max(y_wt(:, 4)) > KCIR_o(1)%has to be lysogeny
        if max(y_wt(:, 4))/maxKCI < 1
            maxKCI = max(y_wt(:, 4));
        end
    elseif max(y_wt(:, 4)) > KCIR_o(1) && max(y_wt(:, 5)) > KCroR_o(1) ...
            && max(y_wt(:, 4)) < maxKCI && max(y_wt(:, 5)) < maxKCro
        %both are larger than minimum thresholds
        if max(y_wt(:, 4))/maxKCI > max(y_wt(:, 5))/maxKCro
            maxKCI = max(y_wt(:, 4));
        else
            maxKCro = max(y_wt(:, 5));
        end
    end
end

KCIR = [KCIR_o(1), maxKCI];
KCroR = [KCroR_o(1), maxKCro];
KCI = mean(KCIR);
KCro = mean(KCroR);

end