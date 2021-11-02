function PD = getPD(sol, KCI, KCro)
%Computes the phase diagram for a single parameter set, given KCI, KCro, 
%no distinction in mixed outcome types. Note that generally sol is already
%in units of X/KX, so KCI, KCro can be passed as rescaling constants

%Phase Diagram stuff
tauFail = 1e3;
rlam_PD = sol.rlam_PD;
MOI_PD = sol.MOI_PD;
PD_CI = zeros(length(rlam_PD), length(MOI_PD));
PD_Cro = zeros(length(rlam_PD), length(MOI_PD));
PD = zeros(length(rlam_PD), length(MOI_PD));
for j = 1:length(rlam_PD)
    for k = 1:length(MOI_PD)
        yPD = zeros(1, 4); %[CI]/KCI, [Cro]/KCro, tauCI, tauCro
        yPD(1, 1) = sol.PDMat(j, k, 1)/KCI;
        yPD(1, 2) = sol.PDMat(j, k, 2)/KCro;
        yPD(1, end-1) = sol.tauDec_PD(j, k, 1);
        yPD(1, end) = sol.tauDec_PD(j, k, 2);
        PD_CI(j, k) = yPD(:, 1); 
        PD_Cro(j, k) = yPD(:, 2); 
        %CI Cro
        if PD_CI(j, k) < 1 && PD_Cro(j, k) < 1
            PD(j, k) = 3; %fail
        elseif PD_CI(j, k) >= 1 && PD_Cro(j, k) < 1
            PD(j, k) = 1; %lys.
        elseif PD_CI(j, k) < 1 && PD_Cro(j, k) >= 1
            PD(j, k) = 2; %lyt.
        elseif PD_CI(j, k) >= 1 && PD_Cro(j, k) >= 1 
            PD(j, k) = 4; %mixed outcome
        end
    end
end



end