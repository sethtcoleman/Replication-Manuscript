function [f, g] = getFlux_S(t, y, n, prod, degr, K, tau, V0, convFac)
%Returns the positive and negative fluxes for a given set of trajectories. y includes viral copy
%number.

%Map of ODEs:
%dydt = [   
%           [cI];
%           [cro];
%           [cII];
%           [CI];
%           [Cro];
%           [CII];
%           [lambda];
%       ];

[rcI_PRM, acI_PRM, r_PRE, rCI, rcro, rCro, rcII, rCII, rM] = ...
    deal(prod(1), prod(2), prod(3), prod(4), prod(5), prod(6), prod(7), ...
    prod(8), prod(9));

[kdil, kcI, kCI, kcro, kCro, kcII, kCII, kM] = deal(degr(1), degr(2), ...
    degr(3), degr(4), degr(5), degr(6), degr(7), degr(8));

[nPRM_CIu, nPRM_CId, nPRM_Cro, nPRE, nCro_Cro, nCro_CI, nCII_Cro, nCII_CI, ...
    nM_Cro, nM_CI, nDeg_CII] = ...
    deal(n(1), n(2), n(3), n(4), n(5), n(6), n(7), n(8), n(9), n(10), n(11));

[KPRM_CIu, KPRM_CId, KPRM_Cro, KPRE, KCro_Cro, KCro_CI, KCII_Cro, KCII_CI, ...
    KM_Cro, KM_CI, KDeg_CII] = ...
    deal(K(1), K(2), K(3), K(4), K(5), K(6), K(7), K(8), K(9), K(10), K(11));

%Prefactors----------------------------------------------------------------
PRM_prod_norm1 = 1./...
    (1 + (y(:, 4)./KPRM_CIu).^nPRM_CIu + (y(:, 4)./KPRM_CId).^nPRM_CId + ...
    (y(:, 5)./KPRM_Cro).^nPRM_Cro);
PRM_prod_norm2 = (y(:, 4)./KPRM_CIu).^nPRM_CIu./...
    (1 + (y(:, 4)./KPRM_CIu).^nPRM_CIu + (y(:, 4)./KPRM_CId).^nPRM_CId + ...
    (y(:, 5)./KPRM_Cro).^nPRM_Cro);
PRE_prod_norm = (y(:, 6).^nPRE)./(KPRE.^nPRE + y(:, 6).^nPRE);
Cro_prod_norm = 1./(1 + (y(:, 5)./KCro_Cro).^nCro_Cro + (y(:, 4)./KCro_CI).^nCro_CI);
CII_prod_norm = 1./(1 + (y(:, 5)./KCII_Cro).^nCII_Cro + (y(:, 4)./KCII_CI).^nCII_CI);
CII_deg_norm = (KDeg_CII.^nDeg_CII)./(KDeg_CII.^nDeg_CII + y(:, 6).^nDeg_CII);
rep_prod_norm = heaviSideTrue(t-tau)./(1 + (y(:, 5)./KM_Cro).^nM_Cro + (y(:, 4)./KM_CI).^nM_CI);

f = zeros(length(t), 9);
g = zeros(length(t), 13);
 
%[1] cI
f(:, 1) = y(:, 7).*rcI_PRM.*PRM_prod_norm1; %PRM,basal
f(:, 2) = y(:, 7).*rcI_PRM.*acI_PRM.*PRM_prod_norm2; %PRM,active
f(:, 3) = y(:, 7).*r_PRE.*PRE_prod_norm; %PRE
g(:, 1) = kcI.*y(:, 1); %RNA degr
g(:, 2) = kdil.*y(:, 1); %dil

%[2] cro
f(:, 4) = y(:, 7).*rcro.*Cro_prod_norm; %cro prod
g(:, 3) = kcro.*y(:, 2); %RNA degr
g(:, 4) = kdil.*y(:, 2); %dil

%[3] cII
f(:, 5) = y(:, 7).*rcII.*CII_prod_norm; %cII prod
g(:, 5) = kcII.*y(:, 3); %RNA degr
g(:, 6) = kdil.*y(:, 3); %dil

%[4] CI
f(:, 6) = rCI.*y(:, 1); %CI prod
g(:, 7) = kCI.*y(:, 4); %Cro degr
g(:, 8) = kdil.*y(:, 4); %dil

%[5] Cro
f(:, 7) = rCro.*y(:, 2); %Cro prod
g(:, 9) = kCro.*y(:, 5); %Cro degr
g(:, 10) = kdil.*y(:, 5); %dil

%[6] CII
f(:, 8) = rCII.*y(:, 3); %CII prod
g(:, 11) = kCII.*CII_deg_norm.*y(:, 6); %CII degr
g(:, 12) = kdil.*y(:, 6);

%[7] lambda
f(:, 9) = rM.*y(:, 7).*rep_prod_norm;
g(:, 13) = kdil.*y(:, 7); %dil



end