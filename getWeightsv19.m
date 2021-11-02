function w = getWeightsv19(t, y, n, prod, degr, K, tau, V0, convFac)
%Returns the weights for the different regulatory terms

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

%PRM-----------------------------------------------------------------------
PRM_norm = (1 + (y(:, 4)./KPRM_CIu).^nPRM_CIu + (y(:, 4)./KPRM_CId).^nPRM_CId + ...
    (y(:, 5)./KPRM_Cro).^nPRM_Cro);
w(:, 1) = 1./PRM_norm; %basal
w(:, 2) = (y(:, 4)./KPRM_CIu).^nPRM_CIu./PRM_norm; %active
w(:, 3) = (y(:, 4)./KPRM_CId).^nPRM_CId./PRM_norm; %CI repr.
w(:, 4) = (y(:, 5)./KPRM_Cro).^nPRM_Cro./PRM_norm; %Cro repr.

%PRE-----------------------------------------------------------------------
PRE_norm = (1 + (y(:, 6)./KPRE).^nPRE);
w(:, 5) = 1./PRE_norm; %basal
w(:, 6) = (y(:, 6)./KPRE).^nPRE./PRE_norm; %active

%PCro----------------------------------------------------------------------
Cro_norm = (1 + (y(:, 5)./KCro_Cro).^nCro_Cro + (y(:, 4)./KCro_CI).^nCro_CI);
w(:, 7) = 1./Cro_norm; %basal
w(:, 8) = (y(:, 5)./KCro_Cro).^nCro_Cro./Cro_norm; %Cro repr.
w(:, 9) = (y(:, 4)./KCro_CI).^nCro_CI./Cro_norm; %CI repr.

%PCII----------------------------------------------------------------------
CII_norm = (1 + (y(:, 5)./KCII_Cro).^nCII_Cro + (y(:, 4)./KCII_CI).^nCII_CI);
w(:, 10) = 1./CII_norm; %basal
w(:, 11) = (y(:, 5)./KCII_Cro).^nCII_Cro./CII_norm; %Cro repr.
w(:, 12) = (y(:, 4)./KCII_CI).^nCII_CI./CII_norm; %CI repr.

%CII degr------------------------------------------------------------------
CII_deg_norm = (1 + (y(:, 6)./KDeg_CII).^nDeg_CII);
w(:, 13) = 1./CII_deg_norm; %basal
w(:, 14) = (y(:, 6)./KDeg_CII).^nDeg_CII./CII_deg_norm; %CII repr.

%Repl----------------------------------------------------------------------
rep_norm = (1 + (y(:, 5)./KM_Cro).^nM_Cro + (y(:, 4)./KM_CI).^nM_CI);
w(:, 15) = heaviSideTrue(t - tau)./rep_norm; %unrepr.
w(:, 16) = heaviSideTrue(t - tau).*(y(:, 5)./KM_Cro).^nM_Cro./rep_norm; %Cro repr.
w(:, 17) = heaviSideTrue(t - tau).*(y(:, 4)./KM_CI).^nM_CI./rep_norm; %CI repr.
w(:, 18) = 1./rep_norm; %unrepr.
w(:, 19) = (y(:, 5)./KM_Cro).^nM_Cro./rep_norm; %Cro repr.
w(:, 20) = (y(:, 4)./KM_CI).^nM_CI./rep_norm; %CI repr.


end