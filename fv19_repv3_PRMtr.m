function dydt = fv19_repv3_PRMtr(t, y, n, prod, degr, K, tau, V0, convFac)
%Returns differential equations at point (t, y), WT. Splits cI mRNA into two
%equations with different translation rates

%Initialize ODEs
dydt = zeros(8, 1);

%Map of ODEs:
%dydt = [   
%           [cI_PRM];
%           [cI_PRE];
%           [cro];
%           [cII];
%           [CI];
%           [Cro];
%           [CII];
%           [lambda];
%       ];

[rcI_PRM, acI_PRM, acI_PRE_PRM, rCI_PRM, rCI_PRE, rcro, rCro, rcII, rCII, rM] = ...
    deal(prod(1), prod(2), prod(3), prod(4), prod(5), prod(6), prod(7), ...
    prod(8), prod(9), prod(10));

[kdil, kcI, kCI, kcro, kCro, kcII, kCII, kM] = deal(degr(1), degr(2), ...
    degr(3), degr(4), degr(5), degr(6), degr(7), degr(8));

[nPRM_CIu, nPRM_CId, nPRM_Cro, nPRE, nCro_Cro, nCro_CI, nCII_Cro, nCII_CI, ...
    nM_Cro, nM_CI, nDeg_CII] = ...
    deal(n(1), n(2), n(3), n(4), n(5), n(6), n(7), n(8), n(9), n(10), n(11));

[KPRM_CIu, KPRM_CId, KPRM_Cro, KPRE, KCro_Cro, KCro_CI, KCII_Cro, KCII_CI, ...
    KM_Cro, KM_CI, KDeg_CII] = ...
    deal(K(1), K(2), K(3), K(4), K(5), K(6), K(7), K(8), K(9), K(10), K(11));

%Prefactors----------------------------------------------------------------
PRM_prod_norm1 = 1/...
    (1 + (y(5)/KPRM_CIu)^nPRM_CIu + (y(5)/KPRM_CId)^nPRM_CId + ...
    (y(6)/KPRM_Cro)^nPRM_Cro);
PRM_prod_norm2 = (y(5)/KPRM_CIu)^nPRM_CIu/...
    (1 + (y(5)/KPRM_CIu)^nPRM_CIu + (y(5)/KPRM_CId)^nPRM_CId + ...
    (y(6)/KPRM_Cro)^nPRM_Cro);
PRE_prod_norm = (y(7)^nPRE)/(KPRE^nPRE + y(7)^nPRE);
Cro_prod_norm = 1/(1 + (y(6)/KCro_Cro)^nCro_Cro + (y(5)/KCro_CI)^nCro_CI);
CII_prod_norm = 1/(1 + (y(6)/KCII_Cro)^nCII_Cro + (y(5)/KCII_CI)^nCII_CI);
CII_deg_norm = (KDeg_CII^nDeg_CII)/(KDeg_CII^nDeg_CII + y(7)^nDeg_CII);
rep_prod_norm = heaviSideTrue(t-tau)/(1 + (y(6)/KM_Cro)^nM_Cro + (y(5)/KM_CI)^nM_CI);
 
%[1] cI PRM
dydt(1) = y(8)*rcI_PRM*PRM_prod_norm1 + y(7)*rcI_PRM*acI_PRM*PRM_prod_norm2 - ...
    kcI*y(1) - kdil*y(1);

%[2] cI PRE
dydt(2) = y(8)*acI_PRE_PRM*rcI_PRM*PRE_prod_norm - ...
    kcI*y(2) - kdil*y(2);

%[3] cro
dydt(3) = y(8)*rcro*Cro_prod_norm - kcro*y(3) - kdil*y(3);

%[4] cII
dydt(4) = y(8)*rcII*CII_prod_norm - kcII*y(4) - kdil*y(4);

%[5] CI
dydt(5) = rCI_PRM*y(1) + rCI_PRE*y(2) - kCI*y(5) - kdil*y(5);

%[5] Cro
dydt(6) = rCro*y(3) - kCro*y(6) - kdil*y(6);

%[6] CII
dydt(7) = rCII*y(4) - kCII*CII_deg_norm*y(7) - kdil*y(7);

%[7] lambda
dydt(8) = rM*y(8)*rep_prod_norm - kdil*y(8) - kM*y(8);

end