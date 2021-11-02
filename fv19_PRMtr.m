function dydt = fv19_PRMtr(t, y, n, prod, degr, K, MOI, V0, convFac, ss)
%Returns differential equations at point y. Splits cI mRNA into two
%species, PRM first, with PRM having a different translational rate

if nargin < 10
    ss = 0;
end

%Initialize ODEs
dydt = zeros(7, 1);

%Map of ODEs:
%dydt = [   
%           [cI_PRM];
%           [cI_PRE];
%           [cro];
%           [cII];
%           [CI];
%           [Cro];
%           [CII];
%       ];

[rcI_PRM, acI_PRM, acI_PRE_PRM, rCI_PRM, rCI_PRE, rcro, rCro, rcII, rCII] = deal(prod(1),...
    prod(2), prod(3), prod(4), prod(5), prod(6), prod(7), prod(8), prod(9));

[kdil, kcI, kCI, kcro, kCro, kcII, kCII, kM] = deal(degr(1), degr(2), ...
    degr(3), degr(4), degr(5), degr(6), degr(7), degr(8));

[nPRM_CIu, nPRM_CId, nPRM_Cro, nPRE, nCro_Cro, nCro_CI, nCII_Cro, nCII_CI, ...
    nDeg_CII] = ...
    deal(n(1), n(2), n(3), n(4), n(5), n(6), n(7), n(8), n(9));

[KPRM_CIu, KPRM_CId, KPRM_Cro, KPRE, KCro_Cro, KCro_CI, KCII_Cro, KCII_CI, ...
    KDeg_CII] = ...
    deal(K(1), K(2), K(3), K(4), K(5), K(6), K(7), K(8), K(9));

%Viral copy number
V = V0*exp(kdil*t);
if ss ~= 1
    lambda = MOI*convFac/V;
elseif ss == 1
    lambda = MOI*convFac/V0;
end

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
 
%[1] cI (PRM)
dydt(1) = lambda*rcI_PRM*PRM_prod_norm1 + lambda*rcI_PRM*acI_PRM*PRM_prod_norm2 - ...
    kcI*y(1) - kdil*y(1);

%[2] cI (PRE)
dydt(2) = lambda*acI_PRE_PRM*rcI_PRM*PRE_prod_norm - ...
    kcI*y(2) - kdil*y(2);

%[3] cro
dydt(3) = lambda*rcro*Cro_prod_norm - kcro*y(3) - kdil*y(3);

%[4] cII
dydt(4) = lambda*rcII*CII_prod_norm - kcII*y(4) - kdil*y(4);

%[5] CI
dydt(5) = rCI_PRM*y(1) + rCI_PRE*y(2) - kCI*y(5) - kdil*y(5);

%[6] Cro
dydt(6) = rCro*y(3) - kCro*y(6) - kdil*y(6);

%[7] CII
dydt(7) = rCII*y(4) - kCII*CII_deg_norm*y(7) - kdil*y(7);

if lambda < 0
   error('\lambda < 0!'); 
end

end