function [cost, err_mat] = getCostThuv2(tm1, ym1, PCR, V, convFac, ...
    i_taus, normWeight)
%Calculates the error in Thu's <MOI> ~ 1 qPCR experiment. Uses log space,
%assumes data is already transformed.

err_mat = zeros(1, length(i_taus));

%Find time indices that correspond to data.
taus = zeros(length(i_taus), 1);
taus_index = 1;
for j = PCR(1).time(i_taus)
    [~, taus(taus_index)] = min(abs(tm1-j));
    taus_index = taus_index +1;
end

scaled_index = 1;
for j = transpose(taus)
    err_mat(scaled_index) = ...
        normWeight(scaled_index)*( log(ym1(j, end)*V(j)/convFac) - ...
        PCR(1).mean(i_taus(scaled_index)) );
 
    scaled_index = scaled_index + 1;
end

% %Rescale by number of observations
% err_mat = err_mat./length(taus);

cost = sum(sum(err_mat.^2));

end