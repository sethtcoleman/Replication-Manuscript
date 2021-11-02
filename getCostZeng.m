function [cost, err_mat] = getCostZeng(tm1, ym1, data, V, convFac, ...
    i_taus, normWeight)
%Calculates the error in Shao et al., iScience 2019 data

err_mat = zeros(1, length(i_taus));

%Find time indices that correspond to data.
taus = zeros(length(i_taus), 1);
taus_index = 1;
for j = data(i_taus, 1)'
    [~, taus(taus_index)] = min(abs(tm1-j));
    taus_index = taus_index +1;
end

scaled_index = 1;
for j = transpose(taus)
    err_mat(scaled_index) = ...
        normWeight(scaled_index)*( ym1(j, 1)*V(j)/convFac - ...
        data(i_taus(scaled_index), 2) );
 
    scaled_index = scaled_index + 1;
end

% %Rescale by number of observations
% err_mat = err_mat./length(taus);

cost = sum(sum(err_mat.^2));

end