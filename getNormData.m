function ty_data_norm = getNormData(ty_data)
%Takes TY's RNA data and normalizes, dividing by the maximum for each
%RNA (cI, cro, cII). 

ty_data_norm = ty_data;
ty_data_norm(:, [2, 5, 8]) = ty_data_norm(:, [2, 5, 8])./max(ty_data_norm(:, 2));
ty_data_norm(:, [3, 6, 9]) = ty_data_norm(:, [3, 6, 9])./max(ty_data_norm(:, 3));
ty_data_norm(:, [4, 7, 10]) = ty_data_norm(:, [4, 7, 10])./max(ty_data_norm(:, 4));

end