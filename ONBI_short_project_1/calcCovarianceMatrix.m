function[dia_endo_covariance_matrix, dia_endo_mean_shape] = calcCovarianceMatrix(data)
% reshape shape coordinate matrices (to vector form) and find average shape
%JA 27/04/2015

%% find the mean shape
dia_endo_sum_of_shapes = zeros(size(data(1).diastolic.endo.xyz(:),1),1);
for i = 1:400
    dia_endo_tmp = data(1).diastolic.endo.xyz(:);
    dia_endo_sum_of_shapes = dia_endo_tmp + dia_endo_sum_of_shapes;
    dia_endo_mean_shape =  dia_endo_sum_of_shapes/400;
end
%% find the Covariance matrix
dia_endo_sum = zeros(size(data(1).diastolic.endo.xyz(:),1));
for n = 1:400
    dia_endo_part1 = (data(1).diastolic.endo.xyz(:) - dia_endo_mean_shape);
    dia_endo_part2 = dia_endo_part1.';
    dia_endo_sum = (dia_endo_part1*dia_endo_part2) + dia_endo_sum;
end
factor = (1/(400-1));
dia_endo_covariance_matrix = factor * dia_endo_sum;
end