function[covariance_matrix] = calcCovarianceMatrix(data)
% reshape shape coordinate matrices (to vector form) and find average shape
%JA 27/04/2015
sum_of_shapes = zeros(size(data(1).diastolic.endo.xyz(:),1),1);
for i = 1:400
    tmp = data(1).diastolic.endo.xyz(:);
    sum_of_shapes = tmp + sum_of_shapes;
    mean_shape =  sum_of_shapes/400;
end
% Covariance matrix
sum = zeros(size(data(1).diastolic.endo.xyz(:),1));
for n = 1:400
    part1 = (data(1).diastolic.endo.xyz(:) - mean_shape);
    part2 = part1.';
    sum = (part1*part2) + sum;
end
factor = (1/(400-1));
covariance_matrix = factor * sum;
end