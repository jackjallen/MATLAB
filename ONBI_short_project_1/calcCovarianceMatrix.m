% function[dia_endo_covariance_matrix, dia_endo_mean_shape, sys_endo_covariance_matrix, sys_endo_mean_shape, sys_myo_covariance_matrix] = calcCovarianceMatrix(data, indices, sys_myo_mean)
function[sys_myo_covariance_matrix] = calcCovarianceMatrix(data, indices, sys_myo_mean)

% reshape shape coordinate matrices (to vector form) and find average shape
%JA 27/04/2015

%% find the mean shape
% dia_endo_sum_of_shapes = zeros(size(data(2).diastolic.endo.xyz(:),1),1);
% sys_endo_sum_of_shapes = zeros(size(data(2).systolic.endo.xyz(:),1),1);
% sys_myo_sum_of_shapes = zeros(size(data(2).systolic.myo.xyz(:),1),1);
% for i = indices'
%     dia_endo_tmp = data(i).diastolic.endo.xyz(:);
%     dia_endo_sum_of_shapes = dia_endo_tmp + dia_endo_sum_of_shapes;
%     dia_endo_mean_shape =  dia_endo_sum_of_shapes/400;
%     
%     sys_endo_tmp = data(i).systolic.endo.xyz(:);
%     sys_endo_sum_of_shapes = sys_endo_tmp + sys_endo_sum_of_shapes;
%     sys_endo_mean_shape =  sys_endo_sum_of_shapes/400;
%     
%     i
% end


%% find the Covariance matrix
% dia_endo_sum = zeros(size(data(2).diastolic.endo.xyz(:),1));
% sys_endo_sum = zeros(size(data(2).systolic.endo.xyz(:),1));
sys_myo_sum = zeros(size(data(2).systolic.myo.xyz(:),1));
for n = sort(indices')
    n
    
%     dia_endo_part1 = (data(n).diastolic.endo.xyz - dia_endo_mean_shape);
%     dia_endo_part2 = dia_endo_part1.';
%     dia_endo_sum = (dia_endo_part1*dia_endo_part2) + dia_endo_sum;
%     
%     sys_endo_part1 = (data(n).systolic.endo.xyz - sys_endo_mean_shape);
%     sys_endo_part2 = sys_endo_part1.';
%     sys_endo_sum = (sys_endo_part1*sys_endo_part2) + sys_endo_sum;
%     
    sys_myo_part1 = (data(n).systolic.myo.xyz(:) - sys_myo_mean(:));
    sys_myo_part2 = sys_myo_part1.';
    sys_myo_sum = (sys_myo_part1*sys_myo_part2) + sys_myo_sum;
    
end
factor = 1/( (size(indices,1)-1) );
% dia_endo_covariance_matrix = factor * dia_endo_sum;
% sys_endo_covariance_matrix = factor *sys_endo_sum;
sys_myo_covariance_matrix = factor * sys_myo_sum;
end




