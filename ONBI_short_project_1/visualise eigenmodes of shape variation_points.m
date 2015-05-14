%% visualise eigenmodes of shape variation - points
% dia_endo_new_shape_min = reshape(dia_endo_new_shape_min(:,eigenmode), [1089 3]);
% dia_endo_new_shape_max = reshape(dia_endo_new_shape_max(:,eigenmode), [1089 3]);
% dia_endo_mean_shape = reshape(dia_endo_mean_shape(:,eigenmode), [1089 3]);
% sys_endo_new_shape_min = reshape(sys_endo_new_shape_min(:,eigenmode), [1089 3]);
% sys_endo_new_shape_max = reshape(sys_endo_new_shape_max(:,eigenmode), [1089 3]);
% sys_endo_mean_shape = reshape(sys_endo_mean_shape(:,eigenmode), [1089 3]);

sys_myo_new_shape_1_min = reshape(sys_myo_new_shape_min(:,1), [2178 3]);
sys_myo_new_shape_1_max = reshape(sys_myo_new_shape_max(:,1), [2178 3]);
sys_myo_new_shape_2_min = reshape(sys_myo_new_shape_min(:,2), [2178 3]);
sys_myo_new_shape_2_max = reshape(sys_myo_new_shape_max(:,2), [2178 3]);

figure
hold on
title 'PCA - eigenmode visualisation - diastolic endocardium'
plot3(dia_endo_mean_shape(:,1),dia_endo_mean_shape(:,2),dia_endo_mean_shape(:,3),'g.')
plot3(dia_endo_new_shape_max(:,1),dia_endo_new_shape_max(:,2),dia_endo_new_shape_max(:,3),'ro')
plot3(dia_endo_new_shape_min(:,1),dia_endo_new_shape_min(:,2),dia_endo_new_shape_min(:,3),'bo')
legend 'mean' 'max b' 'min b'

figure
title 'PCA - eigenmode visualisation - systolic endocardium'
hold on
plot3(sys_endo_mean_shape(:,1),sys_endo_mean_shape(:,2),sys_endo_mean_shape(:,3),'g.')
plot3(sys_endo_new_shape_max(:,1),sys_endo_new_shape_max(:,2),sys_endo_new_shape_max(:,3),'ro')
plot3(sys_endo_new_shape_min(:,1),sys_endo_new_shape_min(:,2),sys_endo_new_shape_min(:,3),'bo')
legend 'mean' 'max b' 'min b'

figure
title 'PCA - eigenmode visualisation - systolic myocardium'
hold on
plot3(sys_myo_mean(:,1),sys_myo_mean(:,2),sys_myo_mean(:,3),'g.')
plot3(sys_myo_new_shape_2_max(:,1),sys_myo_new_shape_2_max(:,2),sys_myo_new_shape_2_max(:,3),'ro')
plot3(sys_myo_new_shape_2_min(:,1),sys_myo_new_shape_2_min(:,2),sys_myo_new_shape_2_min(:,3),'bo')
legend 'mean' 'max b' 'min b'

figure
for c = -1:0.1:1
plot3D(sys_myo_mean)
hold on
 c   

sys_myo_new_shape(:,1) = sys_myo_mean(:) + principle_sys_myo_eigenvectors(:,1)*c*sys_myo_max_b(1,1);
sys_myo_new_shape(:,2) = sys_myo_mean(:) + principle_sys_myo_eigenvectors(:,2)*c*sys_myo_max_b(2,1);
% 
 plot3D(reshape(sys_myo_new_shape(:,1), [2178 3]))
plot3D(reshape(sys_myo_new_shape(:,2), [2178 3]))


pause
axis equal
hold off
end
