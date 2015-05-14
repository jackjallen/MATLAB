%% visualise eigenmodes of shape variation - meshes

dia_endo_new_shape_min = reshape(dia_endo_mean_shape + dia_endo_principle_eigenvector.*dia_endo_min_b, [1089 3]);
dia_endo_new_shape_max = reshape(dia_endo_mean_shape + dia_endo_principle_eigenvector.*dia_endo_max_b, [1089 3]);
dia_endo_mean_shape = reshape(dia_endo_mean_shape, [1089 3]);
figure 
hold on
plot3(dia_endo_mean_shape(:,1),dia_endo_mean_shape(:,2),dia_endo_mean_shape(:,3),'g.')
plot3(dia_endo_new_shape_max(:,1),dia_endo_new_shape_max(:,2),dia_endo_new_shape_max(:,3),'ro')
plot3(dia_endo_new_shape_min(:,1),dia_endo_new_shape_min(:,2),dia_endo_new_shape_min(:,3),'bo')
title 'dia endo eigenmode variation'
legend 'mean' 'max b' 'min b'

