%% Bullseye plots
diaPts = [1:2178];
sysPts = [2179:4356];


% figure('name','myocardium thickness bullseye plot: mean , DETERMINE, diastole')
figure ('name','myocardium thickness bullseye plots')
subplot 221
colormap jet
axis tight
axis square
tmpShape = reshape(mean(all_training_dia_sys_myo_shapes(data(1).DETERMINE_indices,:),1), [2*2178 , 3]);
shape = tmpShape(diaPts,:);
run('preProcessData_BullsEyePlots_NewShape.m');
title 'mean , DETERMINE, diastole'
% figure('name','myocardium thickness bullseye plot: mean , DETERMINE, systole')
subplot 222
colormap jet
axis tight
axis square
tmpShape = reshape(mean(all_training_dia_sys_myo_shapes(data(1).DETERMINE_indices,:),1), [2*2178 , 3]);
shape = tmpShape(sysPts,:);
run('preProcessData_BullsEyePlots_NewShape.m');
title 'mean , DETERMINE, systole'
%
% figure('name','myocardium thickness bullseye plot: mean , MESA, diastole')
subplot 223
colormap jet
axis tight
axis square
tmpShape = reshape(mean(all_training_dia_sys_myo_shapes(data(1).MESA_indices,:),1), [2*2178 , 3]);
shape = tmpShape(diaPts,:);
% run('preProcessData_BullsEyePlots_NewShape.m');
run('preProcessData_BullsEyePlots_NewShape__three_sections.m')
title 'mean , MESA, diastole'
% figure('name','myocardium thickness bullseye plot: mean , MESA, systole')
subplot 224
colormap jet
axis tight
axis square
tmpShape = reshape(mean(all_training_dia_sys_myo_shapes(data(1).MESA_indices,:),1), [2*2178 , 3]);
shape = tmpShape(sysPts,:);
% run('preProcessData_BullsEyePlots_NewShape.m');
run('preProcessData_BullsEyePlots_NewShape__three_sections.m')
title 'mean , MESA, systole'



figure('name','myocardium thickness bullseye plot: mode of variation , -/+*c*std')
colormap jet
subplot 121
axis tight
axis square
tmpShape = reshape(dia_sys_myo_new_shape_minus(:,visualisingMode1), [2*2178 , 3]);
shape = tmpShape(1:2178,:);
run('preProcessData_BullsEyePlots_NewShape.m');
subplot 122
axis tight
axis square
tmpShape = reshape(dia_sys_myo_new_shape_plus(:,visualisingMode1), [2*2178 , 3]);
shape = tmpShape(1:2178,:);
run('preProcessData_BullsEyePlots_NewShape.m')

figure('Name','Myocardium thickness bullseye plot - Mean') %,'NumberTitle','off')
subplot 221
colormap jet
axis tight
axis square
shape = reshape(dia_sys_myo_mean(:), [2*2178 , 3]);
shape = shape(1:2178,:);
run('preProcessData_BullsEyePlots_NewShape.m');
title 'diastole mean'
subplot 222
shape = reshape(dia_sys_myo_mean(:), [2*2178 , 3]);
shape = shape(2179:end,:);
run('preProcessData_BullsEyePlots_NewShape__three_sections.m')
title 'systole, mean'
subplot 223
shape = reshape(dia_sys_myo_mean(:), [2*2178 , 3]);
shape = shape(1:2178,:);
run('preProcessData_BullsEyePlots_NewShape__three_sections.m')
title 'diastole, mean'