%% visualise modes of variation - static
% close all
visualisingMode1 = 1
visualisingMode2 = 2
nStd = 3 %number of standard deviations
x = [-60 ; 40];
y = [-65 ;50];
z = [-110 ; 20];
figure('name','diastolic modes of variation')
run('visualisemodes_diastole.m')
run('visualisemodes_diastole__combined_modes.m')
figure('name','systolic modes of variation')
run('visualisemodes_systole.m')
run('visualisemodes_systole__combined_modes.m')

% visualiseModes(data, dia_sys_myo_mean, principle_dia_sys_myo_eigenvectors, dia_sys_myo_max_b, 'sys', eMode, c,x,y,z)
