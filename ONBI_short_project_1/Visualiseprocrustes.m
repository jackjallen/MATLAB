%% Visualise procrustes

disp('visualise procrustes output')

plot3D(dia_myo_reference)
plot3D(sys_myo_reference)

% plotMESH = @(M,varargin)patch('vertices',M.xyz,'faces',M.tri,'edgecolor','k','facecolor','b',varargin{:});

%reference endocardium
%diastolic
figure
plot3(dia_endo_reference(:,1),dia_endo_reference(:,2),dia_endo_reference(:,3),'+')
axis equal; axis tight;
title 'reference endocardium (diastolic)'
%systolic
figure
plot3(sys_endo_reference(:,1),sys_endo_reference(:,2),sys_endo_reference(:,3),'+')
axis equal; axis tight;
title 'reference endocardium (diastolic)'

%reference epicardium
%distolic
figure
hold on
plot3(dia_epi_reference(:,1),dia_epi_reference(:,2),dia_epi_reference(:,3),'+')
axis equal; axis tight;
title 'reference epicardium (diastolic)'
%systolic
figure
plot3(sys_epi_reference(:,1),sys_epi_reference(:,2),sys_epi_reference(:,3),'+')
axis equal; axis tight;
title 'reference epicardium (systolic)'

%reference (whole LV shape from patient 1)
figure
plot3(dia_myo_reference(:,1), dia_myo_reference(:,2), dia_myo_reference(:,3),'o');
title 'reference LV'

%patient i endocardium and epicardium
%diastolic
figure
hold on
plot3(data(i).diastolic.endo.xyz(:,1), data(i).diastolic.endo.xyz(:,2), data(i).diastolic.endo.xyz(:,3),'o'); 
plot3(data(i).diastolic.epi.xyz(:,1), data(i).diastolic.epi.xyz(:,2), data(i).diastolic.epi.xyz(:,3),'o'); 
title 'patient i endo and epi, diastolic'
%systolic
figure
hold on
plot3(data(i).systolic.endo.xyz(:,1), data(i).systolic.endo.xyz(:,2), data(i).systolic.endo.xyz(:,3),'o'); 
plot3(data(i).systolic.epi.xyz(:,1), data(i).systolic.epi.xyz(:,2), data(i).systolic.epi.xyz(:,3),'o');
title 'patient i endo and epi, systolic'

%patient i whole LV shape
%diastolic
figure
hold on
plot3(data(i).diastolic.myo.xyz(:,1), data(i).diastolic.myo.xyz(:,2), data(i).diastolic.myo.xyz(:,3),'o'); 
title 'patient i, whole LV, diastolic'
%systolic

figure
for i = sort(data(1).MESA_indices')
    i
subplot 121
hold on
% plot3(dia_myo_reference(:,1), dia_myo_reference(:,2), dia_myo_reference(:,3),'go');
plot3(data(i).systolic.endo.xyz(:,1), data(i).systolic.endo.xyz(:,2), data(i).systolic.endo.xyz(:,3),'.'); 
plot3(data(i).systolic.epi.xyz(:,1), data(i).systolic.epi.xyz(:,2), data(i).systolic.epi.xyz(:,3),'.'); 
% plot3(data(i).systolic.myo.xyz(:,1), data(i).systolic.myo.xyz(:,2), data(i).systolic.myo.xyz(:,3),'.'); 

% plot3D(systolic_myo_reshaped(i).xyz)

title 'patient i, whole LV, systolic'
subplot 122
hold on
% plot3(dia_myo_reference(:,1), dia_myo_reference(:,2), dia_myo_reference(:,3),'go');
 plot3(data(i).diastolic.endo.xyz(:,1), data(i).diastolic.endo.xyz(:,2), data(i).diastolic.endo.xyz(:,3),'.'); 
 plot3(data(i).diastolic.epi.xyz(:,1), data(i).diastolic.epi.xyz(:,2), data(i).diastolic.epi.xyz(:,3),'.');
% plot3(data(i).diastolic.myo.xyz(:,1), data(i).diastolic.myo.xyz(:,2), data(i).diastolic.myo.xyz(:,3),'.'); 

% plot3D(diastolic_myo_reshaped(i).xyz)

title 'patient i, whole LV, systolic'
pause;

end
