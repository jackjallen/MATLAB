%% plot surface areas
figure;
hold on
nbins = 25;
histogram(DETERMINE_dia_endo_areas,nbins)
histogram(MESA_dia_endo_areas,nbins)
legend 'DETERMINE' ' MESA'
title 'diastolic endocardium surface areas'
xlabel 'surface areas (mm^2)'
ylabel 'frequency'
print('compare diastolic endocardium surface areas','-dpng')

figure;
hold on
nbins = 25;
histogram(DETERMINE_dia_epi_areas,nbins)
histogram(MESA_dia_epi_areas,nbins)
legend 'DETERMINE' ' MESA'
title 'diastolic epicardium surface areas'
xlabel 'surface areas (mm^2)'
ylabel 'frequency'
print('compare diastolic epicardium surface areas','-dpng')

figure;
hold on
nbins = 25;
histogram(DETERMINE_sys_endo_areas,nbins)
histogram(MESA_sys_endo_areas,nbins)
legend 'DETERMINE' ' MESA'
title 'systolic endocardium surface areas'
xlabel 'surface areas (mm^2)'
ylabel 'frequency'
print('compare systolic endocardium surface areas','-dpng')


figure;
hold on
nbins = 25;
histogram(DETERMINE_sys_epi_areas,nbins)
histogram(MESA_sys_epi_areas,nbins)
legend 'DETERMINE' ' MESA'
title 'systolic epicardium surface areas'
xlabel 'surface areas (mm^2)'
ylabel 'frequency'
print('compare systolic epicardium surface areas','-dpng')

figure;
hold on
nbins = 25;
histogram(DETERMINE_dia_myo_areas,nbins)
histogram(MESA_dia_myo_areas,nbins)
legend 'DETERMINE' ' MESA'
title 'diastolic myocardium surface areas'
xlabel 'surface areas (mm^2)'
ylabel 'frequency'
print('compare diastolic myocardium surface areas','-dpng')

figure;
hold on
nbins = 25;
histogram(DETERMINE_sys_myo_areas,nbins)
histogram(MESA_sys_myo_areas,nbins)
legend 'DETERMINE' ' MESA'
title 'systolic myocardium surface areas'
xlabel 'surface areas (mm^2)'
ylabel 'frequency'
print('compare systolic myocardium surface areas','-dpng')

%% calculate surface area to volume ratios
for i = 1:100 %100 cases in each class
    DETERMINE_dia_endo_AVratio(i,1) = DETERMINE_dia_endo_areas(i,1)/DETERMINE_diastolic_endoVolumes(i,1);
    DETERMINE_dia_epi_AVratio(i,1) = DETERMINE_dia_epi_areas(i,1)/DETERMINE_diastolic_epiVolumes(i,1);
    DETERMINE_sys_endo_AVratio(i,1) = DETERMINE_sys_endo_areas(i,1)/DETERMINE_systolic_endoVolumes(i,1);
    DETERMINE_sys_epi_AVratio(i,1) = DETERMINE_sys_epi_areas(i,1)/DETERMINE_systolic_epiVolumes(i,1);
    
    MESA_dia_endo_AVratio(i,1) = MESA_dia_endo_areas(i,1)/MESA_diastolic_endoVolumes(i,1);
    MESA_dia_epi_AVratio(i,1) = MESA_dia_epi_areas(i,1)/MESA_diastolic_epiVolumes(i,1);
    MESA_sys_endo_AVratio(i,1) = MESA_sys_endo_areas(i,1)/MESA_systolic_endoVolumes(i,1);
    MESA_sys_epi_AVratio(i,1) = MESA_sys_epi_areas(i,1)/MESA_systolic_epiVolumes(i,1);
    
    MESA_dia_myo_AVratio(i,1) = MESA_dia_myo_areas(i,1)/MESA_diastolic_myovolumes(i,1);
    MESA_sys_myo_AVratio(i,1) = MESA_sys_myo_areas(i,1)/MESA_systolic_myovolumes(i,1);
    
    DETERMINE_dia_myo_AVratio(i,1) = DETERMINE_dia_myo_areas(i,1)/DETERMINE_diastolic_myovolumes(i,1);
    DETERMINE_sys_myo_AVratio(i,1) = DETERMINE_sys_myo_areas(i,1)/DETERMINE_systolic_myovolumes(i,1);
    
    
end

figure;
hold on
nbins = 25;
histogram(DETERMINE_dia_endo_AVratio,nbins)
histogram(MESA_dia_endo_AVratio,nbins)
legend 'DETERMINE' ' MESA'
title 'diastolic endocardium surface area to volume ratio'
xlabel 'area/volume (mm^-1)'
ylabel 'frequency'
print('compare diastolic endocardium surface areas to volume ratios','-dpng')

figure;
hold on
nbins = 25;
histogram(DETERMINE_sys_epi_AVratio,nbins)
histogram(MESA_sys_epi_AVratio,nbins)
legend 'DETERMINE' ' MESA'
title 'systolic endocardium surface area to volume ratio'
xlabel 'area/volume (mm^-1)'
ylabel 'frequency'
print('compare systolic endocardium surface area to volume ratios','-dpng')

figure;
hold on
nbins = 25;
histogram(DETERMINE_sys_endo_AVratio,nbins)
histogram(MESA_sys_endo_AVratio,nbins)
legend 'DETERMINE' ' MESA'
title 'systolic endocardium surface area to volume ratio'
xlabel 'area/volume (mm^-1)'
ylabel 'frequency'
print('compare systolic endocardium surface area to volume ratios','-dpng')

figure;
hold on
nbins = 25;
histogram(DETERMINE_sys_epi_AVratio,nbins)
histogram(MESA_sys_epi_AVratio,nbins)
legend 'DETERMINE' ' MESA'
title 'systolic epicardium surface area to volume ratio'
xlabel 'area/volume (mm^-1)'
ylabel 'frequency'
print('compare systolic epicardium surface area to volume ratios','-dpng')

figure;
hold on
nbins = 25;
histogram(DETERMINE_dia_myo_AVratio,nbins)
histogram(MESA_dia_myo_AVratio,nbins)
legend 'DETERMINE' ' MESA'
title 'diastolic myocardium surface area to volume ratio'
xlabel 'area/volume (mm^-1)'
ylabel 'frequency'
print('compare diastolic myocardium surface areas to volume ratios','-dpng')
