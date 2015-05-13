function plotMyoVolumes(DETERMINE_systolic_myovolumes, MESA_systolic_myovolumes, DETERMINE_diastolic_myovolumes, MESA_diastolic_myovolumes)

figure
hold on
nbins = 25;
histogram(DETERMINE_systolic_myovolumes*0.001,nbins)
histogram(MESA_systolic_myovolumes*0.001,nbins)
legend 'DETERMINE' 'MESA'
title 'systolic myocardium volumes'
xlabel 'myocardium volume (ml)'
ylabel 'frequency'
print('compare systolic myocardium volumes','-dpng')

% plot diastolic myocardium volumes
figure
hold on
nbins = 25;
histogram(DETERMINE_diastolic_myovolumes*0.001,nbins)
histogram(MESA_diastolic_myovolumes*0.001,nbins)
legend ' DETERMINE ' ' MESA '
title 'diastolic myocardium volumes'
xlabel 'myocardium volume (ml)'
ylabel 'frequency'
print('compare diastolic myocardium volumes','-dpng')

end