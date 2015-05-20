function plotMyoVolumes(data)

figure
hold on
nbins = 25;

histogram(cell2mat({data(data(1).DETERMINE_indices').sys_myovolumes})*0.001,nbins)
histogram(cell2mat({data(data(1).MESA_indices').sys_myovolumes})*0.001,nbins)
legend 'DETERMINE' 'MESA'
title 'systolic myocardium volumes'
xlabel 'myocardium volume (ml)'
ylabel 'frequency'
print('compare systolic myocardium volumes','-dpng')

% plot diastolic myocardium volumes
figure
hold on
nbins = 25;
histogram(cell2mat({data(data(1).DETERMINE_indices').dia_myovolumes})*0.001,nbins)
histogram(cell2mat({data(data(1).MESA_indices').dia_myovolumes})*0.001,nbins)
legend ' DETERMINE ' ' MESA '
title 'diastolic myocardium volumes'
xlabel 'myocardium volume (ml)'
ylabel 'frequency'
print('compare diastolic myocardium volumes','-dpng')

end