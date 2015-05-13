
figure;
hold on
nbins = 25;
histogram(DETERMINE_myoSV*0.001,nbins)
histogram(MESA_myoSV*0.001,nbins)
legend 'DETERMINE' ' MESA'
title 'stroke volumes calculated using myocardium volumes'
xlabel ' "Stroke volume" ml '
ylabel 'frequency'
% print('compare stroke volumes calculated from myocardium volumes','-dpng')

%plot myocardium "ejection fractions"
figure;
hold on
nbins = 25;
histogram(DETERMINE_myoEF,nbins)
histogram(MESA_myoEF,nbins)
legend 'DETERMINE' ' MESA'
title 'ejection fractions calculated using myocardium volumes'
xlabel '"Ejection fraction" (%)'
ylabel 'frequency'
% print('compare ejection fractions calculated from myocardium volumes','-dpng')