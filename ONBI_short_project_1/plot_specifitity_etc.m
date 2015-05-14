function plotROC (specificity, sensitivity)
%plot specificity-sensitivity ROC curves
figure
plot(specificity, sensitivity)
title 'ejection fraction (EF) ROC'
xlabel ' specificity'
ylabel ' sensitivity'
end


%% plot accuracy versus threshold
figure
plot(accuracy);
%find threshold that gives greatest accuracy
threshold = find(accuracy == max(accuracy));
hold on
plot([threshold threshold],[0.48 0.66], 'g')
plot([0 100],[max(accuracy)  max(accuracy)], 'g')
trial = 55;
plot([trial trial],[0.48 0.66], 'r')
ylabel 'accuracy'
xlabel 'ejection fraction threshold %'
 
%% plot ejection fraction 
figure
hold on
nbins =30;
histogram(DETERMINE_EF, nbins)
histogram(MESA_EF, nbins)
title 'Ejection fractions'
xlabel 'Ejection fraction (%)'
ylabel 'frequency'
hold on
% print('compare ejection fractions','-dpng')
plot([threshold threshold],[0 15], 'g')
legend 'DETERMINE' 'MESA' 'Threshold'

%% plot stroke volume (SV)
figure
hold on
nbins =30;
histogram(DETERMINE_SV*0.001, nbins)
histogram(MESA_SV*0.001, nbins)
legend 'DETERMINE' 'MESA'
title 'Stroke volumes'
xlabel 'Stroke volume (ml)'
ylabel 'frequency'
% print('compare stroke volumes','-dpng')
