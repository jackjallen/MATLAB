function plotAccuracyAV(data, accuracy, positive, negative, positiveName, negativeName, AccXMax)
% plot accuracy versus threshold
positiveName
negativeName
figure
subplot 121
plot(accuracy(1,1:AccXMax));
max_accuracy = max(accuracy)
%find threshold that gives greatest accuracy
threshold = find(accuracy == max(accuracy));
threshold = threshold(1) %incase there are two maximum accuracies
hold on
plot([threshold threshold],[0.5 1], 'g')
plot([0 AccXMax],[max(accuracy)  max(accuracy)], 'g')
ylabel 'accuracy'
% xlabel 'surface area to volume ratio (AVR) '

subplot 122
nbins = 30;
hold on
% histogram(cell2mat({data(data(1).MESA_indices).MESA_sys_endo_AVratio})*100,nbins)
% histogram(cell2mat({data(data(1).DETERMINE_indices).DETERMINE_sys_endo_AVratio})*100,nbins)
histogram(positive*100,nbins)
histogram(negative*100,nbins)
% xlabel 'surface area to volume ratio (AVR) (%)'
% ylabel 'frequency'
% legend 'positive' 'negative'
legend positive negative
% subplot 122
% 
% hold on
% nBins = 30;
% histogram(cell2mat({data((data(1).MESA_indices)).ejectionFraction}),nBins)
% histogram(cell2mat({data((data(1).DETERMINE_indices)).ejectionFraction}),nBins)
% plot([threshold threshold], [0 14], 'r')

% xlabel 'Ejection Fraction (%)'
% ylabel 'frequency'
% legend ('MESA', 'DETERMINE', 'threshold', 'Location', 'NorthWest')

end