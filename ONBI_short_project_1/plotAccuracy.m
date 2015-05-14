function plotAccuracy(accuracy)
% plot accuracy versus threshold
figure
plot(accuracy);
%find threshold that gives greatest accuracy
threshold = find(accuracy == max(accuracy));
threshold = mean(threshold)
hold on
plot([threshold threshold],[0.5 1], 'g')
plot([0 100],[max(accuracy)  max(accuracy)], 'g')
ylabel 'accuracy'
xlabel 'ejection fraction threshold (%)'
 

end