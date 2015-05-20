function plotROC(specificity, sensitivity)
%plot specificity-sensitivity ROC curves
% figure
plot(specificity, sensitivity)
% title 'ejection fraction (EF) ROC'
xlabel ' specificity'
ylabel ' sensitivity'
end
