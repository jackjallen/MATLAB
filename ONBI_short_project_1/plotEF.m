function plotEF(data)

figure
hold on
nBins = 30;
histogram(cell2mat({data((data(1).MESA_indices)).ejectionFraction}),nBins)
histogram(cell2mat({data((data(1).DETERMINE_indices)).ejectionFraction}),nBins)

xlabel 'Ejection Fraction (%)'
ylabel 'frequency'
legend 'MESA' 'DETERMINE'

end