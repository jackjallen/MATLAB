% LDA classification - 2 parameters

LDAcross_validated_classification_rates =  zeros(size(trainingdata,2),size(trainingdata,2));
tic
for param1 = 1:size(trainingdata,2)
    tic
    param1
    for param2 = 1
        param2
        obj = fitcdiscr(trainingdata(:,[param1;param2]), names);
        % obj = fitcdiscr(trainingdata(:,:), names); %using all 12 gives 94% classified
        resuberror = resubLoss(obj); %proportion of misclassified observations
        LDAclassification_rates(param1,param2) = 1 - resuberror;
        
        cvmodel = crossval(obj,'leaveout','on');
        cverror = kfoldLoss(cvmodel);
        LDAcross_validated_classification_rates(param1, param2) = 1 - cverror;
        
        
    end
    toc
end
toc
max_rate = max(max(LDAcross_validated_classification_rates))

size(LDAcross_validated_classification_rates,1)
figure
% subplot 121
plot([0 size(LDAcross_validated_classification_rates,1)],[max_rate max_rate])
legend (['max classification = ' num2str(max_rate)])
hold on
plot([0 size(LDAcross_validated_classification_rates,1)],[0.87 0.87])
bar( LDAcross_validated_classification_rates(:,1))
ylabel 'classification success'
ylim ([0.8 1])
xlim ([1 size(LDAcross_validated_classification_rates,1)])
set(gca,'XTick',[1:1:size(LDAcross_validated_classification_rates,1)])
% set(gca,'XTickLabel',trainingdataLabels)
xlabel 'parameter'
set(gca,'yMinorTick','on')
title 'EF + ...'
% subplot 122
% histogram(LDAclassification_rates)
% xlabel

[bestParam1,bestParam2] = find(LDAcross_validated_classification_rates==max(max(LDAcross_validated_classification_rates)))
[best1] = trainingdataLabels(bestParam1,:)
[best2] = trainingdataLabels(bestParam2,:)