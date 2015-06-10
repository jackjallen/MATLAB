%lda - classification - 4 params

LDAcross_validated_classification_rates =  zeros(size(param1,2),size(param2,2));
for param4 = 1:34
    for param1 = 33
    param4
    for param2 = 27
        for param3 = 1 % param3 is EF
            tic
            obj = fitcdiscr(trainingdata(:,[param4, param1, param2 param3]), names);
            
            resuberror = resubLoss(obj); %proportion of misclassified observations
            LDAclassification_rates(param4) = 1 - resuberror;
            cvmodel = crossval(obj,'leaveout','on');
            cverror = kfoldLoss(cvmodel);
            LDAcross_validated_classification_rates(param4) = 1 - cverror;
            
        end
    end
end
end
hold on
% subplot 121
plot([0 size(trainingdata,2)],[max(max(LDAcross_validated_classification_rates)) max(max(LDAcross_validated_classification_rates))])
legend (['max classification = ' num2str(max(max(LDAcross_validated_classification_rates)))])
plot([0 size(trainingdata,2)],[0.87 0.87])
plot(LDAcross_validated_classification_rates(1:param4),'o')
ylabel 'classification success'
ylim ([0.8 1])
xlim ([1 size(trainingdata,2)])
xlabel 'training data column (parameter index)'
set(gca,'yMinorTick','on')
title 'EF + systolic endocardium area to volume ratio + ...'