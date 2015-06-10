LDAcross_validated_classification_rates =  zeros(size(param1,2),size(param2,2));
for param1 = 1:size(trainingdata,2)
    param1
    for param2 = 32
        for param3 = 1 % param3 is EF
            tic
            obj = fitcdiscr(trainingdata(:,[param1, param2 param3]), names);
            
            resuberror = resubLoss(obj); %proportion of misclassified observations
            LDAclassification_rates(param1, param2) = 1 - resuberror;
            cvmodel = crossval(obj,'leaveout','on');
            cverror = kfoldLoss(cvmodel);
            LDAcross_validated_classification_rates(param1, param2) = 1 - cverror;
            
        end
    end
end

hold on
% subplot 121
plot([0 size(trainingdata,2)],[max(max(LDAcross_validated_classification_rates)) max(max(LDAcross_validated_classification_rates))])
legend (['max classification = ' num2str(max(max(LDAcross_validated_classification_rates)))])
plot([0 size(trainingdata,2)],[0.87 0.87])
plot(LDAcross_validated_classification_rates(:,param2),'o')
ylabel 'classification success'
ylim ([0.8 1])
xlim ([1 size(trainingdata,2)])
xlabel 'training data column (parameter index)'
set(gca,'yMinorTick','on')
title 'EF + systolic endocardium area to volume ratio + ...'
% subplot 122

%parameters that give best classification when used with EF
max(max(LDAcross_validated_classification_rates))
[bestParam1,bestParam2] = find(LDAcross_validated_classification_rates==max(max(LDAcross_validated_classification_rates)))
[best1] = trainingdataLabels(bestParam1,:)
[best2] = trainingdataLabels(bestParam2,:)
% figure
% subplot 121
% bar3(LDAclassification_rates(:,:))
% zlim ([0.8 1])
% subplot 122
% imagesc((LDAclassification_rates(:,:)))
% title 'LDA - 3 parameters'
% xlabel 'training data column'
% ylabel 'training data column'
% colorbar

% LDAclassification_rates = table(LDAclassification_rates(:,:))

