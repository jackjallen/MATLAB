
%SVM - 3 parameters
tic
SVMclassification_rates = zeros((size(param1,2)),(size(param2,2)));
for param1 = 2:33
    for param2 = 29
        for   param3 = 1 % ejection fractions
            tic
            param1
            param2
            % SVMModel = fitcsvm(trainingdata, group, 'Standardize', true)
            SVMModel = fitcsvm(trainingdata(:,[param1;param2]), names, 'KernelFunction', 'linear' );
            %SVMModel = fitcsvm(trainingdata(:,:), names, 'KernelFunction', 'linear' );
            
            
            % cross-validation of the SVM
            %             CVSVMModel = crossval(SVMModel);
            %             misclassification_rate = kfoldLoss(CVSVMModel);
            %             %         classification_rates(param1,param2) = 1 - misclassification_rate
            %             SVMclassification_rates(param1, param2) = 1 - misclassification_rate
            %
        end
    end
end

toc