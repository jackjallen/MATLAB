function[data, accuracy, sensitivity, specificity] = calcAccuracy(data)

for i = 1:401
    for n = 1:100
        
        if data(1).DETERMINE_indices(n)==i
            positive(n) = data(i).ejectionFraction;
        else
        end
        if data(1).MESA_indices(n)==i
%             n;
%             i;
            negative(n) = data(i).ejectionFraction;
        else
        end
    end
end


for i = 1:100

% if we class over threshold as MESA ('negative') and under threshold as
% DETERMINE ('positive')
false_negative = sum(positive>i);
true_negative = sum(negative>i);
true_positive = sum(positive<i);
false_positive = sum(negative<i);
%fraction, not percentage
accuracy(i) = (true_positive + true_negative)/(true_positive + false_positive + false_negative + true_negative);

sensitivity(i) = true_positive / (true_positive + false_negative );
specificity(i) = true_negative / (true_negative + false_positive );
end

