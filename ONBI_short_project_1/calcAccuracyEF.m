function[data, accuracy, sensitivity, specificity] = calcAccuracy(positive, negative)

% for n = 1:100
%     n
%     for d = data(1).DETERMINE_indices(n)
%        d
%         positive(n) = data(d).ejectionFraction;       
%     end
%     for m = data(1).MESA_indices(n)
%       m
%         negative(n) = data(m).ejectionFraction;       
%     end
%     
% end
positive;
negative;

for i = 1:100
    % if we class over threshold as MESA ('negative') and under threshold as
    % DETERMINE ('positive')
    i
    
    nfalse_negative = sum(positive>i);
    ntrue_negative = sum(negative>i);
    ntrue_positive = sum(positive<i);
    nfalse_positive = sum(negative<i);
    %fraction, not percentage
    accuracy(i) = (ntrue_positive + ntrue_negative)/(ntrue_positive + nfalse_positive + nfalse_negative + ntrue_negative);
    
    sensitivity(i) = ntrue_positive / (ntrue_positive + nfalse_negative );
    specificity(i) = ntrue_negative / (ntrue_negative + nfalse_positive );
end

end

