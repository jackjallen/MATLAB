function[data, accuracy, sensitivity, specificity] = calcAccuracyAVratio(data)

for n = 1:100
    n;
    for d = data(1).DETERMINE_indices(n)
       d;
        positive(n) = cell2mat({data(d).DETERMINE_sys_endo_AVratio});     
    end
    for m = data(1).MESA_indices(n)
      m;
        negative(n) = cell2mat({data(m).MESA_sys_endo_AVratio});
    end
    
end
positive =  positive;
negative = negative;

for i = 1:50
    % if we class over threshold as MESA ('negative') and under threshold as
    % DETERMINE ('positive')
    
    
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

