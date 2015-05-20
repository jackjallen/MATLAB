function[data] = extractEndoEpi(data, shape_nRows)

for i = 1:400
% d = find(concatIndices'==i);

% i
% reshape myo shapes from vector to matrix
% data(i).diastolic.myo.xyz = reshape(data(i).diastolic.myo.xyz,[2*shape_nRows, shape_nCols]);
% data(i).systolic.myo.xyz = reshape(data(i).systolic.myo.xyz,[2*shape_nRows, shape_nCols]);

% diastolic_myo_reshaped(i).xyz = reshape(diastolic_myo_shapes(:,i),size(dia_myo_reference));
% systolic_myo_reshaped(i).xyz = reshape(systolic_myo_shapes(:,i),size(sys_myo_reference));

% extract endo and epi shapes from full shape matrices (data(i).diastolic.myo.xyz)
data(i).diastolic.endo.xyz = data(i).diastolic.myo.xyz(1:shape_nRows, :);
data(i).diastolic.epi.xyz = data(i).diastolic.myo.xyz(shape_nRows+1:2*shape_nRows, :);
data(i).systolic.endo.xyz = data(i).systolic.myo.xyz(1:shape_nRows, :);
data(i).systolic.epi.xyz = data(i).systolic.myo.xyz(shape_nRows+1:2*shape_nRows, :);

end