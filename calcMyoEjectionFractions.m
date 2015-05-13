function [data] = calcMyoEjectionFractions(data)

for i = 1:401
    %find all the myocardium stroke volumes (SVs) and ejection fractions (EFs)
    data(i).myoSV = data(i).diastolic.myoVolume - data(i).systolic.myoVolume;
    data(i).myoEF = (data(i).myoSV/data(i).diastolic.myoVolume)*100;
    
    %from the full set of EFs and SVs, extract those of the two classes
    for n = 1:100
        if data(1).DETERMINE_indices(n)==i
            data(i).DETERMINE_myoSV = data(i).diastolic.myoVolume - data(i).systolic.myoVolume;
            data(i).DETERMINE_myoEF = (data(i).DETERMINE_myoSV / data(i).diastolic.myoVolume )*100;
        end
        
        if data(1).MESA_indices(n)==i
            data(i).MESA_myoSV = data(i).diastolic.myoVolume - data(i).systolic.myoVolume;
            data(i).MESA_myoEF = (data(i).MESA_myoSV / data(i).diastolic.myoVolume )*100;
        end
    end
    
end