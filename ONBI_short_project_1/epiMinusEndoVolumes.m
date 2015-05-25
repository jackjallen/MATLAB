function [data] = epiMinusEndoVolumes(data)

for i = 1:401
    data(i).diastolic.epiMinusEndoVolume = data(i).diastolic.epi.volume - data(i).diastolic.endo.volume ;
    data(i).systolic.epiMinusEndoVolume = data(i).systolic.epi.volume - data(i).systolic.endo.volume ;
    
    for d = 1:100
        if data(1).DETERMINE_indices(d)==i
            data(i).DETERMINE.diastolic.epiMinusEndoVolume = data(i).diastolic.epiMinusEndoVolume
            data(i).DETERMINE.systolic.epiMinusEndoVolume = data(i).systolic.epiMinusEndoVolume
            
            
        else
        end
    end
    
    for d = 1:100
        if data(1).MESA_indices(d)==i
            %         i;
            %         m;
            data(i).MESA.diastolic.epiMinusEndoVolume = data(i).diastolic.epiMinusEndoVolume
            data(i).MESA.systolic.epiMinusEndoVolume = data(i).systolic.epiMinusEndoVolume
            
            
            
        else
        end
    end
    
    
end