function[data] = storeMyoVolumes(data)


for i = 1:401
%     data(i).dia_myovolumes = data(i).diastolic.myoVolume;
%     data(i).sys_myovolumes = data(i).systolic.myoVolume;
%     
    data(i).dia_epiMinusEndoVolumes = data(i).diastolic.epiMinusEndoVolume;
    data(i).sys_epiMinusEndoVolumes = data(i).systolic.epiMinusEndoVolume;
    
end

% for n = data(1).DETERMINE_indices'
%     
%     data(n).DETERMINE.diastolic.myovolume =  data(n).dia_myovolumes;
%     data(n).DETERMINE.systolic.myovolume =  data(n).sys_myovolumes;    
% end
% 
% for n = data(1).MESA_indices'
%     data(n).MESA.diastolic.myovolume = data(n).dia_myovolumes;
%     data(n).MESA.systolic.myovolume = data(n).sys_myovolumes;    
% end


end

