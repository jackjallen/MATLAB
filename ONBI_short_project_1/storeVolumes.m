function [data] = storeVolumes(data, sys_epi_volumes, sys_endo_volumes, dia_epi_volumes, dia_endo_volumes)

for i = 1:401
    data(i).diastolic.endo.volume = dia_endo_volumes(i,1);
    data(i).diastolic.epi.volume = dia_epi_volumes(i,1);
    data(i).systolic.endo.volume = sys_endo_volumes(i,1);
    data(i).systolic.epi.volume = sys_epi_volumes(i,1);
    
    for d = 1:100
    if data(1).DETERMINE_indices(d)==i
        data(i).DETERMINE.diastolic.endo.volume = data(i).diastolic.endo.volume;
        data(i).DETERMINE.diastolic.epi.volume = data(i).diastolic.epi.volume;
        data(i).DETERMINE.systolic.endo.volume = data(i).systolic.endo.volume;
        data(i).DETERMINE.systolic.epi.volume = data(i).systolic.epi.volume;
      
    else
    end
    end
    
    for d = 1:100
    if data(1).MESA_indices(d)==i
%         i;
%         m;
    data(i).MESA.diastolic.endo.volume = data(i).diastolic.endo.volume;
    data(i).MESA.diastolic.epi.volume = data(i).diastolic.epi.volume;
    data(i).MESA.systolic.endo.volume = data(i).systolic.endo.volume;
    data(i).MESA.systolic.epi.volume = data(i).systolic.epi.volume;
    
    else
    end
    end
    
    
end
