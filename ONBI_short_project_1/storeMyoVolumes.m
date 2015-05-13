function[data] = storeMyoVolumes(data)


for i = 1:401
    data(i).dia_myovolumes = data(i).diastolic.myoVolume;
    data(i).sys_myovolumes = data(i).systolic.myoVolume;
    
    for n = 1:100
        
        if data(1).DETERMINE_indices(n)==i
            
            data(n).DETERMINE_diastolic_myovolumes =  data(i).dia_myovolumes;
            data(n).DETERMINE_systolic_myovolumes =  data(i).sys_myovolumes;
            
        end
        
        if data(1).MESA_indices(n) ==i
            data(n).MESA_diastolic_myovolumes = data(i).dia_myovolumes;
            data(n).MESA_systolic_myovolumes = data(i).sys_myovolumes;
        end
    end
end

end

