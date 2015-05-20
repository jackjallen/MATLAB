function[data] = calcAVR(data)

%% calculate surface area (A) to volume (V) ratios (R)
%100 cases in each class

for d = data(1).DETERMINE_indices'

    data(d).DETERMINE_dia_endo_AVratio = data(d).DETERMINE_dia_endo_total_area/((data(d).diastolic.endo.volume)^(2/3));
    data(d).DETERMINE_dia_epi_AVratio= data(d).DETERMINE_dia_epi_total_area/((data(d).diastolic.epi.volume)^(2/3));
    data(d).DETERMINE_sys_endo_AVratio = data(d).DETERMINE_sys_endo_total_area/((data(d).systolic.endo.volume)^(2/3));
    data(d).DETERMINE_sys_epi_AVratio = data(d).DETERMINE_sys_epi_total_area/((data(d).systolic.epi.volume)^(2/3));
    data(d).DETERMINE_dia_myo_AVratio = data(d).DETERMINE_dia_myo_total_area/((data(d).diastolic.myoVolume)^(2/3));
    data(d).DETERMINE_sys_myo_AVratio = data(d).DETERMINE_sys_myo_total_area/((data(d).systolic.myoVolume)^(2/3));
end

for m = data(1).MESA_indices'
    data(m).MESA_dia_endo_AVratio = data(m).MESA_dia_endo_total_area/((data(m).diastolic.endo.volume)^(2/3));
    data(m).MESA_dia_epi_AVratio = data(m).MESA_dia_epi_total_area/((data(m).diastolic.epi.volume)^(2/3));
    data(m).MESA_sys_endo_AVratio = data(m).MESA_sys_endo_total_area/((data(m).systolic.endo.volume)^(2/3));
    data(m).MESA_sys_epi_AVratio = data(m).MESA_sys_epi_total_area/((data(m).systolic.epi.volume)^(2/3));
    data(m).MESA_dia_myo_AVratio = data(m).MESA_dia_myo_total_area/((data(m).diastolic.myoVolume)^(2/3));
    data(m).MESA_sys_myo_AVratio = data(m).MESA_sys_myo_total_area/((data(m).systolic.myoVolume)^(2/3));
end


end