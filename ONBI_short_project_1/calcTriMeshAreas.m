function [data] = calcTriMeshAreas(data)

for i = 1:401 
    for d = find(data(1).MESA_indices==i)
       
       data(i).MESA_dia_endo_sides = calcTriSides(data(i).diastolic.endo.tri, data(i).diastolic.endo.xyz); 
       data(i).MESA_dia_endo_areas(d,1) = calcTriMeshArea(data(i).MESA_dia_endo_sides);
        
       data(i).MESA_dia_epi_sides = calcTriSides(data(i).diastolic.epi.tri, data(i).diastolic.epi.xyz);
       data(i).MESA_dia_epi_areas(d,1) = calcTriMeshArea(data(i).MESA_dia_epi_sides);
       
       data(i).MESA_sys_endo_sides = calcTriSides(data(i).systolic.endo.tri, data(i).systolic.endo.xyz);
       data(i).MESA_sys_endo_areas(d,1) = calcTriMeshArea(data(i).MESA_sys_endo_sides);
     
       data(i).MESA_sys_epi_sides = calcTriSides(data(i).systolic.epi.tri, data(i).systolic.epi.xyz);
       data(i).MESA_sys_epi_areas(d,1) = calcTriMeshArea(data(i).MESA_sys_epi_sides);
      
       data(i).MESA_dia_myo_sides = calcTriSides(data(i).diastolic.myo.tri, data(i).diastolic.myo.xyz);
       data(i).MESA_dia_myo_areas(d,1) = calcTriMeshArea(data(i).MESA_dia_myo_sides);
    
       data(i).MESA_sys_myo_sides = calcTriSides(data(i).systolic.myo.tri, data(i).systolic.myo.xyz);
       data(i).MESA_sys_myo_areas(d,1) = calcTriMeshArea(data(i).MESA_sys_myo_sides);
    end

    for d = find(data(1).DETERMINE_indices==i)
       data(i).DETERMINE_dia_endo_sides = calcTriSides(data(i).diastolic.endo.tri, data(i).diastolic.endo.xyz); 
       data(i).ETERMINE_dia_endo_areas(d,1) = calcTriMeshArea(data(i).DETERMINE_dia_endo_sides);
          
       data(i).DETERMINE_dia_epi_sides = calcTriSides(data(i).diastolic.epi.tri, data(i).diastolic.epi.xyz);
       data(i).DETERMINE_dia_epi_areas(d,1) = calcTriMeshArea(data(i).DETERMINE_dia_epi_sides);
       
       data(i).DETERMINE_sys_endo_sides = calcTriSides(data(i).systolic.endo.tri, data(i).systolic.endo.xyz);
       data(i).DETERMINE_sys_endo_areas(d,1) = calcTriMeshArea(data(i).DETERMINE_sys_endo_sides);
       
       data(i).DETERMINE_sys_epi_sides = calcTriSides(data(i).systolic.epi.tri, data(i).systolic.epi.xyz);
       data(i).DETERMINE_sys_epi_areas(d,1) = calcTriMeshArea(data(i).DETERMINE_sys_epi_sides);
       
       data(i).DETERMINE_dia_myo_sides = calcTriSides(data(i).diastolic.myo.tri, data(i).diastolic.myo.xyz);
       data(i).DETERMINE_dia_myo_areas(d,1) = calcTriMeshArea(data(i).DETERMINE_dia_myo_sides);
       
       data(i).DETERMINE_sys_myo_sides = calcTriSides(data(i).systolic.myo.tri, data(i).systolic.myo.xyz);
       data(i).DETERMINE_sys_myo_areas(d,1) = calcTriMeshArea(data(i).DETERMINE_sys_myo_sides);
    end
end
