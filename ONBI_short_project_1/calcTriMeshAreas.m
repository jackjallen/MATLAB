function [data] = calcTriMeshAreas(data)

for i = 1:401
    
    %add area calc for all cases?
end

for m =data(1).MESA_indices'
    
   
   [data(m).MESA_dia_endo_total_area, ~, ~] = trifacet_area3D(data(m).diastolic.endo.tri, data(m).diastolic.endo.xyz);    
   [data(m).MESA_dia_epi_total_area, ~, ~] = trifacet_area3D(data(m).diastolic.epi.tri, data(m).diastolic.epi.xyz);
    
   [data(m).MESA_sys_endo_total_area, ~,~] = trifacet_area3D(data(m).systolic.endo.tri, data(m).systolic.endo.xyz);
   [data(m).MESA_sys_epi_total_area, ~, ~] = trifacet_area3D(data(m).systolic.epi.tri, data(m).systolic.epi.xyz);

   [data(m).MESA_dia_myo_total_area, ~, ~] = trifacet_area3D(data(m).diastolic.myo.tri, data(m).diastolic.myo.xyz);
   [data(m).MESA_sys_myo_total_area, ~, ~] = trifacet_area3D(data(m).systolic.myo.tri, data(m).systolic.myo.xyz);
end

for d = data(1).DETERMINE_indices'
 [data(d).DETERMINE_dia_endo_total_area, ~, ~] = trifacet_area3D(data(d).diastolic.endo.tri, data(d).diastolic.endo.xyz);    
   [data(d).DETERMINE_dia_epi_total_area, ~, ~] = trifacet_area3D(data(d).diastolic.epi.tri, data(d).diastolic.epi.xyz);
    
   [data(d).DETERMINE_sys_endo_total_area, ~,~] = trifacet_area3D(data(d).systolic.endo.tri, data(d).systolic.endo.xyz);
   [data(d).DETERMINE_sys_epi_total_area, ~, ~] = trifacet_area3D(data(d).systolic.epi.tri, data(d).systolic.epi.xyz);

   [data(d).DETERMINE_dia_myo_total_area, ~, ~] = trifacet_area3D(data(d).diastolic.myo.tri, data(d).diastolic.myo.xyz);
   [data(d).DETERMINE_sys_myo_total_area, ~, ~] = trifacet_area3D(data(d).systolic.myo.tri, data(d).systolic.myo.xyz);
end
end
