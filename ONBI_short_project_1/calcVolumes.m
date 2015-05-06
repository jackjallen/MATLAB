function [sys_epi_volumes, sys_endo_volumes, dia_epi_volumes, dia_endo_volumes] = calcVolumes(data)

% Function to calculate the volumes and centers of mass of the different
% areas of the LV.
%
% Receives a structure containing the 3D data pointc coordinates.
% Returns the volumes of the four different areas.
%
% JA 27/04/15 (scaled up from Ernesto's code from 24/04/2015 so as to include all 400 patients and all 4 area descriptions)

%%
% vtkCleanPolyData(EPI_ED) (JA:function provided by Ernesto)fix the possible replicated nodes and spurious
% edges
for i = 1:400 
%calculate the hole edge line coordinates xyz B for all patients
% %diastolic, endo
% data(i).diastolic.endo.xyz = reshape(data(i).diastolic.endo.xyz, [1089, 3] );
% data(i).diastolic.epi.xyz = reshape(data(i).diastolic.epi.xyz, [1089, 3] );
% data(i).systolic.endo.xyz = reshape(data(i).systolic.endo.xyz, [1089, 3] );
% data(i).systolic.epi.xyz = reshape(data(i).systolic.epi.xyz, [1089, 3] );

data(i).diastolic.endo.B = vtkFeatureEdges( vtkCleanPolyData(data(i).diastolic.endo) , 'BoundaryEdgesOn' , [] , 'FeatureEdgesOff' , [] );  %extracting the boundary
data(i).diastolic.endo.B.xyz = data(i).diastolic.endo.B.xyz( [2 1 3:end], : );  %fixing the connectivity.
%diastolic, epi
data(i).diastolic.epi.B = vtkFeatureEdges( vtkCleanPolyData(data(i).diastolic.epi) , 'BoundaryEdgesOn' , [] , 'FeatureEdgesOff' , [] );  %extracting the boundary
data(i).diastolic.epi.B.xyz = data(i).diastolic.epi.B.xyz( [2 1 3:end], : );  %fixing the connectivity.
%systolic, endo
data(i).systolic.endo.B = vtkFeatureEdges( vtkCleanPolyData(data(i).systolic.endo) , 'BoundaryEdgesOn' , [] , 'FeatureEdgesOff' , [] );  %extracting the boundary
data(i).systolic.endo.B.xyz = data(i).systolic.endo.B.xyz( [2 1 3:end], : );  %fixing the connectivity.
%systolic, epi
data(i).systolic.epi.B = vtkFeatureEdges( vtkCleanPolyData(data(i).systolic.epi) , 'BoundaryEdgesOn' , [] , 'FeatureEdgesOff' , [] );  %extracting the boundary
data(i).systolic.epi.B.xyz = data(i).systolic.epi.B.xyz( [2 1 3:end], : );  %fixing the connectivity.

%% find the centers of the upper holes (e.g diastolic_endo_C) and add to hole contour
% diastolic, endo
diastolic_endo_C = mean( data(i).diastolic.endo.B.xyz , 1 );
data(i).diastolic.endo.B.xyz = [ data(i).diastolic.endo.B.xyz ; diastolic_endo_C ];
% diastolic, epi
diastolic_epi_C = mean( data(i).diastolic.epi.B.xyz , 1 );
data(i).diastolic.epi.B.xyz = [ data(i).diastolic.epi.B.xyz ; diastolic_epi_C ];
% systolic, endo
systolic_endo_C = mean( data(i).systolic.endo.B.xyz , 1 );
data(i).systolic.endo.B.xyz = [ data(i).systolic.endo.B.xyz ; systolic_endo_C ];
% systolic, epi
systolic_epi_C = mean( data(i).systolic.epi.B.xyz , 1 );
data(i).systolic.epi.B.xyz = [ data(i).systolic.epi.B.xyz ; systolic_epi_C ];

data(i).diastolic.endo.B.tri = [];
data(i).diastolic.epi.B.tri = [];
data(i).systolic.endo.B.tri = [];
data(i).systolic.epi.B.tri = [];
%% make list of triangle face indices for the lids
%diastolic, endo
for t = 1:size(data(i).diastolic.endo.B.xyz,1)-2
    data(i).diastolic.endo.B.tri(end+1,:) = [ t , t+1 , size(data(i).diastolic.endo.B.xyz,1) ];
end
data(i).diastolic.endo.B.tri(end+1,:) = [ t+1 , 1 , size(data(i).diastolic.endo.B.xyz,1) ];

%diastolic, epi
for t = 1:size(data(i).diastolic.epi.B.xyz,1)-2
    data(i).diastolic.epi.B.tri(end+1,:) = [ t , t+1 , size(data(i).diastolic.epi.B.xyz,1) ];
end
data(i).diastolic.epi.B.tri(end+1,:) = [ t+1 , 1 , size(data(i).diastolic.epi.B.xyz,1) ];

%systolic, endo
for t = 1:size(data(i).systolic.endo.B.xyz,1)-2
    data(i).systolic.endo.B.tri(end+1,:) = [ t , t+1 , size(data(i).systolic.endo.B.xyz,1) ];
end
data(i).systolic.endo.B.tri(end+1,:) = [ t+1 , 1 , size(data(i).systolic.endo.B.xyz,1) ];

%systolic, epi
for t = 1:size(data(i).systolic.epi.B.xyz,1)-2
    data(i).systolic.epi.B.tri(end+1,:) = [ t , t+1 , size(data(i).systolic.epi.B.xyz,1) ];
end
data(i).systolic.epi.B.tri(end+1,:) = [ t+1 , 1 , size(data(i).systolic.epi.B.xyz,1) ];

%% append both meshes
sys_epi_Last_prevID  = size( data(i).systolic.epi.xyz , 1 );
sys_endo_Last_prevID  = size( data(i).systolic.endo.xyz , 1 );
dia_epi_Last_prevID  = size( data(i).diastolic.epi.xyz , 1 );
dia_endo_Last_prevID  = size( data(i).diastolic.endo.xyz , 1 );

data(i).systolic.epi.xyz = [ data(i).systolic.epi.xyz ; data(i).systolic.epi.B.xyz ];
data(i).systolic.endo.xyz = [ data(i).systolic.endo.xyz ; data(i).systolic.endo.B.xyz ];
data(i).diastolic.epi.xyz = [ data(i).diastolic.epi.xyz ; data(i).diastolic.epi.B.xyz ];
data(i).diastolic.endo.xyz = [ data(i).diastolic.endo.xyz ; data(i).diastolic.endo.B.xyz ];

data(i).systolic.epi.tri = [ data(i).systolic.epi.tri ; data(i).systolic.epi.B.tri + sys_epi_Last_prevID ];
data(i).systolic.endo.tri = [ data(i).systolic.endo.tri ; data(i).systolic.endo.B.tri + sys_endo_Last_prevID ];
data(i).diastolic.epi.tri = [ data(i).diastolic.epi.tri ; data(i).diastolic.epi.B.tri + dia_epi_Last_prevID ];
data(i).diastolic.endo.tri = [ data(i).diastolic.endo.tri ; data(i).diastolic.endo.B.tri + dia_endo_Last_prevID ];

%% make sure that every triangle in the new, closed mesh points outwards
data(i).systolic.epi = FixNormals( data(i).systolic.epi );
data(i).systolic.endo = FixNormals( data(i).systolic.endo );
data(i).diastolic.epi = FixNormals( data(i).diastolic.epi );
data(i).diastolic.endo = FixNormals( data(i).diastolic.endo );

%% calculate LV volume and center of mass (CoM)
[data(i).systolic.epi.volume, data(i).systolic.epi.centerofmass] = MeshVolume(data(i).systolic.epi);
[data(i).systolic.endo.volume, data(i).systolic.endo.centerofmass] = MeshVolume(data(i).systolic.endo);
[data(i).diastolic.epi.volume, data(i).diastolic.epi.centerofmass] = MeshVolume(data(i).diastolic.epi);
[data(i).diastolic.endo.volume, data(i).diastolic.endo.centerofmass] = MeshVolume(data(i).diastolic.endo);

%is the calculated volume in the right ball park?
data(i).systolic.epi.differencevolume = prod(diff( BBMesh( data(i).systolic.epi ) , 1  , 1 ) ) - data(i).systolic.epi.volume;  %%it should be positive!!
data(i).systolic.endo.differencevolume = prod(diff( BBMesh( data(i).systolic.endo ) , 1  , 1 ) ) - data(i).systolic.endo.volume;  %%it should be positive!!
data(i).diastolic.epi.differencevolume = prod(diff( BBMesh( data(i).diastolic.epi ) , 1  , 1 ) ) - data(i).diastolic.epi.volume; %%it should be positive!!
data(i).diastolic.endo.differencevolume = prod(diff( BBMesh( data(i).diastolic.endo ) , 1  , 1 ) ) - data(i).diastolic.endo.volume; %%it should be positive!!

%store the difference volumes, in order plot and check are all positive
% sys_epi_diff_volumes(i,1) = data(i).systolic.epi.differencevolume;
% sys_endo_diff_volumes(i,1) = data(i).systolic.endo.differencevolume;
% dia_epi_diff_volumes(i,1) = data(i).diastolic.epi.differencevolume;
% dia_endo_diff_volumes(i,1) = data(i).diastolic.endo.differencevolume;

%store the volumes in vectors
sys_epi_volumes(i,1) = data(i).systolic.epi.volume;
sys_endo_volumes(i,1) = data(i).systolic.endo.volume;
dia_epi_volumes(i,1) = data(i).diastolic.epi.volume;
dia_endo_volumes(i,1) = data(i).diastolic.endo.volume;


end
% hold on
% plot(sys_epi_diff_volumes)
% plot(sys_endo_diff_volumes)
% plot(dia_epi_diff_volumes)
% plot(dia_endo_diff_volumes)
% legend('sys epi diff volumes', 'sys endo diff volumes','dia epi diff volumes','dia endo diff volumes')




end
   