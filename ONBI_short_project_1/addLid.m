function [data.systolic.epi.tri, data.systolic.endo.tri, data.diastolic.epi.tri, data.diastolic.endo.tri] = addLid(data)
% this function was written at an attempt at splitting calcVolumes into
% multiple functions.
%JA 28/04/2015
for i = 1:400 
%calculate the hole edge line coordinates xyz B for all patients
%diastolic, endo
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

% appended_sys_epi = data(i).systolic.epi.tri;
% appended_sys_endo = data(i).systolic.endo.tri;
% appended_dia_epi = data(i).diastolic.epi.tri;
% appended_dia_endo = data(i).diastolic.endo.tri;

end

end