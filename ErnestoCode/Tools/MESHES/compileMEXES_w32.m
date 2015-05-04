function compileMEXES_w32(filename)

filenames=dir('*.cpp')
nFiles = length(filenames);

for i = 1:nFiles
    if( (~(strcmp('vtkPlantilla.cpp', filenames(i).name))) && ...
            (~(strcmp('carveIntersect.cpp', filenames(i).name))) && ...
            (~(strcmp('carveMesh2poly2Mesh.cpp', filenames(i).name))) && ...
            (~(strcmp('carvePlantilla.cpp', filenames(i).name))) )
        if(nargin == 0)
            [pathstr, name, ext] = fileparts(filenames(i).name);
            eval(['mex -v ' filenames(i).name ' -output Z:\MESHES_VTK\' name '.mexw32']);
        else
            if(strcmp(filename, filenames(i).name))
                mex -v filenames(i).name;
            end
        end
    end
end

% eval(['!chmod u+rwx .mexw32']);
% eval(['!chmod g+rwx .mexw32']);
% eval(['!chmod o-rwx .mexw32']);

% mex -v vtkFeatureEdges.cpp
% mex -v vtkPolyDataPointSampler.cpp 
% mex -v vtkBooleanOperationPolyDataFilter.cpp
% mex -v vtkFillHolesFilter.cpp
% mex -v vtkPolyDataReader.cpp
% mex -v vtkButterflySubdivisionFilter.cpp
% mex -v vtkHull.cpp
% mex -v vtkPolyDataWriter.cpp                  
% mex -v vtkCleanPolyData.cpp
% mex -v vtkIntersectionPolyDataFilter.cpp
% mex -v vtkQuadricDecimation.cpp
% mex -v vtkClipPolyData.cpp
% mex -v vtkLinearSubdivisionFilter.cpp
% mex -v vtkShrinkPolyData.cpp                  
% mex -v vtkClosestElement.cpp
% mex -v vtkMassProperties.cpp
% mex -v vtkSmoothPolyDataFilter.cpp
% mex -v vtkClosestPoint.cpp
% mex -v vtkWindowedSincPolyDataFilter.cpp
% mex -v vtkDecimatePro.cpp
% mex -v  vtkPolyDataConnectivityFilter.cpp
% mex -v vtkDensifyPolyData.cpp
% mex -v vtkPolyDataNormals.cpp