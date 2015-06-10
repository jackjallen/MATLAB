function make( file )

if (nargin == 0)
    file = 'default';
end

    cwd = pwd;  CLEANUP = onCleanup( @() cd(cwd) );

    allProjects = {
        'IntersectMeshes'
        'IntersectMeshes_OK'
        'carveIntersect'
        'IntersectMeshes_OK2'
        'IntersectMeshes_OK2_single'
        'IntersectMeshes_OK2_reduced'
        'triangle_lib'
        'triangle_lib_single'
        'triangle_batch'
        'vtkFeatureEdges'
        'vtkPolyDataPointSampler'
        'vtkBooleanOperationPolyDataFilter'
        'vtkFillHolesFilter'
        'vtkPolyDataReader'
        'vtkButterflySubdivisionFilter'
        'vtkHull'
        'vtkPolyDataWriter'
        'vtkCleanPolyData'
        'vtkIntersectionPolyDataFilter'
        'vtkQuadricDecimation'
        'vtkClipPolyData'
        'vtkLinearSubdivisionFilter'
        'vtkShrinkPolyData'
        'vtkClosestElement'
        'vtkMassProperties'
        'vtkSmoothPolyDataFilter'
        'vtkClosestPoint'
        'vtkWindowedSincPolyDataFilter'
        'vtkDecimatePro'
        ' vtkPolyDataConnectivityFilter'
        'vtkDensifyPolyData'
        'vtkPolyDataNormals'
        'vtkContour'
        'boundary3d'
        'DeletePoints'
        'IsInsidePolyData'
        'all_vtk'
        };
    
    if(ispc)
        %Defines
        DEFINES_TRI = [' -DNO_TIMER -DANSI_DECLARATORS -DTRILIBRARY '];
        
        %Flags
        COMP_FLAGS=' -v -O ';
        
        %miatools
        MIATOOLS_MEXOPTS=[' -f ' ['"C:\Documents and Settings\mmonla\Datos de programa\MathWorks\MATLAB\R2012b\mexopts_clean.bat" '] ];
        LIBS_MIATOOLS_PATH={'Z:\GridData\lib\', 'Z:\Tools\parseargs\libParseargs_st'};
        DIRS_MIATOOLS_PATH=['Z:\Tools'];
        INCLUDE_DIRS_MIATOOLS= [' -I' LIBS_MIATOOLS_PATH{1} ' -I' LIBS_MIATOOLS_PATH{2} ' -I' DIRS_MIATOOLS_PATH ];
        LIBS_MIATOOLS=['-L' DIRS_MIATOOLS_PATH ' -ldl'];

        
        %VTK
        INCLUDE_DIRS_VTK=[' -IZ:\MESHES_VTK -IC:\Projects\VTK5.10.1 -IC:\Projects\VTK5.10.1\Common ' ...
            '-IC:\Projects\VTK5.10.1\Utilities -IC:\Projects\VTK5.10.1\VolumeRendering -IC:\Projects\VTK5.10.1\Rendering '...
            '-IC:\Projects\VTK5.10.1\Charts -IC:\Projects\VTK5.10.1\Chemistry -IC:\Projects\VTK5.10.1_src\VolumeRendering '...
            '-IC:\Projects\VTK5.10.1_src\Hybrid -IC:\Projects\VTK5.10.1_src\Widgets -IC:\Projects\VTK5.10.1_src\Rendering '...
            '-IC:\Projects\VTK5.10.1_src\Charts -IC:\Projects\VTK5.10.1_src\Chemistry -IC:\Projects\VTK5.10.1_src\Rendering\Testing\Cxx '...
            '-IC:\Projects\VTK5.10.1_src\IO -IC:\Projects\VTK5.10.1_src\Imaging -IC:\Projects\VTK5.10.1_src\Graphics '...
            '-IC:\Projects\VTK5.10.1_src\GenericFiltering -IC:\Projects\VTK5.10.1_src\Filtering -IC:\Projects\VTK5.10.1_src\Common '...
            '-IC:\Projects\VTK5.10.1_src\Utilities -IC:\Projects\VTK5.10.1_src\Common\Testing\Cxx -IC:\Projects\VTK5.10.1\Utilities\vtknetcdf '...
            '-IC:\Projects\VTK5.10.1_src\Utilities\vtknetcdf -IC:\Projects\VTK5.10.1_src\Utilities\vtknetcdf\include '...
            '-IC:\Projects\VTK5.10.1\Utilities\vtklibproj4 -IC:\Projects\VTK5.10.1_src\Utilities\vtklibproj4 '...
            '-IC:\Projects\VTK5.10.1\Utilities\DICOMParser -IC:\Projects\VTK5.10.1_src\Utilities\DICOMParser '...
            '-IC:\Projects\VTK5.10.1\Utilities\vtkfreetype\include -IC:\Projects\VTK5.10.1_src\Utilities\vtkfreetype\include '...
            '-IC:\Projects\VTK5.10.1\Utilities\LSDyna -IC:\Projects\VTK5.10.1_src\Utilities\LSDyna -IC:\Projects\VTK5.10.1\Utilities\MaterialLibrary '...
            '-IC:\Projects\VTK5.10.1_src\Utilities\MaterialLibrary -IC:\Projects\VTK5.10.1\Utilities\vtkmetaio '...
            '-IC:\Projects\VTK5.10.1_src\Utilities\vtkmetaio -IC:\Projects\VTK5.10.1\Utilities\verdict -IC:\Projects\VTK5.10.1_src\Utilities\verdict '...
            '-IC:\Projects\VTK5.10.1\Utilities\vtkhdf5 -IC:\Projects\VTK5.10.1_src\Utilities\vtkhdf5 -IC:\Projects\VTK5.10.1\Utilities\vtkhdf5\src '...
            '-IC:\Projects\VTK5.10.1_src\Utilities\vtkhdf5\src -IC:\Projects\VTK5.10.1\Utilities\vtkhdf5\hl\src '...
            '-IC:\Projects\VTK5.10.1_src\Utilities\vtkhdf5\hl\src -IC:\Projects\VTK5.10.1_src\Utilities\utf8\source '...
            '-IC:\Projects\VTK5.10.1_src\Utilities\ftgl\src -IC:\Projects\VTK5.10.1\Utilities\ftgl '];
        
        LIBS_VTK_PATH='Z:\MESHES_VTK\vtk_libs\w32';
        INCLUDE_LIBS_VTK=[' -L' LIBS_VTK_PATH];
        LIBS_VTK=[' -lvtkCommon -lvtkDICOMParser -lvtkexoIIc -lvtkexpat '...
            '-lvtkFiltering -lvtkfreetype -lvtkftgl -lvtkGenericFiltering ' ...
            '-lvtkGraphics -lvtkHybrid -lvtkImaging -lvtkIO -lvtkjpeg -lvtkNetCDF -lvtkpng ' ...
            '-lvtkRendering -lvtksys -lvtktiff -lvtkVolumeRendering -lvtkWidgets -lvtkzlib '];
        %eval(['!set LIBS_VTK=' LIBS_VTK]);
        
        %LIBPATH_VTK = ['/LIBPATH:"' LIBS_VTK_PATH ' " ' LIBS_VTK ];
        

        %triangle
        DIRS_TRI_PATH='Z:\MESHES_VTK\triangle_libs\w32\';
        INCLUDE_DIRS_TRI=[' -I' DIRS_TRI_PATH ' '];
        LIBS_TRI=' -LZ:\MESHES_VTK\triangle_libs ';
        LIBS_TRIANGLE='triangle.lib ';
        
        %carve
        INCLUDE_DIRS_CARVE=' -IC:\Projects\carve1.4.0_src\include\carve -IC:\Projects\carve1.4.0_src\include ';
        LIBS_CARVE_PATH='Z:\MESHES_VTK\carve_libs\w32\';
        LIBS_CARVE=[' -L' LIBS_CARVE_PATH ' -lcarve.lib -lcarve_misc.lib -lcarve_fileformats.lib -lcarve_ui.lib '];
        
        %Parseargs, mzLib
        INCLUDE_DIRS_MZLIBS=[' -I' LIBS_MIATOOLS_PATH{1}];
        PARSEARGS_DEPENDENCY=[LIBS_MIATOOLS_PATH{1} 'libParseargs.lib '];
    end
    
    if(isunix)
        COMP_FLAGS=' -v -O ';
        %miatools
        MIATOOLS_MEXOPTS=' -f /extra/disco1/miaTools/Tools/mexopts_clean.sh ';
        INCLUDE_DIRS_MIATOOLS=' -I/extra/disco1/miaTools/MESHES_VTK -I/extra/disco1/miaTools/Tools -I/extra/disco1/miaTools/Tools/mztimer ';
        LIBS_MIATOOLS=' -L/extra/disco1/miaTools/Tools -ldl ';
        
        %VTK
        INCLUDE_DIRS_VTK=' -I/usr/local/include/vtk-5.10 ';
        
        LIBS_VTK_PATH='/extra/disco1/miaTools/MESHES_VTK/vtk_libs';
        INCLUDE_LIBS_VTK=[' -L' LIBS_VTK_PATH];
        LIBS_VTK=[' -lvtkCommon -lvtkDICOMParser ' ...
            '-lvtkexoIIc -lvtkexpat -lvtkFiltering -lvtkfreetype -lvtkftgl -lvtkGenericFiltering ' ...
            '-lvtkGraphics -lvtkHybrid -lvtkImaging -lvtkIO -lvtkjpeg -lvtkNetCDF -lvtkpng ' ...
            '-lvtkRendering -lvtksys -lvtktiff -lvtkVolumeRendering -lvtkWidgets -lvtkzlib '];
        
        %triangle
        INCLUDE_DIRS_TRI=' -I/extra/disco1/miaTools/MESHES_VTK/triangle_libs ';
        DIRS_TRI_PATH='/extra/disco1/miaTools/MESHES_VTK/triangle_libs/';
        LIBS_TRI=' -L/extra/disco1/miaTools/MESHES_VTK/triangle_libs ';
        
        %carve
        INCLUDE_DIRS_CARVE=' -I/usr/local/include -I/usr/local/src/carve-1.4.0/include -I/usr/local/CARVE1.4.0/include/ ';
        LIBS_CARVE=' -L/extra/disco1/miaTools/MESHES_VTK/carve_libs -lcarve -lcarve_misc -lcarve_fileformats -lcarve_ui ';
        
        %Parseargs, mzLib
        INCLUDE_DIRS_MZLIBS=' -I/extra/disco1/miaTools/Tools/parseargs/libParseargs_st -I/extra/disco1/miaTools/GridData/include ';
        DIRS_MIATOOLS_PATH='/extra/disco1/miaTools/Tools/';
        PARSEARGS_DEPENDENCY='libParseargs.a ';
    end
    
    switch file
        
        case {'IntersectMeshes.cpp','IntersectMeshes'}
            %para Windows
            if(ispc)
                if(~isarch64)
                    %-O -IZ:\GridData\ -IZ:\Tools -IZ:\MESHES_VTK\triangle_libs\  -IZ:\MESHES_VTK -IC:\Projects\VTK5.10.1 -IC:\Projects\VTK5.10.1\Common -IC:\Projects\VTK5.10.1\Utilities -IC:\Projects\VTK5.10.1\VolumeRendering -IC:\Projects\VTK5.10.1\Rendering -IC:\Projects\VTK5.10.1\Charts -IC:\Projects\VTK5.10.1\Chemistry -IC:\Projects\VTK5.10.1_src\VolumeRendering -IC:\Projects\VTK5.10.1_src\Hybrid -IC:\Projects\VTK5.10.1_src\Widgets -IC:\Projects\VTK5.10.1_src\Rendering -IC:\Projects\VTK5.10.1_src\Charts -IC:\Projects\VTK5.10.1_src\Chemistry -IC:\Projects\VTK5.10.1_src\Rendering\Testing\Cxx -IC:\Projects\VTK5.10.1_src\IO -IC:\Projects\VTK5.10.1_src\Imaging -IC:\Projects\VTK5.10.1_src\Graphics -IC:\Projects\VTK5.10.1_src\GenericFiltering -IC:\Projects\VTK5.10.1_src\Filtering -IC:\Projects\VTK5.10.1_src\Common -IC:\Projects\VTK5.10.1_src\Utilities -IC:\Projects\VTK5.10.1_src\Common\Testing\Cxx -IC:\Projects\VTK5.10.1\Utilities\vtknetcdf -IC:\Projects\VTK5.10.1_src\Utilities\vtknetcdf -IC:\Projects\VTK5.10.1_src\Utilities\vtknetcdf\include -IC:\Projects\VTK5.10.1\Utilities\vtklibproj4 -IC:\Projects\VTK5.10.1_src\Utilities\vtklibproj4 -IC:\Projects\VTK5.10.1\Utilities\DICOMParser -IC:\Projects\VTK5.10.1_src\Utilities\DICOMParser -IC:\Projects\VTK5.10.1\Utilities\vtkfreetype\include -IC:\Projects\VTK5.10.1_src\Utilities\vtkfreetype\include -IC:\Projects\VTK5.10.1\Utilities\LSDyna -IC:\Projects\VTK5.10.1_src\Utilities\LSDyna -IC:\Projects\VTK5.10.1\Utilities\MaterialLibrary -IC:\Projects\VTK5.10.1_src\Utilities\MaterialLibrary -IC:\Projects\VTK5.10.1\Utilities\vtkmetaio -IC:\Projects\VTK5.10.1_src\Utilities\vtkmetaio -IC:\Projects\VTK5.10.1\Utilities\verdict -IC:\Projects\VTK5.10.1_src\Utilities\verdict -IC:\Projects\VTK5.10.1\Utilities\vtkhdf5 -IC:\Projects\VTK5.10.1_src\Utilities\vtkhdf5 -IC:\Projects\VTK5.10.1\Utilities\vtkhdf5\src -IC:\Projects\VTK5.10.1_src\Utilities\vtkhdf5\src -IC:\Projects\VTK5.10.1\Utilities\vtkhdf5\hl\src -IC:\Projects\VTK5.10.1_src\Utilities\vtkhdf5\hl\src -IC:\Projects\VTK5.10.1_src\Utilities\utf8\source -IC:\Projects\VTK5.10.1_src\Utilities\ftgl\src -IC:\Projects\VTK5.10.1\Utilities\ftgl IntersectMeshes_OK2.cpp Z:\MESHES_VTK\triangle_libs\w32\triangle.lib -LZ:\MESHES_VTK\vtk_libs\w32 -lvtkCommon -lvtkDICOMParser -lvtkexoIIc -lvtkexpat -lvtkFiltering -lvtkfreetype -lvtkftgl -lvtkGenericFiltering -lvtkGraphics -lvtkHybrid -lvtkImaging -lvtkIO -lvtkjpeg -lvtkNetCDF -lvtkpng. -lvtkRendering -lvtksys -lvtktiff -lvtkVolumeRendering -lvtkWidgets -lvtkzlib -output IntersectMeshes_OK2.mexa64
                    eval(['mex' [COMP_FLAGS MIATOOLS_MEXOPTS INCLUDE_DIRS_MIATOOLS INCLUDE_DIRS_TRI INCLUDE_DIRS_VTK 'IntersectMeshes.cpp ' ...
                        DIRS_TRI_PATH 'triangle.lib ' INCLUDE_LIBS_VTK LIBS_VTK]]);
                else
                end
            end
            
            %Para linux           
            if(isunix)
                %-O -f /home/mia/mmonla/.matlab/R2012a/mexopts_carve.sh IntersectMeshes_OK.cpp -output IntersectMeshes_OK.mexa64
                eval(['mex' [COMP_FLAGS MIATOOLS_MEXOPTS INCLUDE_DIRS_MIATOOLS INCLUDE_DIRS_TRI INCLUDE_DIRS_VTK 'IntersectMeshes.cpp ' ...
                    DIRS_TRI_PATH 'libtriangle.a' INCLUDE_LIBS_VTK LIBS_VTK]]);
                
            end
        case {'IntersectMeshes_OK.cpp','IntersectMeshes_OK'}
            %para Windows
            if(ispc)
                if(~isarch64)
                    %-O -IZ:\GridData\ -IZ:\Tools -IZ:\MESHES_VTK\triangle_libs\  -IZ:\MESHES_VTK -IC:\Projects\VTK5.10.1 -IC:\Projects\VTK5.10.1\Common -IC:\Projects\VTK5.10.1\Utilities -IC:\Projects\VTK5.10.1\VolumeRendering -IC:\Projects\VTK5.10.1\Rendering -IC:\Projects\VTK5.10.1\Charts -IC:\Projects\VTK5.10.1\Chemistry -IC:\Projects\VTK5.10.1_src\VolumeRendering -IC:\Projects\VTK5.10.1_src\Hybrid -IC:\Projects\VTK5.10.1_src\Widgets -IC:\Projects\VTK5.10.1_src\Rendering -IC:\Projects\VTK5.10.1_src\Charts -IC:\Projects\VTK5.10.1_src\Chemistry -IC:\Projects\VTK5.10.1_src\Rendering\Testing\Cxx -IC:\Projects\VTK5.10.1_src\IO -IC:\Projects\VTK5.10.1_src\Imaging -IC:\Projects\VTK5.10.1_src\Graphics -IC:\Projects\VTK5.10.1_src\GenericFiltering -IC:\Projects\VTK5.10.1_src\Filtering -IC:\Projects\VTK5.10.1_src\Common -IC:\Projects\VTK5.10.1_src\Utilities -IC:\Projects\VTK5.10.1_src\Common\Testing\Cxx -IC:\Projects\VTK5.10.1\Utilities\vtknetcdf -IC:\Projects\VTK5.10.1_src\Utilities\vtknetcdf -IC:\Projects\VTK5.10.1_src\Utilities\vtknetcdf\include -IC:\Projects\VTK5.10.1\Utilities\vtklibproj4 -IC:\Projects\VTK5.10.1_src\Utilities\vtklibproj4 -IC:\Projects\VTK5.10.1\Utilities\DICOMParser -IC:\Projects\VTK5.10.1_src\Utilities\DICOMParser -IC:\Projects\VTK5.10.1\Utilities\vtkfreetype\include -IC:\Projects\VTK5.10.1_src\Utilities\vtkfreetype\include -IC:\Projects\VTK5.10.1\Utilities\LSDyna -IC:\Projects\VTK5.10.1_src\Utilities\LSDyna -IC:\Projects\VTK5.10.1\Utilities\MaterialLibrary -IC:\Projects\VTK5.10.1_src\Utilities\MaterialLibrary -IC:\Projects\VTK5.10.1\Utilities\vtkmetaio -IC:\Projects\VTK5.10.1_src\Utilities\vtkmetaio -IC:\Projects\VTK5.10.1\Utilities\verdict -IC:\Projects\VTK5.10.1_src\Utilities\verdict -IC:\Projects\VTK5.10.1\Utilities\vtkhdf5 -IC:\Projects\VTK5.10.1_src\Utilities\vtkhdf5 -IC:\Projects\VTK5.10.1\Utilities\vtkhdf5\src -IC:\Projects\VTK5.10.1_src\Utilities\vtkhdf5\src -IC:\Projects\VTK5.10.1\Utilities\vtkhdf5\hl\src -IC:\Projects\VTK5.10.1_src\Utilities\vtkhdf5\hl\src -IC:\Projects\VTK5.10.1_src\Utilities\utf8\source -IC:\Projects\VTK5.10.1_src\Utilities\ftgl\src -IC:\Projects\VTK5.10.1\Utilities\ftgl IntersectMeshes_OK2.cpp Z:\MESHES_VTK\triangle_libs\w32\triangle.lib -LZ:\MESHES_VTK\vtk_libs\w32 -lvtkCommon -lvtkDICOMParser -lvtkexoIIc -lvtkexpat -lvtkFiltering -lvtkfreetype -lvtkftgl -lvtkGenericFiltering -lvtkGraphics -lvtkHybrid -lvtkImaging -lvtkIO -lvtkjpeg -lvtkNetCDF -lvtkpng. -lvtkRendering -lvtksys -lvtktiff -lvtkVolumeRendering -lvtkWidgets -lvtkzlib -output IntersectMeshes_OK2.mexa64
%                     eval(['mex' [COMP_FLAGS MIATOOLS_MEXOPTS INCLUDE_DIRS_MIATOOLS INCLUDE_DIRS_TRI INCLUDE_DIRS_VTK 'IntersectMeshes_OK.cpp ' ...
%                         DIRS_TRI_PATH 'triangle.lib ' INCLUDE_LIBS_VTK LIBS_VTK]]);
                    
                    disp('Compiling Release.....');
                    
                    disp('Compiling IntersectMeshes:Release.....');
                    [status,result] = system('msbuild.exe C:\Projects\carve1.4.0\carve.sln /t:Clean;IntersectMeshes /p:Configuration=Release');
                    if(status ~=0), error('Error compilando IntersectMeshes: Release. Error: \n\n%s', result); end;
                    
                    disp('Compiling Debug.....');
                    
                    disp('Compiling IntersectMeshes:Debug.....');
                    [status,result] = system('msbuild.exe C:\Projects\carve1.4.0\carve.sln /t:Clean;IntersectMeshes /p:Configuration=Debug');
                    if(status ~=0), error('Error compilando IntersectMeshes: Debug. Error: \n\n%s', result); end;                    
                else
                end
                disp('EL PROYECTO SE COMPILA EN MODO RELEASE.');
                disp('LOS FICHEROS IntersectMeshes.mexw32 Y IntersectMeshesD.mexw32 SE COPIAN EN Z:\MESHES_VTK, C:\Projects\carve1.4.0\libs Y Z:\MESHES_VTK\carve_libs\w32.');
                disp('ESTE PROYECTO SE COMPILA CON VISUAL STUDIO QUE ESTA ES: C:\Projects\carve1.4.0\Mex\IntersectMeshes');
            end
            
            %Para linux
            if(isunix)
                %-O -f /home/mia/mmonla/.matlab/R2012a/mexopts_carve.sh IntersectMeshes_OK.cpp -output IntersectMeshes_OK.mexa64
                eval(['mex' [COMP_FLAGS MIATOOLS_MEXOPTS INCLUDE_DIRS_MIATOOLS INCLUDE_DIRS_MZLIBS INCLUDE_DIRS_TRI INCLUDE_DIRS_VTK 'IntersectMeshes_OK.cpp ' ...
                    DIRS_TRI_PATH 'libtriangle.a' INCLUDE_LIBS_VTK LIBS_VTK LIBS_MIATOOLS]]);
                
            end
            
        case {'IntersectMeshes_OK2.cpp','IntersectMeshes_OK2'}
            %para Windows
            if(ispc)
                if(~isarch64)
                    disp('Compiling Release.....');
                    
                    disp('Compiling IntersectMeshes2:Release.....');
                    [status,result] = system('msbuild.exe C:\Projects\carve1.4.0\carve.sln /t:Clean;IntersectMeshes2 /p:Configuration=Release');
                    if(status ~=0), error('Error compilando IntersectMeshes2: Release. Error: \n\n%s', result); end;
                    
                    disp('Compiling Debug.....');
                    
                    disp('Compiling IntersectMeshes2:Debug.....');
                    [status,result] = system('msbuild.exe C:\Projects\carve1.4.0\carve.sln /t:Clean;IntersectMeshes2 /p:Configuration=Debug');
                    if(status ~=0), error('Error compilando IntersectMeshes2: Debug. Error: \n\n%s', result); end;                    
                else
                end
                disp('EL PROYECTO SE COMPILA EN MODO RELEASE.');
                disp('LOS FICHEROS IntersectMeshes2.mexw32 Y IntersectMeshes2D.mexw32 SE COPIAN EN Z:\MESHES_VTK, C:\Projects\carve1.4.0\libs Y Z:\MESHES_VTK\carve_libs\w32.');
                disp('ESTE PROYECTO SE COMPILA CON VISUAL STUDIO QUE ESTA ES: C:\Projects\carve1.4.0\Mex\IntersectMeshes2');
            end
            
            %Para linux
            if(isunix)
                make triangle_lib
                %-O -f /home/mia/mmonla/.matlab/R2012a/mexopts_carve.sh IntersectMeshes_OK.cpp -output IntersectMeshes_OK.mexa64
                eval(['mex' [COMP_FLAGS MIATOOLS_MEXOPTS INCLUDE_DIRS_MIATOOLS INCLUDE_DIRS_MZLIBS INCLUDE_DIRS_TRI INCLUDE_DIRS_VTK 'IntersectMeshes_OK2.cpp ' ...
                    DIRS_TRI_PATH 'libtriangle.a' INCLUDE_LIBS_VTK LIBS_VTK LIBS_MIATOOLS]]);
                
                eval(['mex ' [' -v -g ' MIATOOLS_MEXOPTS INCLUDE_DIRS_MIATOOLS INCLUDE_DIRS_MZLIBS INCLUDE_DIRS_TRI INCLUDE_DIRS_VTK 'IntersectMeshes_OK2.cpp ' ...
                    DIRS_TRI_PATH 'libtriangleD.a' INCLUDE_LIBS_VTK LIBS_VTK LIBS_MIATOOLS ' -output IntersectMeshes_OK2D.mexa64']]);
                
            end
        case {'IntersectMeshes_OK2_single'}
            %para Windows
            if(ispc)
                if(~isarch64)
%                     disp('Compiling Release.....');
%                     
%                     disp('Compiling IntersectMeshes2:Release.....');
%                     [status,result] = system('msbuild.exe C:\Projects\carve1.4.0\carve.sln /t:Clean;IntersectMeshes2 /p:Configuration=Release');
%                     if(status ~=0), error('Error compilando IntersectMeshes2: Release. Error: \n\n%s', result); end;
%                     
%                     disp('Compiling Debug.....');
%                     
%                     disp('Compiling IntersectMeshes2:Debug.....');
%                     [status,result] = system('msbuild.exe C:\Projects\carve1.4.0\carve.sln /t:Clean;IntersectMeshes2 /p:Configuration=Debug');
%                     if(status ~=0), error('Error compilando IntersectMeshes2: Debug. Error: \n\n%s', result); end;
                else
                end
                disp('EL PROYECTO SE COMPILA EN MODO RELEASE.');
                disp('LOS FICHEROS IntersectMeshes2.mexw32 Y IntersectMeshes2D.mexw32 SE COPIAN EN Z:\MESHES_VTK, C:\Projects\carve1.4.0\libs Y Z:\MESHES_VTK\carve_libs\w32.');
                disp('ESTE PROYECTO SE COMPILA CON VISUAL STUDIO QUE ESTA ES: C:\Projects\carve1.4.0\Mex\IntersectMeshes2');
            end
            
            %Para linux
            if(isunix)
                make triangle_lib_single
                !cp IntersectMeshes_OK2.cp IntersectMeshes_OK2_single.cp
                %-O -f /home/mia/mmonla/.matlab/R2012a/mexopts_carve.sh IntersectMeshes_OK.cpp -output IntersectMeshes_OK.mexa64
                eval(['mex' [COMP_FLAGS MIATOOLS_MEXOPTS INCLUDE_DIRS_MIATOOLS INCLUDE_DIRS_MZLIBS INCLUDE_DIRS_TRI INCLUDE_DIRS_VTK 'IntersectMeshes_OK2_single.cpp ' ...
                    DIRS_TRI_PATH 'libtriangle_single.a' INCLUDE_LIBS_VTK LIBS_VTK LIBS_MIATOOLS]]);
                
                eval(['mex ' [' -v -g ' MIATOOLS_MEXOPTS INCLUDE_DIRS_MIATOOLS INCLUDE_DIRS_MZLIBS INCLUDE_DIRS_TRI INCLUDE_DIRS_VTK 'IntersectMeshes_OK2_single.cpp ' ...
                    DIRS_TRI_PATH 'libtriangleD_single.a' INCLUDE_LIBS_VTK LIBS_VTK LIBS_MIATOOLS ' -output IntersectMeshes_OK2D_single.mexa64']]);
                
            end
        case {'IntersectMeshes_OK2_reduced.cpp','IntersectMeshes_OK2_reduced'}
            %para Windows
            if(ispc)
                if(~isarch64)
%                     disp('Compiling Release.....');
%                     
%                     disp('Compiling IntersectMeshes2:Release.....');
%                     [status,result] = system('msbuild.exe C:\Projects\carve1.4.0\carve.sln /t:Clean;IntersectMeshes2 /p:Configuration=Release');
%                     if(status ~=0), error('Error compilando IntersectMeshes2: Release. Error: \n\n%s', result); end;
%                     
%                     disp('Compiling Debug.....');
%                     
%                     disp('Compiling IntersectMeshes2:Debug.....');
%                     [status,result] = system('msbuild.exe C:\Projects\carve1.4.0\carve.sln /t:Clean;IntersectMeshes2 /p:Configuration=Debug');
%                     if(status ~=0), error('Error compilando IntersectMeshes2: Debug. Error: \n\n%s', result); end;
                else
                end
                disp('EL PROYECTO SE COMPILA EN MODO RELEASE.');
                disp('LOS FICHEROS IntersectMeshes2.mexw32 Y IntersectMeshes2D.mexw32 SE COPIAN EN Z:\MESHES_VTK, C:\Projects\carve1.4.0\libs Y Z:\MESHES_VTK\carve_libs\w32.');
                disp('ESTE PROYECTO SE COMPILA CON VISUAL STUDIO QUE ESTA ES: C:\Projects\carve1.4.0\Mex\IntersectMeshes2');
            end
            
            %Para linux
            if(isunix)
                %-O -f /home/mia/mmonla/.matlab/R2012a/mexopts_carve.sh IntersectMeshes_OK.cpp -output IntersectMeshes_OK.mexa64
                eval(['mex' [COMP_FLAGS MIATOOLS_MEXOPTS INCLUDE_DIRS_MIATOOLS INCLUDE_DIRS_MZLIBS INCLUDE_DIRS_TRI INCLUDE_DIRS_VTK 'IntersectMeshes_OK2_reduced.cpp ' ...
                    INCLUDE_LIBS_VTK LIBS_VTK LIBS_MIATOOLS]]);
                
                eval(['mex ' [' -v -g ' MIATOOLS_MEXOPTS INCLUDE_DIRS_MIATOOLS INCLUDE_DIRS_MZLIBS INCLUDE_DIRS_TRI INCLUDE_DIRS_VTK 'IntersectMeshes_OK2_reduced.cpp ' ...
                    INCLUDE_LIBS_VTK LIBS_VTK LIBS_MIATOOLS ' -output IntersectMeshes_OK2D_reduced.mexa64']]);
                
            end
            
        case {'carveIntersect.cpp','carveIntersect'}
            if(ispc)
                if(~isarch64)
                    disp('Compiling Release.....');
                    
                    disp('Compiling carveIntersect:Release.....');
                    [status,result] = system('msbuild.exe C:\Projects\carve1.4.0\carve.sln /t:Clean;carveIntersect /p:Configuration=Release');
                    if(status ~=0), error('Error compilando carveIntersect: Release. Error: \n\n%s', result); end;
                    
                    disp('Compiling Debug.....');
                    
                    disp('Compiling carveIntersect:Debug.....');
                    [status,result] = system('msbuild.exe C:\Projects\carve1.4.0\carve.sln /t:Clean;carveIntersect /p:Configuration=Debug');
                    if(status ~=0), error('Error compilando carveIntersect: Debug. Error: \n\n%s', result); end;
                else
                end
%                 omitir = ['C:\ARCHIV~1\MICROS~2\VC\lib\libcmt.lib'];
%                 LIBS_CARVE = [LIBS_CARVE_PATH 'carve.lib ' LIBS_CARVE_PATH 'carve_misc.lib ' LIBS_CARVE_PATH 'carve_fileformats.lib ' LIBS_CARVE_PATH 'carve_ui.lib '];
%                 eval(['mex ' ['-v -O -DUNICODE -D_UNICODE -D_IMPORT_AS_ST_LIB '] ['LINKFLAGS="$LINKFLAGS /NODEFAULTLIB:"'] omitir ['" " ']  ...
%                     [INCLUDE_DIRS_MIATOOLS INCLUDE_DIRS_CARVE ...
%                     'carveIntersect.cpp ' PARSEARGS_DEPENDENCY LIBS_CARVE LIBS_MIATOOLS]]);

%                 LIBS_CARVE = [LIBS_CARVE_PATH 'carve.lib ' LIBS_CARVE_PATH 'carve_misc.lib ' LIBS_CARVE_PATH 'carve_fileformats.lib ' LIBS_CARVE_PATH 'carve_ui.lib '];
%                 eval(['mex ' ['-v -O -DUNICODE -D_UNICODE -D_IMPORT_AS_ST_LIB ']  ...
%                     [INCLUDE_DIRS_MIATOOLS INCLUDE_DIRS_CARVE ...
%                     'carveIntersect.cpp ' PARSEARGS_DEPENDENCY LIBS_CARVE LIBS_MIATOOLS]]);

                disp('EL PROYECTO SE COMPILA EN MODO RELEASE.');
                disp('LOS FICHEROS carveIntersect.mexw32 Y carveIntersectD.mexw32 SE COPIAN EN Z:\MESHES_VTK, C:\Projects\carve1.4.0\libs Y Z:\MESHES_VTK\carve_libs\w32.');
                disp('ESTE PROYECTO SE COMPILA CON VISUAL STUDIO QUE ESTA ES: C:\Projects\carve1.4.0\Mex\carveIntersect');
            end
            if(isunix)
                %Primero nos aseguramos de que esta el parserags en la
                %version necesaria
                cd /extra/disco1/miaTools
                make parseargs_st_mex
                cd(cwd)
                
                %-O -D_IMPORT_AS_ST_LIB -ldl -f /home/mia/mmonla/.matlab/R2012a/mexopts_carve.sh carveIntersect.cpp -output carveIntersect.mexa64
                eval(['mex' [COMP_FLAGS ' -D_IMPORT_AS_ST_LIB' ...
                    MIATOOLS_MEXOPTS INCLUDE_DIRS_MIATOOLS INCLUDE_DIRS_CARVE INCLUDE_DIRS_MZLIBS ...
                    'carveIntersect.cpp ' ...
                    DIRS_MIATOOLS_PATH PARSEARGS_DEPENDENCY LIBS_CARVE LIBS_MIATOOLS ]]);
                
                disp('EL PROYECTO SE COMPILA EN MODO RELEASE.');
                disp('TAMBIEN PUEDE COMPILARSE DESDE nsight CON EL PROYECTO carveIntersect EN EL WORKSPACE InGrid-workspace.');
            end
            
%         case {'carveIntersectD.cpp','carveIntersectD'}
%             if(ispc)
%                 if(~isarch64)
%                 else
%                 end
%                 
%                 disp('ESTE PROYECTO SE COMPILA CON VISUAL STUDIO QUE ESTA ES: C:\Projects\carve1.4.0\Mex\carveIntersect');
%             end
%             if(isunix)
%                 %Primero nos aseguramos de que esta el parserags en la
%                 %version necesaria
%                 cd /extra/disco1/miaTools
%                 make parseargs_st_mex
%                 cd(cwd)
%                 
%                 %-O -D_IMPORT_AS_ST_LIB -ldl -f /home/mia/mmonla/.matlab/R2012a/mexopts_carve.sh carveIntersect.cpp -output carveIntersect.mexa64
%                 eval(['mex' [COMP_FLAGS ' -D_IMPORT_AS_ST_LIB' ...
%                     MIATOOLS_MEXOPTS INCLUDE_DIRS_MIATOOLS INCLUDE_DIRS_CARVE INCLUDE_DIRS_MZLIBS ...
%                     'carveIntersect.cpp ' ...
%                     DIRS_MIATOOLS_PATH PARSEARGS_DEPENDENCY LIBS_CARVE LIBS_MIATOOLS ]]);
%                 
%                 disp('EL PROYECTO SE COMPILA EN MODO DEBUG.');
%                 disp('TAMBIEN PUEDE COMPILARSE DESDE nsight CON EL PROYECTO carveIntersect EN EL WORKSPACE InGrid-workspace.');
%             end
            
        case {'triangle_lib'}
            if(ispc)
                cd Y:\triangle
                make triangle_lib
            end
            if(isunix)
                cd /extra/disco1/miaTools/MESHES_VTK/triangle_libs
                make triangle_lib
            end
        case {'triangle_lib_single'}
            if(ispc)
                cd Y:\triangle
                make triangle_lib_single
            end
            if(isunix)
                cd /extra/disco1/miaTools/MESHES_VTK/triangle_libs
                make triangle_lib_single
            end
            
        case {'triangle_batch'}
            if(ispc)
                cd Y:\triangle
                make triangle_batch
            end
            if(isunix)
                cd /extra/disco1/miaTools/MESHES_VTK/triangle_libs
                make triangle_batch
            end
        case {'vtkFeatureEdges'}
            if(ispc)
                if(~isarch64)
                    eval(['mex ' [COMP_FLAGS MIATOOLS_MEXOPTS INCLUDE_DIRS_MIATOOLS INCLUDE_DIRS_VTK ' vtkFeatureEdges.cpp ' INCLUDE_LIBS_VTK LIBS_VTK]]);
                else
                end
            end
            if(isunix)
                eval(['mex ' [COMP_FLAGS MIATOOLS_MEXOPTS INCLUDE_DIRS_MIATOOLS INCLUDE_DIRS_VTK ' vtkFeatureEdges.cpp ' INCLUDE_LIBS_VTK LIBS_VTK]]);
                !chmod ug+rwx vtkFeatureEdges.mexa64
            end
        case {'vtkPolyDataPointSampler'}
            if(ispc)
                if(~isarch64)
                    eval(['mex ' [COMP_FLAGS MIATOOLS_MEXOPTS INCLUDE_DIRS_MIATOOLS INCLUDE_DIRS_VTK ' vtkPolyDataPointSampler.cpp ' INCLUDE_LIBS_VTK LIBS_VTK]]);
                else
                end
            end
            if(isunix)
                eval(['mex ' [COMP_FLAGS MIATOOLS_MEXOPTS INCLUDE_DIRS_MIATOOLS INCLUDE_DIRS_VTK ' vtkPolyDataPointSampler.cpp ' INCLUDE_LIBS_VTK LIBS_VTK]]);
                !chmod ug+rwx vtkPolyDataPointSampler.mexa64
            end
        case {'vtkBooleanOperationPolyDataFilter'}
            if(ispc)
                if(~isarch64)
                    eval(['mex ' [COMP_FLAGS MIATOOLS_MEXOPTS INCLUDE_DIRS_MIATOOLS INCLUDE_DIRS_VTK ' vtkBooleanOperationPolyDataFilter.cpp ' INCLUDE_LIBS_VTK LIBS_VTK]]);
                else
                end
            end
            if(isunix)
                eval(['mex ' [COMP_FLAGS MIATOOLS_MEXOPTS INCLUDE_DIRS_MIATOOLS INCLUDE_DIRS_VTK ' vtkBooleanOperationPolyDataFilter.cpp ' INCLUDE_LIBS_VTK LIBS_VTK]]);
                !chmod ug+rwx vtkBooleanOperationPolyDataFilter.mexa64
            end
        case {'vtkFillHolesFilter'}
            if(ispc)
                if(~isarch64)
                    eval(['mex ' [COMP_FLAGS MIATOOLS_MEXOPTS INCLUDE_DIRS_MIATOOLS INCLUDE_DIRS_VTK ' vtkFillHolesFilter.cpp ' INCLUDE_LIBS_VTK LIBS_VTK]]);
                else
                end
            end
            if(isunix)
                eval(['mex ' [COMP_FLAGS MIATOOLS_MEXOPTS INCLUDE_DIRS_MIATOOLS INCLUDE_DIRS_VTK ' vtkFillHolesFilter.cpp ' INCLUDE_LIBS_VTK LIBS_VTK]]);
                !chmod ug+rwx vtkFillHolesFilter.mexa64
            end
        case {'vtkPolyDataReader'}
            if(ispc)
                if(~isarch64)
                    eval(['mex ' [COMP_FLAGS MIATOOLS_MEXOPTS INCLUDE_DIRS_MIATOOLS INCLUDE_DIRS_VTK ' vtkPolyDataReader.cpp ' INCLUDE_LIBS_VTK LIBS_VTK]]);
                else
                end
            end
            if(isunix)
                eval(['mex ' [COMP_FLAGS MIATOOLS_MEXOPTS INCLUDE_DIRS_MIATOOLS INCLUDE_DIRS_VTK ' vtkPolyDataReader.cpp ' INCLUDE_LIBS_VTK LIBS_VTK]]);
                !chmod ug+rwx vtkPolyDataReader.mexa64
            end
        case {'vtkButterflySubdivisionFilter'}
            if(ispc)
                if(~isarch64)
                    eval(['mex ' [COMP_FLAGS MIATOOLS_MEXOPTS INCLUDE_DIRS_MIATOOLS INCLUDE_DIRS_VTK ' vtkButterflySubdivisionFilter.cpp ' INCLUDE_LIBS_VTK LIBS_VTK]]);
                else
                end
            end
            if(isunix)
                eval(['mex ' [COMP_FLAGS MIATOOLS_MEXOPTS INCLUDE_DIRS_MIATOOLS INCLUDE_DIRS_VTK ' vtkButterflySubdivisionFilter.cpp ' INCLUDE_LIBS_VTK LIBS_VTK]]);
                !chmod ug+rwx vtkButterflySubdivisionFilter.mexa64
            end

        case {'vtkHull'}
            if(ispc)
                if(~isarch64)
                    eval(['mex ' [COMP_FLAGS MIATOOLS_MEXOPTS INCLUDE_DIRS_MIATOOLS INCLUDE_DIRS_VTK ' vtkHull.cpp ' INCLUDE_LIBS_VTK LIBS_VTK]]);
                else
                end
            end
            if(isunix)
                eval(['mex ' [COMP_FLAGS MIATOOLS_MEXOPTS INCLUDE_DIRS_MIATOOLS INCLUDE_DIRS_VTK ' vtkHull.cpp ' INCLUDE_LIBS_VTK LIBS_VTK]]);
                !chmod ug+rwx vtkHull.mexa64
            end
        case {'vtkPolyDataWriter'}
            if(ispc)
                if(~isarch64)
                    eval(['mex ' [COMP_FLAGS MIATOOLS_MEXOPTS INCLUDE_DIRS_MIATOOLS INCLUDE_DIRS_VTK ' vtkPolyDataWriter.cpp ' INCLUDE_LIBS_VTK LIBS_VTK]]);
                else
                end
            end
            if(isunix)
                eval(['mex ' [COMP_FLAGS MIATOOLS_MEXOPTS INCLUDE_DIRS_MIATOOLS INCLUDE_DIRS_VTK ' vtkPolyDataWriter.cpp ' INCLUDE_LIBS_VTK LIBS_VTK]]);
                !chmod ug+rwx vtkPolyDataWriter.mexa64
            end
        case {'vtkCleanPolyData'}
            if(ispc)
                if(~isarch64)
                    eval(['mex ' [COMP_FLAGS MIATOOLS_MEXOPTS INCLUDE_DIRS_MIATOOLS INCLUDE_DIRS_VTK ' vtkCleanPolyData.cpp ' INCLUDE_LIBS_VTK LIBS_VTK]]);
                else
                end
            end
            if(isunix)
                eval(['mex ' [COMP_FLAGS MIATOOLS_MEXOPTS INCLUDE_DIRS_MIATOOLS INCLUDE_DIRS_VTK ' vtkCleanPolyData.cpp ' INCLUDE_LIBS_VTK LIBS_VTK]]);
                !chmod ug+rwx vtkCleanPolyData.mexa64
            end
        case {'vtkIntersectionPolyDataFilter'}
            if(ispc)
                if(~isarch64)
                    eval(['mex ' [COMP_FLAGS MIATOOLS_MEXOPTS INCLUDE_DIRS_MIATOOLS INCLUDE_DIRS_VTK ' vtkIntersectionPolyDataFilter.cpp ' INCLUDE_LIBS_VTK LIBS_VTK]]);
                else
                end
            end
            if(isunix)
                eval(['mex ' [COMP_FLAGS MIATOOLS_MEXOPTS INCLUDE_DIRS_MIATOOLS INCLUDE_DIRS_VTK ' vtkIntersectionPolyDataFilter.cpp ' INCLUDE_LIBS_VTK LIBS_VTK]]);
                !chmod ug+rwx vtkIntersectionPolyDataFilter.mexa64
            end
        case {'vtkQuadricDecimation'}
            if(ispc)
                if(~isarch64)
                    eval(['mex ' [COMP_FLAGS MIATOOLS_MEXOPTS INCLUDE_DIRS_MIATOOLS INCLUDE_DIRS_VTK ' vtkIntersectionPolyDataFilter.cpp ' INCLUDE_LIBS_VTK LIBS_VTK]]);
                else
                end
            end
            if(isunix)
                eval(['mex ' [COMP_FLAGS MIATOOLS_MEXOPTS INCLUDE_DIRS_MIATOOLS INCLUDE_DIRS_VTK ' vtkIntersectionPolyDataFilter.cpp ' INCLUDE_LIBS_VTK LIBS_VTK]]);
                !chmod ug+rwx vtkQuadricDecimation.mexa64
            end
        case {'vtkClipPolyData'}
            if(ispc)
                if(~isarch64)
                    eval(['mex ' [COMP_FLAGS MIATOOLS_MEXOPTS INCLUDE_DIRS_MIATOOLS INCLUDE_DIRS_VTK ' vtkClipPolyData.cpp ' INCLUDE_LIBS_VTK LIBS_VTK]]);
                else
                end
            end
            if(isunix)
                eval(['mex ' [COMP_FLAGS MIATOOLS_MEXOPTS INCLUDE_DIRS_MIATOOLS INCLUDE_DIRS_VTK ' vtkClipPolyData.cpp ' INCLUDE_LIBS_VTK LIBS_VTK]]);
                !chmod ug+rwx vtkClipPolyData.mexa64
            end
        case {'vtkLinearSubdivisionFilter'}
            if(ispc)
                if(~isarch64)
                    eval(['mex ' [COMP_FLAGS MIATOOLS_MEXOPTS INCLUDE_DIRS_MIATOOLS INCLUDE_DIRS_VTK ' vtkLinearSubdivisionFilter.cpp ' INCLUDE_LIBS_VTK LIBS_VTK]]);
                else
                end
            end
            if(isunix)
                eval(['mex ' [COMP_FLAGS MIATOOLS_MEXOPTS INCLUDE_DIRS_MIATOOLS INCLUDE_DIRS_VTK ' vtkLinearSubdivisionFilter.cpp ' INCLUDE_LIBS_VTK LIBS_VTK]]);
                !chmod ug+rwx vtkLinearSubdivisionFilter.mexa64
            end
        case {'vtkShrinkPolyData'}
            if(ispc)
                if(~isarch64)
                    eval(['mex ' [COMP_FLAGS MIATOOLS_MEXOPTS INCLUDE_DIRS_MIATOOLS INCLUDE_DIRS_VTK ' vtkShrinkPolyData.cpp ' INCLUDE_LIBS_VTK LIBS_VTK]]);
                else
                end
            end
            if(isunix)
                eval(['mex ' [COMP_FLAGS MIATOOLS_MEXOPTS INCLUDE_DIRS_MIATOOLS INCLUDE_DIRS_VTK ' vtkShrinkPolyData.cpp ' INCLUDE_LIBS_VTK LIBS_VTK]]);
                !chmod ug+rwx vtkShrinkPolyData.mexa64
            end
        case {'vtkClosestElement'}
            if(ispc)
                if(~isarch64)
                    eval(['mex ' [COMP_FLAGS MIATOOLS_MEXOPTS INCLUDE_DIRS_MIATOOLS INCLUDE_DIRS_VTK ' vtkClosestElement.cpp ' INCLUDE_LIBS_VTK LIBS_VTK]]);
                else
                end
            end
            if(isunix)
                eval(['mex ' [COMP_FLAGS MIATOOLS_MEXOPTS INCLUDE_DIRS_MIATOOLS INCLUDE_DIRS_VTK ' vtkClosestElement.cpp ' INCLUDE_LIBS_VTK LIBS_VTK]]);
                !chmod ug+rwx vtkClosestElement.mexa64
            end
        case {'vtkMassProperties'}
            if(ispc)
                if(~isarch64)
                    eval(['mex ' [COMP_FLAGS MIATOOLS_MEXOPTS INCLUDE_DIRS_MIATOOLS INCLUDE_DIRS_VTK ' vtkMassProperties.cpp ' INCLUDE_LIBS_VTK LIBS_VTK]]);
                else
                end
            end
            if(isunix)
                eval(['mex ' [COMP_FLAGS MIATOOLS_MEXOPTS INCLUDE_DIRS_MIATOOLS INCLUDE_DIRS_VTK ' vtkMassProperties.cpp ' INCLUDE_LIBS_VTK LIBS_VTK]]);
                !chmod ug+rwx vtkMassProperties.mexa64
            end
        case {'vtkSmoothPolyDataFilter'}
            if(ispc)
                if(~isarch64)
                    eval(['mex ' [COMP_FLAGS MIATOOLS_MEXOPTS INCLUDE_DIRS_MIATOOLS INCLUDE_DIRS_VTK ' vtkSmoothPolyDataFilter.cpp ' INCLUDE_LIBS_VTK LIBS_VTK]]);
                else
                end
            end
            if(isunix)
                eval(['mex ' [COMP_FLAGS MIATOOLS_MEXOPTS INCLUDE_DIRS_MIATOOLS INCLUDE_DIRS_VTK ' vtkSmoothPolyDataFilter.cpp ' INCLUDE_LIBS_VTK LIBS_VTK]]);
                !chmod ug+rwx vtkSmoothPolyDataFilter.mexa64
            end
        case {'vtkClosestPoint'}
            if(ispc)
                if(~isarch64)
                    eval(['mex ' [COMP_FLAGS MIATOOLS_MEXOPTS INCLUDE_DIRS_MIATOOLS INCLUDE_DIRS_VTK ' vtkClosestPoint.cpp ' INCLUDE_LIBS_VTK LIBS_VTK]]);
                else
                end
            end
            if(isunix)
                eval(['mex ' [COMP_FLAGS MIATOOLS_MEXOPTS INCLUDE_DIRS_MIATOOLS INCLUDE_DIRS_VTK ' vtkClosestPoint.cpp ' INCLUDE_LIBS_VTK LIBS_VTK]]);
                !chmod ug+rwx vtkClosestPoint.mexa64
            end
        case {'vtkWindowedSincPolyDataFilter'}
            if(ispc)
                if(~isarch64)
                    eval(['mex ' [COMP_FLAGS MIATOOLS_MEXOPTS INCLUDE_DIRS_MIATOOLS INCLUDE_DIRS_VTK ' vtkWindowedSincPolyDataFilter.cpp ' INCLUDE_LIBS_VTK LIBS_VTK]]);
                else
                end
            end
            if(isunix)
                eval(['mex ' [COMP_FLAGS MIATOOLS_MEXOPTS INCLUDE_DIRS_MIATOOLS INCLUDE_DIRS_VTK ' vtkWindowedSincPolyDataFilter.cpp ' INCLUDE_LIBS_VTK LIBS_VTK]]);
                !chmod ug+rwx vtkWindowedSincPolyDataFilter.mexa64
            end
        case {'vtkDecimatePro'}
            if(ispc)
                if(~isarch64)
                    eval(['mex ' [COMP_FLAGS MIATOOLS_MEXOPTS INCLUDE_DIRS_MIATOOLS INCLUDE_DIRS_VTK ' vtkDecimatePro.cpp ' INCLUDE_LIBS_VTK LIBS_VTK]]);
                else
                end
            end
            if(isunix)
                eval(['mex ' [COMP_FLAGS MIATOOLS_MEXOPTS INCLUDE_DIRS_MIATOOLS INCLUDE_DIRS_VTK ' vtkDecimatePro.cpp ' INCLUDE_LIBS_VTK LIBS_VTK]]);
                !chmod ug+rwx vtkDecimatePro.mexa64
            end
        case {'vtkPolyDataConnectivityFilter'}
            if(ispc)
                if(~isarch64)
                    eval(['mex ' [COMP_FLAGS MIATOOLS_MEXOPTS INCLUDE_DIRS_MIATOOLS INCLUDE_DIRS_VTK ' vtkPolyDataConnectivityFilter.cpp ' INCLUDE_LIBS_VTK LIBS_VTK]]);
                else
                end
            end
            if(isunix)
                eval(['mex ' [COMP_FLAGS MIATOOLS_MEXOPTS INCLUDE_DIRS_MIATOOLS INCLUDE_DIRS_VTK ' vtkPolyDataConnectivityFilter.cpp ' INCLUDE_LIBS_VTK LIBS_VTK]]);
                !chmod ug+rwx vtkPolyDataConnectivityFilter.mexa64
            end
            case {'vtkDensifyPolyData'}
            if(ispc)
                if(~isarch64)
                    eval(['mex ' [COMP_FLAGS MIATOOLS_MEXOPTS INCLUDE_DIRS_MIATOOLS INCLUDE_DIRS_VTK ' vtkDensifyPolyData.cpp ' INCLUDE_LIBS_VTK LIBS_VTK]]);
                else
                end
            end
            if(isunix)
                eval(['mex ' [COMP_FLAGS MIATOOLS_MEXOPTS INCLUDE_DIRS_MIATOOLS INCLUDE_DIRS_VTK ' vtkDensifyPolyData.cpp ' INCLUDE_LIBS_VTK LIBS_VTK]]);
                !chmod ug+rwx vtkDensifyPolyData.mexa64
            end
        case {'vtkPolyDataNormals'}
            if(ispc)
                if(~isarch64)
                    eval(['mex ' [COMP_FLAGS MIATOOLS_MEXOPTS INCLUDE_DIRS_MIATOOLS INCLUDE_DIRS_VTK ' vtkPolyDataNormals.cpp ' INCLUDE_LIBS_VTK LIBS_VTK]]);
                else
                end
            end
            if(isunix)
                eval(['mex ' [COMP_FLAGS MIATOOLS_MEXOPTS INCLUDE_DIRS_MIATOOLS INCLUDE_DIRS_VTK ' vtkPolyDataNormals.cpp ' INCLUDE_LIBS_VTK LIBS_VTK]]);
                !chmod ug+rwx vtkPolyDataNormals.mexa64
            end
        case {'vtkContour'}
            if(ispc)
                if(~isarch64)
                    eval(['mex ' [COMP_FLAGS MIATOOLS_MEXOPTS INCLUDE_DIRS_MIATOOLS INCLUDE_DIRS_VTK ' vtkContour.cpp ' INCLUDE_LIBS_VTK LIBS_VTK]]);
                else
                end
            end
            if(isunix)
                eval(['mex ' [COMP_FLAGS MIATOOLS_MEXOPTS INCLUDE_DIRS_MIATOOLS INCLUDE_DIRS_VTK ' vtkContour.cpp ' INCLUDE_LIBS_VTK LIBS_VTK]]);
                !chmod ug+rwx vtkContour.mexa64
            end
        case {'boundary3d'}
            if(ispc)
                if(~isarch64)
                    eval(['mex ' [COMP_FLAGS MIATOOLS_MEXOPTS INCLUDE_DIRS_MIATOOLS INCLUDE_DIRS_VTK ' boundary3d.cpp ' INCLUDE_LIBS_VTK LIBS_VTK]]);
                else
                end
            end
            if(isunix)
                eval(['mex ' [COMP_FLAGS MIATOOLS_MEXOPTS INCLUDE_DIRS_MIATOOLS INCLUDE_DIRS_VTK ' boundary3d.cpp ' INCLUDE_LIBS_VTK LIBS_VTK]]);
                !chmod ug+rwx boundary3d.mexa64
            end
        case {'DeletePoints'}
            if(ispc)
                if(~isarch64)
                    eval(['mex ' [COMP_FLAGS MIATOOLS_MEXOPTS INCLUDE_DIRS_MIATOOLS INCLUDE_DIRS_VTK ' DeletePoints.c ' INCLUDE_LIBS_VTK LIBS_VTK]]);
                else
                end
            end
            if(isunix)
                eval(['mex ' [COMP_FLAGS MIATOOLS_MEXOPTS INCLUDE_DIRS_MIATOOLS INCLUDE_DIRS_VTK ' DeletePoints.c ' INCLUDE_LIBS_VTK LIBS_VTK]]);
                !chmod ug+rwx DeletePoints.mexa64
            end
        case {'IsInsidePolyData'}
            if(ispc)
                if(~isarch64)
                    eval(['mex ' [COMP_FLAGS MIATOOLS_MEXOPTS INCLUDE_DIRS_MIATOOLS INCLUDE_DIRS_VTK ' IsInsidePolyData.cpp ' INCLUDE_LIBS_VTK LIBS_VTK]]);
                else
                end
            end
            if(isunix)
                eval(['mex ' [COMP_FLAGS MIATOOLS_MEXOPTS INCLUDE_DIRS_MIATOOLS INCLUDE_DIRS_VTK ' IsInsidePolyData.cpp ' INCLUDE_LIBS_VTK LIBS_VTK]]);
                !chmod ug+rwx IsInsidePolyData.mexa64
            end
            
        case {'all_vtk'}
            make vtkFeatureEdges;
            make vtkPolyDataPointSampler
            make vtkBooleanOperationPolyDataFilter
            make vtkFillHolesFilter
            make vtkPolyDataReader
            make vtkButterflySubdivisionFilter
            make vtkHull
            make vtkPolyDataWriter
            make vtkCleanPolyData
            make vtkIntersectionPolyDataFilter
            make vtkQuadricDecimation
            make vtkClipPolyData
            make vtkLinearSubdivisionFilter
            make vtkShrinkPolyData
            make vtkClosestElement
            make vtkMassProperties
            make vtkSmoothPolyDataFilter
            make vtkClosestPoint
            make vtkWindowedSincPolyDataFilter
            make vtkDecimatePro
            make vtkPolyDataConnectivityFilter
            make vtkDensifyPolyData
            make vtkPolyDataNormals
            make vtkContour
            make boundary3d
            make DeletePoints
            make IsInsidePolyData
            
                    
        case {'all'}
            for p = 1:numel( allProjects)
               make(  allProjects{p} );
            end
                    
        otherwise
            fprintf('Valid projects are:\n\n');
            cellfun( @(p)fprintf('make %s\n',p) , allProjects );
            fprintf('\n\n');
            error('invalid projectName.');
            
            
    end

end