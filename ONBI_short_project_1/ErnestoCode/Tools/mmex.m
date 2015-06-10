function OK_ = mmex( mname , varargin )
% 
% mmex file.mmex
%         -> compile the project defined in file.mmex
% mmex file.mmex [options]
%         -> compile with options
% mmex file1.mmex file2.mmex ...
%         -> compile projects in files file1.mmex, file2.mmex, ...
% mmex dirname/file.mmex
%         -> go to dirname, compile project file.mmex, and came back to the current dir
% mmex dirname
%         -> list mmex files in dirname
% mmex dirname -R
%         -> list recursively mmex files in dirname
% mmex dirname -all
%         -> make all mmex files in dirname
% mmex .
%         -> list mmexs in current dir
% mmex . -R
%         -> list recursively mmexs in current dir
% mmex . -all
%         -> make all mmex in current dir
% mmex . -R -all
%         -> make all mmex in current dir and below
% mmex . -ALL
%         -> same as below
%

  VERBOSE = true;

  if isfile( mname )
    vprintf( '\n')
    vprintf( '================================================================================================\n')
    vprintf( '**MEXING  ''%s''\n', mname );
    vprintf( '   ---------------------------------------------------------------------------------------------\n')

    
    OK = false;
    try
      
      OK = true;
    end
    
    
    vprintf( '\n')
    vprintf( '\n')
    vprintf( '   ---------------------------------------------------------------------------------------------\n')
    if OK
    vprintf( '**MEXING ''%s'' OK......\n', mname);
    else
    vprintf( '**SOME ERROR in MEXING  ''%s''......\n', mname );
    end
    vprintf( '================================================================================================\n')
    vprintf( '\n')
    vprintf( '\n')
    
    if nargout, OK_ = OK; end
    return;
  end

  if isdir( mname )
    if isequal( mname , '.' ), mname = ''; end
    
    RECURSIVE  = false;
    compileALL = false;
    while numel( varargin )
      if isequal( varargin{1} , '-R' )
        RECURSIVE = true;
        varargin(1) = []; continue;
      end
      if isequal( varargin{1} , '-all' )
        compileALL = true;
        varargin(1) = []; continue;
      end
      if isequal( varargin{1} , '-ALL' )
        RECURSIVE = true;
        compileALL = true;
        varargin(1) = []; continue;
      end
      break;
    end
    
    if RECURSIVE
      files = rdir( mname , '.*\.mmex$', Inf );
    else
      files = rdir( mname , '.*\.mmex$', 0   );
    end
    errorsOnCompile = {};
    for f = 1:numel(files)
      if compileALL
        OK = false;
        try
          OK = mmex(files(f).name,varargin{:});
        end
        if ~OK
          errorsOnCompile{end+1} = files(f).name;
        end
      else
        fprintf('mmex ''%s''\n',files(f).name);
      end
    end
    if numel( errorsOnCompile )
      fprintf('**********************************************\n')
      fprintf('**ERRORS in:\n')
      for i = 1:numel( errorsOnCompile )
      fprintf('%s\n',errorsOnCompile{i});
      end
      fprintf('**********************************************\n')
    end
    
    return;
  end
  
  [p,f,e] = fileparts( mname );
  if isempty(e) && isfile( [ mname , '.mmex' ] )
    try
      mmex( [ mname , '.mmex' ] , varargin{:} );
    catch LE
      throwAsCalles(LE)
    end
    return;
  end
  
  
  function vprintf(varargin)
    if VERBOSE, fprintf(varargin{:}); end
  end
  
end


% 
% 
% 
% 
% ( fname , varargin )
% 
%   [p,n,e] = fileparts(fname);
%   if isempty(e)
%     fname = [ fname '.c' ];
%   end
%   
%   fullfname = which(fname);
%   if isempty(fullfname)
%     fprintf('no file  %s\n', fname );
%     return;
%   end
% 
%   fprintf( 'mexing  %s\n', fullfname );
%   [p,fn,e] = fileparts(fullfname);
% %   cwd = pwd;
% 
%   switch computer
%     case {'PCWIN','PCWIN64'}
%       myMEXdir = [ fileparts(fileparts(which('mmex')))  filesep 'Tools' ];
%       outdir   = [ fileparts(fileparts(which('mmex')))  filesep 'MEXs' ];
% %       myMEXdir = 'g:\Work\MEXs';
% %       outdir   = 'g:\Work\MEXs';
% 
%       copyfile( fullfname , [outdir filesep fn e] , 'f' );
%       
%       AT_END_DELETE = [outdir filesep fn e];
%       
%     case {'GLNX86','GLNXA64'}
%       myMEXdir = [ fileparts(fileparts(which('mmex')))  filesep 'Tools' ];
% 
%       [a,hn] = system('echo $HOSTNAME' );
%       hn( find(hn=='.',1,'first'):end ) = [];
%       if any( strcmp( hn , {'apolo','selene','arcadia'}) )
%         outdir   = [ fileparts(fileparts(which('mmex')))  filesep 'MEXs/SAA' ];
%       else
%         outdir   = [ fileparts(fileparts(which('mmex')))  filesep 'MEXs' ];
%       end
%       
%       fullfname = removeComments( fullfname );
%       AT_END_DELETE = fullfname;
%       
%       varargin = [ varargin , '-output' , fn ];
% 
%   end
%   
%   cwd = pwd;
%   cd(outdir);
%   try
%   switch computer
%     case {'PCWIN','PCWIN64'}
%       mex( '-O' , ['-I' myMEXdir] , varargin{:} , fullfname );
%     case {'GLNX86','GLNXA64'}
%       mex( '-O' , ['-I' myMEXdir] , varargin{:} , fullfname , ' -lut');
%   end
%   catch
%     le = lasterror;
%     disp( le.message );
%   end
%   cd( cwd );
% 
%   if ~isempty( AT_END_DELETE ), delete( AT_END_DELETE ); end
% 
%   function fnn = removeComments( fn )
%     [p,fnn,e] = fileparts( fn );
%     fnn = [p filesep 'no_comments.' fnn e ];
%     
%     fidr = fopen( fn , 'r' );
%     fidw = fopen( fnn , 'w' );
%     while 1
%       l = fgetl( fidr );
%       if ~ischar( l ), break; end
%       
%       l = regexprep( l , '^\s*\<_inline\>' , '__inline__' );
%       
%       pos = findstr(l,'//');
%       if ~isempty(pos)
%         l = l(1:pos(1)-1);
%       end
%       
%       fprintf(fidw,'%s\n',l);
%     end
%     fclose(fidr);
%     fclose(fidw);
%   end
% 
% end
% 
% %{
% mmex jacobi3x3.c
% mmex val2ind.c
% mmex Interp1D.c
% mmex getInterval.c
% mmex inpoly.c
% mmex boundary.c
% 
% mmex MOS_GrayImage.c
% 
% mmex toboggan_image.c
% mmex follow_image_dijkstra.c
% mmex image_dijkstra.c
% mmex ApplyContrastFunction.c
% mmex Interp3DGridOn3DGrid.c
% mmex Interp3DGridOn3DGrid_single.c
% mmex InterpPointsOn3DGrid.c
% mmex InterpPointsOn3DGrid_single.c
% 
% mmex Euler4D_constant.c
% mmex Euler4D_constant_single.c
% mmex EulerBilinear_single.c
% mmex EulerTrilinear_single.c
% mmex Euler4D.c
% mmex Euler4D_single.c
% mmex EulerBilinear.c
% mmex EulerTrilinear.c
% mmex ExponentialLinear.c
% 
% mmex Solve3x3Sys.c
% mmex Solve3x3Sys_single.c
% mmex Solve2x2Sys.c
% mmex Solve2x2Sys_single.c
% mmex Solve2x2SymSys.c
% mmex Solve2x2SymSys_single.c
% mmex Solve3x3SymSys.c
% mmex Solve3x3SymSys_single.c
% %}
% 
% 
% 
% %{
% @echo off
% rem MSVC80OPTS.BAT
% rem
% rem    Compile and link options used for building MEX-files
% rem    using the Microsoft Visual C++ compiler version 8.0
% rem
% rem StorageVersion: 1.0
% rem C++keyFileName: MSVC80OPTS.BAT
% rem C++keyName: Microsoft Visual C++ 2005
% rem C++keyManufacturer: Microsoft
% rem C++keyVersion: 8.0
% rem C++keyLanguage: C++
% rem
% rem    $Revision: 1.1.10.6 $  $Date: 2007/11/07 17:44:06 $
% rem    Copyright 1984-2007 The MathWorks, Inc.
% rem
% rem ********************************************************************
% rem General parameters
% rem ********************************************************************
% 
% set MATLAB=%MATLAB%
% set VS80COMNTOOLS=%VS80COMNTOOLS%
% set VSINSTALLDIR=C:\Win\Visual Studio
% set VCINSTALLDIR=%VSINSTALLDIR%\VC
% set PATH=%VCINSTALLDIR%\BIN\;%VCINSTALLDIR%\PlatformSDK\bin;%VSINSTALLDIR%\Common7\IDE;%VSINSTALLDIR%\SDK\v2.0\bin;%VSINSTALLDIR%\Common7\Tools;%VSINSTALLDIR%\Common7\Tools\bin;%VCINSTALLDIR%\VCPackages;%MATLAB_BIN%;%PATH%
% set INCLUDE=%VCINSTALLDIR%\ATLMFC\INCLUDE;%VCINSTALLDIR%\INCLUDE;%VCINSTALLDIR%\PlatformSDK\INCLUDE;%VSINSTALLDIR%\SDK\v2.0\include;%INCLUDE%
% set LIB=%VCINSTALLDIR%\ATLMFC\LIB;%VCINSTALLDIR%\LIB;%VCINSTALLDIR%\PlatformSDK\lib;%VSINSTALLDIR%\SDK\v2.0\lib;%MATLAB%\extern\lib\win32;%LIB%
% set MW_TARGET_ARCH=win32
% 
% rem ********************************************************************
% rem Compiler parameters
% rem ********************************************************************
% set COMPILER=cl
% set COMPFLAGS= /c /Zp8 /GR /W3 /EHa- /D_CRT_SECURE_NO_DEPRECATE /D_SCL_SECURE_NO_DEPRECATE /D_SECURE_SCL=0 /DMATLAB_MEX_FILE /nologo /MD
% rem set OPTIMFLAGS=/O2  /DNDEBUG
% set OPTIMFLAGS= /O2 /DNDEBUG  /Oy /GL /Oi /Ot /Ox /fp:fast
% set DEBUGFLAGS=/Zi /Fd"%OUTDIR%%MEX_NAME%%MEX_EXT%.pdb"
% set NAME_OBJECT=/Fo
% 
% rem ********************************************************************
% rem Linker parameters
% rem ********************************************************************
% set LIBLOC=%MATLAB%\extern\lib\win32\microsoft
% set LINKER=link
% set LINKFLAGS=/dll /export:%ENTRYPOINT% /MAP /LIBPATH:"%LIBLOC%" libmx.lib libmex.lib libmat.lib /implib:%LIB_NAME%.x /MACHINE:X86 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib
% set LINKOPTIMFLAGS=
% set LINKDEBUGFLAGS=/DEBUG /PDB:"%OUTDIR%%MEX_NAME%%MEX_EXT%.pdb"
% set LINK_FILE=
% set LINK_LIB=
% set NAME_OUTPUT=/out:"%OUTDIR%%MEX_NAME%%MEX_EXT%"
% set RSP_FILE_INDICATOR=@
% 
% rem ********************************************************************
% rem Resource compiler parameters
% rem ********************************************************************
% set RC_COMPILER=rc /fo "%OUTDIR%mexversion.res"
% set RC_LINKER=
% 
% set POSTLINK_CMDS=del "%OUTDIR%%MEX_NAME%.map"
% set POSTLINK_CMDS1=del %LIB_NAME%.x
% set POSTLINK_CMDS2=mt -outputresource:"%OUTDIR%%MEX_NAME%%MEX_EXT%";2 -manifest "%OUTDIR%%MEX_NAME%%MEX_EXT%.manifest"
% set POSTLINK_CMDS3=del "%OUTDIR%%MEX_NAME%%MEX_EXT%.manifest" 
% %}
% 
% 
% 
% 
% 
% 
% %{
% 
% 
% 
% 
% 
% 
% %}
