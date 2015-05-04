function psfrag( inFile , outFile , frags , margins )
% left bottom right top

  if nargin < 2, error('inFile outFile were expected'); end
  if nargin < 3 || isempty( frags ), frags = {}; end
  if nargin < 4, margins = [-1 -1 1 1]; end
  if numel(margins) ~= 4, error('margins se esperaba de tamano 4.'); end

  if ~isfile( inFile )
    error('no existe inFile:  %s', inFile );
  end
  if isempty(fileparts(inFile))
    inFile = which(inFile);
  end
  
  DIR = tmpname('matlabfrag_');
  mkdir(DIR);
  
  cwd = pwd;
  cd(DIR);
  
  CLEANUPc = onCleanup( @()cd(cwd) ); 
  CLEANUPa = onCleanup( @()delete(fullfile(DIR,'*.*')) ); 
  CLEANUPb = onCleanup( @()rmdir(DIR) ); 

  
  [inP,inF,inE] = fileparts( inFile );
  if strcmp( inE , '.plt');
    inE = '.ps';
  end
  copyfile( inFile , fullfile( DIR , [ 'inFile' , inE ] ) );
  
fileprintf(fullfile(DIR,'a.tex'),'\documentclass[a3]{article}');
% fileprintf(fullfile(DIR,'a.tex'),'\usepackage[a3paper]{geometry}');
fileprintf(fullfile(DIR,'a.tex'),'\usepackage[dvips]{graphicx}');
fileprintf(fullfile(DIR,'a.tex'),'\usepackage{psfrag}');
fileprintf(fullfile(DIR,'a.tex'),'\usepackage{amscd,amsmath}');
fileprintf(fullfile(DIR,'a.tex'),'\usepackage{amssymb}');
fileprintf(fullfile(DIR,'a.tex'),'\usepackage{mathabx}');
fileprintf(fullfile(DIR,'a.tex'),'\usepackage{color}');
fileprintf(fullfile(DIR,'a.tex'),'\usepackage[T1]{fontenc}');
fileprintf(fullfile(DIR,'a.tex'),'\usepackage{anyfontsize}');

fileprintf(fullfile(DIR,'a.tex'),'\usepackage{xparse}');
fileprintf(fullfile(DIR,'a.tex'),'\NewDocumentCommand{\safecommand}{ m O{0} m }{');
fileprintf(fullfile(DIR,'a.tex'),'\providecommand{#1}{}');
fileprintf(fullfile(DIR,'a.tex'),'\renewcommand{#1}[#2]{#3}');
fileprintf(fullfile(DIR,'a.tex'),'}');


fileprintf(fullfile(DIR,'a.tex'),'\DeclareFontFamily{OT1}{pzc}{}');
fileprintf(fullfile(DIR,'a.tex'),'\DeclareFontShape{OT1}{pzc}{m}{it}{<-> s * [1.10] pzcmi7t}{}');
fileprintf(fullfile(DIR,'a.tex'),'\DeclareMathAlphabet{\mathpzc}{OT1}{pzc}{m}{it}');

fileprintf(fullfile(DIR,'a.tex'),'\setlength{\hoffset}{-1in}');
fileprintf(fullfile(DIR,'a.tex'),'\setlength{\voffset}{-1in}');
fileprintf(fullfile(DIR,'a.tex'),'\setlength{\oddsidemargin}{5mm}');
fileprintf(fullfile(DIR,'a.tex'),'\setlength{\textwidth}{200mm}');

fileprintf(fullfile(DIR,'a.tex'),'\pagestyle{empty}');

fileprintf(fullfile(DIR,'a.tex'),'\makeatletter');
fileprintf(fullfile(DIR,'a.tex'),'\def\mathoperator#1{\mathop{\operator@font #1}\nolimits}');
fileprintf(fullfile(DIR,'a.tex'),'\makeatother');


if iscell( frags ) && numel( frags ) == 2 && iscell( frags{1} ) && iscell( frags{2} )

fileprintf(fullfile(DIR,'a.tex'),' ');
for i = 1:numel(frags{1})
fileprintf(fullfile(DIR,'a.tex'),frags{1}{i});
end
fileprintf(fullfile(DIR,'a.tex'),' ');
  
frags = frags{2};
end


fileprintf(fullfile(DIR,'a.tex'),'\begin{document}');
fileprintf(fullfile(DIR,'a.tex'),'\begin{figure}');

fileprintf(fullfile(DIR,'a.tex'),' ');
for i = 1:numel(frags)
fileprintf(fullfile(DIR,'a.tex'),frags{i});
end
fileprintf(fullfile(DIR,'a.tex'),' ');

% fileprintf(fullfile(DIR,'a.tex'),'\\includegraphics[clip=true, scale=1]{%s}', [ 'inFile' , inE ] );
%fileprintf(fullfile(DIR,'a.tex'),'\\includegraphics[clip=true, width=\\textwidth]{%s}', [ 'inFile' , inE ] );

fileprintf(fullfile(DIR,'a.tex'),'\\includegraphics[clip=true,width=\\textwidth,height=\\textheight,keepaspectratio]{%s}', [ 'inFile' , inE ] );

fileprintf(fullfile(DIR,'a.tex'),'\end{figure}');
fileprintf(fullfile(DIR,'a.tex'),'\end{document}');


fprintf('\nLATEXING .');

  [status,result] = system( 'C:\\Win\\TEX\\MiKTeX\\miktex\\bin\\latex.exe -quiet a.tex' , '-echo');
fprintf('.');
  if status, keyboard; end
  [status,result] = system( 'C:\Win\TEX\MiKTeX\miktex\bin\dvips.exe a.dvi' );
fprintf('.');

%   %crop the PS
%   [status,result] = system( '"C:\Archivos de programa\Accesorios\Postscript\gs9.07\bin\gswin32c.exe" -dBATCH -dNOPAUSE -sDEVICE=bbox -q -r1440x1440 a.ps' );
%   result = sscanf( result( find(result==':',1)+1:find(result==10,1) ) , '%d' );
%   bb     = result(:) + [-1;-1;1;1];

  
  %crop the PS
  [status,result] = system( '"C:\Archivos de programa\Accesorios\Postscript\gs9.07\bin\gswin32c.exe" -dBATCH -dNOPAUSE -sDEVICE=pcxmono -q -r1440x1440 -sOutputFile="a.pcx" "a.ps"' );
fprintf('.');
  I = imread('a.pcx');
  bb = [0 0 0 0];
  bb([1 3]) = range( find(         ~all(I,1)   ) );
  bb([2 4]) = range( find( flipud( ~all(I,2) ) ) );
  bb = bb/20;
  bb = [ floor(bb(1:2)) , ceil(bb(3:4)) ] + margins(:).';
fprintf('.');
  
  POSTSCRIPT = readFile( 'a.ps' );
fprintf('.');
  
%   POSTSCRIPT = POSTSCRIPT( cellfun('isempty' , regexp( POSTSCRIPT , '^\%', 'once' ) ) );
%   POSTSCRIPT = POSTSCRIPT( cellfun('isempty' , regexp( POSTSCRIPT , '^\%DVIPS', 'once' ) ) );
  
  fid = fopen('a.ps','w');
  fprintf(fid,'%%!PS-Adobe-3.0 EPSF-3.0\n');
  fprintf(fid,'%%%%BoundingBox: %d %d %d %d\n', bb );
  cellfun( @(s) fprintf(fid,'%s\n',s) , POSTSCRIPT );
  fclose(fid);
fprintf('.');
  %end crop the PS
  
  [outP,outF,outE] = fileparts( outFile );
fprintf('.');
  
  switch lower( outE )
    case '.pdf'
%       [status,result] = system( '"C:\Archivos de programa\Accesorios\Postscript\gs9.07\bin\gswin32c.exe" -sDEVICE=pdfwrite -sOutputFile="a.pdf" -r720 -dBATCH -dNOPAUSE  -dEPSCrop -dCompatibilityLevel=1.0 "a.ps"' );
%       [status,result] = system( '"C:\Archivos de programa\Accesorios\Postscript\gs9.07\bin\gswin32c.exe"         -q -dNOPAUSE -dBATCH -dEPSCrop -dEmbedAllFonts=true -dUseFlateCompression=true -dHaveTrueTypes -dSubsetFonts=true -dNOPLATFONTS -dUseCIEColor=true -dColorConversionStrategy=/UseDeviceIndependentColor -dProcessColorModel=/DeviceRGB -dAntiAliasColorImages=false -dUseFlateCompression=true -dLZWEncodePages=true -dCompatibilityLevel=1.6 -dAutoFilterColorImages=false -dAutoFilterGrayImages=false  -dColorImageFilter=/FlateEncode -dGrayImageFilter=/FlateEncode     -sDEVICE=pdfwrite      -sOutputFile="a.pdf" -f "a.eps"' );

      [status,result] = system( '"C:\Archivos de programa\Accesorios\Postscript\gs9.07\bin\gswin32c.exe" -q -dNOPAUSE -dBATCH -dEPSCrop -dEmbedAllFonts=true -dUseFlateCompression=true -dAutoRotatePages=/None -dHaveTrueTypes -dSubsetFonts=true -dNOPLATFONTS -dUseCIEColor=true -dColorConversionStrategy=/UseDeviceIndependentColor -dProcessColorModel=/DeviceRGB -dAntiAliasColorImages=false -dAntiAliasGrayImages=false -dUseFlateCompression=true -dLZWEncodePages=true -dAutoFilterColorImages=false -dAutoFilterGrayImages=false  -dColorImageFilter=/FlateEncode -dGrayImageFilter=/FlateEncode -dCompatibilityLevel=1.2 -sDEVICE=pdfwrite -sOutputFile="a.pdf" -f "a.ps"');
fprintf('.');

      
      movefile( 'a.pdf' , 'a.out' );
fprintf('.');
  end


  cd(cwd)
  
  while 1
fprintf('.');
    try
      copyfile( fullfile(DIR,'a.out')  , outFile );
      break
    end
    [p,f,e] = fileparts(outFile);
    outFile = fullfile( p , [f,'_',e] );
  end


fprintf('  DONE\n');

end
