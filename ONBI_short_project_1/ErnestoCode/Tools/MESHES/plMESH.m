function plMESH( varargin )
%
% plMESH( mesh )
%

interp= 'E:\MoreWork\VTKViewer\bin\wish';
script= 'E:\MoreWork\VTKViewer\VTKViewer.tcl';

commando = sprintf( '"%s" "%s"', interp , script );

for i=1:nargin
  write_VTP( varargin{i} , sprintf( '%s\\%s.plMESH.vtp' , pwd , inputname(i) ) );
  commando = sprintf( '%s \"%s\\%s.plMESH.vtp\"', commando,  pwd , inputname(i) );

  if isfield( varargin{i} , 'lmk' )
    write_landmarks( sprintf( '%s\\%s.plMESH.lmk' , pwd , inputname(i) ) , varargin{i}.lmk );
  end

end

commando = sprintf( '%s &', commando);

system( commando ) ;
