function [pV , vname ] = createWorkSpaceVariable( V , prefix )

  if nargin < 2
    prefix = 'var';
  end

  vars = evalin('base','who');
  id = 1;
  while 1
    vname = sprintf('%s__%04d', prefix , id );
    if ~any( strcmp( vars , vname ) )
      break;
    end
    id = id+1;
  end

  assignin( 'base' , vname , V );
  pV = pVAR( vname );

end
