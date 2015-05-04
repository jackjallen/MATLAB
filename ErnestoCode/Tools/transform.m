function [ p , d_tp ] = transform( p , varargin )

  [ varargin,i,formatpoints ]= parseargs( varargin,'Format','formatpoints','pointsformat','fp','pf','$DEFS$',false );
  
  if ~formatpoints, [ varargin,formatpoints ]= parseargs(varargin,'Rows'                        ,'$FORCE$','rows'       ); end
  if ~formatpoints, [ varargin,formatpoints ]= parseargs(varargin,'RowsH'                       ,'$FORCE$','rowsh'      ); end
  if ~formatpoints, [ varargin,formatpoints ]= parseargs(varargin,'r2d','rows2d'                ,'$FORCE$','rows2d'     ); end
  if ~formatpoints, [ varargin,formatpoints ]= parseargs(varargin,'r2dh','rows2dh'              ,'$FORCE$','rows2dh'    ); end
  if ~formatpoints, [ varargin,formatpoints ]= parseargs(varargin,'col','Columns'               ,'$FORCE$','columns'    ); end
  if ~formatpoints, [ varargin,formatpoints ]= parseargs(varargin,'colsh','ColumnsH'            ,'$FORCE$','columnsh'   ); end
  if ~formatpoints, [ varargin,formatpoints ]= parseargs(varargin,'c2d','col2d','columns2d'     ,'$FORCE$','columns2d'  ); end
  if ~formatpoints, [ varargin,formatpoints ]= parseargs(varargin,'c2dh','cols2dh','columns2dh' ,'$FORCE$','columns2dh' ); end

  if ~formatpoints, formatpoints= 'rows'; end
  
  if isempty(p)
      if nargout>1, d_tp=[]; end;
      return;
  end;
  
%   dH = zeros(16,0);
  dH = eye(16);
  
  if numel( varargin ) == 1
    
    H = varargin{1};
    if iscell( H )
      dH = H{2};
      H  = H{1};
    end

  else

    if nargout < 2
       H     = maketransform( varargin{:} );
    else
      [H,dH] = maketransform( varargin{:} );
    end

  end

  if isequal( H , eye(4) )   &&  nargout < 2
    return;
  end

  
  if size(H,1) < 4 || size(H,2) < 4
    H(4,4)= 1;
  end

  
  if strcmpi( formatpoints , 'rows' )
    if isequal( H(4,:) , [0 0 0 1] )

      if nargout > 1
%         d_tp = kron( [ p , ones(size(p,1),1) ] , eye(3,3) )*dH([1 2 3 5 6 7 9 10 11 13 14 15],:);
%         d_tp = d_tp( [ 1:3:end 2:3:end 3:3:end ] , : );
        
        d_tp = [ bsxfun( @plus , p * dH([ 1  5  9],:) , dH(13,:) ) ; ...
                 bsxfun( @plus , p * dH([ 2  6 10],:) , dH(14,:) ) ; ...
                 bsxfun( @plus , p * dH([ 3  7 11],:) , dH(15,:) ) ];
        
      end
      
      if ~isequal( H(1:3,1:3) , eye(3) )
        p = p*H(1:3,1:3).';
      end

      if ~isequal( H(1:3,4) , [0;0;0] )
        p = bsxfun( @plus , p , H(1:3,4).' );
      end

      return;

    end
  end
  
  
  if nargout > 1
    error('falta por implementar');
  end
  
  switch lower(formatpoints)
    case {'r','rows'}
      p= p(:,[1 2 3]).';
    case {'rh','rowsh'};
      p= p(:,[1 2 3 4]).';
      p= bsxfun( @rdivide , p ,  p(4,:) );
    case {'r2d','rows2d'};
      p= p(:,[1 2]).';
    case {'r2dh','rows2dh'};
      p([1 2 4],:)= p(:,[1 2 3]).';
      p= tp./repmat( p(4,:),4,1 );

    case {'c','col','columns'}
      p= p([1 2 3],:);
    case {'ch','colsh','columnsh'};
      p= p([1 2 3 4],:);
      p= tp./repmat( p(4,:),4,1 );
    case {'c2d','col2d','columns2d'};
      p= p([1 2],:);
    case {'c2dh','cols2dh','columns2dh'};
      p([1 2 4],:)= p([1 2 3],:);
      p= tp./repmat( p(4,:),4,1 );
  end
  p(4,:)=1;

  p= H*p;
  try
    p= bsxfun( @rdivide , p ,  p(4,:) );
  catch
    p= p(1:3,:) ./ repmat(p(4,:),3,1);
  end
  p(4,:)=1;


  switch lower(formatpoints)
    case {'r','rows'}
      p= p([1 2 3],:).';
    case {'rh','rowsh'};
      p= p.';
    case {'r2d','rows2d'};
      try
        if abs( p(3,:) ) < 1e-10
          p= p([1 2],:).';
        else
          p= p([1 2 3],:).';
        end
      catch
          p= p([1 2],:).';
      end
    case {'r2dh','rows2dh'};
      if abs( p(3,:) ) < 1e-10
        p= p([1 2 4],:).';
      else
        p= p.';
      end

    case {'c','col','columns'}
      p= p([1 2 3],:);
    case {'ch','colsh','columnsh'};

    case {'c2d','col2d','columns2d'};
      if abs( p(3,:) ) < 1e-10
        p= p([1 2],:);
      else
        p= p([1 2 3],:);
      end
    case {'c2dh','cols2dh','columns2dh'};
      if abs( p(3,:) ) < 1e-10
        p= p([1 2 4],:);
      else
        
      end
  end

end
