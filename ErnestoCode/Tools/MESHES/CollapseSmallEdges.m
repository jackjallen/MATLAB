function M = CollapseSmallEdges( M , elen )

  if isfield( M , 'vertices' ) && ~isfield( M , 'xyz' )
    M.xyz = M.vertices;
    M = rmfield( M , 'vertices' );
  end
  if isfield( M , 'faces' ) && ~isfield( M , 'tri' )
    M.tri = M.faces;
    M = rmfield( M , 'faces' );
  end

  M.tri = sort( M.tri , 2 );
  M.tri( M.tri(:,1) == M.tri(:,2) , : ) = [];
  M.tri( M.tri(:,2) == M.tri(:,3) , : ) = [];


  Es = tri2edges( M.tri );
  
  fprintf('%d edges encontrados\n',size(Es,1));

  ALL_Es = Es;

  t0 = clock;
  B = vtkFeatureEdges( M , 'BoundaryEdgesOn',[],'FeatureEdgesOff',[],'NonManifoldEdgesOff',[],'ManifoldEdgesOff',[],'ColoringOff',[] );
  if  ~isempty( B )  &&  isfield( B, 'xyz' )  &&  numel( B.xyz ) > 0
    B = unique( vtkClosestPoint( M , B.xyz ) );
    B = ismembc( Es , B );
    Es = Es( ~xor( B(:,1) , B(:,2) ) , : );
  end

  fprintf('quitados los del borde, %d edges  ( %g secs )\n',size(Es,1),etime(clock,t0) );


  t0 = clock;

  Ls2 = sum( ( M.xyz( Es(:,1) , : ) - M.xyz( Es(:,2) , : ) ).^2 , 2 );
  
  Es( Ls2 > elen^2 , : ) = [];
  
  fprintf('despues de evaluar su longitud, quedan, %d edges  ( %g secs )\n',size(Es,1),etime(clock,t0) );
  

  ord = 1:size(Es,1);
  [ignore,ord] = sort( ord );
  ord = ord( 1:min(end,200000) );
  
  Es = Es( ord , : );

  fprintf('me quedo con unos pocos , %d\n'  , size(Es,1) );
  
  
  


  fprintf('empiezo con la seleccion_mx\n' );

  t0 = clock;

  Es = selec_edges_mx( Es );
  
  fprintf('termino ( %g secs )\n', etime(clock,t0) );
  
  
  Es( any( ~Es , 2 ) , : ) = [];
  
  fprintf('quedan %d edges\n', size(Es,1) );

  

  Es = Es(1:min(end,50000),:);
  fprintf('elegi  , %d\n'  , size(Es,1) );
  
  
  ALL_Es( ~any( ismembc( ALL_Es , unique(Es) ) , 2 ) , : ) = [];
  
  t0 = clock;
  fprintf('empiezo a evaluar los tripletes\n' );
  
  Ts = get_tripletes_mx( ALL_Es , Es ).';
  
  fprintf('ok...( %g secs )\n', etime(clock,t0) )
  
  t0 = clock;
  fprintf('ahora hago la diff de triangles\n');
  
  Ts = complex2tri( setdiff( tri2complex( Ts ) , tri2complex( M.tri ) ) );
  
  if ~isempty( Ts )
    Ts = tri2edges( Ts );
    Es = complex2edges( setdiff( edges2complex( Es ) , edges2complex( Ts ) ) );
  end

  
  fprintf('termino de evaluar tripletes ( %g secs ) y quedan ( %d ) edges \n', etime(clock,t0) , size( Es,1 ) );
  
  
  
  fprintf('finalmente limpiar mesh\n');
  
  if isempty( Es ), return; end
  
  centers = ( M.xyz( Es(:,1) , : ) + M.xyz( Es(:,2) , : ) )/2;

  M.xyz( Es(:,1) , : ) = centers;
%   M.xyz( Es(:,2) , : ) = centers;

  idx = 1:max(M.tri(:));
  idx( Es(:,2) ) = Es(:,1);
  
  M.tri = idx( M.tri );
  M.tri = sort( M.tri , 2 );

  M.tri( M.tri(:,1) == M.tri(:,2) , : ) = [];
  M.tri( M.tri(:,2) == M.tri(:,3) , : ) = [];
  M.tri( M.tri(:,1) == M.tri(:,3) , : ) = [];
  
%   M = vtkCleanPolyData( M , 'SetAbsoluteTolerance', 0 , 'PointMergingOn',[],'ConvertStripsToPolysOff',[],'ConvertPolysToLinesOff',[],'ConvertLinesToPointsOff',[],'ToleranceIsAbsoluteOn',[] );
%   if any( M.tri(:) == 0 )
%     M.tri( any( ~M.tri ,2) , : ) = [];
%     M = vtkCleanPolyData( M , 'SetAbsoluteTolerance', 0 , 'PointMergingOn',[],'ConvertStripsToPolysOn',[],'ConvertPolysToLinesOn',[],'ConvertLinesToPointsOn',[],'ToleranceIsAbsoluteOn',[] );
%   end  


  M = DeletePoints( M , setdiff( 1:size(M.xyz,1) , unique( M.tri ) ) );

  fprintf('done\n\n');
  
  
  
  function e = tri2edges( t )
    t = sort( t , 2 );
    e = [ t(:,[1 2]) ; t(:,[2 3]) ; t(:,[1 3]) ];

    ndx = 1:size(e,1);
    [ignore,ind] = sort( e( ndx , 2 ),'ascend'); ndx = ndx(ind);
    [ignore,ind] = sort( e( ndx , 1 ),'ascend'); ndx = ndx(ind);
    e = e( ndx , : );

    e = e( [ true ; ~all( ~diff( e , 1 , 1 ) , 2 ) ] , : );
  end

  function c = tri2complex( t )
    c = t(:,1) + t(:,2)/1e7 + 1i* t(:,3);
  end

  function t = complex2tri( c )
    t = [ round(real(c))  ,  ( real(c) - round(real(c)) )*1e7 , imag( c ) ];
  end

  function c = edges2complex( e )
    c = e(:,1) + 1i*e(:,2);
  end

  function e = complex2edges( c )
    e = [ real(c) , imag(c) ];
  end

end
