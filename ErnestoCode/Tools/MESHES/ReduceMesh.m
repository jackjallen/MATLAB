function m= ReduceMesh( m , r )

  if nargin < 2
    r= 0.5;
  end

  [m.tri,m.xyz]= reducepatch( m.tri, m.xyz, r );

end
