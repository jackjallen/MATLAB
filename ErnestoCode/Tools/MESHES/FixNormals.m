function M = FixNormals( M )
%
%  M = FixNormals( M )
%

  N    = ComputeNormals( M );
  Nvtk = vtkComputeNormals( M , 'SetComputeCellNormals',true,'SetComputePointNormals',false,'SetAutoOrientNormals',true);

  change = find( sum( N .* Nvtk , 2 ) < 0 );
  M.tri(change,:) = M.tri(change,[2 1 3]);

end
