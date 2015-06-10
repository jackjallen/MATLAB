function xyz = Barycentric2XYZ( MESH , tri , bary_coords , newtri )


    
    for f = fieldnames( MESH ).'
      if numel( f{1} ) < 3, continue; end
      if strncmp( f{1} , 'xyz' , 3 )
        
        xyz.(f{1}) = bsxfun( @times , MESH.(f{1})( MESH.tri(tri,1) , : ) , bary_coords(:,1) ) + ...
                     bsxfun( @times , MESH.(f{1})( MESH.tri(tri,2) , : ) , bary_coords(:,2) ) + ...
                     bsxfun( @times , MESH.(f{1})( MESH.tri(tri,3) , : ) , bary_coords(:,3) ) ;
            
      end
    end
    
    if nargin > 3
      
      if isscalar( newtri ) && islogical( newtri ) && isfalse( newtri )
        xyz = xyz.xyz;
      else
        xyz.tri = newtri;
      end
      
    end
    
    
end
