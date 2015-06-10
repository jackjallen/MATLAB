function M = MercatorMesh( M )


  for f=fieldnames(M).'
    if strcmp( f{1} , 'xyz' ) || strncmp( f{1} , 'lmk' , 3 )
      M.(f{1}) =  [ atan2( M.(f{1})(:,1) , M.(f{1})(:,2) ) , atan( M.(f{1})(:,3)./sqrt( M.(f{1})(:,1).^2 + M.(f{1})(:,2).^2 ) ) ];
      M.(f{1})(end,3) = 0;
    end
  end
    
  M.tri( any( reshape( M.xyz( M.tri , 1 ) , [] , 3 ) >  pi/2 , 2 ) & ...
         any( reshape( M.xyz( M.tri , 1 ) , [] , 3 ) < -pi/2 , 2 ) ,:) = [];
  M.tri( any( reshape( M.xyz( M.tri , 2 ) , [] , 3 ) >  pi/4 , 2 ) & ...
         any( reshape( M.xyz( M.tri , 2 ) , [] , 3 ) < -pi/4 , 2 ) ,:) = [];
  

end
