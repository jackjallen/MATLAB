function M = closeMyocardium( M )

  if numel( unique( M.triRegion ) ) ~= 2
    error('Tow pieces must be provided');
  end
  
  Bids = find( M.xyzIsBoundary > 0 );
%   Bids = Bids([1:9,26:32,18:25,10:17]);

  if 0
  newT = bsxfun(@plus, [1,2,33,33,2,34] , [1:8,26:31,18:24,10:16].' ) - 1;
  newT = reshape( newT , [] , 3 );
  newT = [ newT ; 9 41 26 ; 26 58 41 ; 32 18 64 ; 64 18 50 ; 42 57 25 ; 25 10 42 ; 1 33 49; 49 17 1 ];

  plot3d( M.xyz(Bids,:) , '.' );
  text( M.xyz(Bids,1),M.xyz(Bids,2),M.xyz(Bids,3),arrayfun(@(n)sprintf('%d',n),1:numel(Bids),'un',0),'verticalalignment','bottom','horizontalalignment','left');
  patch( 'vertices', M.xyz( Bids , : ) , 'faces' , newT ,'facecolor','c');
  view([-1 4])
  end

  newT = [1,2,33;1,17,49;1,33,49;2,3,34;2,33,34;3,4,35;3,34,35;4,5,36;4,35,36;5,6,37;5,36,37;6,7,38;6,37,38;7,8,39;7,38,39;8,9,40;8,39,40;9,26,41;9,40,41;10,11,42;10,25,42;11,12,43;11,42,43;12,13,44;12,43,44;13,14,45;13,44,45;14,15,46;14,45,46;15,16,47;15,46,47;16,17,48;16,47,48;17,48,49;18,19,50;18,32,64;18,50,64;19,20,51;19,50,51;20,21,52;20,51,52;21,22,53;21,52,53;22,23,54;22,53,54;23,24,55;23,54,55;24,25,56;24,55,56;25,42,57;25,56,57;26,27,58;26,41,58;27,28,59;27,58,59;28,29,60;28,59,60;29,30,61;29,60,61;30,31,62;30,61,62;31,32,63;31,62,63;32,63,64];
  
  M.tri = [ M.tri ; Bids( newT ) ];
  M.triRegion(end+1:size(M.tri,1)) = -M.triRegion(1);

  M = FixNormals( M );
%   plotMESH( M ,'td',ComputeNormals(M)); camlight

end
