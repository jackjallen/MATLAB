function [L,S] = breakPolyline( L1 , L2 )

if 0

  t = linspace(0,2*pi-1,500).';
  r1 = sin(4*t)+2;  L1 = [ r1.*cos(t) , r1.*sin(t) ];
  t = linspace(0.2,2*pi-0.5,10000).';
  r2 = sin(17*t)+2.5;  L2 = [ r2.*cos(t) , r2.*sin(t) ];
  
  plot( L1(:,1),L1(:,2) , 'r' , L2(:,1),L2(:,2) , 'k' );
 
  L = breakPolyline( L1 , L2 );
  for s = 1:numel(L)
    hplot( L{s}(:,1) , L{s}(:,2) , 'o-','color' , colorith( rem(s,5) + 1 ),'linewidth',2 );
  end
  
end

  Xs = (1:size(L1,1)).';

  AB = diff( L1 , 1 , 1 );
  CD = diff( L2 , 1 , 1 ).';
  DEN = 1./( AB(:,1)*CD(2,:) - AB(:,2)*CD(1,:) );

  t1 = bsxfun(@plus , L1(1:end-1,:)*[0 -1;1 0]*CD   , sum( L2(1:end-1,:) .* ( L2(2:end,:)*[0 -1;1 0] ) , 2 ).' ).* DEN;
  t2 = bsxfun(@plus , AB*[0 -1;1 0]*L2(1:end-1,:).' , sum( L1(1:end-1,:) .* ( L1(2:end,:)*[0 1;-1 0] ) , 2 )   ).* DEN;

  w = ( 0 <= t1 & t1 <= 1 & 0 <= t2 & t2 <= 1 );
  
  if ~any(w(:))
    L = {L1};
    S = { [ (1:size(L1,1)).' , zeros( size(L1,1) , 1 ) ] };
    return;
    
  end
  
  [S1,S2] = find( w );
  
  t1 = bsxfun( @plus , Xs(1:end-1) , t1 );
  t1 = t1(w);
  [t1,idx] = unique( t1 );
  S2 = S2(idx);
  
  if t1(1) ~= 1
    t1 = [ 1 ; t1(:) ];
    S2 = [ 0 ; S2(:) ];
  end
  if t1(end) ~= size(L1,1)
    t1 = [ t1(:) ; size(L1,1) ];
  end
  
  
%   t1 = unique( [1 ; t1(w) ; size(L1,1) ] );
  
  L = cell(numel(t1)-1,1); S = cell(numel(t1)-1,1);
  
  for s = 1:(numel(t1)-1)
    i0 = t1(s  );
    i1 = t1(s+1);
    switch ( ~rem(i0,1)*2 + ~rem(i1,1) )
      case 3        %[ free  frees  free ]
        ps = i0:i1;
        L{s} = L1(ps,:);
        S{s} = [ ps(:) , zeros(numel(ps),1) ];
        
      case 1        %[ constrained  frees  free ]
        ps = ceil(i0):i1;
        ii0 = floor(i0);
        L{s} = [ L1(ii0,:) + (t1(s)-ii0)*( L1(ii0+1,:) - L1(ii0,:) )  ;...
                 L1(ps,:)                                         ];
        S{s} = [ ii0   , S2(s)                                    ;...
                 ps(:) , zeros(numel(ps),1)                       ];
               
               
      case 2        %[ free  frees  constrained ]
        ps = i0:floor(i1);
        ii1 = ceil(i1);
        L{s} = [ L1(ps,:)                                          ;...
                 L1(ii1-1,:) + (t1(s+1)-ii1+1)*( L1(ii1,:) - L1(ii1-1,:) ) ];
        S{s} = [ ps(:) , zeros(numel(ps),1)                       ;...
                 ii1-1 , S2(s+1)                                    ];
               
      case 0        %[ constrained  frees  constrained ]
        ps = ceil(i0):floor(i1);
        ii0 = floor(i0);
        ii1 = ceil(i1);
        L{s} = [ L1(ii0,:) + (t1(s)-ii0)*( L1(ii0+1,:) - L1(ii0,:) )  ;...
                 L1(ps,:)                                             ;...
                 L1(ii1-1,:) + (t1(s+1)-ii1+1)*( L1(ii1,:) - L1(ii1-1,:) ) ];
        S{s} = [ ii0   , S2(s)                                    ;...
                 ps(:) , zeros(numel(ps),1)                       ;...
                 ii1-1 , S2(s+1)                                  ];
        
    end
    
  end
  
  
  
%   p1 = ceil(  S1(1) + t1( S1(1) , S2(1) ) );
%   ps = 1:p1-1;
%   
%   
%   L{1} = [ L1(ps,:)                                                         ;...
%              L1(p1,:) + t1( S1(1) , S2(1) )* ( L1(p1  ,:) - L1(p1-1,:) )  ];
%   S{1} = [ ps(:) , zeros(numel(ps),1)   ;...
%              S1(1) S2(1)              ];
%   for s = 1:numel(S1)-1
%     p0 = floor( S1(s  ) + t1( S1(s  ) , S2(s  ) ) );
%     p1 = ceil(  S1(s+1) + t1( S1(s+1) , S2(s+1) ) );
%     
%     ps = p0+1:p1-1;
%     L{s+1} = [ L1(p0,:) + t1( S1(s  ) , S2(s  ) )* ( L1(p0+1,:) - L1(p0  ,:) )  ;...
%                L1(ps,:)                                                         ;...
%                L1(p1,:) + t1( S1(s+1) , S2(s+1) )* ( L1(p1  ,:) - L1(p1-1,:) )  ];
%     S{s+1} = [ S1(s) S2(s)                  ;...
%                ps(:) , zeros(numel(ps),1)   ;...
%                S1(s+1) S2(s+1)              ];
%   end
%   p0 = floor( S1(end) + t1( S1(end) , S2(end) ) );
% 
%   ps = p0+1:size(L1,1);
%   L{end+1} = [ L1(p0,:) + t1( S1(end) , S2(end) )* ( L1(p0+1,:) - L1(p0  ,:) )  ;...
%                L1(ps,:)                                                         ];
%   S{end+1} = [ S1(end) S2(end)              ;...
%                ps(:) , zeros(numel(ps),1)   ];
%   

             
%   t1 = bsxfun(@plus, Xs(1:end-1) , t1 );
%   newPoints = unique( [ 1 ; vec( t1(w(:)) ) ; size(L1,1) ] );
%   
%   L = {};
%   for c = 1:numel(newPoints)-1
%     ids = unique( [ newPoints(c) , ceil(newPoints(c)):newPoints(c+1) , newPoints(c+1) ] );
%     if numel(ids) > 1
%       L{end+1} = Interp1D( L1 , Xs , ids );
%     end
%   end

end
