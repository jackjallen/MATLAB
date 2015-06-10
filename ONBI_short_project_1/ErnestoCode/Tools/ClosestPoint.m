function [id,dist] = ClosestPoint( P , Q )

%%

  %NPs = [ 10 20 50 100 200 500 1000 2000 5000 ];
  NPs = unique( round( geospace( 60 , 2000 , 35 ) ));

  %NQs = [ 5  10 20 50 100 200 ];
  NQs = unique( round( geospace( 10 , 150 , 15 ) ));
  
  ppqq = ndmat( 1:numel(NPs) , 1:numel(NQs) );
  T1 = zeros(numel(NPs),numel(NQs));  T2 = T1;  T3 = T2;
  
  for it = 1:1
    for r = shuffle( 1:size(ppqq,1) )

        pp = ppqq(r,1);
        qq = ppqq(r,2);

        disp( [ NPs( pp ) , NQs( qq ), sum( T1(:) == 0 )] );

        P = ndmat( linspace(-10,10,ceil(cbrt(NPs(pp)))) , linspace(-10,10,ceil(cbrt(NPs(pp)))) , linspace(-10,10,ceil(cbrt(NPs(pp)))) );
%         P = randn( NPs(pp) ,3);

        Q = ndmat( linspace(-10,10,NQs(qq)) , linspace(-10,10,NQs(qq)) , linspace(-10,10,NQs(qq)) );
%         Q = randn( NQs(qq)^3 , 3 );

        tic;
        [id,dist1] = ClosestPoint( P , Q ); 
        t = toc;
        T1(pp,qq) = T1(pp,qq) + t;

        tic;
        delau = delaunayn( P );
        [id,dist2] = dsearchn( P , delau , Q );
        t = toc;
        T2(pp,qq) = T2(pp,qq) + t;

        tic;
        GLtree = BuildGLTree3D( P.' );
        [id,dist3] = KNNSearch3D( P.' , Q.' , GLtree , 1 );
        DeleteGLTree3D( GLtree );
        t = toc;
        T3(pp,qq) = T3(pp,qq) + t;

    end

%   range( dist1 - dist2 )
%   range( dist1 - dist3 )
  
%     cla
%     surf( NQs , NPs , T1 , 'facecolor','r' )
%     hold on;
%     surf( NQs , NPs , T2 , 'facecolor','g' )
%     surf( NQs , NPs , T3 , 'facecolor','b' )
%     hold off
% 
%     set(gca,'XScale','log','yscale','log','ZScale','log');
%     drawnow
    save g:\closestPoint_Ts_randn


  end

%%


% end



%%

    cla
    surf( log10(NQs) , log10(NPs) , log10(T1) , 'facecolor','m' )
    hold on; 
    surf( log10(NQs) , log10(NPs) , log10(T2) , 'facecolor','g' )
    surf( log10(NQs) , log10(NPs) , log10(T3) , 'facecolor','b' )
    hold off

    set(gca,'XTick',log10(NQs),'xtickLabel',arrayfun(@(x)num2str(x),NQs,'un',0) );
    set(gca,'YTick',log10(NPs),'ytickLabel',arrayfun(@(x)num2str(x),NPs,'un',0) );
    
    xlabel( 'Q' )
    ylabel( 'P' )    
    
    view(0,-90)
    
 