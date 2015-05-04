function N = nonii( I , classtype )
%
% Number Of Numbers In Interval
%
% N = nonii( [x1 x2] );
%

  if numel(I) ~= 2 || I(1) > I(2), error('an increasing 2 vector is expected.'); end
  x1 = I(1); x2 = I(2);
  
  if nargin < 2
     classtpye = class( I );
  end
  switch lower( classtpye )
    case 'double'
      x1 = double(x1);
      x2 = double(x2);
      if      x1 >= 0 && x2 > 0

        b1 = sprintf( '%08x.%08x' , typecast( x1 , 'uint32' ) );
        l1 = cell2mat( arrayfun( @(h) dec2bin(hex2dec(h),4) , b1( 1:8  )  , 'un',0) );
        u1 = cell2mat( arrayfun( @(h) dec2bin(hex2dec(h),4) , b1(10:end)  , 'un',0) );
        n1 = bin2dec( [ u1(13:end)  l1 ] );
        e1 = bin2dec( u1(2:12) );
        % s1 = bin2dec( u1(1) );
        % realpow(-1,s1)* ( 1+realpow(2,-52)*n1 ) *  realpow(2,e1-1023)

        b2 = sprintf( '%08x.%08x' , typecast( x2 , 'uint32' ) );
        l2 = cell2mat( arrayfun( @(h) dec2bin(hex2dec(h),4) , b2( 1:8  )  , 'un',0) );
        u2 = cell2mat( arrayfun( @(h) dec2bin(hex2dec(h),4) , b2(10:end)  , 'un',0) );
        n2 = bin2dec( [ u2(13:end)  l2 ] );
        e2 = bin2dec( u2(2:12) );

        N = realpow(2,52)*( e2 - e1 - 1 ) + n2 - n1 + hex2dec('fffffffffffff') + 1 + 1;

      elseif  x1 < 0 && x2 < 0
        
        N = nonii( -[ x2 x1 ] );
        
      elseif  x1 < 0 && x2 >= 0
        
        N = nonii( [eps(0) -x1] ) + nonii( [0 x2] );
        
      end
      
      
      if 0
        next = @(x) x+eps(x);
        x=x1; i=1;
        while ~isequal(x,x2)
          x=next(x);
          i=i+1; 
        end
        i
      end
        
    case {'single','float'}
      x1 = single(x1);
      x2 = single(x2);
      if      x1 >= 0 && x2 > 0
        b1 = sprintf( '%04x.%04x' , typecast( x1 , 'uint16' ) );
        l1 = cell2mat( arrayfun( @(h) dec2bin(hex2dec(h),4) , b1( 1:4   )  , 'un',0) );
        u1 = cell2mat( arrayfun( @(h) dec2bin(hex2dec(h),4) , b1( 6:end )  , 'un',0) );
        n1 = bin2dec( [ u1(10:end)  l1 ] );
        e1 = bin2dec( u1(2:9) );
        % s1 = bin2dec( u1(1) );
        % realpow(-1,s1)* ( 1+realpow(2,-23)*n1 ) *  realpow(2,e1-127)

        b2 = sprintf( '%04x.%04x' , typecast( x2 , 'uint16' ) );
        l2 = cell2mat( arrayfun( @(h) dec2bin(hex2dec(h),4) , b2( 1:4   )  , 'un',0) );
        u2 = cell2mat( arrayfun( @(h) dec2bin(hex2dec(h),4) , b2( 6:end )  , 'un',0) );
        n2 = bin2dec( [ u2(10:end)  l2 ] );
        e2 = bin2dec( u2(2:9) );

        N = realpow(2,23)*( e2 - e1 - 1 ) + n2 - n1 + hex2dec('7fffff') + 1 + 1;

      elseif  x1 < 0 && x2 < 0
        
        N = nonii( -[ x2 x1 ] );
        
      elseif  x1 < 0 && x2 >= 0
        
        N = nonii( [eps(single(0)) -x1] ) + nonii( [single(0) x2] );
        
      end
      
      
      if 0
        next = @(x) x+eps(x);
        x=x1; i=1;
        while ~isequal(x,x2)
          x=next(x);
          i=i+1; 
        end
        i
      end      
      
      
  end
    



end