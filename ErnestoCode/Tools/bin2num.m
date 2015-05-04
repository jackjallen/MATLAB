function B = bin2num( B , type )

  if ischar( B ) , B = B == '1'; end
  B = reshape( B , 8 , numel(B)/8 );

  B = uint8(  single([128,64,32,16,8,4,2,1])  * B );
%   B = uint8(  fliplr(single([128,64,32,16,8,4,2,1]))  * B );

  if nargin > 1
    B = typecast( B(:) , type );
  end
  
end
