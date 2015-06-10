function rs= LocalCoord( xyz , p1,p2,p3)
% 
% rst= LocalCoord( xyz ,triangle )
% 
% convert from parametrics global coordinates (xyz) 
%  to coordinates (rst)
% 
% triangle= [ p1x p1y p1z
%             p2x p2y p2z
%             p3x p3y p3z ]
%

rs= (xyz - p1)*inv( [ p2-p1; p3-p1]' )';

