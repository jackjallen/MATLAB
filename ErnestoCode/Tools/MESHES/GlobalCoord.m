function xyz= GlobalCoord( rs , p1,p2,p3)
% 
% xyz= GlobalCoord( rst , triangle )
% 
% convert from parametrics coordinates (rst) to
% global coordinates (xyz) 
% 
% triangle= [ p1x p1y p1z
%             p2x p2y p2z
%             p3x p3y p3z ]
%

xyz= p1 + rs(1)*(p2-p1) + rs(2)*(p3-p1);

