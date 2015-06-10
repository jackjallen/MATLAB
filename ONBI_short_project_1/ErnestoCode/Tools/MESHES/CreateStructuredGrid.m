function mesh= CreateStructuredGrid(x,y)
% 
% mesh= create_structured_grid(x,y)
% 
% example:
% 
% x= 0.5:.5:2.4
% y=-1:2
% 
%   .     .     .     .           . --  . --  . --  .
%   (13)  (14)  (15)  (16)        |  /  |  /  |  /  |
%   .     .     .     .           . --  . --  . --  .
%   (9)   (10)  (11)  (12)        |  /  |  /  |  /  |
%   .     .     .     .           . --  . --  . --  .
%   (5)   (6)   (7)   (8)         |  /  |  /  |  /  |
%   .     .     .     .           . --  . --  . --  .
%   (1)   (2)   (3)   (4)
%   
%   (1) = ( 0.5 , -1 )
%   (4) = ( 2.5 , -1 )
%   (13)= ( 0.5 , 2  )
%   (16)= ( 2.5 , 2  )
%   
%   

N= size(x(:),1);
M= size(y(:),1);

[X,Y]= ndgrid( x,y );
xy= [X(:) Y(:)];

tri= [];
for i=1:M-1
  for j=1:N-1
    k= (i-1)*N + j;
    tri= [ tri ; k k+1 k+N+1 ];
    tri= [ tri ; k k+N+1 k+N ];
  end
end

mesh= MakeMesh('p',xy, 't', tri);


