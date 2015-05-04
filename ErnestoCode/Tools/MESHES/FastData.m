function [V,A,N]=FastData(M);

if isempty(M.tri)
    V=0;
    A=0;
    N=[];
else
N = cross( ( M.xyz(M.tri(:,2),:) - M.xyz(M.tri(:,1),:) ) , ( M.xyz(M.tri(:,3),:) - M.xyz(M.tri(:,1),:) ) , 2 );
A = sqrt( sum( N.^2 , 2 ) );
N(A==0,:)=0;
N(A~=0,:) = bsxfun( @rdivide , N(A~=0,:) , A(A~=0) );
A = A/2;
V=  sum(  ( M.xyz( M.tri(:,1),: ) + M.xyz( M.tri(:,2),: ) + M.xyz( M.tri(:,3),: ) ) .* N , 2 );
V = 2*sum( V .* A )/( 3 * factorial(3) );
end;