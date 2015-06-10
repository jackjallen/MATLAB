function xx = GenRandDistribution(x, y, N)
  
  xx = [];
  while numel(xx) < N
    xxx = rand(N,1)*(x(end)-x(1)) +  x(1);
    xxx = xxx( ( rand(N,1)*max(y) ) < ( Interp1D(y(:), x(:), xxx(:)) ) );
    
    xx = [ xx ; xxx ];
  end
  
  xx = xx(1:N);

end


%{

xy = [ 0  1  2 3 3.1 3.2 3.3 3.5 3.6 3.7 4 5 7; ...
       1 .2 .5 0  0   2   0   0   4   0  1 2 0 ];
xx = GenRandDistribution( xy(1,:) , xy(2,:) , 500000 );
ed = linspace( -1 , 8 , 50001 );
id = getInterval( xx , ed );
c = accumarray( id , 1 )/numel(xx);
c = [ 0 ; c ]; c(numel(ed)+1) = 0;
plot( dualVector(ed) , c , '.-')

%}