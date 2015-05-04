function x= delimit( x , minn, maxx )
% 
% delimit ( x )  = delimint( x, 0 ,1 );
% delimit ( x , range )= delimit( x, range(:,1) , range(:,2) );
% delimit ( x , min , max )= delimit( x, max, min);
% 

  switch nargin
    case 1
      range= [0 1];
    case 2
      range= minn;
      if size(range,2) ~= 2
        error('range has to be a 2-column argument.');
      end
    case 3
      range(:,1)= minn(:);
      range(:,2)= maxx(:);
  end

  minimo= min(range,[],2);
  maximo= max(range,[],2);

  if numel( minimo )==1
    minimo= repmat( minimo, numel(x),1);
  end
  if numel( maximo )==1
    maximo= repmat( maximo, numel(x),1);
  end
  
  i= x(:)<minimo;
  x( i )= minimo(i);

  j= x(:)>maximo;
  x( j )= maximo(j);

end
