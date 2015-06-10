function u_ = getUSER
  
  u = getenv('USER');

%   [state, u ] = system( 'who am i' );
%   u = u( 1 : ( find( u==' ' , 1)-1 ) );
%   
%   u = strrep( strrep( u , char(10) , '' ) , char(13) , '' );

  if nargout == 0
    fprintf('%s\n',u);
    return;
  end
  
  u_ = u;

end