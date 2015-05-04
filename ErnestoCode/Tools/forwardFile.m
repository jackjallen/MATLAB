function forwardFile( fid , N , varargin )
%{

fid=fopen('fichero.txt','w');
fprintf(fid,'%f\n',1:1e7);
fclose(fid);


fid=fopen('fichero.txt','r');

data = [];
while ~feof(fid)
  d = textscan(fid,'%f',1);
  data(end+1,1) = d{1};
  forwardFile( fid , 1050000 );
end
data

fclose(fid);

%}


  try
    positions = ftell( fid );
  catch
    error( 'no parece ser un fid' );
  end

  [varargin,i,buffersize] = parseargs( varargin ,'BufferSize','buffer' ,'$DEFS$', 1e5 );
  [varargin,i,delimiters] = parseargs( varargin ,'DELimiters'          ,'$DEFS$', uint8( sprintf('\n') ) );
  delimiters = uint8( delimiters );

  while ~feof( fid )
    data = fread( fid , buffersize , '*uint8' );
    data = ismembc( data , delimiters );

    S = sum( data );
    if S >= N
      fseek( fid , -buffersize  , 'cof' );
      
      positions = find(data,N);
      fseek( fid , positions(end) , 'cof' );
      return;
    end
    
    N = N - S;
  end

end
