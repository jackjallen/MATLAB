function T = OneLineCode( T )

  if ischar( T ) || isnumeric( T )
    T = readFile( T );
  end
  T = T(:);

  T = removeComments( T );

  fname = [ tempname '.m'];
  
  fid = fopen( fname , 'w' );
  cellfun( @(s) fprintf( fid , '%s\n' , s ) , T );
  fclose(fid);
  try
    mlintTXTname = tempname;
    fid2 = fopen( mlintTXTname , 'w' );
%     fprintf(fid2,'# Copyright 2009 The MathWorks, Inc.\n\n');
    fprintf(fid2,'1MCC\n0ACCUM\n0AGROW\n0ALIGN\n0AND2\n0ASGSL\n0ATCA\n0ATNCA\n0ATNPP\n0ATPPP\n0ATTF\n0ATUNK\n0BDLGI\n');
    fprintf(fid2,'0BDLOG\n0BDSCA\n0BDSCI\n0CFERR\n0CLTWO\n0COMNL\n0CPROP\n0CTCH\n0CTPCT\n0DEFNU\n0DUALC\n0EFIND\n0EQEFF\n0ERTXT\n0EVLC\n0EXIST\n');
    fprintf(fid2,'0FAFD\n0FALFD\n0FDEPR\n0FLPST\n0FNAN\n0FNDEF\n0FNDSB\n0FOBS\n0FREAD\n0FXSET\n0FXUP\n0GFLD\n0GTARG\n0INUSD\n0INUSL\n0ISCEL\n');
    fprintf(fid2,'0ISCHR\n0ISLOG\n0ISMT\n0ISSTR\n0LERR\n0LNGNM\n0LOGL\n0LSTER\n0LTARG\n0M3COL\n0MCABF\n0MCAP\n0MCCD\n0MCDIR\n0MCDP\n');
    fprintf(fid2,'0MCFIL\n0MCHLP\n0MCKBD\n0MCLL\n0MCMFL\n0MCMLR\n0MCPRD\n0MCPRT\n0MCSVP\n0MCTBX\n0MCUOA\n0MGONE\n0MIPC1\n0MNEFF\n0MXFND\n0N2UNI\n');
    fprintf(fid2,'0NASGU\n0NBRAK\n0NO4LP\n0NOANS\n0NOCOM\n0NODEF\n0NOEFF\n0NOPRT\n0NOPRV\n0NOPTS\n0NOSEM\n0NOV6\n0NUSED\n0OR2\n0PFBFN\n0PFBNS\n');
    fprintf(fid2,'0PFBR\n0PFDF\n0PFEVB\n0PFEVC\n0PFGP\n0PFGV\n0PFIFN\n0PFIIN\n0PFLD\n0PFNF\n0PFNST\n0PFOUS\n0PFPIE\n0PFPRT\n0PFRFH\n0PFRIN\n');
    fprintf(fid2,'0PFRNC\n0PFRNG\n0PFRUS\n0PFSAME\n0PFSV\n0PFTIN\n0PFTUS\n0PFUIX\n0PFUNK\n0PFXST\n0PROP\n0PSIZE\n0PUSE\n0RAND\n0RANDN\n0RGXP1\n');
    fprintf(fid2,'0RHSFN\n0RMFLD\n0SEPEX\n0SETNU\n0SFLD\n0SPERR\n0SPWRN\n0ST2NM\n0STCI\n0STCMP\n0STFLD\n0STLOW\n0STNCI\n0STOUT\n0STRSZ\n0STTOK\n');
    fprintf(fid2,'0TLEV\n0TRIM1\n0TRIM2\n0TRYNC\n0UDIM\n0UNRCH\n0USENS\n0V6ON\n0VCAT\n0VUNUS\n0WLAST\n0WNOFF\n0WNON\n0WNTAG\n0WTXT\n0XFERR\n');
    fclose(fid2);
    
    ml = mlintmex( fname , '-struct', [ '-config=' mlintTXTname ] , '-id' );
    ml = ml{1};
  catch
    try, delete( fname );        end
    try, delete( mlintTXTname ); end
    error('error en mlint!!!');
  end
  try, delete( fname ); end

  
  for n = 1:numel(ml)
    if strcmp( ml(n).id , 'EOFER' )
      fprintf( 'error :  %s\n\n' , ml(n).message );
    else
      fprintf( 'error :  %s\n' , ml(n).message );
      fprintf( '   en :    %s\n' , T{ ml(n).loc(1) } );
      fprintf( '%*s^%*s^\n\n' , ml(n).loc(2) + 9 , '' , ml(n).loc(3) - ml(n).loc(2) , '' );
    end
  end
  if numel(ml)
    try, delete( mlintTXTname ); end
    error('check the code...');
  end
  

  for l = 1:numel(T)
    str = strtrim( T{l} );
    if numel(str)>3 && strcmp( str(end-3:end) , ' ...' )
      T{l} = str(1:end-3);
    elseif str(end) ~= ',' && str(end) ~= ';'
      T{l} = [ str , ' ,' ];
    end
    T{l} = [ T{l} ,  ' ' ];
  end


  T = strtrim( cell2mat( T(:).' ) );
  if strcmp( T(end) , ',' ),  T(end) = [];  end
  

  
  fid = fopen( fname , 'w' );
  fprintf( fid , '%s\n' , T );
  fclose(fid);
  try
    ml = mlintmex( fname , '-struct', [ '-config=' mlintTXTname ] , '-id' );
    ml = ml{1};
  catch
    try, delete( fname );        end
    try, delete( mlintTXTname ); end
    error('error en mlint 2!!!');
  end
  try, delete( fname );        end
  try, delete( mlintTXTname ); end

  
  for n = 1:numel(ml)
    disp( ml(n) );
    fprintf( '\n\n' );
  end
  
end
