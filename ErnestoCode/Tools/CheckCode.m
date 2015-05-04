function codeHasErrors = CheckCode( T , messages )

  if nargin < 2
    messages = {'1MCC';'0ACCUM';'0AGROW';'0ALIGN';'0AND2';'0ASGSL';'0ATCA';'0ATPPP';'0ATTF';'0BDLGI';'0BDLOG';'0BDSCA';'0BDSCI';'0CFERR';'0COMNL';'0CPROP';'0CTCH';'0DEFNU';'0DUALC';'0EFIND';'0EQEFF';'0EVLC';'0EXIST';'0FDEPR';'0FLPST';'0FNAN';'0FNDEF';'0FNDSB';'0FOBS';'0FREAD';'0FXSET';'0FXUP';'0GFLD';'0INUSD';'0INUSL';'0ISCEL';'0ISCHR';'0ISLOG';'0ISMT';'0ISSTR';'0LERR';'0LNGNM';'0LOGL';'0LSTER';'0M3COL';'0MCABF';'0MCAP';'0MCCD';'0MCDP';'0MCHLP';'0MCKBD';'0MCLL';'0MCMFL';'0MCMLR';'0MCPRD';'0MCPRT';'0MCSVP';'0MCTBX';'0MCUOA';'0MIPC1';'0MNEFF';'0MXFND';'0N2UNI';'0NASGU';'0NBRAK';'0NO4LP';'0NOANS';'0NOCOM';'0NODEF';'0NOEFF';'0NOPRT';'0NOPTS';'0NOSEM';'0NUSED';'0OR2';'0PFBNS';'0PFGV';'0PFIIN';'0PFOUS';'0PFPRT';'0PFRIN';'0PFRNG';'0PFRUS';'0PROP';'0PSIZE';'0RAND';'0RANDN';'0RGXP1';'0RMFLD';'0SEPEX';'0SETNU';'0SFLD';'0SPERR';'0SPWRN';'0ST2NM';'0STCI';'0STCMP';'0STFLD';'0STLOW';'0STNCI';'0STOUT';'0STRSZ';'0STTOK';'0TLEV';'0TRIM1';'0TRIM2';'0TRYNC';'0UDIM';'0UNRCH';'0USENS';'0VCAT';'0VUNUS';'0WLAST';'0WNOFF';'0WNON';'0WNTAG';'0XFERR'};
%     messages = {'0ACCUM';'0AGROW';'0ALIGN';'0AND2';'11ASGLU';'0ATCA';'0ATNCA';'0ATNPP';'0ATPPP';'0ATTF';'0ATUNK';'0BDLGI';'0BDLOG';'0BDOPT';'0BDSCA';'0BDSCI';'0CFERR';'0CFIG';'0CLTWO';'0COLND';'0COMNL';'0CONT';'0CPROP';'0CTCH';'0CTPCT';'0DEFNU';'0DSARG';'0DUALC';'0EFIND';'0EQEFF';'0ERTXT';'0EVLC';'0EXIST';'0FAFD';'0FALFD';'0FDEPR';'0FLPST';'0FNAN';'0FNDEF';'0FNDSB';'0FOBS';'0FREAD';'0FXSET';'0FXUP';'0GFLD';'0GPFST';'0GTARG';'0INUSD';'0INUSL';'0ISCEL';'0ISCHR';'0ISLOG';'0ISMT';'0ISSTR';'0LERR';'0LNGNM';'0LOGL';'0LSTER';'0LTARG';'0M3COL';'0MCABF';'0MCAP';'0MCCD';'0MCDIR';'0MCDP';'0MCFIL';'0MCHLP';'0MCKBD';'0MCLL';'0MCMFL';'0MCMLR';'0MCPRD';'0MCPRT';'0MCSVP';'0MCTBX';'0MCUOA';'0MGONE';'0MIPC1';'0MNEFF';'0MXFND';'0N2UNI';'0NBRAK';'0NO4LP';'0NOANS';'0NOBRK';'0NOCLS';'0NOCOM';'0NOEFF';'0NOPRT';'0NOPRV';'0NOPTS';'0NOSEM';'0NOV6';'0NPERS';'0NUSED';'0OR2';'0PFBFN';'0PFBNS';'0PFBR';'0PFDF';'0PFEVB';'0PFEVC';'0PFGP';'0PFGV';'0PFIFN';'0PFIIN';'0PFLD';'0PFNF';'0PFNST';'0PFOUS';'0PFPIE';'0PFPRT';'0PFRFH';'0PFRIN';'0PFRNC';'0PFRNG';'0PFRUS';'0PFSAME';'0PFSV';'0PFTIN';'0PFTUS';'0PFUIX';'0PFUNK';'0PFXST';'0PROP';'0PSIZE';'0PUSE';'0RAND';'0RANDN';'0REDEF';'0RGXP1';'0RHSFN';'0RMFLD';'0ROWLN';'0SBTMP';'0SEPEX';'0SFLD';'0SPERR';'0SPWRN';'0ST2NM';'0STCI';'0STCMP';'0STFLD';'0STLOW';'0STNCI';'0STOUT';'0STRSZ';'0STTOK';'0TLEV';'0TRIM1';'0TRIM2';'0TRYNC';'0TWOCM';'0UDIM';'0UNRCH';'0V6ON';'0VCAT';'0VUNUS';'0WLAST';'0WNOFF';'0WNON';'0WNTAG';'0WTXT';'0XFERR';'0XFIG'}
  end

  fname = [ tempname '.m'];
  
  fid = fopen( fname , 'w' );
  cellfun( @(s) fprintf( fid , '%s\n' , s ) , T );
  fclose(fid);

  try
    mlintTXTname = [ tempname '.txt' ];
    fid = fopen( mlintTXTname , 'w' );
    cellfun( @(s) fprintf( fid , '%s\n' , s ) , messages );
    fclose(fid);
    
    ml = mlintmex( fname , '-struct', [ '-config=' mlintTXTname ] , '-id' );
    ml = ml{1};
  catch
    try, delete( fname );         end
    try, delete( mlintTXTname );  end
    error('error en mlint!!!');
  end
  try, delete( fname );         end
  try, delete( mlintTXTname );  end


  for ee = 1:numel(ml)
    if strcmp( ml(ee).id , 'ASGLU' )
      ml(ee).id = -1;
    end
  end
  ml = ml( arrayfun( @(ee) ischar( ml(ee).id ) , 1:numel(ml) ) );
  
  
  if numel( ml ) > 0
    
    intEditor = handle(awtcreate('com.mathworks.mlwidgets.interactivecallbacks.InteractiveCallbackEditor', ...
      'Ljava.awt.Rectangle;Ljava.lang.String;Ljava.lang.String;', ...
      java.awt.Rectangle(200,400,600,400),...
      [] , sprintf('%s\n',T{:}) ) );
    intEditor.setVisible(true);
    
    codeHasErrors = true;

  else

    codeHasErrors = false;
    
  end
  
end
