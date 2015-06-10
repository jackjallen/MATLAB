function [ hn_ , hnl ] = getHOSTNAME()

  [state,hn] = system('hostname -s');
  hn = strrep( strrep( hn , char(10) , '' ) , char(13) , '' );
  
  if nargout == 0
    fprintf('%s\n',hn);
    return;
  end
  
  hn_ = hn;
  
  if nargout > 1
    [state,hnl] = system('hostname');
    hnl = strrep( strrep( hnl , char(10) , '' ) , char(13) , '' );
  end


end