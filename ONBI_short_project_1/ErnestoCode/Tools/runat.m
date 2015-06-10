function runat( secs , F )

  persistent runat______

  if ~isempty( runat______ )
    try, stop( runat______ );                                 end
    if isvalid( runat______ )  &&  get( runat______ , 'TasksExecuted' ) == 0
      try, feval( get( runat______ , 'TimerFcn' ) , 0 , 0 );    end
    end
  end
  try
    delete( runat______ );
  end
      
    
  
  runat______ = timer( 'TimerFcn', @(varargin) exec_and_delete( F ) , 'StartDelay', secs ,'ExecutionMode','singleShot');
  start( runat______ );


  function exec_and_delete( F )
    feval( F );
    stop( runat______ );
    delete( runat______ );
  end
  
end
