function [ isC , clusterID ] = isCONDOR()

  isC       = false;
  clusterID = [];
  if ~isunix, return; end
  
  try
    CONDOR_prefix = condor('bin');
  catch
    CONDOR_prefix = 'ssh -o ConnectTimeout=1 -o NumberOfPasswordPrompts=0 -o BatchMode=yes hermes.cps.unizar.es  /opt/condor-7.4.1/bin/';
  end
  
    
  mpid = feature('getpid');
  
  pid = mpid;    allPIDS = pid;
  cmd = '';
  while pid
    prev_cmd = cmd;
    [status,cmd] = system( sprintf('export COLUMNS=10000; ps -o cmd=  %d' , pid ) );
    
    %fprintf('%g -->  %s\n',pid,cmd);
    
    if ~isempty( regexp( cmd ,'condor_starter','ONCE') )
      isC = true;
      break;
    end
    
    [status,pid] = system( sprintf('ps -o ppid=  %d' , pid ) );
    pid = str2double(pid);   allPIDS(end+1)= pid;
  end
  
  allPIDS = sort( allPIDS );
  
  if isC  && nargout > 1

    prev_cmd = regexp( prev_cmd , '(\S)+\s+(\S)+\s+([^\n]*)' ,'tokens' );
    
    
    condorQ = [ ...
      CONDOR_prefix ...
      'condor_q ' ...
      sprintf('-submitter %s ' , getUSER ) ...
      '-format \"%d.\" ClusterId -format \"%d\\n\" ProcId ' ...
      '-constraint ' ...
      '\" '...
      sprintf('Cmd  == \\\\\\\"%s\\\\\\\" ' , prev_cmd{1}{2} )...
      '\&\& '...
      sprintf('Args == \\\\\\\"%s\\\\\\\" ' , strrep(strrep(strrep(strrep(prev_cmd{1}{3},'''','\'''),'(','\('),')','\)'),';','\;') )  ...
      '\&\& '...
      'JobStatus == 2 ' ...
      '\&\& '...
      '\( ' ...
      sprintf('substr\\(RemoteHost,6,8\\) == \\\\\\\"%s\\\\\\\" ', getHOSTNAME ) ...
      '\|\| ' ...
      sprintf('substr\\(RemoteHost,7,8\\) == \\\\\\\"%s\\\\\\\" ', getHOSTNAME ) ...
      '\) ' ...
      '\" '];

    for loop_id = 1:15
      [status,clusters] = system( condorQ );

      %disp(['status   = ',uneval(status)  ]);
      %disp(['clusters = ',uneval(clusters)]);
      
      
      if ~status && ~isempty(regexp(clusters,'ssh', 'once'))
        status = 1024;
      end
      
      if ~status, break; end
      pause(1);
    end
    
    clusters = textscan(clusters,'%s','delimiter','\n');
    clusters = clusters{1};
    
    if numel( clusters ) == 1

      clusterID = clusters{1};
      
    else
      
      %fprintf('-----------\n\n'); disp(['clusters   = ' , uneval(clusters)]); fprintf('-----------\n\n');

      for loop_id = 1:15
        for i = 1:numel( clusters )
          if isempty( clusters{i} ), continue; end
          
          %fprintf('\n\n%s\n',clusters{i});

          condor_SSH = [...
            CONDOR_prefix ...
            'condor_ssh_to_job ' ...
            clusters{i} ' ' ...
            '''echo \$_CONDOR_JOB_PIDS''' ];

          [status,job_pid] = system( condor_SSH );

          %disp(['job_pid = ' uneval(job_pid)]);

          job_pid = sscanf( job_pid , '%g' );
          if isempty(job_pid), continue; end

          %disp(['job_pid = ' uneval(job_pid)]);
          %disp(['allPIDS = ' uneval(allPIDS)]);

          if any( ismembc( job_pid , allPIDS ) )
            clusterID = clusters{i};
            break;
          end
          clusters{i} = [];
        end
        
        if ~isempty(clusterID), break; end
        pause(1);
      end
      
    end

  end
  
  
end

