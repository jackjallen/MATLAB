function make( file )

if (nargin == 0)
    file = 'default';
end
    cwd = pwd;  CLEANUP = onCleanup( @() cd(cwd) );
    
    allProjects = {
        'parseargs_st_mex'
        'parseargs_st_batch'
        'libParseargs_st_mex'
        'libParseargs_st_batch'
        'inv4x4'
        };

    %Defines
        
    if(ispc)
        
    end
    
    if(isunix)
        
    end
    
    switch file
        
        case {'parseargs_st_mex'}
            %para Windows
            if(ispc)
               cd Z:\Tools\parseargs\
               make parseargs_st_mex
            end
            
            %Para linux           
            if(isunix)
                cd /extra/disco1/miaTools/Tools/parseargs/parseargs_st
                make parseargs_st_mex
            end
            
        case {'parseargs_st_batch'}
            %para Windows
            if(ispc)
               cd Z:\Tools\parseargs\
               make parseargs_st_batch
            end
            
            %Para linux
            if(isunix)
                cd /extra/disco1/miaTools/Tools/parseargs/parseargs_st
                make parseargs_st_batch
            end
            
        case {'libParseargs_st_mex'}
            %para Windows
            if(ispc)
               cd Z:\Tools\parseargs\
               make libParseargs_st_mex
            end
            
            %Para linux
            if(isunix)
                cd /extra/disco1/miaTools/Tools/parseargs/libParseargs_st
                make libParseargs_st_mex
            end
            
        case {'libParseargs_st_batch'}
            %para Windows
            if(ispc)
               cd Z:\Tools\parseargs\
               make libParseargs_st_batch
            end
            
            %Para linux
            if(isunix)
                cd /extra/disco1/miaTools/Tools/parseargs/libParseargs_st
                make libParseargs_st_batch
            end
            
        case {'inv4x4'}
            %para Windows
            if(ispc)
                cd Z:\Tools\inv4x4
                make inv4x4
            end
            
            %Para linux
            if(isunix)
                cd /extra/disco1/miaTools/Tools/inv4x4
                make inv4x4
            end
            
        case {'all'}
            for p = 1:numel( allProjects)
               make(  allProjects{p} );
            end
                    
        otherwise
            fprintf('Valid projects are:\n\n');
            cellfun( @(p)fprintf('make %s\n',p) , allProjects );
            fprintf('\n\n');
            error('invalid projectName.');
            
            
    end

end