function isarch = isarch64()

    [C] = computer;
    if(ispc)
        if( strcmp(C ,'PCWIN') == 1 )
            isarch = false;
        else if( strcmp(C ,'PCWIN64') == 1 )
                isarch = true;
            else
                error('Arch unknown.');
            end
        end
    end

    if(isunix)
        if( strcmp(C ,'GLNX86') == 1 )
            isarch = false;
        else if( strcmp(C ,'GLNXA64') == 1 )
                isarch = true;
            else
                error('Arch unknown.');
            end
        end
    end

end