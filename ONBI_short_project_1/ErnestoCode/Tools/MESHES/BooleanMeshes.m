function [varargout]=BooleanMeshes(M,op,N,varargin)

% function [varargout]=IntersectBooleanMeshes(M,op,N,varargin)
%           op:
%           'carve':
%           'bool':
%           'NoFix':
%           'RemoveDup':
%           'Clean':
%           'PreReord':
%           'Tol':

% save('/extra/disco1/users/david/dataerror','M','N');
k=0.001;

method='normal';
[varargin,method]       = parseargs(varargin,'Carve','$FORCE$',{'carve',method});
[varargin,method]       = parseargs(varargin,'Triangle','$FORCE$',{'triang',method});
saved_varargin=varargin;
[varargin,fix]          = parseargs(varargin,'NoFix','$FORCE$',{false,true});
[varargin,removedup]    = parseargs(varargin,'RemoveDup','$FORCE$');
[varargin,clean]        = parseargs(varargin,'Clean','$FORCE$');
[varargin,istol,tol]    = parseargs(varargin,'tol','$FORCE$');
[varargin,prereord]     = parseargs(varargin,'PreReord','$FORCE$');

if ~istol, tol=0.0; end;

switch op
    case {'intersect','intersection','i','and','&'}
        op_carve = 'i';
        op_bool = 1;
    case {'union','join','u','or','|','+'}
        op_carve = 'u';
        op_bool = 0;
    case {'difference','-','dif'}
        op_carve = '-';
        op_bool = 2;
    otherwise
        error('Unknown boolean operation');
end;


if prereord
    
    M=vtkPolyDataNormals(M,'ComputePointNormalsOff', [], 'ComputeCellNormalsOn'     , [], 'ConsistencyOn'            , [], 'AutoOrientNormalsOn',[]);
    M=rmfield(M,'triNormals');
    
end;

% save dataerror M N;

if nargout<=1
    switch method
        case 'normal'
            try
                [varargout{1}]=carveIntersect(M,N,op_carve);
            catch
                disp(['error in carveIntersect ( op = ' op ' ). Trying with triangle']);
                try
                   [M_out,M_in,M_on,M_no,N_out,N_in]=IntersectMeshes_OK(Mnew,N);
                    switch op_bool
                        case 0, [varargout{1}]=AppendMeshes(M_out,N_out); 
                                [varargout{1}]=AppendMeshes(varargout{1},M_on); % union
                        case 1, [varargout{1}]=AppendMeshes(M_in,N_in); 
                                [varargout{1}]=AppendMeshes(varargout{1},M_on); % interseccion
                        case 2, N_in.tri=N_in.tri(:,[1 3 2]);
                                [varargout{1}]=AppendMeshes(M_out,N_in); 
                                [varargout{1}]=AppendMeshes(varargout{1},M_no); % diferencia
                    end;
                catch
                err=true;
                while err
                    disp('error in IntersectMeshes_OK ');
                    err=false;
                    Mnew = AddNoise2Mesh(k,M);
                    try
                        [M_out,M_in,M_on,M_no,N_out,N_in]=IntersectMeshes_OK(Mnew,N);
                        switch op_bool
                            case 0, [varargout{1}]=AppendMeshes(M_out,N_out); 
                                    [varargout{1}]=AppendMeshes(varargout{1},M_on); % union
                            case 1, [varargout{1}]=AppendMeshes(M_in,N_in); 
                                    [varargout{1}]=AppendMeshes(varargout{1},M_on); % interseccion
                            case 2, N_in.tri=N_in.tri(:,[1 3 2]);
                                    [varargout{1}]=AppendMeshes(M_out,N_in); 
                                    [varargout{1}]=AppendMeshes(varargout{1},M_no);% diferencia
                        end;
                    catch
                        err=true;
                    end;
                end;
            end;
            end;
       case 'carve'
            try
                [varargout{1}]=carveIntersect(M,N,op_carve);
            catch
                err=true;
                while err
                    disp(['error in carveIntersect ( op = ' op ' ). ']);
                    err=false;
                    Mnew = AddNoise2Mesh(k,M);
                    try
                        [varargout{1}]=carveIntersect(Mnew,N,op_carve);
                    catch
                        err=true;
                    end;
                end;
            end;
            
            
        case 'triang'
            try
               [M_out,M_in,M_on,M_no,N_out,N_in]=IntersectMeshes_OK(M,N);
                switch op_bool
                    case 0, [varargout{1}]=AppendMeshes(M_out,N_out); 
                            [varargout{1}]=AppendMeshes(varargout{1},M_on); % union
                    case 1, [varargout{1}]=AppendMeshes(M_in,N_in); 
                            [varargout{1}]=AppendMeshes(varargout{1},M_on); % interseccion
                    case 2, N_in.tri=N_in.tri(:,[1 3 2]);
                            [varargout{1}]=AppendMeshes(M_out,N_in);
                            [varargout{1}]=AppendMeshes(varargout{1},M_no);% diferencia
                end;
            catch
                err=true;
                while err
                    disp('error in IntersectMeshes_OK ');
                    err=false;
                    Mnew = AddNoise2Mesh(k,M);
                    try
                        [M_out,M_in,M_on,M_no,N_out,N_in]=IntersectMeshes_OK(Mnew,N);
                        
                        switch op_bool
                            case 0, [varargout{1}]=AppendMeshes(M_out,N_out); 
                                    [varargout{1}]=AppendMeshes(varargout{1},M_on); % union
                            case 1, [varargout{1}]=AppendMeshes(M_in,N_in); 
                                    [varargout{1}]=AppendMeshes(varargout{1},M_on); % interseccion
                            case 2, N_in.tri=N_in.tri(:,[1 3 2]);
                                    [varargout{1}]=AppendMeshes(M_out,N_in); 
                                    [varargout{1}]=AppendMeshes(varargout{1},M_no);% diferencia
                        end;
                    catch
                        err=true;
                    end;
                end;
            end;

                   
    end; 
    
else
    switch method
        case 'normal'
            try
                if nargout<=3
                    [varargout{1:nargout}]=carveIntersect(M,N,'split');
                else
                    [varargout{1:nargout}]=carveIntersect(M,N,'splitFull');
                end;
            catch
            
                disp('error in carveIntersect ( op = split ). Trying with triangle');
                try
                    [varargout{1:nargout}]=IntersectMeshes_OK(M,N);

                catch
                    err=true;
                    while err
                        disp('error in IntersectMeshes_OK ');
                        err=false;
                        Mnew = AddNoise2Mesh(k,M);
                        try
                            [varargout{1:nargout}]=IntersectMeshes_OK(Mnew,N);
                        catch
                            err=true;
                        end;
                    end;
                end;
            end;
        case 'carve'
            try
            if nargout<=3
                [varargout{1:nargout}]=carveIntersect(M,N,'split');
            else
                [varargout{1:nargout}]=carveIntersect(M,N,'splitFull');
            end;
            catch
                err=true;
                while err
                    disp('error in carveIntersect ( op = split ). ');
                    err=false;
                    Mnew = AddNoise2Mesh(k,M);
                    try
                        if nargout<=3
                            [varargout{1:nargout}]=carveIntersect(Mnew,N,'split');
                        else
                            [varargout{1:nargout}]=carveIntersect(Mnew,N,'splitFull');
                        end;
                    catch
                        err=true;
                    end;
                end;
            end;
        case 'triang'
            try
                [varargout{1:nargout}]=IntersectMeshes_OK22(M,N);

            catch
                err=true;
                while err
                    disp('error in IntersectMeshes_OK ');
                    err=false;
                    Mnew = AddNoise2Mesh(k,M);
                    try
                        [varargout{1:nargout}]=IntersectMeshes_OK(Mnew,N);
                    catch
                        err=true;
                    end;
                end;
            end;
            
            
            
    end;
            
end;


for i=1:nargout
    if (~isfield(varargout{i},'tri'))||isempty(varargout{i}.tri)
        varargout{i}=struct('xyz',[],'tri',[]);
    end;
    if ~isempty(varargout{i}.tri)
        if clean
        [varargout{i}]=vtkCleanPolyData(varargout{i},'SetAbsoluteTolerance',tol,'ToleranceIsAbsoluteOn',[],'PointMergingOn',[] ,'ConvertPolysToLinesOff',[], 'ConvertLinesToPointsOff',[] );
        end;
        if fix
        varargout{i}.tri = varargout{i}.tri( ~ ( any( ~varargout{i}.tri , 2 ) | varargout{i}.tri(:,1) == varargout{i}.tri(:,2) | varargout{i}.tri(:,2) == varargout{i}.tri(:,3) | varargout{i}.tri(:,1) == varargout{i}.tri(:,3) ), : );
        end;
        if removedup
        [varargout{i}.tri]=(unique(varargout{i}.tri,'rows'));
        end;
    end;
end;

