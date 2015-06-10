function [varargout]=IntersectBooleanMeshes(M,op,N,varargin)
% function [varargout]=IntersectBooleanMeshes(M,op,N,varargin)
%           op:
%           'carve':
%           'bool':
%           'NoFix':
%           'RemoveDup':
%           'Clean':
%           'PreReord':
%           'Tol':

method='normal';
[varargin,method]       = parseargs(varargin,'Carve','$FORCE$',{'carve',method});
[varargin,method]       = parseargs(varargin,'Bool','$FORCE$',{'bool',method});
[varargin,fix]          = parseargs(varargin,'NoFix','$FORCE$',{false,true});
[varargin,removedup]    = parseargs(varargin,'RemoveDup','$FORCE$');
[varargin,clean]        = parseargs(varargin,'Clean','$FORCE$');
[varargin,istol,tol]    = parseargs(varargin,'Tol','$FORCE$');
[varargin,prereord]     = parseargs(varargin,'PreReord','$FORCE$');
if ~istol, tol=1e-6; end;

if prereord
    M=vtkPolyDataNormals(M,'ComputePointNormalsOff', [], 'ComputeCellNormalsOn'     , [], 'ConsistencyOn'            , []);
    N=vtkPolyDataNormals(N,'ComputePointNormalsOff', [], 'ComputeCellNormalsOn'     , [], 'ConsistencyOn'            , []);
    M=rmfield(M,'triNormals');
    N=rmfield(N,'triNormals');
end;
if nargout<=1
    switch method
        case 'normal'
            try
                [varargout{1}]=carveIntersect(M,N,'i');
            catch
                [varargout{1}]=vtkBooleanOperationPolyDataFilter(M,N,'SetOperationToIntersection',[]);
            end;
        case 'carve'
            try
                [varargout{1}]=carveIntersect(M,N,'i');
            catch
                [varargout{1}]=[];
            end;
        case 'bool'
            [varargout{1}]=vtkBooleanOperationPolyDataFilter(M,N,'SetOperationToIntersection',[],'SetTolerance',tol);
    end;    
else
    [varargout{1:nargout}]=IntersectMeshes(M,N);
end;

if clean
    for i=1:nargout
        [varargout{i}]=vtkCleanPolyData(varargout{i},'SetAbsoluteTolerance',tol,'SetToleranceIsAbsolute',true,'SetPointMerging',true ,'SetConvertPolysToLines',true, 'SetConvertLinesToPoints',true );
    end;
end;


if fix
    for i=1:nargout
       
        varargout{i}.tri = varargout{i}.tri( ~ ( any( ~varargout{i}.tri , 2 ) | varargout{i}.tri(:,1) == varargout{i}.tri(:,2) | varargout{i}.tri(:,2) == varargout{i}.tri(:,3) | varargout{i}.tri(:,1) == varargout{i}.tri(:,3) ), : );
        
    end;
end;

if removedup
    for i=1:nargout
        [varargout{i}.tri]=(unique(varargout{i}.tri,'rows'));
    end;
end;