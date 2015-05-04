function [E dE]=Meshes_Energy_Function_Rigid(M1,M2,transformModel,p0,func,met)


np=numel(p0);

[M dM]=maketransform(transformModel,p0);


TpM2 = transformMesh(M2,M);


a=FastData(BooleanMeshes(TpM2,'i',M1,met));
b=FastData(M1)-a;
c=FastData(TpM2)-a;



switch func
    case {'int' 'inter' 'intersection' 'intersect' 'interseccion'}, 
        E=a;
    case {'union' 'junction' 'joint'}, 
        E=a+b+c;
    case {'jaccard' 'jac' 'jacard'}
        E=1-a/(a+b+c);
    case {'-int' '-inter' '-intersection' '-intersect' '-interseccion'}
        E=-a;
    case {'hamming','ham'}
        E=b+c;
end;


if nargout>1
    
    iM = maketransform(transformModel,p0,'inv');
    
    [M_OUT,M_IN]=BooleanMeshes(TpM2,'i',M1,'clean');
    
        
    d_a=zeros(1,np);
    d_b=zeros(1,np);
    d_c=zeros(1,np);
    d_d=zeros(1,np);
    
    % OUT
    
    ip_OUT = transform(M_OUT.xyz,iM);
    [p_OUT dp_OUT]  = transform(ip_OUT,{M,dM});
    
    [kk,A_OUT,N_OUT]=FastData(M_OUT);
    N=size(p_OUT,1);
    dTp=dp_OUT(vec([1:N;N+1:2*N;2*N+1:3*N]),:);
    dTp=mat2cell(dTp, 3*ones(N,1), np);
    for j=1:size(M_OUT.tri,1)
        inc=A_OUT(j)*N_OUT(j,:)*(dTp{M_OUT.tri(j,1)}+dTp{M_OUT.tri(j,2)}+dTp{M_OUT.tri(j,3)})/3;
        d_c=d_c+inc;
    end;
    
    % IN
    ip_IN = transform(M_IN.xyz,iM);
    [p_IN dp_IN]  = transform(ip_IN,{M,dM});
    
    [kk,A_IN,N_IN]=FastData(M_IN);
    N=size(p_IN,1);
    dTp=dp_IN(vec([1:N;N+1:2*N;2*N+1:3*N]),:);
    dTp=mat2cell(dTp, 3*ones(N,1), np);
    for j=1:size(M_IN.tri,1)
        inc=A_IN(j)*N_IN(j,:)*(dTp{M_IN.tri(j,1)}+dTp{M_IN.tri(j,2)}+dTp{M_IN.tri(j,3)})/3;
        d_a=d_a+inc;
        d_b=d_b-inc;
    end;
    
switch func
    case {'int' 'inter' 'intersection' 'intersect' 'interseccion'}, 
         dE=d_a;
    case {'union' 'junction' 'joint'}, E=area(union(M1,TpM2));
        dE=d_a+d_b+d_c;
    case {'jaccard' 'jac' 'jacard'}
        dE=-d_a/(a+b+c)+a*(d_a+d_b+d_c)/(a+b+c)^2;
    case {'-int' '-inter' '-intersection' '-intersect' '-interseccion'}
        dE=-d_a;
    case {'hamming','ham'}
        dE=d_b+d_c;
end;

end;