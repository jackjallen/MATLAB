function exploretensorvariability( T , Vdims , fcn , names )

%   Tdims = 1:ndims(T);
%   Tdims = setdiff( Tdims, Vdims );
  dims  = 1:ndims(T);
  
  for ii= 1:numel(Vdims)
    i= Vdims(ii);
    try
      name = names{ii};
      var{i} = eEntry('range', [1 size(T,i)],'step',1,'ReturnFcn', @(x) round(x),'slider2edit', @(x) sprintf('%s: %d',name, round(x)),'callback',@(x) call ,'ivalue',1 );
    catch
      var{i} = eEntry('range', [1 size(T,i)],'step',1,'ReturnFcn', @(x) round(x),'slider2edit', @(x) sprintf('%d',round(x)),'callback',@(x) call ,'ivalue',1 );
    end
    p = get( var{i}.panel , 'Position' );
    set( var{i}.panel , 'Position' , [ 0 (numel(Vdims)-ii)*p(4) p(3) p(4) ] );
    drawnow;
    var{i}.continuous  = 1;
  end
  
  call;
  
  function call
    s={};
    for d= dims
      if ismember(d,Vdims)
        s{d}= var{d}.value;
      else
        s{d}= ':';
      end
    end
    feval( fcn , T( s{:} ) );
  end


end
