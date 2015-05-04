function L = MeshToboggan( M , D )

  if ~isfield( M , 'ESUP' ) || ~isfield( M , 'PSUP' )
    M= CreateESUP( CreatePSUP( M ));
  end

  if numel(D) ~= size(M.xyz,1) , error('D tiene que tener otrp tamanio'); end

  L = D(:)+NaN;
  lastlabel = 0;

  
  nolabeled = find(isnan(L));
  while ~isempty(nolabeled)
    labelnode(nolabeled(1));
    nolabeled = find(isnan(L));
  end
  
  function label = labelnode(node)
    neig = PSUP( node , M );
    d = D(node);
    dd = D(neig);
    if ~any(dd<d)
      lastlabel = lastlabel+1;
      label = lastlabel;
    else
      for n = neig
        if D(n) < d
          l=L(n);
          if isnan(l)
            label = labelnode(n);
          else
            label = l;
          end
        end
      end
    end
    L(node)= label;
  end

end
