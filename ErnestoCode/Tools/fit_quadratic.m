function a = fit_quadratic( aa , ee )

%       while numel(aa) > 3
%         [mm,mmid] = max(ee);
%         ee(mmid)=[];
%         aa(mmid)=[];
%       end
  
  polyW = diag( diff( dualVector(aa) ) );
  polyM = [ ones(numel(aa),1) aa(:) aa(:).^2 ];

  setWarning('off','MATLAB:rankDeficientMatrix');
  setWarning('off','MATLAB:nearlySingularMatrix');
  setWarning('off','MATLAB:singularMatrix');
  ABC = linsolve( polyM.' * polyW * polyM ,  polyM.' * polyW * ee(:) );
  restoreWarning(  'MATLAB:rankDeficientMatrix');
  restoreWarning(  'MATLAB:nearlySingularMatrix');
  restoreWarning(  'MATLAB:singularMatrix');

  ABC(4) = 0;
  
  if ABC(3) <= 0
    a = Inf;
    return;
  end

  a = -0.5 * ABC(2)/ABC(3);
end
