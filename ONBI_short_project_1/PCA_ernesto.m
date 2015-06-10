DATA = training_dia_sys_myo_shapes';
m = mean(DATA,2);
RES = bsxfun(@minus,DATA,m);
[u,s,v] = svd(RES, 'econ');
diag(s)

Bcoefs = u.' * RES;

max( abs( m + u*Bcoefs(:,2) - DATA(:,2) ))

trainingdata = [Bcoefs(1:100,1:100)' ; Bcoefs(1:100,101:200)']