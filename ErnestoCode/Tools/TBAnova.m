function [pval,t,stats] = TBAnova(y,group,varargin)
% function to replace MATLAB anovan for multidimensional data vector y
% see anovan help function for sintaxis


% First:  argument checking and selection of defaults
error(nargchk(2,Inf,nargin,'struct'));
dostats = (nargout>=3);      % true to compute output stats structure


if ~iscell(group) && ~isa(group,'categorical')
    if isnumeric(group) && size(group,1)==size(y,1)
        group = num2cell(group,1);
    else
        error('stats:anovan:BadGroup',...
             'GROUP must be a cell array or matrix of grouping variables.')
    end
end
group = group(:);
ng = length(group);

% Which style of calling sequence are we using?
newstyle = true;
if ~isempty(varargin) && nargin<=6
   va1 = varargin{1};
   if ~ischar(va1) || isempty(va1)
      newstyle = false;
   elseif ~isempty(strmatch(va1,{'l' 'i' 'f' 'linear' 'interaction' 'full'},...
                            'exact'))
      newstyle = false;
   end
end

if newstyle
   % Calling sequence with named arguments
   okargs =   {'model' 'sstype' 'varnames' 'display' 'alpha' ...
               'random'     'continuous' 'nested'};
   defaults = {1       3        ''         'on'      0.05    ...
               false(ng,1)  []           []};
   [eid emsg model sstype varnames displayopt alpha ...
             randomvar continuous vnested] = ...
                internal.stats.getargs(okargs,defaults,varargin{:});
   if ~isempty(eid)
      error(sprintf('stats:anovan:%s',eid),emsg);
   end
else
   % Old-style calling sequence with fixed argument positions
   model = 1;
   sstype = 3;
   varnames = '';
   displayopt = 'on';
   continuous = [];
   vnested = [];
   if nargin>=3, model = varargin{1}; end
   if nargin>=4, sstype = varargin{2}; end
   if nargin>=5, varnames = varargin{3}; end
   if nargin>=6, displayopt = varargin{4}; end
   alpha = 0.05;
   randomvar = false(ng,1);
end

% Check optional arguments
[varnames,termlist,randomvar,sstype,continuous,vnested] = argcheck(ng,...
                    varnames,model,randomvar,sstype,alpha,continuous,vnested);
doems = any(randomvar);        % true to compute expected mean squares
if isequal(lower(sstype),'h')
    sstype = 2;
    dohier = true;
else
    dohier = false;
end
if any(vnested(:))
    varnames = makenestednames(varnames,vnested);
end



% STEP 1:  Remove NaN's and prepare grouping variables.
% Also, make sure all groups are rows, and create group index and name arrays.
[y,allgrps,nanrow,varinfo] = removenans(y,group,continuous);
[n,NN] = size(y);

% STEP 2:  Create dummy variable arrays for each grouping variable.
% STEP 2a: For type 3 ss, create constraints for each grouping variable.
varinfo.varnames = varnames;
varinfo = makedummyvars(continuous,allgrps,vnested,varinfo);

% STEP 3:  Create dummy variable arrays for each term in the model.
nterms = size(termlist,1);
fulltermlist = showtermnesting(termlist,vnested);
tnested = gettermnesting(fulltermlist,vnested,continuous);

[~,sindex] = sortrows(fulltermlist);
terminfo = makedummyterms(sindex, termlist, varinfo, ...
                          continuous, tnested);

% STEP 4:  Create the full design matrix
[dmat,cmat,termname] = makedesignmat(terminfo,n);

% STEP 5:  Fit the full model and compute the residual sum of squares
% % % % Find the null space of the constraints matrix
% % % [Qc,Rc,Ec] = qr(cmat'); pc = Rrank(Rc); Qc0 = Qc(:,pc+1:end);
% % % 
% % % % Do qr decomposition on design matrix projected to null space
% % % Dproj = dmat*Qc0; [Qd,Rd,Ed] = qr(Dproj,0);
% % % dfx = Rrank(Rd); Qd = Qd(:,1:dfx); Rd = Rd(1:dfx,1:dfx);

% STEP 5:  Fit the full model and compute the residual sum of squares
mu =  mean(y); y = bsxfun(@minus, y, mu); % two passes to improve accuracy
mu2 =  mean(y); mu = mu + mu2;  y = bsxfun(@minus, y, mu2); 
sst = sum(y.^2);
if ~dostats
   [ssx,dfx,dmat2,y] = dofit(dmat, y, cmat, sst, mu);
else
   % Do initial fit requesting stats structure, then convert from full-rank
   % to overdetermined form, then find a label for each coefficient
   [ssx,dfx,dmat2,y,stats] = dofit(dmat,y,cmat,sst,mu);
   [codes,names] = convertstats(stats,terminfo,varinfo,continuous);
end
   
sse = sst - ssx;
dfe = n - dfx;

% % % 
% % % y2 = Qd' * y;            % rotate y into that space
% % % % zq = Rd \ y2;            % fit rotated y to projected design matrix
% % % 
% % % z = zeros(length(Ed),NN); % coefficient vector extend to full size ...
% % % z(Ed(1:dfx),:) = Rd \ y2;       % zq ... and in correct order for null space
% % % 
% % % % b = ;             % coefficients back in full space
% % % ssx = sum(y2.^2); % norm(y2)^2;        % sum of squares explained by fit
% % % 
% % % dmat2 = Qd' * dmat;
% % % 
% % % if (dostats)
% % % 
% % % 
% % %    stats.source = 'anovan';
% % %    
% % %    % Get residuals
% % %    yhat = Dproj * z;
% % %    stats.resid = y - yhat;
% % %    
% % %    % Calculate coefficients, then adjust for previously removed mean
% % %    stats.coeffs = Qc0 * z; % coefficients back in full space
% % %    if ~isempty(stats.coeffs)
% % %       stats.coeffs(1,:) = stats.coeffs(1,:) + mu;            % undo mean adjustment before storing
% % %    end
% % % %    stats.coeffs = t;
% % % 
% % %    RR = zeros(dfx,size(Dproj,2));
% % %    RR(:,Ed(1:dfx)) = Rd;
% % %    stats.Rtr = RR';
% % %    [qq,rr,ee] = qr(dmat,0);
% % %    pp = Rrank(rr);
% % %    rowbasis = zeros(pp,size(rr,2));
% % %    rowbasis(:,ee) = rr(1:pp,:);
% % %    stats.rowbasis = rowbasis;
% % %    stats.dfe = n - dfx;
% % %    stats.mse = (sst - ssx) / max(1, stats.dfe);
% % %    stats.nullproject = Qc0;
% % % 
% % % end
% % % 
% % % y = y2;
% % % 
% % % 
% % % sse = sst - ssx;
% % % dfe = n - dfx;

% STEP 6:  Determine which models to compare for testing each term
ssw  = -ones(nterms, NN);      % sum of squares with this term
sswo = ssw;                   % sum of squares without this term
dfw  = -ones(nterms, 1);                   % residual d.f. with this term
dfwo = dfw;                   % residual d.f. without this term


switch sstype
    case 1
        modw = tril(ones(nterms)); % get model with this term
        k = nterms;                % locations of model with all terms
    case 2
        modw = termsnotcontained(fulltermlist,continuous,dohier);
        k = (sum(modw,2) == nterms);
        TnotC = modw;
    case 3
        modw = ones(nterms);
        k = 1:nterms;
        TnotC = termsnotcontained(fulltermlist,continuous,dohier);
end
modw = logical(modw);
modwo = logical(modw - eye(nterms));   % get model without this term

% STEP 7:  Fit each model, get its residual SS and d.f.
ssw(k,:) = repmat(ssx,nterms,1);        % for full model we already know the results
dfw(k) =dfx;

% Fit each model separately
dfboth = [dfw; dfwo];
if doems
   Llist = cell(nterms,1);
end

% Consider interactions before their components for type 3 ss, so
% examine the terms in decreasing order of the number factors in the term
if sstype==3
   termsize = sum(fulltermlist,2)';
   [~,sindices] = sort(termsize);
   sindices = [sindices(:); sindices(:)];
else
   sindices = [(1:size(termlist,1)), (1:size(termlist,1))]';
end

% Fit all required subsets of the full model
fitsubmodels;   % nested function

% STEP 7A:  Compute expected means squares if requested
if doems
   emsMat = makeemsmat(dmat, cmat, nterms, Llist, termname);
end

% STEP 7A:  Compute expected means squares if requested
if doems
   emsMat = makeemsmat(dmat, cmat, nterms, Llist, termname);
end

% STEP 8:  Compute the sum of squares attributed to each term
ssterm = ssw - sswo;
dfterm = dfw - dfwo;
ssterm(dfterm==0) = 0;    % make this exact

% STEP 9:  Compute the mean square for each term
msterm = bsxfun(@rdivide,ssterm, max(1, dfterm)); %ssterm ./ max(1, dfterm);
mse = sse * (dfe>0) / max(1, dfe);

% STEP 9a:  Equate computed and expected mean squares, then solve for
%           variance component estimates
randomterm = find(termlist*randomvar > 0);  % indices of random terms
if isempty(randomterm)
   msdenom = repmat(mse, size(msterm,1),1);
   dfdenom = repmat(dfe, size(msterm,1),size(msterm,2));
else
   [msdenom,dfdenom,varest,varci,txtdenom,txtems,denommat,rtnames] = ...
                getRandomInfo(msterm,dfterm,mse,dfe,emsMat,randomterm,...
                              terminfo.tnames,alpha);
end


% STEP 10:  Compute an F statistic for each term
fstat = Inf(size(msterm));
fstat(msdenom>0) = msterm(msdenom>0) ./ msdenom(msdenom>0);
fstat(~(msdenom>0 | msterm>0)) = NaN;

% STEP 11:  Compute the p-value for each term
pval = nan(size(fstat));
t = bsxfun(@and, dfdenom>0 , repmat(dfterm,1,size(dfdenom,2))>0);
%pval(t) = fpval(fstat(t), dfterm(t), dfdenom(t));
pval(t) = fcdf(1./fstat,dfdenom, repmat(dfterm,1,size(t,2)));
% if isempty(randomterm)
%     t = bsxfun(@and, dfdenom>0 , dfterm>0);
%     for ii=find(t)'
%         pval(ii,:) = 1 - fcdf(fstat(ii,:), repmat(dfterm(ii),1,size(dfdenom,2)), dfdenom(ii,:));
%     end
% else
%     pval = (1-fcdf(fstat, repmat(dfterm,1,size(dfdenom,2)), dfdenom));
% end
% p = pval;

% Create anova table as a cell array
sing = (dfterm < terminfo.dfterm0);  % is non-nested term singular?

[tbl, t] = maketable;   % nested function
if ~isempty(tnested)
    sing(any(tnested,2)) = false;
end

if (nargout > 3), terms = termlist; end


% Store additional information for multiple comparisons
if dostats
   stats = makestats(stats);  % nested function
end


% ----- nested function
function stats = makestats(stats)
   stats.terms = termlist;
   stats.nlevels = varinfo.vlevels;
   stats.nlevels(continuous) = 1;
   stats.continuous = continuous;
   stats.vmeans = varinfo.vmeans;
   stats.termcols = [1; terminfo.termlength];
   stats.coeffnames = names;
   stats.vars = codes;
   stats.varnames = varinfo.varnames; % names of grouping variables
   stats.grpnames = varinfo.allnames; % names of levels of these variables
%    if length(stats.resid)<length(nanrow)
%       r(~nanrow) = stats.resid;
%       r(nanrow) = NaN;
%       stats.resid = r;
%    end
   stats.vnested = [];

   % Extra info for models with random effects
   if doems
      stats.ems = emsMat;
      stats.denom = denommat;
      stats.dfdenom = dfdenom;
      stats.msdenom = msdenom;
      stats.varest = varest;
      stats.varci = varci;
      stats.txtdenom = txtdenom;
      stats.txtems = txtems;
      stats.rtnames = rtnames;
   else
      stats.ems = [];
      stats.denom = [];
      stats.dfdenom = [];
      stats.msdenom = [];
      stats.varest = [];
      stats.varci = [];
      stats.txtdenom = [];
      stats.txtems = [];
      stats.rtnames = [];
   end
end


% ----- nested function
function [tbl,t] = maketable

   tbl = [ssterm(:,1) dfterm sing msterm(:,1) fstat(:,1) pval(:,1); sse(:,1) dfe 0 mse(:,1) NaN NaN];
   tbl = num2cell(tbl);
   tbl(end,end-1:end) = {[]};
   tbl = [{'Source' 'Sum Sq.' 'd.f.' 'Singular?' 'Mean Sq.' 'F' 'Prob>F'}; ...
          terminfo.tnames tbl];

   % Add random effects information if required
   if ~isempty(randomterm)
     rterms = [randomterm; nterms+1];
      C = size(tbl,2);
      tbl(1,end+1:end+8) = {'Type' 'Expected MS' 'MS denom' 'd.f. denom' 'Denom. defn.' ...
                            'Var. est.' 'Var. lower bnd' 'Var. upper bnd'};
      tbl(1+(1:nterms),C+1)          = {'fixed'};
      tbl(1+rterms,C+1)              = {'random'};
      tbl(1+(1:length(txtems)),C+2)  = txtems;
      tbl(1+(1:length(msdenom)),C+3) = num2cell(msdenom(:));
      tbl(1+(1:length(dfdenom)),C+4) = num2cell(dfdenom(:));
      tbl(1+(1:length(txtdenom)),C+5)= txtdenom;
      tbl(1+rterms,C+6)              = num2cell(varest);
      tbl(1+rterms,C+7)              = num2cell(varci(:,1));
      tbl(1+rterms,C+8)              = num2cell(varci(:,2));
   end

   % So far this is the table for the output argument
   t = tbl;
   
   for i=1:size(ssterm,1)
       t{i+1,1+1} = ssterm(i,:);
       t{i+1,4+1} = msterm(i,:);
       t{i+1,5+1} = fstat(i,:);
       t{i+1,6+1} = pval(i,:);
   end
   t{end,2}=sse; t{end,4+1}=mse;
   
   % Add total information
   tbl{end+1,1} = 'Total';
   tbl{end,2} = sst(1);
   tbl{end,3} = n-1;
   tbl{end,4} = 0;
   

   t{end+1,1} = 'Total';
   t{end,2} = sst;
   t{end,3} = n-1;
   t{end,4} = 0;
      
 
   % For display, mark terms that are rank deficient
   t{1,1} = ['  ' t{1,1}];
   for r=2:size(t,1)
      x = t{r,1};
      if (t{r,4})
         x = ['# ' x];
      else
         x = ['  ' x];
      end
      t{r,1} = x;
   end
   t = t(:,[1 2 3 4 5 6 7]);


end % of nested function maketable



% ----- nested function
function fitsubmodels
   for j=length(sindices):-1:1
      % Find the next model index to fit
      k = sindices(j);
      
      % Look in unsorted arrays to see if we have already fit this model
      if j>nterms
         k0 = k+nterms;
      else
         k0 = k;
      end
      if dfboth(k0)~=-1
         continue
      end
      
      % Find the model with this index
      if (j > nterms)
         thismod = modwo(k, :);
      else
         thismod = modw(k, :);
      end
      
      % Get the design matrix for this model
      keepterms = find(thismod);
      clist = ismember(termname, [0 keepterms]);
      X = dmat2(:,clist);
      C = cmat(:,clist);

      % Fit this model
      [ssx0,dfx0] = dofit(X, y, C);
      
      % If this model is the "without" model for a term, compute L
      if doems && j>nterms
         if sstype==1
            prevterms = (1:nterms)<k;
         else
            prevterms = TnotC(k,:);
         end
         termscontained = ~prevterms;
         prevterms(k) = false;
         prevcols = ismember(termname, [0 find(prevterms)]);
         curcols = (termname == k);
         if sstype==3
            prevL = cat(1,Llist{termscontained});
         else
            prevL = [];
         end
         Llist{k} = getL(dmat, curcols, prevcols, prevL);
      end
      
      % Use these results for each term that requires them
      mod0 = repmat(thismod, nterms, 1);
      k = find(all(modw == mod0, 2));
      if ~isempty(k) 
          ssw(k,:) = ssx0; 
          dfw(k) = dfx0;
          dfboth(k) = 0;
      end
      k = find(all(modwo == mod0, 2));
      if ~isempty(k)
          sswo(k,:) = ssx0;
          dfwo(k) = dfx0;
          dfboth(nterms+k) = 0;
      end
   end
end









end % of main function

% -----------------------------
function [msDenom,dfDenom,varest,varci,txtdenom,txtems,denommat,rtnames] = ...
                      getRandomInfo(msTerm,dfTerm,msErr,dfErr,ems,...
                                    randomterms,tnames,alpha)
%GETRANDOMINFO Get info for random effects, such as denominator for F tests

NN = size(msTerm,2);
nterms = size(msTerm,1);
msDenom = repmat(NaN,nterms,NN);
dfDenom = repmat(NaN,nterms,NN);
txtdenom = cell(nterms,1);
txtems = cell(nterms+1,1);
txtems{end} = 'V(Error)';

% For convenience, combine error info with term info
msTerm = [msTerm; msErr];
dfTerm = [dfTerm; dfErr];
randomterms = [randomterms; size(msTerm,1)];

% Find the matrix giving expected mean squares for random rows, and
% compute the coefficients of these expected values.
randrows = ems(randomterms,:);
emsmat = ems(1:end-1,:);
emsmat(1:size(emsmat,1)+1:end) = 0;
[q,r,e] = qr(randrows',0);
p = Rrank(r);
B = r(1:p,:)\(q(:,1:p)'*emsmat');

% Loop over all terms
for j=1:nterms
   % First compute F statistic denominator information
   % See which terms are involved in this combination
   b = B(:,j);
   tol = sqrt(eps(class(b))) * max(abs(b));
   nonzeroterms = find(abs(b) >= tol);
   indices = e(nonzeroterms);
   tnums = randomterms(indices);
   if length(nonzeroterms) == 1
      % If it's a single term, use it exactly
      msDenom(j,:) = msTerm(tnums,:);
      dfDenom(j,:) = dfTerm(tnums);
      txtdenom{j} = sprintf('MS(%s)',tnames{tnums});
   else
      linprod = zeros(numel(tnums), NN);
      % Otherwise use the linear combination and Satterthwaite's approximation
      for ii = 1:numel(tnums)
          linprod(ii,:) = b(nonzeroterms(ii)) .* msTerm(tnums(ii),:);
      end
      msDenom(j,:) = sum(linprod,1);
      if any(dfTerm(tnums) == 0)
          dfDenom(j,:) = NaN;
      else
          dfDenom(j,:) = msDenom(j,:).^2 ./ sum(bsxfun(@rdivide,linprod.^2, dfTerm(tnums)),1);
%           dfDenom(j) = msDenom(j)^2 / sum((linprod.^2) ./ dfTerm(tnums));
      end
      txtdenom{j} = linprodtxt(b(nonzeroterms),tnames(tnums),'MS');
   end
   
   % Next get an expression for the expected mean square
   prefix = repmat('Q',nterms+1,1); % quadratic form for fixed effects
   prefix(randomterms,1) = 'V';     % variance for random effects
   txtems{j} = linprodtxt(ems(j,:),tnames,prefix);
end

% Get the variance component estimates
%  [varest,varci] = varcompest(ems,randomterms,msTerm,dfTerm,alpha);
varest = []; % nan(2, NN);
varci = []; % nan(2, 2, NN);

% Package up the results
denommat = B';
rtnames = tnames(randomterms);
end

% ---------------------
function p = Rrank(R,tol)
%RRANK Compute rank of R factor after qr decomposition
if (min(size(R))==1)
   d = abs(R(1,1));
else
   d = abs(diag(R));
end
if nargin<2
   tol = 100 * eps(class(d)) * max(size(R));
end
p = sum(d > tol*max(d));
end

% --------------------------
function m = termsnotcontained(terms,vnested,continuous,dohier)
%TERMSNOTCONTAINED Creates a logical matrix indicating term containment
%   m(i,j)==1  iff  t(i) is not contained by t(j)

% For starters, omit continuous variables from this test
contnum = find(continuous);
t = terms;
t(:,contnum) = 0;

% The main test:  no overlap between terms included in t(i) and not in t(j)
m = (t*~t') > 0;

% Now consider continuous variables (X)
nterms = size(terms,1);
for j=1:length(contnum)
   v = terms(:,contnum(j));
   vv = repmat(v,1,nterms);
   if dohier
      % A cannot be contained in B if A has a higher power of X
      m = m | vv>vv';
   else
      % A cannot be contained in B if they don't match on X
      m = m | vv~=vv';
   end
end

% Set diagonals to 1 because we want proper containment
m(1:(nterms+1):end) = 1;
end

% ----------------------
function ems = getems(L,R,termnums)
%GETEMS Final step in expected mean square computation
%   The L array is an array of constraints on the coefficients,
%   R is the R factor in a QR decomposition X, and termnums is a
%   vector mapping R column numbers to term numbers.  The output is
%   a vector of expected mean square coefficients with one element
%   per term.  The first element is for the constant term, and
%   the last is for the error term.  See "Computing Expected Mean
%   Squares" by Goodnight and Speed, 1978.

% In G&S notation, form the matrix L*pinv(X'*X)*L'.
dfL = size(L,1);
temp = L/R;
LML = temp*temp';

% Create the Cholesky decomposition U and then compute U\L.
U = cholcov(LML);
C = U'\L;

% Sum the squared values of the columns other than the column for
% the constant term, and assign the sums to terms.
sumsq = sum(C(:,2:end).^2,1);
cols = termnums(2:end);
rows = ones(size(cols));
ems = accumarray([rows(:),cols(:)], sumsq) / dfL;
ems(:,end+1) = 1;  % add entry for error term
end

% --------------------------------
function terminfo = makedummyterms(sindex,termlist,varinfo,continuous,tnested)

ncols = 1;
[nterms,nfactors] = size(termlist);

termdum = cell(nterms, 1);      % dummy variable = design matrix cols
termconstr = cell(nterms,1);    % constraints to make each term well defined
levelcodes = cell(nterms, 1);   % codes for levels of each M row
tnames = cell(nterms, 1);       % name for entire term, e.g. A*B
dfterm0 = zeros(nterms, 1);     % nominal d.f. for each term
termvars = cell(nterms, 1);     % list of vars in each term
termlength = zeros(size(sindex));% length of each term (number of columns)

auxtermlist = zeros(0,nfactors);  % may need to create temp terms for nesting
auxtermdum = cell(0,1);           % may need their dummy variables

isnested = false;

% For each term,
for j=1:nterms
   % Get dummy columns, var list, term name, etc
   sj = sindex(j);
   tm = termlist(sj,:);
   if ~isempty(tnested)
      isnested = any(tnested(sj,:));
   end
   [tdum,vars,tn,df0,tconstr] = maketerm(tm,isnested,varinfo,j,sindex,...
                                         termlist,termdum,termconstr,tnames,dfterm0);

   % Store this term's information
   k = size(tdum, 2);
   termlength(sindex(j),1) = k;
   ncols = ncols + k;
   termdum{sj} = tdum;
   termvars{sj} = vars;
   levelcodes{sj} = fliplr(fullfact(varinfo.vlevels(vars(end:-1:1))));
   tnames{sj,1} = tn;
   dfterm0(sj) = df0;

   % For a nested term, figure out the constraint now
   if isnested
      % First get dummy variables for this term and its nesters
      if any(tm(continuous))
         % Ignore any continuous variable contributions to this term
         tmcat = tm;                % this term, cat predictors only
         tmcat(continuous) = 0;
         ntmcat = termlist(tnested(sj,:),:);
         ntmcat(:,continuous) = 0;  % nesting terms, cat predictors only
         [Ydum,auxtermlist,auxtermdum] =findtermdum(tmcat,termlist,termdum,...
                              auxtermlist,auxtermdum,varinfo);
         Ydum = Ydum{1};
         [Xdum,auxtermlist,auxtermdum] =findtermdum(ntmcat,termlist,termdum,...
                              auxtermlist,auxtermdum,varinfo);
      else
         % No continuous variables, so work with the whole term
         Ydum = tdum;
         Xdum = termdum(tnested(sj,:)>0);
      end

      % Constrain the part of this term in its nesters to be 0
      tconstr = gettermconstr(Ydum, Xdum);
   end
   termconstr{sj} = tconstr;
end
tnames{length(tnames)+1,1} = 'Error';

% Package up term info into a structure
terminfo.termdum = termdum;
terminfo.termconstr = termconstr;
terminfo.levelcodes = levelcodes;
terminfo.tnames = tnames;
terminfo.dfterm0 = dfterm0;
terminfo.termvars = termvars;
terminfo.termlength = termlength;
end


% ----------------------------
function [tdum,vars,tn,df0,tconstr] = maketerm(tm,isnested,varinfo,j,sindex,termlist,...
                                             termdum,termconstr,tnames,dfterm0)
% Make term info such as dummy vars, name, constraints

% Loop over elements of the term
df0 = 1;
tdum = [];         % empty term so far
tconstr = 1;       % empty constraints so far
tn = '';           % blank name so far
vars = find(tm);   % list of variables making up terms
pwrs = tm(vars);   % and their exponents
tm = (tm>0);       % and a boolean mask for them
for varidx = 1:length(vars)
   % Process each variable participating in this term
   varnum = vars(varidx);          % variable name
   thispwr = pwrs(1);              % power of this variable
   tm(varnum) = 0;                 % term without this variable
   pwrs(1) = [];                   % powers without this variable
   df0 = df0 * varinfo.df(varnum); % d.f. so far

   % Combine its dummy variable with the part computed so far
   G = varinfo.vdum{varnum};       % dummy vars for this grouping var
   thisname = varinfo.varnames{varnum};
   if thispwr>1
       G = G .^ thispwr;
       thisname = sprintf('%s^%d',thisname,thispwr);
   end
   tdum = termcross(G,tdum);   % combine G into term dummy vars

   % Construct the term name and constraints matrix
   if nargout>1
      if (isempty(tn))
         tn = thisname;
         if ~isnested
            tconstr = varinfo.vconstr{varnum};
         end
      else
         tn = [tn '*' thisname];
         if ~isnested
            tconstr = [kron(varinfo.vconstr{varnum},eye(size(tconstr,2)));
                       kron(eye(length(varinfo.vconstr{varnum})),tconstr)];
         end
      end
   end

   if varidx<length(vars) && j>1
      % If the rest of this term is computed, take advantage of that
      prevterms = termlist(sindex(1:j-1),:);
      oldtm = find(all(prevterms(:,tm) == repmat(pwrs,size(prevterms,1),1),2));
      oldtm = oldtm(~any(prevterms(oldtm,~tm),2));
      if (length(oldtm) > 0)
         k = sindex(oldtm(1));
         tdum = termcross(termdum{k}, tdum);
         if nargout>1
            if ~isnested
               oconstr = termconstr{k};
               tconstr = [kron(tconstr,              eye(size(oconstr,2)));
                          kron(eye(size(tconstr,2)), oconstr)];
            end
            tn = [tn '*' tnames{k}];
            df0 = df0 * dfterm0(k);
         end
         break;
      end
   end
end
if (isempty(tn)), tn = 'Constant'; end
end

% -----------------------------
function varinfo = makedummyvars(continuous,allgrps,vnested,varinfo)
% Create main effect term values for each variable

ng = size(allgrps,2);
vdum = cell(ng,1);      % dummy terms for variables
vconstr = cell(ng,1);   % constraints
vlevels = ones(ng,1);   % levels corresponding to each dummy term
vmeans = zeros(ng,1);   % variable means, used for continuous factors
% varnames = varinfo.varnames;
allnames = varinfo.allnames;

for j=1:ng
   if continuous(j)
      % For continuous variables, use the variable itself w/o constraints
      vdum{j} = allgrps(:,j);
      vconstr{j} = ones(0,1);
      vmeans(j) = mean(allgrps(:,j));
   elseif isempty(vnested) || ~any(vnested(j,:))
      % Create dummy variable arrays for each grouping variable
      % using a sum-to-zero constraint
      vdum{j} = idummy(allgrps(:,j), 3);
      nlevgrp  = size(vdum{j},2);
      vconstr{j} = ones(1,nlevgrp);  % sum-to-zero constraint for categorical
      vlevels(j) = nlevgrp;
   else
      % Create dummy variable arrays for a nested variable
      nesternums = find(vnested(j,:));
      allvars = allgrps(:,[nesternums, j]);
      [ugrps,ignore2,grpnum] = unique(allvars,'rows');
      vdum{j} = idummy(grpnum,3);
      vconstr{j} = []; % to be computed later
      vlevels(j) = size(vdum{j},2);
      
      % Fix up names for nested variables
      levelnames = varinfo.allnames{j};
      namesj = cell(vlevels(j),1);
      nester1names = allnames{nesternums(1)};
      for rnum = 1:vlevels(j)
          nesterlist = nester1names{ugrps(rnum,1)};
          for k = 2:length(nesternums)
              nesterknames = allnames{nesternums(k)};
              nesterlist = sprintf('%s,%s',nesterlist,...
                                   nesterknames{ugrps(rnum,k)});
          end
          namesj{rnum} = sprintf('%s(%s)',levelnames{ugrps(rnum,end)},...
                                          nesterlist);
      end
      varinfo.allnames{j} = namesj;
   end
end

varinfo.vdum = vdum;
varinfo.vconstr = vconstr;
varinfo.vlevels = vlevels;
varinfo.vmeans = vmeans;
end

% -----------------------------
function d = idummy(x, method)
%DUMMY  Creates a matrix of dummy variables for a discrete variable
%   D=IDUMMY(X,METHOD) creates an array D of dummy variables for the
%   grouping variable I (integers 1,...,g), using the method specified:
%
%   method = 1:   0/-1/1 coding, full rank
%   method = 2:   0/1 coding, full rank
%   method = 3:   0/1 coding, overdetermined

%   Copyright 1993-2005 The MathWorks, Inc. 
%   $Revision: 1.4.2.2 $  $Date: 2005/11/18 14:28:47 $

if (nargin < 2)
   method = 1;
end

n = length(x);
g = max(x);
ncols = g - (method ~= 3);
d = repmat(0, n, ncols);

if (g > 1 || method==3)
   % Fill in -1 for the first level
   if (method == 1)
      d((x == 1),:) = -1;
   end
   
   % Fill in 1 in the appropriate column for other levels
   m3 = (method == 3);
   for j=(2-m3):g
      d((x == j),j-1+m3) = 1;
   end
end
end

% -----------------------------
function ab = termcross(a,b)
%TERMCROSS Multiply dummy variables for two terms to get interaction

%   Copyright 1993-2004 The MathWorks, Inc. 
%   $Revision: 1.3.2.1 $  $Date: 2004/01/24 09:36:40 $
if (isempty(a)), ab = b; return, end
if (isempty(b)), ab = a; return, end

na = size(a,2);
nb = size(b,2);
acols = repmat((1:na), 1, nb);
bcols = reshape(repmat((1:nb), na, 1), 1, na*nb);
ab = a(:,acols) .* b(:,bcols);

end

% ------------------------------------------------
function txt=linprodtxt(coeffs,names,prefix)
%LINPRODTXT Create a text representing a linear combination of some things

txt = '';
plustxt = '';

% If things will display as 0 or 1, make them exact
tol = sqrt(eps(class(coeffs))) * max(abs(coeffs));
coeffs(abs(coeffs)<tol) = 0;
coeffs(abs(coeffs-1)<tol) = 1;
coeffs(abs(coeffs+1)<tol) = -1;

% Add each component of the linear combination
preftxt = prefix;
for i=1:length(coeffs)
   bi = coeffs(i);
   if bi~=0
      % Get the sign
      if bi<0
         signtxt = '-';
      else
         signtxt = plustxt;
      end
      
      % Get the coefficient and multiplication symbol, unless it's one
      if abs(bi)==1
         coefftxt = '';
      else
         coefftxt = sprintf('%g*',abs(bi));
      end
      
      % If each row has its own prefix, get this one
      if size(prefix,1)>1
         preftxt = deblank(prefix(i,:));
      end
      
      txt = sprintf('%s%s%s%s(%s)',txt,signtxt,coefftxt,preftxt,names{i});
      plustxt = '+';
   end
end
end
% -----------------------------------------
function [varest,varci] = varcompest(ems,Rterms,msTerm,dfTerm,alpha)
%VARCOMPEST Estimate variance components from obs. & exp. mean squares

nterms = size(msTerm,1) - 1;          % not including error term

Fterms = find(~ismember(1:nterms, Rterms));
Fterms = Fterms(:);

% Get A and B so that expected means square is Av+Bq
% where v is the variances and q is quadratic forms in the fixed effects
B = ems(Rterms, Fterms);
A = ems(Rterms, Rterms);
const = msTerm(Rterms,:) - B*msTerm(Fterms,:);
Ainv = pinv(A);          % should be well-conditioned and easily invertible
varest = Ainv * const;

% Get confidence bounds for the variance components
varci = repmat(NaN,size(varest,1),2);
if all(dfTerm>0)
   L = bsxfun(@times, msTerm , dfTerm ./ chi2inv(1-alpha/2,dfTerm));
   U = bsxfun(@times, msTerm , dfTerm ./ chi2inv(alpha/2,dfTerm));
   AinvB = Ainv * B;
   for j=find(varest'>0)
      Rcoeffs = Ainv(j,:);
      Fcoeffs = AinvB(j,:);
      t = (Rcoeffs>0);
      s = (Fcoeffs>0);
      cL = 0;
      cU = 0;
      if any(t)
         cL = cL + dot(Rcoeffs(t),  L(Rterms(t)));
         cU = cU + dot(Rcoeffs(t),  U(Rterms(t)));
      end
      if any(~t)
         cL = cL + dot(Rcoeffs(~t), U(Rterms(~t)));
         cU = cU + dot(Rcoeffs(~t), L(Rterms(~t)));
      end
      if any(s)
         cL = cL + dot(Fcoeffs(s),  L(Fterms(s)));
         cU = cU + dot(Fcoeffs(s),  U(Fterms(s)));
      end
      if any(~s)
         cL = cL + dot(Fcoeffs(~s), U(Fterms(~s)));
         cU = cU + dot(Fcoeffs(~s), L(Fterms(~s)));
      end
      varci(j,1) = max(0,cL);
      varci(j,2) = max(0,cU);
   end   
end
end

% -----------------------------
function [dmat,cmat,termname] = makedesignmat(terminfo,n)

termlength = terminfo.termlength;
ncols = sum(termlength);

nconstr = sum(cellfun('size',terminfo.termconstr,1));
dmat = ones(n, ncols+1);      % to hold design matrix
cmat = zeros(nconstr,ncols);  % to hold constraints matrix
cbase = 0;                    % base from which to fill in cmat
termname = zeros(ncols,1);
termstart = cumsum([2; termlength(1:end-1)]);
termend = termstart + termlength - 1;
for j=1:length(termlength)
   clist = termstart(j):termend(j);
   dmat(:, clist) = terminfo.termdum{j};
   C = terminfo.termconstr{j};
   nC = size(C,1);
   cmat(cbase+1:cbase+nC,clist) = C;
   termname(clist) = j;
   cbase = cbase + nC;
end
end

% -----------------------------
function L = getL(X, curcols, incols, L2)
%GETL Get L hypothesis matrix for the current term after adjusting for
%     terms that remain "in" for this test, making the whole thing
%     orthogonal to the matrix L2 containing hypotheses for terms that
%     contain the current term.
%
%        Type   incols                  L2
%        ----   --------------------    ----------------
%        1      lower-numbered terms    empty
%        2      non-containing terms    empty
%        3      non-containing terms    containing terms

x1 = X(:,curcols);     % the current term
x0 = X(:,incols);      % terms that remain "in" when testing x1

% Find hypothesis for the current term.  First remove effects of terms it does
% not contain (these terms remain in the model when testing the current term).
[Q0,R,E] = qr(x0,0);
tol = 100 * max(size(x0)) * eps(R(1,1));
if min(size(R))>1
   p = sum(abs(diag(R)) > tol);
else
   p = 1;
end

Q0 = Q0(:,1:p);
adjx1 = x1 - Q0 * (Q0' * x1);

% Now find a full-rank basis for the adjusted effect.
[Q2,R,E] = qr(adjx1,0);
tol = 100 * max(size(R)) * eps(R(1,1));
if min(size(R))>1
   p = sum(abs(diag(R)) > tol);
else
   p = 1;
end
Q2 = Q2(:,1:p);
L = Q2' * X;

% Make L rows orthogonal to the L rows for containing terms.
if ~isempty(L2)
   L2tr = L2';
   [Q,R,E] = qr(L2tr,0);
   L = (L' - Q * (Q' * L'))';

   % Make L full rank
   [Q,R,E] = qr(L,0);
   tol = 100 * max(size(L)) * eps(R(1,1));
   if min(size(R))>1
      p = sum(abs(diag(R)) > tol);
   else
      p = 1;
   end
   if p<size(L,1)
      L = L(1:p,:);
   end
end
end



% -------------------
function [varnames,termlist,randomvar,sstype,continuous,vnested] = ...
        argcheck(ng,varnames,model,randomvar,sstype,alpha,continuous,vnested)
% variable names
if isempty(varnames)
   varnames = cellstr([repmat('X',ng,1) strjust(num2str((1:ng)'),'left')]);
else
   if (iscell(varnames)), varnames = varnames(:); end
   if (size(varnames,1) ~= ng)
      error('stats:anovan:InputSizeMismatch','VARNAMES must have %d rows', ng);
   elseif (~(ischar(varnames) || iscellstr(varnames)))
      error('stats:anovan:BadVarNames',...
            'VARNAMES must be a character matrix or cell array of strings.');
   elseif (ischar(varnames))
      varnames = cellstr(varnames);
   end
end

% sum-of-squares type
if isempty(sstype)
   sstype = 3;
elseif ~isscalar(sstype) || (sstype ~= 1 && sstype ~= 2 && sstype ~=  3 && lower(sstype) ~= 'h')
   error('stats:anovan:BadSumSquares','SSTYPE must be 1, 2, or 3.');
end

% continuous must be a logical vector indicating continuous variables
if islogical(continuous) && length(continuous)==ng && numel(continuous)==ng
   continuous = continuous(:);
elseif isnumeric(continuous) && all(ismember(continuous(:),1:ng))
   continuous = ismember(1:ng,continuous);
else
   error('stats:anovan:BadContinuous',...
         'CONTINUOUS must be a vector of indices between 1 and %d.', ng);
end

% nested must be a binary matrix
if ~isempty(vnested)
    if ~isequal(size(vnested),[ng ng]) || ~all(vnested(:)==0 | vnested(:)==1)
        error('stats:anovan:BadNested',...
         'The NESTED parameter must be a %d-by-%d matrix of 0''s and 1''s.',...
         ng,ng);
    end
    
    % find all indirect nesting
    vnested = nest_d2i(vnested);
end

% alpha (confidence level is 100*(1-alpha)%
if ~isnumeric(alpha) || numel(alpha)~=1 || alpha<=0 || alpha>=1
   error('stats:anovan:BadAlpha','ALPHA must be a scalar between 0 and 1.');
end

% Randomvar may be a logical vector indicating the random variables
if islogical(randomvar) && length(randomvar)==ng && numel(randomvar)==ng
   randomvar = randomvar(:);
elseif isnumeric(randomvar) && all(ismember(randomvar(:),1:ng))
   % or it may be an index vector -- convert to logical
   randomvar = ismember((1:ng)', randomvar(:));
else
   error('stats:anovan:BadRandomVar',...
         'RANDOMVAR must be a vector of indices between 1 and %d.', ng);
end
if any(randomvar(continuous))
   error('stats:anovan:BadRandomVar',...
        'You cannot specify a continuous predictor as having a random effect');
end

% model
if isempty(model)
   model = 1;
elseif ~isnumeric(model)
   if strcmp(model, 'linear') || strcmp(model, 'l')
      model = 1;
   elseif strcmp(model, 'interaction') || strcmp(model, 'i')
      model = 2;
   elseif strcmp(model, 'full') || strcmp(model, 'f')
      model = ng;
   else
      error('stats:anovan:BadModel','MODEL is invalid.');
   end
elseif min(size(model))>1 || ...
       (isequal(size(model),[1 ng]) && all(ismember(model,0:1)))
   % Matrix form of input
   if (any(model~=floor(model)|model<0))
      error('stats:anovan:BadModel',...
            'MODEL matrix entries must be non-negative integers.');
   elseif size(model,2)~=ng
      error('stats:anovan:InputSizeMismatch',...
            'MODEL matrix must have %d columns.',ng);
   end
else
   if (any(model~=floor(model)|model<=0))
      % Original vector form, no longer advertised
      error('stats:anovan:BadModel','MODEL is invalid.');
   end
   model = model(:);
end

% convert model to numeric matrix list of terms
if (length(model)==1) || (any(size(model)==1) && size(model,2)~=ng)
   % Convert from scalar or coded integer to matrix form
   termlist = makemodel(model, ng, vnested);
else
   termlist = model;
end

% look for invalid model terms
t = all(termlist==0,2);
termlist(t,:) = [];      % remove constant term -- it's implicit
if any(any(termlist(:,~continuous)>1))
    error('stats:anovan:BadModel',...
         'MODEL cannot have squared or higher terms for categorical factors.');
end
if ~isempty(vnested)
   % Must not have A(B) and B(A)
   if any(any(vnested & vnested'))
       error('stats:anovan:BadNesting',...
             'Circular nesting (such as A(B) and B(A)) is not allowed.');
   end
   
   % Must not have an interaction between B and A(B)
   [nestee,nester] = find(vnested);
   bad = any(termlist(:,nestee) & termlist(:,nester), 2);
   if any(bad)
      error('stats:anovan:BadModel',...
           ['MODEL must not have an interaction between a nested factor'...
            '\nand the factor that nests it.']);
   end
end
end


% ----------------------------
function fulltermlist = showtermnesting(termlist,vnested)

fulltermlist = termlist;

[nestee,nester] = find(vnested);
for j=1:length(nester)
    t = fulltermlist(:,nestee(j))>0;
    fulltermlist(t,nester(j)) = 1;
end
end

% -----------------------------
function tnested = gettermnesting(fulltermlist,vnested,continuous)
% Create matrix with (i,j) indicating if term i is nested in term j

% No nested variables implies no nested terms
if isempty(vnested)
    tnested = [];
    return
end

% Work with categorical and continuous terms separately
ctermlist = fulltermlist(:,continuous);
fulltermlist(:,continuous) = 0;

nterms = size(fulltermlist,1);
tnested = zeros(nterms);

% If A(B), then any term containing A may be nested in one containing B
[nestee,nester] = find(vnested);
for j=1:length(nester)
    hasnestee = (fulltermlist(:,nestee(j))>0);
    hasnester = (fulltermlist(:,nester(j))>0);
    tnested(hasnestee,hasnester) = 1;
end

% But terms are not nested within themselves
tnested(1:nterms+1:end) = 0;

% And there is no nesting relationship if the continuous vars don't match
for j=1:size(ctermlist,2)
    c = repmat(ctermlist(:,j),1,nterms);
    tnested = tnested & (c == c');
end

% And the full representation for a nesting term must be a subset of the
% full representation for a nested term
[nestee,nester] = find(tnested);
for j=1:length(nester)
    if any(fulltermlist(nestee(j),:) < fulltermlist(nester(j),:))
        tnested(nestee(j),nester(j)) = 0;
    end
end

end

% -----------------------------------
function [y,allgrps,nanrow,varinfo] = removenans(y,group,continuous)

% Find NaNs among response and group arrays
[n,NN] = size(y);
nanrow = find(isnan(y));  % !!! SO. change that behaviour more optimal
ng = length(group);
for j=1:ng
   gj = group{j};
   if (size(gj,1) == 1), gj = gj(:); end
   if (size(gj,1) ~= n)
      error('stats:anovan:InputSizeMismatch',...
            'Group variable %d must have %d elements.',j,n);
   end
   if (ischar(gj)), gj = cellstr(gj); end
   if ~isvector(gj)
       error('stats:anovan:BadGroup',...
             'Grouping variable must be a vector or a character array.');
   end

   group{j} = gj;
   if (isnumeric(gj) & any(isnan(gj)))
      nanrow = (nanrow | repmat(find(isnan(gj)),1,NN));
   elseif isa(gj,'categorical')
      nanrow = (nanrow | repmat(isundefined(gj),1,NN));
   elseif (~isnumeric(gj))
      nanrow = (nanrow | repmat(strcmp(gj,''),1,NN));
   end
end

% Remove rows with NaN anywhere
y(nanrow,:) = [];
[n,NN] = size(y);
ng = length(group);
dfvar = zeros(ng,1);
allgrps = zeros(n, ng);
allnames = cell(ng,1);

% Get arrays describing the groups
for j=1:ng
   gj = group{j};
   gj(nanrow,:) = [];
   group{j} = gj;
   if continuous(j)
      dfvar(j) = 1;
      allgrps(:,j) = gj;
   else
      [gij,gnj] = grp2idx(gj);
      nlevels = size(gnj,1);
      dfvar(j) = nlevels - 1;
      allgrps(:,j) = gij;
      allnames{j} = gnj;
   end
end

% The df and allnames information does not yet reflect nesting.  We will
% fix up allnames for nested variables later, and not use df for them.

varinfo.df = dfvar;
varinfo.allnames = allnames;
end


% --------------------------
function v = makemodel(p,ng,vnested)
%MAKEMODEL Helper function to make model matrix
%    P = max model term size, NG = number of grouping variables
% or
%    P = vector of term codes
% VNESTED is a matrix specifying direct or indirect nesting.

% We want a matrix with one row per term.  Each row has a 1 for
% the variables participating in the term.  There will be rows
% summing up to 1, 2, ..., p.

if numel(p)==1
   % Create model matrix from a scalar max order value
   vgen = 1:ng;
   v = eye(ng);                      % linear terms
   for j=2:min(p,ng)
      c = nchoosek(vgen,j);          % generate column #'s with 1's
      nrows = size(c,1);             % generate row #'s
      r = repmat((1:nrows)',1,j);    %    and make it conform with c
      m = zeros(nrows,ng);           % create a matrix to get new rows
      m(r(:)+nrows*(c(:)-1)) = 1;    % fill in 1's
      v = [v; m];                    % append rows
   end
else
   % Create model matrix from terms encoded as bit patterns
   nterms = length(p);
   v = zeros(nterms,ng);
   for j=1:nterms
      tm = p(j);
      while(tm)
         % Get last-numbered effect remaining
         lne = 1 + floor(log2(tm));
         tm = bitset(tm, lne, 0);
         v(j,lne) = 1;
      end
   end
end

% Remove terms forbidden by nesting
[nestee,nester] = find(vnested);
bad = any(v(:,nestee) & v(:,nester), 2);
v(bad,:) = [];

end

% --------------------------
function [ssx,dfx,dmat2,y2,stats] = dofit(dmat,y,cmat,sst,mu)
%DOFIT Do constrained least squares fit and reduce data for subsequent fits

% Find the null space of the constraints matrix
[Qc,Rc,~] = qr(cmat');
pc = Rrank(Rc);
Qc0 = Qc(:,pc+1:end);

% Do qr decomposition on design matrix projected to null space
Dproj = dmat*Qc0;
[Qd,Rd,Ed] = qr(Dproj,0);
Dproj = []; % no longer needed but potentially big
dfx = Rrank(Rd);
Qd = Qd(:,1:dfx);
Rd = Rd(1:dfx,1:dfx);

% Fit y to design matrix in null space
y2 = Qd' * y;            % rotate y into that space
zq = Rd \ y2;            % fit rotated y to projected design matrix

z = zeros(length(Ed),size(y,2)); % coefficient vector extend to full size ...
z(Ed(1:dfx),:) = zq;       % ... and in correct order for null space

b = Qc0 * z;             % coefficients back in full space
 
ssx = sum(y2.^2);        % sum of squares explained by fit

% Return reduced design matrix if requested
if nargout>2
   dmat2 = Qd' * dmat;
end

% Prepare for multiple comparisons if requested
if (nargout >= 3)
   stats.source = 'anovan';
   
   % Get residuals
   % yhat = Dproj * z;
   yhat = Qd * Rd * zq;
   stats.resid = y - yhat;
   
   % Calculate coefficients, then adjust for previously removed mean
   t = b;
   if ~isempty(t)
      t(1,:) = t(1,:) + mu;            % undo mean adjustment before storing
   end
   stats.coeffs = t;

   RR = zeros(dfx,length(Ed));
   RR(:,Ed(1:dfx)) = Rd;
   stats.Rtr = RR';
   [~,rr,ee] = qr(dmat,0);
   pp = Rrank(rr);
   rowbasis = zeros(pp,size(rr,2));
   rowbasis(:,ee) = rr(1:pp,:);
   stats.rowbasis = rowbasis;
   stats.dfe = length(y) - dfx;
   stats.mse = (sst - ssx) / max(1, stats.dfe);
   stats.nullproject = Qc0;
end
end

%---------------------
function [codes,names] = convertstats(stats,terminfo,varinfo,continuous)

allnames = varinfo.allnames;
varnames = varinfo.varnames;
ng = length(varnames);
fullcoef = stats.coeffs;
codes = zeros(size(fullcoef,1), ng);
names = cell(size(fullcoef,1),1);
base = 1;
names{1} = 'Constant';

% compute names
for j=1:length(terminfo.levelcodes)
   M = terminfo.levelcodes{j};
   varlist = terminfo.termvars{j};
   v1 = varlist(1);
   vals1 = allnames{v1};
   nrows = size(M,1);
   codes(base+(1:nrows),varlist) = M;
   for k=1:nrows
      if continuous(v1)
          nm = varnames{v1};
      else
          nm = sprintf('%s=%s',varnames{v1},vals1{M(k,1)});
      end
      for i=2:length(varlist)
         v2 = varlist(i);
         vals2 = allnames{v2};
         if continuous(v2)
             nm2 = varnames{v2};
         else
             nm2 = sprintf('%s=%s',varnames{v2},vals2{M(k,i)});
         end
         nm = sprintf('%s * %s',nm,nm2);
      end
      names{base+k} = nm;
   end
   base = base+nrows;
end
end

