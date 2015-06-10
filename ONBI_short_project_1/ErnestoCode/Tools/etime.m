function [seconds_,str,ymwdhms] = etime( t1 , t0 )
%ETIME  Elapsed time.
%   ETIME(T1,T0) returns the time in seconds that has elapsed between
%   vectors T1 and T0.  The two vectors must be six elements long, in
%   the format returned by CLOCK:
%
%       T = [Year Month Day Hour Minute Second]
%
%   Time differences over many orders of magnitude are computed accurately.
%   The result can be thousands of seconds if T1 and T0 differ in their
%   first five components, or small fractions of seconds if the first five
%   components are equal.
%
%     t0 = clock;
%     operation
%     etime(clock,t0)
%
%   See also TIC, TOC, CLOCK, CPUTIME, DATENUM.

%   Copyright 1984-2002 The MathWorks, Inc. 
%   $Revision: 5.9.4.1 $  $Date: 2002/09/30 12:01:19 $

% Compute time difference accurately to preserve fractions of seconds.

  if nargin < 1
    error('one or two arguments expected.');
  end
  if nargin < 2
    t0 = t1;
    t1 = clock;
  end
  
  if iscell( t1 )
    for j = 1:numel(t1)
      seconds_(j,:) = etime( t0 , t1{j} );
    end
    return;
  end
  if iscell( t0 )
    for i = 1:numel(t0)
      seconds_(:,i) = etime( t0{i} , t1 );
    end
    return;
  end
  
  
  if      isscalar( t0 )
    t0 = datevec( t0 );
  elseif  ischar( t0 )
    t0 = datevec( t0 );
  elseif numel( t0 ) ~= 6
    error('invalid t0');
  end
    
  if      isscalar( t1 )
    t1 = datevec( t1 );
  elseif  ischar( t1 )
    t1 = datevec( t1 );
  elseif numel( t1 ) ~= 6
    error('invalid t1');
  end

  seconds = 86400*(datenummx(t1(:,1:3)) - datenummx(t0(:,1:3))) + ...
      (t1(:,4:6) - t0(:,4:6))*[3600; 60; 1];
  
  if nargout >  0, seconds_ = seconds; end
  if nargout == 1, return; end
  
  secs = abs( seconds );
  
  ymwdhms(1) = floor( secs / ( 365*24*60*60 ) );  secs = secs - ymwdhms(1)*365*24*60*60;
  ymwdhms(2) = floor( secs / (  30*24*60*60 ) );  secs = secs - ymwdhms(2)* 30*24*60*60;
  ymwdhms(3) = floor( secs / (   7*24*60*60 ) );  secs = secs - ymwdhms(3)*  7*24*60*60;
  ymwdhms(4) = floor( secs / (     24*60*60 ) );  secs = secs - ymwdhms(4)*    24*60*60;
  ymwdhms(5) = floor( secs / (        60*60 ) );  secs = secs - ymwdhms(5)*       60*60;
  ymwdhms(6) = floor( secs / (           60 ) );  secs = secs - ymwdhms(6)*          60;
  ymwdhms(7) = secs;
  

  str = [];
  
  if      ymwdhms(1) > 1
    str = sprintf('%s%d years',str,ymwdhms(1));
  elseif  ymwdhms(1) == 1
    str = sprintf('%s1 year',str);
  end
  
  if      ymwdhms(2) > 1
    if ~isempty( str ), str = [ str ', ']; end
    str = sprintf('%s%d months',str,ymwdhms(2));
  elseif  ymwdhms(2) == 1
    if ~isempty( str ), str = [ str ', ']; end
    str = sprintf('%s1 month',str);
  end
  
  if      ymwdhms(3) > 1
    if ~isempty( str ), str = [ str ', ']; end
    str = sprintf('%s%d weeks',str,ymwdhms(3));
  elseif  ymwdhms(3) == 1
    if ~isempty( str ), str = [ str ', ']; end
    str = sprintf('%s1 week',str);
  end
  
  if      ymwdhms(4) > 1
    if ~isempty( str ), str = [ str ', ']; end
    str = sprintf('%s%d days',str,ymwdhms(4));
  elseif  ymwdhms(4) == 1
    if ~isempty( str ), str = [ str ', ']; end
    str = sprintf('%s1 day',str);
  end
  
  if      ymwdhms(5) > 1
    if ~isempty( str ), str = [ str ', ']; end
    str = sprintf('%s%d hours',str,ymwdhms(5));
  elseif  ymwdhms(5) == 1
    if ~isempty( str ), str = [ str ', ']; end
    str = sprintf('%s1 hour',str);
  end
  
  if      ymwdhms(6) > 1
    if ~isempty( str ), str = [ str ', ']; end
    str = sprintf('%s%d minutes',str,ymwdhms(6));
  elseif  ymwdhms(6) == 1
    if ~isempty( str ), str = [ str ', ']; end
    str = sprintf('%s1 minute',str);
  end

  if      ymwdhms(7) ~= 1
    if ~isempty( str ), str = [ str ' and ']; end
    str = sprintf('%s%f seconds',str,ymwdhms(7));
    %str = sprintf('%s%s seconds',str,uneval(ymwdhms(7)));
  elseif  ymwdhms(7) == 1
    if ~isempty( str ), str = [ str ' and ']; end
    str = sprintf('%s1 second',str);
  end

  if seconds < 0
    str = [ 'minus '  str ];
  end
  
  if nargout == 0
    disp( str );
  end
  
end
