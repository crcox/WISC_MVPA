function v = value( x, data )
if nargin == 1,
    global cvx___
    data = cvx___.x;
end
nx = size( data, 1 );
nb = size( x.basis_, 1 );
if nx < nb,
    data( end + 1 : nb, : ) = NaN;
elseif nx > nb,
    data( nb + 1 : end, : ) = [];
end
v = cvx_reshape( data.' * x.basis_, x.size_ );
if any( x.size_ == 1 ), v = full( v ); end

% Copyright 2008 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
