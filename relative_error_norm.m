function e = relative_error_norm(x_hat,x)
% e = relative_error_norm(x_hat,x)
%--------------------------------------------------------------------------
% PURPOSE
%  Computes the L2 relative error norm between x_hat and x.
%
% INPUT: x_hat                      calculated solution
%        x                          exact solution 
%
% OUTPUT: e                         the relative L2 norm
%--------------------------------------------------------------------------
e = norm(x_hat-x,2) / norm(x,2);
%--------------------------end---------------------------------------------
