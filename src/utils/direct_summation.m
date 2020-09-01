function F = direct_summation(f,x,M)
% F = direct_summation(f,x,M)
%--------------------------------------------------------------------------
% PURPOSE
%  Computes the discrete Fourier transform of the uniform or non-uniform
%  positioned data f by direct summation.
%
% INPUT: f = [f_1; f_2; ... ;f_N]   1-dimensional input data
%        x = [x_1; x_2; ... ;x_N]   uniform or non-uniform positions 
%                                   in [0,2*pi]
%        M                          number of frequencies k s.t.
%                                   -M/2 <= k < M/2
%
% OUTPUT: F = [F_1; F_2; ... ;F_M]  Fourier coefficients
%--------------------------------------------------------------------------

%-Dimensions and parameters------------------------------------------------
% Number of data points
N = length(f);

% Set uniform positions if x is not provided
if nargin < 2
    x = (0:N-1)'/N;
    x = 2*pi*x;
end
if nargin < 3
    M = N;
end

%-Compute sum--------------------------------------------------------------
F = zeros(M,1);
idx = 1;
for k = -M/2:M/2-1
    sum = 0;
    for j = 0:N-1
        sum = sum + f(j+1)*exp(-1i*k*x(j+1));
    end
    F(idx) = sum;
    idx = idx + 1;
end

%--------------------------end---------------------------------------------
