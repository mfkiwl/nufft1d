function F = nufft1d(f,x,M,R,M_sp,tau)
% F = nufft1d(f,x,M,R,M_sp,tau)
%--------------------------------------------------------------------------
% PURPOSE
%  Approximates the 1-dimensional non-uniform discrete Fourier transform 
%  using the non-uniform fast Fourier transform (NUFFT) algorithm. See 
%  Greengard and Lee (2004) for more details.
%
% INPUT: f = [f_1; f_2; ... ;f_N]   1-dimensional input data
%        x = [x_1; x_2; ... ;x_N]   non-uniform positions in [0,2*pi]
%        M                          number of frequencies k s.t. 
%                                   -M/2 <= k < M/2
%        R                          oversampling factor
%        M_sp                       number of neighbours to spread data to
%        tau                        spreading factor of Gaussian kernel
%
% OUTPUT: F = [F_1; F_2; ... ;F_M]  Fourier coefficients
%--------------------------------------------------------------------------

%-Dimensions and parameters------------------------------------------------
% Number of data points
N = length(x);

% Set default values
if nargin < 3
    M = N;
end
if nargin < 4
    R = 2;
end
if nargin < 5
    M_sp = 12;
end
if nargin < 6
    tau = (1/M^2)*(pi*M_sp)/(R*(R-0.5));
end

% Number of points in oversampled grid
M_r = R*M;

% Frequencies
k = (-M/2:M/2-1)';

%-Construct oversampled uniform grid---------------------------------------
% Fine upsampled grid in [0,2*pi]
xi = linspace(0,2*pi,M_r)';

% Step size in upsampled grid
h = (2*pi)/M_r;

% Find indices in xi that are the closest to x_j
distances = bsxfun(@(x,y) abs(x-y), x(:), xi');
[~, idx] = min(distances,[],2);

%-Convolve f with g_tau----------------------------------------------------
f_tau = zeros(M_r,1);
for j = 1:N
    for l = -M_sp:M_sp
        f_tau(mod(idx(j)+l,M_r)+1) = f_tau(mod(idx(j)+l,M_r)+1) + ...
            f(j)*exp(-(x(j)-h*idx(j)-h*l)^2/(4*tau));
    end
end

%-Standard FFT-------------------------------------------------------------
% FFT and shift the frequencies -M_r/2 <= k < M_r/2
F_tau = fftshift(fft(f_tau,M_r));

%-Compute Fourier coefficients---------------------------------------------
% Correct for spreading in Fourier domain for -M/2 <= k < M/2
F = sqrt(pi/tau) * exp(k.^2*tau) .* F_tau(M_r/2-M/2+1:M_r/2+M/2)/M_r;

%--------------------------end---------------------------------------------
