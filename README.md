# nufft1d
An implementation of the non-uniform fast Fourier transform (NUFFT) of type 1 in 1D based on the work done by Greengard and Lee (2004) [1]. It achieves an accuracy comparable with the one of the FFT (around `1e-14`) when compared to performing the direct summation. Worth mentioning is however that this is a naive implementation of the algorithm in [1], and does thus not fully utilize all of the computational benefits of the algorithm.

For more information about the NUFFT as well as the mentioned algorithm and its performance, see the included [report](D_Krantz_Non_Uniform_fast_Fourier_Transform_Report.pdf).

## Usage
[Here](example.m) you can see an example of how the main function `nufft1d` can be used (example from the [report](D_Krantz_Non_Uniform_fast_Fourier_Transform_Report.pdf)).

```matlab
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
```

## References
[1] Leslie Greengard and June-Yub Lee. "Accelerating the Nonuniform Fast Fourier Transform". In: *SIAM Review* 46.3 (2004), pp. 443-454. DOI: [10.1137/S003614450343200X](https://doi.org/10.1137/S003614450343200X). eprint: [https://doi.org/10.1137/S003614450343200X](https://doi.org/10.1137/S003614450343200X). URL: [https://doi.org/10.1137/S003614450343200X](https://doi.org/10.1137/S003614450343200X).
