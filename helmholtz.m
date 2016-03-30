clear all;

% Solves the inhomogeneous Helmholtz Equation
% You can use this to see how waves would propagate in a 2D environment
% Or to investigate diffraction, refraction, .etc
% This example program solves for the propagation of wifi signals in a 2D environment -- my apartment!

% Floorplan image location
% Floorplan must be a B&W image. White areas are open space.
floorplan = 'floorplan.png';

% Scale of image in meters/pixel
scale = 0.8/42;

% Size of grid squares for calculation.
h = 0.02;

% Wavelength of source in meters.
lambda = 0.2;

% Position of source in metres.
% Top left of the image is (0,0)
x0 = 1;
y0 = 4;

% Index of refraction of walls.
% Imaginary part accounts for loss.
ior = 1.4+0.05i;

% For profiling / timing
tstart = tic;

% Read in floorplan
fprintf('Reading floorplan . . .\n');
I = imread(floorplan);

% Dimensions of image in meters
D = [0.0 scale*size(I,2) 0.0 scale*size(I,1)];

% N,M are the number of grid squares in the x and y dimensions
N = floor(((D(2)-D(1))/h)+1);
M = floor(((D(4)-D(3))/h)+1);

% Source function
f = @(x,y) exp(-1000*((x-x0).^2+(y-y0).^2));

% Wave number.
% In free space, k=k0=2*pi/lambda
% Otherwise, k=k0*n, where n is the (complex) index of refraction of
% the medium.
k = @(x,y) 2*pi/lambda*((ior-1).*(I(min(floor(y/scale)+1, size(I,1)), ...
    min(floor(x/scale)+1, size(I,2)))==0)+1);

% itp: Convert x,y into an index n in an array
itp = @(func, n) func(D(1)+(mod(n-1,N)+1)*h, D(3)+fix((n-1)/N)*h);

% Versions of F and K that operate on arrays
ff = @(n) itp(f,n);
kk = @(n) itp(k,n);

% Build the (sparse) laplacian matrix for our 2D grid
fprintf('Grid size: %d x %d\n', N, M);
fprintf('Building laplacian matrix . . .\n');
tic;
Dg = repmat([ones(1,N-1) 0], 1,M);
L = spdiags([Dg' ones(N*M,1)], [-1 -N], N*M, N*M) - 2*speye(N*M);
L = (L+L')./h^2;    % Account for spacial dimension / grid size
telapsed = toc;
fprintf('\tElapsed time: %.2f s\n', telapsed);

% Build the wavevector matrix.
% This is just a diagonal matrix with kk(n) as the diagonal elements.
fprintf('Building wavevector matrix . . .\n');
tic;
K = spdiags(arrayfun(kk, 1:N*M)', 0, N*M, N*M);
telapsed = toc;
fprintf('\tElapsed time: %.2f s\n', telapsed);

% Evaluate the function ff(n) over each grid square.
fprintf('Preparing forcing vector . . .\n');
tic;
F = sparse(ff(1:N*M)');
telapsed = toc;
fprintf('\tElapsed time: %.2f s\n', telapsed);

% Solve the linear system! (The hard part)
% Note that we have (L+K^2)*A=F
% So that A = (L+K^2)\F
fprintf('Solving linear system . . .\n');
tic;
% Compute the answer and reshape it into a matrix
A = reshape((L + K.^2)\F,N,M);
telapsed = toc;
fprintf('\tElapsed time: %.2f s\n', telapsed);

% Plot the output!
fprintf('Plotting results . . .\n');
figure;
colormap('jet');
imagesc(abs(A'));
colorbar;
axis equal;
telapsed = toc(tstart);
fprintf('Done. (Elaspsed time: %.2f s)\n', telapsed);