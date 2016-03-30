# helmholtz_solver
Simulates the (time-independent) wave equation in 2D environments using a finite-difference approach.

For example, starting with the floorplan of my apartment:

![Floorplan](floorplan.png?raw=true "Floorplan")

The program can simulate how wifi signals might propagate through the walls from my router:

![Wifi Propagation](img/example_output.png?raw=true "Wifi Propagation")

## How it works
In short, we want to approximate the solution to the inhomogeneous Helmholtz equation:

![Inhomogeneous Helmholtz Equation](img/equation0.png?raw=true "Inhomogeneous Helmholtz Equation")

First, we divide the region over which we would like to solve the equation into a grid.

Next, we represent each cell of that grid with an element in a column vector. I.e., a[n] might represent grid cell (n mod W, floor(n / W)), where W is the width of the grid.

We evaluate the source function f over our grid, and store the result in the column vector F.

Next, we evaluate the [wavenumber](https://en.wikipedia.org/wiki/Wavenumber) k for each grid cell, which is a function of the index of refraction of each cell.
In open space, k = k0 = 2*pi/wavelength. Otherwise, k = k0*n, where n is the index of refraction of the material present. We store the k values in a diagonal matrix K, where the diagonal elements are simply the k values for each cell.

Finally, we create a matrix L which is a discrete approximation of the 2D [laplacian operator](https://en.wikipedia.org/wiki/Laplace_operator#Two_dimensions). I.e., when multiplied on the left, this matrix takes the second spatial derivative of a function defined on our grid using the finite difference approximation.

With these tools in place, solving the problem is as simple as solving the matrix equation M*A=F, where M = (L+K^2) and A is the column vector containing our result.
