# Distributed Navier Stokes solver via grid decomposition - extension submission

In the `./extension/` directory, I present the extended core solver. It has been extended to include a simple square obstacle in the centre of the domain with a variable size configurable via the `half_width` parameter in `main()`. Because of the modular design of the core logic, most of the augmentations to make this possible were fairly trivial. They were largely composed of making the core solver routine aware of the obstacle and imposing boundary conditions on its walls.

The majority of the core functionality is explained in the core [README.md](../core/README.md). Please consult that before this.

## Usage

To run the extension solver, from the root directory run

```bash
mpicxx -std=c++17 -O2 -o solver-extension extension/solver.cpp
mpirun -np nprocs solver-extension --Nx Nx --Ny Nx --Lx Lx --Ly Ly
```

Note the program accepts command line arguments for setting the domain shape and size, but the logic is fairly brittle. This is largely for benchmarking purposes.

Again from the root, to postprocess results (combining them into an animation panel), run

```bash
python3 postprocessing.py --extension
```

This will output a `flow_fields.gif` to `./out`. Note the postprocessing script is somewhat flexible in terms of data location so long as you don't alter the structure/ naming convention of the output data. You can specify a particular file path with

```bash
python3 postprocessing.py --root "./path/to/data"
```

## Logic breakdown 

The only new functionality is the `extension` namespace which encapsulates a `build_square_obstacle` function that, given a half width and the global `Nx` and `Ny`, creates a 'mask' representing the obstacle. For each subdomain, the mask is 1 where the obstacle exists, and 0 elsewhere.

The core routine is also largely the same, but includes the following additions in the setup phase before entering the main solve:

1. The obstacle is built across all processes. Note it is built according to the _global_ domain shape, but divided between each subdomain on the corresponding process.
2. The obstacle halos are eaxchanged once. We need not do this repeatedly, since the obstacle data is fixed, but each proc does need to be aware of wether or not its neighbour contains any portion of the obstacle.

#### Solving

The structure of the main solver loop is similar, but in this extended version we have some other consdierations:

- Within the obstacle bounds, the velocity should be zero.
- At the walls of the obstacle, no slip conditions should be imposed.

This leads to changes in the way stencils should be computed. I handled this using conditonals in each phase of the main solve. For instance, the main loop of the velocity projection became:

```cpp
for (int i = i_start; i <= i_end; i++)
    for (int j = j_start; j <= j_end; j++)
    {
        // zero the velocities inside the obstacle
        if (obstacle(i, j) > 0.5) {
            u(i, j) = 0.0;
            v(i, j) = 0.0;
            continue;
        }

        // check neighbours for obstacle
        // enforce no slip at the boundaries of the obstacle
        const double P_right = (obstacle(i + 1, j) > 0.5) ? P(i, j) : P(i+1, j);
        const double P_left = (obstacle(i - 1, j) > 0.5) ? P(i, j) : P(i-1, j);
        const double P_up = (obstacle(i, j + 1) > 0.5) ? P(i, j) : P(i, j+1);
        const double P_down = (obstacle(i, j - 1) > 0.5) ? P(i, j) : P(i, j-1);

        // use the above in the stencil:
        u(i, j) = u(i, j) - dt * (1. / rho) * (P_right - P_left) / (2. * dx);
        v(i, j) = v(i, j) - dt * (1. / rho) * (P_up - P_down) / (2. * dy);

        double vel = sqrt(u(i, j) * u(i, j) + v(i, j) * v(i, j));

        vmax = max(vmax, vel);
}
```

Here, at each step of the solve I use conditional statements to check if neighbouring cells are within the obstacle domain (indicated by `obstacle() > 0.5`). If the conditonal is true, then in the current cell we should enforce no slip boundary conditions by using the value of the current cell in place of that in the neighbouring cell. These conditionals are evaluated with ternary operators for compactness, and the stencil uses the dynamically set values instead of being computed on the fly by hardcoded neighbour accesses.