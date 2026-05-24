# Distributed Navier Stokes solver via grid decomposition - core submission

In this directory, I present the core of the submission: a distributed Navier Stokes solver for a pressure driven transient simulation which employs gird decomposition. In this document, I break down the important implementation details of this part of the submission. Notes on the directory structure and more general comments can be found in the main ![README.md](./../README.md).

## Usage

To run the core solver, from the root directory run

```bash
mpicxx -std=c++17 -O2 -o solver-core core/solver.cpp
mpirun -np nprocs ./solver-core --Nx Nx --Ny Nx --Lx Lx --Ly Ly
```
Note the program accepts command line arguments for setting the domain shape and size, but the logic is fairly brittle. This is largely for benchmarking purposes.

Again from the root, to postprocess results (combining them into an animation panel), run

```bash
python3 postprocessing.py --core
```

This will output a `flow_fields.gif` to `./out`. Note the postprocessing script is somewhat flexible in terms of data location so long as you don't alter the structure/ naming convention of the output data. You can specify a particular file path with

```bash
python3 postprocessing.py --root "./path/to/data"
```

## Logic breakdown 

In the `./core/` directory, I present the core parallelised Navier Stokes solver. There are a few core components to explain before detailing the full parallel routine. They are summarised in the table below.

|Component|Type|Purpose|Functionality|
|---------|----|-------|-------------|
|`solver.h`|Header file|Encapsulate custom data handling and communications logic|See below.|
|`solver.cpp`|Source file|House the main solver workflow.|Also handles basic arg parsing and some file writing logic.|
|`Domain`|`struct`|Main data container for all of the fields.  Stores data in flattened vector for contiguous storage, exposed only by public getters/ setters.|Stores full data and metadata for an aribtrary grid. Handles grid decomposition into equal aspect ratio sub regions. Handles basic operations: swaps, getting, setting.|
|`SubDomain`|`struct`|Dataclass defining metadata used to instntiate `Domain` objects following grid decomposition, including neighbouring ranks and BC application flags.|`Domain::decompose` returns a vector of `SubDomain` objects which contain information on how to instantiate full `Domain` objects, one for each rank.|
|`MPI_SubDomain`|`MPI_Datatype`|Custom type to communicate decomposition metadata.|Wraps the `SubDomain` type to enable the communication of the mixed/discontiguous types.|
|`MPI_Domainrow`|`MPI_Datatype`|Custom type to communicate the (geometrically) vertical halo buffers of units in the decomposed grid.|`Domain` objects store data in row-major order. Thus columns are contiguous in memory, and so this type is registered as `MPI_Type_contiguous`.|
|`MPI_Domaincol`|`MPI_Datatype`|Custom type to communicate the (geometrically) horizontal halo buffers of units in the decomposed grid.|`Domain` objects store data in row-major order. Thus rows are discontiguous in memory, and so this type is registered as `MPI_Type_vector`.|
|`exchange_halos`|function|Reusable function to exchange halo buffers for the fields of units in the decomposed grid.|Uses non-blocking sends and receives to transfer the above types to/from buffers/interior cells.|
|`CustomTypes`|`namespace`|To encapsulate all custom MPI type logic.|Defines and builds custom MPI types.|
|`utils`|`namespace`|To encapsulate utility functions|Core logic: building local `Domain` objects, exchanging halos. QOL functionality: metadata console logging, IO file checking, metadata writing.|

## The decomposition algorithm

`Domain::decompose` implements a simple grid decomposition that
splits the global `Nx` x `Ny` grid into `p` subdomains with the goal
of producing sub blocks with aspect ratios as close to unity as possible (i.e., as close to square as possible).

Key steps
1. **Choose proc grid (px, py):** iterate over factor pairs of `p` (px rangesfrom 1..p where `p % px == 0`) and compute the block aspect for each pair using

    - Block aspect = (Nx * py) / (Ny * px)

    The algorithm selects the factor pair (px, py) that minimises (Block aspect - 1)^2, i.e., closest to a square block.

2. **Base sizes and remainders:** once px and py are chosen, compute the base sizes and remainders:

    - nx_base = Nx / px, nx_rem = Nx % px
    - ny_base = Ny / py, ny_rem = Ny % py

    The decomposition distributes the remainder cells by giving the first `nx_rem` columns (grid_i = 0..nx_rem-1) one extra cell in x, and the first `ny_rem` rows one extra cell in y. This yields near even load balancing.

3. **Per-rank assignment:** for each rank `r` (0..p-1) compute grid coordinates grid_i = r % px and grid_j = r / px. Then

    - nx_local = nx_base + (grid_i < nx_rem ? 1 : 0)
    - ny_local = ny_base + (grid_j < ny_rem ? 1 : 0)
    - x0_global = grid_i * nx_base + min(grid_i, nx_rem)
    - y0_global = grid_j * ny_base + min(grid_j, ny_rem)

    A `SubDomain` is created with these local sizes and global offsets. The `SubDomain` constructor also computes neighbour ranks (`left_rank`, `right_rank`, `up_rank`, `down_rank`) using `MPI_PROC_NULL` for out-of-bound neighbours (to avoid nested conditionals in later communications), and sets the boundary flags `has_bc_left`, `has_bc_right`, `has_bc_up`, `has_bc_down`.

`Domain::decompose(p)` returns a `std::vector<SubDomain>` with length `p`. This call is performed on rank 0 and the vector is scattered to each process using the custom `MPI_SubDomain` type (see `CustomTypes::build_subdomain_type`). Each process uses its received `SubDomain` to instantiate local `Domain` objects (with halo/buffer), then builds local halo MPI types (`MPI_Domainrow`, `MPI_Domaincol`) and performs halo exchanges as needed.

The algorithm is deterministic and fast (iterates divisors of `p`, O(p) work for small p), and tends to produce compact blocks which reduces the surface-to-volume ratio, improving communication cost. Remainder cells are assigned to the earliest grid columns/rows; this simple scheme preserves contiguity of global indices and keeps the code straightforward.

## Full routine

The routine, in summary, is as follows:

1. MPI is initialised and the number of processes `p` is registered along with the PID `id` of this process.
2. The `MPI_SubDomain` type is built.
3. `setup()` computes global parameters `dx`, `dy`, `dt`.
4. Global simulation metadata is written to `./out/metadata.dat`. If `./out` does not exist, it is created. 
5. Proc 0 sets up a single, full sized `Domain` object `P`, without halo buffers. 
6. `P.decompose(p)` is called. This decomposes the grid into `p` sub units, returning a vector of `SubDomain` metadata objects, one for each proc.
7. On all processes, local `SubDomain` receive buffers are allocated, and proc 0 scatters the decomposition metadata into these buffers using the `MPI_SubDomain` type.
8. All fields (e.g., `u`, `u_old`, `v`) are identical in shape. However, `P` contains the inlet condition, while the other fields are intialised to zero. Thus, the `SubDomain` metadata is used to create a `template_field` that is instantiated once and copied onto all fields other than `P`, which is built seperately, enforcing the inlet condition.
9. The halo exchange types are built. This is done on all processes, since subdomain shape (and hence halo size) may very from process to process.
10. Output files are prepped. Each proc writes its own metadata to its history file. See the directory structure diagram in the main [README.md](../README.md) for notes on how files are stored.
11. The main solver loop is called.

#### Solving

The main solver loop and subroutines can be summarised by the "psuedocode" below:

```plaintxt
while (t < t_final); do
    if vel_max > 0; do
        - recompute dt based on CFL: dt = min(courant * min(dx, dy) / vel_max, dt_min)
    else
        - dt = dt_min

    - reduce dt in place across all processes using MPI_Allreduce with MPI_MIN
    - increment time: t += dt
    - increment iteration counter
    - swap u and v with their previous counterparts (u_old, v_old)

    - compute intermediate velocity
        - compute local start/end indices based on grid location of this proc (use has_bc_* flags)
        - for i = i_start..i_end, j = j_start..j_end:
            - update u(i,j), v(i,j) with diffusion terms
            - apply x-advection using upwind (u_old)
            - apply y-advection using upwind (v_old)

    - exchange u, v halos with neighbours (non-blocking MPI sends/receives under utils::exchange_halos)

    - calculate ppm RHS
        - compute local start/end indices based on grid location of this proc
        - for i = i_start..i_end, j = j_start..j_end:
            - compute local divergence from centered differences of u, v
            - set PPrhs(i,j) = rho / dt * div(u,v)

    - solve pressure Poisson with Jacobi
        - set local iteration counter = 0, local_diff = 0, local_sum = 0, tol = large
        - while tol > rtol and it < max_it; do
            - swap P and P_old
            - reset local_diff, local_sum
            - for i = i_start..i_end, j = j_start..j_end:
                - update P(i,j) via Jacobi stencil using P_old and PPrhs
                - accumulate local_sum += |P(i,j)|
                - accumulate local_diff += |P(i,j) - P_old(i,j)|

            - apply pressure boundary conditions locally (set_pressure_BCs)
            - exchange P halos with neighbours (utils::exchange_halos)

            - reduce local_diff → global_diff with MPI_Allreduce (MPI_SUM)
            - reduce local_sum → global_sum with MPI_Allreduce (MPI_SUM)

            - compute global relative residual: tol = global_diff / max(global_sum, 1e-10)
            - increment Jacobi iteration counter

    - project velocity with pressure gradient
        - compute local start/end indices
        - set local vel_max = 0
        - for i = i_start..i_end, j = j_start..j_end:
            - correct u(i,j), v(i,j) using centered pressure gradient in x and y
            - compute local speed = sqrt(u(i,j)^2 + v(i,j)^2)
            - update local vel_max = max(local vel_max, speed)

    - apply velocity boundary conditions locally (set_velocity_BCs on u)

    - if t >= next_output_time; do
        - increment output counter
        - advance next_output_time by dt_out
        - on rank 0: print iteration, time, Jacobi iterations, vel_max
        - for each field P, u, v:
            - append local subdomain data to rank-specific file
              (each rank writes its own chunk)

end while
```

:wq