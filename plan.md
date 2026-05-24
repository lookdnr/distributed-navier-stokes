# CW scrapbook

## Communications

1. Collectively communicate the appropriate (minimum) time step
2. Transfer domain boundary intermediate velocities to neighbours
3. Transfer domain boundary pressures to neighbours
4. Collectively communicate the residual
5. Transfer domain boundary final velocities to neighbours

### Notes

Note 2, 3, and 5 are the same process, essentially, so it would be good to define generic MPI datatypes so we can write generic functions that communicate boundaries.

We need to change the datatypes to ensure they are contiguous.

```plaintext

  top buffer -> |           |
             -------------------
                |           |
                | Subdomain |   -> proc i + 1
                |           |
             -------------------
                |           |   <- no need to transfer into corners
              ^
              |
          side buffer
```

Probably best to assign each proc an `i`, `j` identifier to indicate grid position to ease communications.

In additon:

- Communications should be non blocking for P2P.
- Decomposition logic should be custom.
- Custom MPI types should be used.

## Deliverables

- [x] Documented source
- [x] Compile + run docs
- [x] Extension docs
- [x] Output animation
- [x] Post proc script for results + anim

## Roadmap

1. ~~Refactor the solver to use contiguous local storage with ghost cells so boundary exchange can use MPI datatypes cleanly.~~
2. ~~Write custom domain decomposition logic that maps ranks onto a 2D process grid.~~
3. ~~Add reusable non-blocking halo exchange functions for the intermediate velocities, pressure field, and final velocities.~~
4. ~~Use custom MPI datatypes for row and column boundaries, then pair them with `MPI_Isend`/`MPI_Irecv` and `MPI_Waitall`.~~
5. ~~Add collective reductions for the timestep and Poisson residual, then finish with per-rank output plus a separate post-processing script.~~

## Working notes:

In Domain.decompose(), we should compute things like

- grid i, j locations
- local nx/ ny
- ranks of neighbouring processes
- whether or not to compute BCs on an edge of the subdomain

This merits a SubDomain class which will store this metadata. Proc 0 can be repsonsible for scattering metadata (as opposed to actual grid data), and each proc can instantiate objects of the Domain class using this metadata.

### Day 1 progress:

- [x] Implemented Domain and SubDomain
- [x] Implemented core decomposition logic
- [x] Implemented custom MPI type for SubDomain
- [x] Implemented scatter 

Next:

- Apply the decomp metadata: create Domain objects on each proc with relevant data
- Create MPI types for communicating boundary problems. Swaps should probably be done with `MPI_Sendrecv`. 

### Day 2:

For communicating halos:

- send left interior: field(1, 1) with MPI_Domainrow
- recv left ghost: field(0, 1) with MPI_Domainrow
- send right interior: field(nx, 1) with MPI_Domainrow
- recv right ghost: field(nx+1, 1) with MPI_Domainrow
- send down interior: field(1, 1) with MPI_Domaincol
- recv down ghost: field(1, 0) with MPI_Domaincol
- send up interior: field(1, ny) with MPI_Domaincol
- recv up ghost: field(1, ny+1) with MPI_Domaincol

#### Solver refactor

Communications are in place. They run, but do still need to be verified.

We need to refactor the solver, also. The serial versions run on global field allocations and the solver internals rely on them. This is completely inefficient. The pattern should be:

- Proc 0 creates some field that is `Nx` by `Ny`, computes the decomposition metadata, and scatters it.
- Each proc instatiates its own subdomains and initialises them (including the inlet pressure condition). 
   - Ideally this would just be creating `P` and another variable, and setting the rest to be a copy of the other. We can't copy out `P` since it contains the inlet condition, but all others are the same shape and initialised to zero. 

So, there is lots we need to do:

- [x] Remove global allocations
- [x] Remove/ gut `setup()`
- [x] Adapt function signatures and implementation
- [ ] (?) ~~implement a copy constructor for easy allocation. Maybe verbosity is better here.~~
- [x] Add to writing logic to enable clear writing to file for each proc

### Running the solver

With everything refactored and file writes working, we just need to make sure the exchange logic is working as intended before getting the solver live.

- [x] Verify halo exchange logic
- [ ] Add collective reductions for the timestep

Note for the last point: the timestep used should be the _global_ minimum.

After that, we can get on with postproc and extensions.

Note: halo exchange has been verified. Change debug_exchange to true to see halo writes

#### Logistics

Halo exchange:

- After computing intermediate velocities excahnge `u`, `v` since PPrhs needs them
- After each Jacobi iteratio excahnge P halos

Reductions:

- Reduce dt to get global minimum
- Inside Jacobi loop, reduce tolerance to ensure global convergence
- (Optional) reduce `vel_max` for console printing

## Day 3: extension

I am most interested in the field obstacle, so let's try that out. Here is a plan:

- We will use a square obstacle at first, as this will be most straigthforwrd. We can create an obstacle mask as a Domain object built from global indices
- In the solve (intermediate v, pressure RHS, pressure Poisson, v projection), we need to check:
   - if the cell is solid, do not update it as a fluid
   - if it is fluid, continue as normal unless one of its neighbours is solid
- In the obstacle, we should for zero velocity. The boundaries of the obstacle are walls, and we will make them no slip. If a fludi cell has a solid neighbour, the obstacle face should have a zero normal pressure gradient.
