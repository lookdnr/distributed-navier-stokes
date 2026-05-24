#include "mpi.h"
#include <cmath>
#include <ctime>
#include <fstream>
#include <iostream>
#include <sstream>

#include "solver.h"

using namespace CustomTypes; // declare namespace for custom MPI types

using namespace std;

int Nx = 401;
int Ny = 101;

double Lx = 0.1, Ly = 0.05;
const double rho = 1000, nu = 1e-6;
const double P_max = 0.5;
const double t_end = 50.0;
const double dt_min = 1.e-3;
const double courant = 0.01;
const double dt_out = 0.5;

int id, p; // MPI variables (rank/ size)

double dx, dy, dt, t;

// write sim data to structured metadata file
void write_solver_parameters(const string root = "./out/core") {
  // create dir if it doesnt exist
  if (!filesystem::exists(root))
    filesystem::create_directories(root);

  std::stringstream fname;
  std::fstream f1;

  fname << root + "/metadata" << ".dat";

  // check is open
  f1.open(fname.str().c_str(), std::ios_base::out);
  if (!f1.is_open()) {
    cerr << "Failed to open " << fname.str() << ": " << strerror(errno) << "\n";
    cout.flush();
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  // get curr data and time
  time_t timestamp;
  time(&timestamp);
  f1 << "# Solver run for " << ctime(&timestamp) << '\n'
     << "# Nx: " << Nx << '\n'
     << "# Ny: " << Ny << '\n'
     << "# Lx: " << Lx << '\n'
     << "# Ly: " << Ly << '\n'
     << "# dx: " << dx << '\n'
     << "# dy: " << dy << '\n'
     << "# dt0: " << dt << '\n'
     << "# rho: " << rho << '\n'
     << "# nu: " << nu << '\n'
     << "# courant: " << courant << '\n'
     << "# P_max: " << P_max << '\n'
     << "# t_end: " << t_end << '\n'
     << "# dt_out: " << dt_out << '\n';
}

void grids_to_file(int timestep, const Domain &P, const Domain &u,
                   const Domain &v, const SubDomain &sd,
                   const bool debug_write = false,
                   const string root = "./out/core") {
  // loop only over the local nx, ny
  int nx = sd.nx_local;
  int ny = sd.ny_local;
  int i_min = 1;
  int j_min = 1;
  int rank = sd.rank;

  // for debugging, write halos
  if (debug_write) {
    nx++;
    ny++;
    i_min--;
    j_min--;
  }

  // Write the output for a single time step to file
  stringstream fname;
  fstream f1;

  fname << root + "/P/P_rank" << rank << ".dat";
  f1.open(fname.str().c_str(), ios_base::out | ios_base::app);

  // write timestep metadata
  f1 << "$ timestep " << timestep << '\n' << "$ t = " << t << '\n';
  for (int i = i_min; i <= nx;
       i++) // note here and all other loops we skip the halo
  {
    for (int j = j_min; j <= ny; j++)
      f1 << P(i, j) << "\t";
    f1 << endl;
  }
  f1.close();

  fname.str("");
  fname << root + "/u/u_rank" << rank << ".dat";
  f1.open(fname.str().c_str(), ios_base::out | ios_base::app);

  f1 << "$ timestep " << timestep << '\n' << "$ t = " << t << '\n';
  for (int i = i_min; i <= nx; i++) {
    for (int j = j_min; j <= ny; j++)
      f1 << u(i, j) << "\t";
    f1 << endl;
  }

  f1.close();

  fname.str("");
  fname << root + "/v/v_rank" << rank << ".dat";
  f1.open(fname.str().c_str(), ios_base::out | ios_base::app);

  f1 << "$ timestep " << timestep << '\n' << "$ t = " << t << '\n';
  for (int i = i_min; i <= nx; i++) {
    for (int j = j_min; j <= ny; j++)
      f1 << v(i, j) << "\t";
    f1 << endl;
  }
  f1.close();
}

// compute global parameters
void setup(void) {
  dx = Lx / (Nx - 1);
  dy = Ly / (Ny - 1);

  t = 0.0;
}

void calculate_ppm_RHS_central(Domain &PPrhs, Domain &u, Domain &v,
                               const SubDomain &sd) {
  const int nx = sd.nx_local;
  const int ny = sd.ny_local;
  const int i_start = sd.has_bc_left ? 2 : 1;
  const int i_end = sd.has_bc_right ? nx - 1 : nx;
  const int j_start = sd.has_bc_down ? 2 : 1;
  const int j_end = sd.has_bc_up ? ny - 1 : ny;

  for (int i = i_start; i <= i_end; i++)
    for (int j = j_start; j <= j_end; j++) {
      PPrhs(i, j) = rho / dt *
                    ((u(i + 1, j) - u(i - 1, j)) / (2. * dx) +
                     (v(i, j + 1) - v(i, j - 1)) / (2. * dy));
    }
}

void set_pressure_BCs(Domain &P, const SubDomain &sd) {
  const int nx = sd.nx_local;
  const int ny = sd.ny_local;

  // fixed pressure at the inlet
  if (sd.has_bc_left) {
    for (int j = 1; j <= ny; ++j)
      P(1, j) = P_max;
  }

  // No slip over top and bottom walls, so zero normal gradient for pressure
  if (sd.has_bc_down) {
    for (int i = 1; i <= nx; ++i)
      P(i, 1) = P(i, 2);
  }

  if (sd.has_bc_up) {
    for (int i = 1; i <= nx; ++i)
      P(i, ny) = P(i, ny - 1);
  }

  // No slip wall over the upper half of the right boundary
  if (sd.has_bc_right) {
    for (int j = 1; j <= ny; ++j) {
      const int global_j = static_cast<int>(sd.y0_global) + (j - 1);
      if (global_j >= Ny / 2)
        P(nx, j) = P(nx - 1, j);
    }
  }
}

int pressure_poisson_jacobi(Domain &P, Domain &P_old, Domain &PPrhs,
                            const SubDomain sd, double rtol = 1.e-5) {
  double tol = 10. * rtol; // global tolerance
  int it = 0;

  const int nx = sd.nx_local;
  const int ny = sd.ny_local;
  const int i_start = sd.has_bc_left ? 2 : 1;
  const int i_end = sd.has_bc_right ? nx - 1 : nx;
  const int j_start = sd.has_bc_down ? 2 : 1;
  const int j_end = sd.has_bc_up ? ny - 1 : ny;

  const int max_it = 10000;

  while (tol > rtol) {
    P.swap(P_old);
    // local accumulators
    double local_diff = 0.0;
    double local_sum = 0.0;
    local_diff = 0.0;
    it++;

    // Jacobi iteration
    for (int i = i_start; i <= i_end; i++)
      for (int j = j_start; j <= j_end; j++) {
        P(i, j) = 1.0 / (2.0 + 2.0 * (dx * dx) / (dy * dy)) *
                  (P_old(i + 1, j) + P_old(i - 1, j) +
                   (P_old(i, j + 1) + P_old(i, j - 1)) * (dx * dx) / (dy * dy) -
                   (dx * dx) * PPrhs(i, j));

        local_sum += fabs(P(i, j));
        local_diff += fabs(P(i, j) - P_old(i, j));
      }

    set_pressure_BCs(P, sd);

    // exchange halos
    utils::exchange_halos(P, sd);

    // reduce the residual and field magnitude across all ranks
    double global_diff = 0.0;
    double global_sum = 0.0;
    MPI_Allreduce(&local_diff, &global_diff, 1, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
    MPI_Allreduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);

    // compute the global relative residual.
    tol = global_diff / std::max(1.e-10, global_sum);

    if (it >= max_it) {
      if (id == 0) {
        std::cerr << "Warning: Jacobi reached max iterations (" << max_it
                  << ") with residual " << tol << std::endl;
      }
      break;
    }
  }
  return it;
}

void calculate_intermediate_velocity(Domain &u, Domain &u_old, Domain &v,
                                     Domain &v_old, const SubDomain &sd) {
  const int nx = sd.nx_local;
  const int ny = sd.ny_local;
  // compoute local start indices based on grid pos
  const int i_start = sd.has_bc_left ? 2 : 1;
  const int i_end = sd.has_bc_right ? nx - 1 : nx;
  const int j_start = sd.has_bc_down ? 2 : 1;
  const int j_end = sd.has_bc_up ? ny - 1 : ny;

  for (int i = i_start; i <= i_end; i++)
    for (int j = j_start; j <= j_end; j++) {
      // viscous diffusion
      u(i, j) = u_old(i, j) +
                dt * nu *
                    ((u_old(i + 1, j) + u_old(i - 1, j) - 2.0 * u_old(i, j)) /
                         (dx * dx) +
                     (u_old(i, j + 1) + u_old(i, j - 1) - 2.0 * u_old(i, j)) /
                         (dy * dy));
      v(i, j) = v_old(i, j) +
                dt * nu *
                    ((v_old(i + 1, j) + v_old(i - 1, j) - 2.0 * v_old(i, j)) /
                         (dx * dx) +
                     (v_old(i, j + 1) + v_old(i, j - 1) - 2.0 * v_old(i, j)) /
                         (dy * dy));
      // advection - upwinding
      if (u_old(i, j) > 0.0) {
        u(i, j) -= dt * u_old(i, j) * (u_old(i, j) - u_old(i - 1, j)) / dx;
        v(i, j) -= dt * u_old(i, j) * (v_old(i, j) - v_old(i - 1, j)) / dx;
      } else {
        u(i, j) -= dt * u_old(i, j) * (u_old(i + 1, j) - u_old(i, j)) / dx;
        v(i, j) -= dt * u_old(i, j) * (v_old(i + 1, j) - v_old(i, j)) / dx;
      }
      if (v_old(i, j) > 0.0) {
        u(i, j) -= dt * v_old(i, j) * (u_old(i, j) - u_old(i, j - 1)) / dy;
        v(i, j) -= dt * v_old(i, j) * (v_old(i, j) - v_old(i, j - 1)) / dy;
      } else {
        u(i, j) -= dt * v_old(i, j) * (u_old(i, j + 1) - u_old(i, j)) / dy;
        v(i, j) -= dt * v_old(i, j) * (v_old(i, j + 1) - v_old(i, j)) / dy;
      }
    }
}

void set_velocity_BCs(Domain &u, const SubDomain &sd) {
  const int nx = sd.nx_local;
  const int ny = sd.ny_local;

  // Zero velocity gradients at pressure inlet and outlet
  if (sd.has_bc_left) {
    for (int j = 1; j <= ny; j++)
      u(1, j) = u(2, j);
  }

  if (sd.has_bc_right) {
    for (int j = 1; j <= ny; j++) {
      const int global_j = static_cast<int>(sd.y0_global) + (j - 1);
      if (global_j < Ny / 2)
        u(nx, j) = u(nx - 1, j);
    }
  }
}

double project_velocity(Domain &u, Domain &v, Domain &P, const SubDomain &sd) {
  const int nx = sd.nx_local;
  const int ny = sd.ny_local;
  // compoute local start indices based on grid pos
  const int i_start = sd.has_bc_left ? 2 : 1;
  const int i_end = sd.has_bc_right ? nx - 1 : nx;
  const int j_start = sd.has_bc_down ? 2 : 1;
  const int j_end = sd.has_bc_up ? ny - 1 : ny;

  double vmax = 0.0;
  for (int i = i_start; i <= i_end; i++)
    for (int j = j_start; j <= j_end; j++) {
      u(i, j) =
          u(i, j) - dt * (1. / rho) * (P(i + 1, j) - P(i - 1, j)) / (2. * dx);
      v(i, j) =
          v(i, j) - dt * (1. / rho) * (P(i, j + 1) - P(i, j - 1)) / (2. * dy);

      double vel = sqrt(u(i, j) * u(i, j) + v(i, j) * v(i, j));

      vmax = max(vmax, vel);
    }

  return vmax;
}

void solve_NS(Domain &P, Domain &P_old, Domain &PPrhs, Domain &u, Domain &u_old,
              Domain &v, Domain &v_old, const SubDomain sd) {
  double vel_max = 0.0;
  int time_it = 0;
  int its;
  int out_it = 0;
  double t_out = dt_out;

  grids_to_file(out_it, P, u, v, sd);

  while (t < t_end) {
    if (vel_max > 0.0) {
      dt = min(courant * min(dx, dy) / vel_max, dt_min);
    } else
      dt = dt_min;

    // reduce to find global minimum and set in place
    MPI_Allreduce(MPI_IN_PLACE, &dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

    t += dt;
    time_it++;
    u.swap(u_old);
    v.swap(v_old);

    calculate_intermediate_velocity(u, u_old, v, v_old, sd);

    // exchange halos since updated u and v required for ppm RHS
    utils::exchange_halos(u, sd);
    utils::exchange_halos(v, sd);

    calculate_ppm_RHS_central(PPrhs, u, v, sd);
    its = pressure_poisson_jacobi(P, P_old, PPrhs, sd, 1.e-5);
    vel_max = project_velocity(u, v, P, sd);
    set_velocity_BCs(u, sd);

    if (t >= t_out) {
      out_it++;
      t_out += dt_out;
      if (id == 0) {
        cout << "Iteration " << time_it << ": t = " << t
             << "\n\tJacobi iterations: " << its << " vel_max: " << vel_max
             << endl;
        cout.flush();
      }
      grids_to_file(out_it, P, u, v, sd);
    }
  }
}

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &p);
  MPI_Comm_rank(MPI_COMM_WORLD, &id);

  // debug flag. set to try to verify halo exchanges
  // if true, fills subdomains with rank number and fill writes contain halos
  const bool debug_exchange = false;

  // parse optional command-line overrides for grid and domain size
  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    if ((arg == "--Nx" || arg == "-nx") && i + 1 < argc) {
      Nx = std::stoi(argv[++i]);
    } else if ((arg == "--Ny" || arg == "-ny") && i + 1 < argc) {
      Ny = std::stoi(argv[++i]);
    } else if ((arg == "--Lx" || arg == "-lx") && i + 1 < argc) {
      Lx = std::stod(argv[++i]);
    } else if ((arg == "--Ly" || arg == "-ly") && i + 1 < argc) {
      Ly = std::stod(argv[++i]);
    } else if (arg == "-h" || arg == "--help") {
      if (id == 0)
        std::cout << "Usage: [--Nx N] [--Ny N] [--Lx val] [--Ly val]"
                  << std::endl;
      MPI_Finalize();
      return 0;
    }
  }

  // build subdomain type
  build_subdomain_type();

  setup(); // set dx, dy, t=0
  write_solver_parameters();

  vector<SubDomain>
      decomp_metadata; // to contain decomp instructions for scatter

  // decompose on proc 0
  if (id == 0) {
    bool use_buffer = false; // no need to use buffers on setup
    Domain P = Domain(Nx, Ny, use_buffer);
    decomp_metadata = P.decompose(p); // compute the decomposition metadata once

    if (p > 1)
      utils::print_metadata(decomp_metadata); // print computed decomp metadata
  }

  // scatter metadata
  SubDomain subdomain;
  MPI_Scatter(decomp_metadata.data(), 1, MPI_SubDomain, &subdomain, 1,
              MPI_SubDomain, 0, MPI_COMM_WORLD);

  // give each proc its own local domain
  // Note all but P are identical, hence the template field
  Domain template_field = utils::build_local_domain(subdomain, true);
  Domain u = template_field;
  Domain u_old = template_field;
  Domain v = template_field;
  Domain v_old = template_field;
  Domain P_old = template_field;
  Domain PPrhs = template_field;

  Domain P = utils::build_local_domain(subdomain, true, P_max, debug_exchange);

  // build halo type once per proc
  // this is required since each proc could have differently shaped halos
  build_halo_types(subdomain);

  utils::prep_files(subdomain);
  solve_NS(P, P_old, PPrhs, u, u_old, v, v_old, subdomain);

  if (id == 0)
    cout << endl << "Solve complete." << endl;

  MPI_Type_free(&MPI_SubDomain);
  MPI_Type_free(&MPI_Domaincol);
  MPI_Type_free(&MPI_Domainrow);
  MPI_Finalize();
  return 0;
}
