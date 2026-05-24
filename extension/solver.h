#pragma once
#include <omp.h>
#include <assert.h>
#include <iostream>
#include <utility> // for std::swap
#include <vector>
#include <math.h>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <filesystem>
#include <string>

struct SubDomain;

struct Domain
{
    public:
        size_t Nx, Ny; // exposed shape descriptors
        bool use_buffer; // 
        
        // constructors
        Domain() { // default
            Nx = 0;
            Ny = 0;
            use_buffer = false;
            buff_size = 0;
            nx_interior = 0;
            ny_interior = 0;
            data.clear();
        }

        Domain(size_t nx, size_t ny, bool use_buffer, const double fill_value=0.0)
        {
            assert(nx >= 0);
            assert(ny >= 0);

            // set shape descriptors
            Nx = nx; Ny = ny;
            nx_interior = nx; ny_interior = ny;

            // set size of halo at domain edges
            buff_size = use_buffer ? 1 : 0;

            // allocate and set to zero
            data.assign((nx_interior + 2*buff_size) * (ny_interior + 2*buff_size), fill_value);
        }

        // overloaded operators: setter and getter.
        double& operator()(size_t i, size_t j) { return data[idx(i, j)]; }
        const double& operator()(size_t i, size_t j) const { return data[idx(i, j)]; }

        // swap to avoid copying arrays
        void swap(Domain& other)
        {
            std::swap(data, other.data);
            std::swap(buff_size, other.buff_size);
            std::swap(nx_interior, other.nx_interior);
            std::swap(ny_interior, other.ny_interior);
            std::swap(Nx, other.Nx);
            std::swap(Ny, other.Ny);
        }

        std::vector<SubDomain> decompose(const int num_procs); // domain decomposition

    private:
        // note data is stored in a flattened, row-major form
        std::vector<double> data; // contiguous flattened storage

        size_t buff_size; // buffer size, 0 by default
        size_t nx_interior, ny_interior; // num interior points: for loop clarity

        // indexer for flattened array
        // note because of the row-major order, COLUMNS are contiguous
        // consider the i = 0 case. the row is fixed, and we sweep the column
        size_t idx(const size_t i, const size_t j) const { return i * (ny_interior + 2*buff_size) + j; }
};

// metadata class
// domain.decompose creates an array of these objects which will be scattered out
// they should be used to inform each proc how to create objects of type Domain
// representing subsects of full Domain objects.
struct SubDomain
{
    // default constructro
    SubDomain()
        : nx_local(0), ny_local(0), x0_global(0), y0_global(0), buff_size(1),
        px(1), py(1), grid_i(0), grid_j(0), rank(-1),
        left_rank(MPI_PROC_NULL), right_rank(MPI_PROC_NULL), up_rank(MPI_PROC_NULL), down_rank(MPI_PROC_NULL),
        has_bc_left(false), has_bc_right(false), has_bc_up(false), has_bc_down(false)
        { }

    SubDomain(int Rank, size_t grid_I, size_t grid_J, size_t pX, size_t pY, size_t Nx_local, size_t Ny_local, size_t X0_global, size_t Y0_global) 
    : rank(Rank), grid_i(grid_I), grid_j(grid_J), px(pX), py(pY), nx_local(Nx_local), ny_local(Ny_local), x0_global(X0_global), y0_global(Y0_global) // assign some variables
    {
        assert(px > 0);
        assert(py > 0);
        assert(rank >= 0);
        assert(rank < px * py);

        // compute neighbour ranks using -1 as a sentinel for OOB
        // these conds say (is this on a boundary) ? no: neighbour = pm 1, yes: neighbour = -1
        // use MPI_PROC_NULL if no neighbours to avoid later nested ifs
        left_rank = (grid_i > 0) ? rank - 1 : MPI_PROC_NULL;
        right_rank = (grid_i < px - 1) ? rank + 1 : MPI_PROC_NULL;
        down_rank = (grid_j > 0) ? rank - px : MPI_PROC_NULL;
        up_rank = (grid_j < py - 1) ? rank + px : MPI_PROC_NULL;

        // boundary flags for edges of global domain
        has_bc_left = (grid_i == 0);
        has_bc_right = (grid_i == (px - 1));
        has_bc_down = (grid_j == 0);
        has_bc_up = (grid_j == (py - 1));
    }

    // data accessing logic
    size_t nx_local, ny_local; // local num points
    size_t x0_global, y0_global; // start indices in terms of the global Domain object

    size_t buff_size = 1; // buffer size: always buffer subdomains

    // grid logic
    size_t px, py; // number of procs in each grid direction
    size_t grid_i, grid_j; // grid coords

    // comms logic
    int rank; // store id of this proc
    int left_rank, right_rank, up_rank, down_rank; // rank of neighbouring processes
    
    // BC flags
    bool has_bc_left = false;
    bool has_bc_right = false;
    bool has_bc_up = false; 
    bool has_bc_down = false;

    // create MPI type to communicate
    void buildMPIType();

};

// namespace to encapsulate custom MPI types and building utils
namespace CustomTypes
{
    // declare types
    MPI_Datatype MPI_SubDomain;
    MPI_Datatype MPI_Domainrow;
    MPI_Datatype MPI_Domaincol;

    // utility to build subdomain type
    void build_subdomain_type() {
        const int n_fields = 18; // number of fields

        // define block lengths
        int blocklengths[n_fields] = {
        1, // rank
        1, 1, // grid_i, grid_j
        1, 1, // px, py
        1, 1, // nx_local, ny_local
        1, 1, // x0_global, y0_global
        1, // buff_size
        1, 1, 1, 1, // left, right, up, down ranks
        1, 1, 1, 1 // has_bc_*
        };

        // define datatypes
        MPI_Datatype types[n_fields] = {
            MPI_INT,
            MPI_UNSIGNED_LONG, MPI_UNSIGNED_LONG,
            MPI_UNSIGNED_LONG, MPI_UNSIGNED_LONG,
            MPI_UNSIGNED_LONG, MPI_UNSIGNED_LONG,
            MPI_UNSIGNED_LONG, MPI_UNSIGNED_LONG,
            MPI_UNSIGNED_LONG,
            MPI_INT, MPI_INT, MPI_INT, MPI_INT,
            MPI_C_BOOL, MPI_C_BOOL, MPI_C_BOOL, MPI_C_BOOL 
        };

        // write displacements
        MPI_Aint displacements[n_fields];
        displacements[0] = offsetof(SubDomain, rank);
        displacements[1] = offsetof(SubDomain, grid_i);
        displacements[2] = offsetof(SubDomain, grid_j);
        displacements[3] = offsetof(SubDomain, px);
        displacements[4] = offsetof(SubDomain, py);
        displacements[5] = offsetof(SubDomain, nx_local);
        displacements[6] = offsetof(SubDomain, ny_local);
        displacements[7] = offsetof(SubDomain, x0_global);
        displacements[8] = offsetof(SubDomain, y0_global);
        displacements[9] = offsetof(SubDomain, buff_size);
        displacements[10] = offsetof(SubDomain, left_rank);
        displacements[11] = offsetof(SubDomain, right_rank);
        displacements[12] = offsetof(SubDomain, up_rank);
        displacements[13] = offsetof(SubDomain, down_rank);
        displacements[14] = offsetof(SubDomain, has_bc_left);
        displacements[15] = offsetof(SubDomain, has_bc_right);
        displacements[16] = offsetof(SubDomain, has_bc_up);
        displacements[17] = offsetof(SubDomain, has_bc_down);

        // make and commit type
        MPI_Type_create_struct(n_fields, blocklengths, displacements, types, &MPI_SubDomain);
        MPI_Type_commit(&MPI_SubDomain);
    }

    void build_halo_types(const SubDomain& sd) {
        const int nx = sd.nx_local;
        const int ny = sd.ny_local;
        const int halo = 1;
        const int stride = ny + 2*halo; // for indexing

        // recall that due to the row-major storage,columns are contiguous
        MPI_Type_contiguous(ny, MPI_DOUBLE, &MPI_Domaincol);
        MPI_Type_commit(&MPI_Domaincol);

        MPI_Type_vector(nx, 1, stride, MPI_DOUBLE, &MPI_Domainrow);
        MPI_Type_commit(&MPI_Domainrow);
    }
}

// domain decomposition implementation
// decompose into num_procs approximately square grid units
// works by computing px, py pairs and computing the aspect ratio of the blocks
// accepts the pair that produces the best (closest to unity) aspect ratio
std::vector<SubDomain> Domain::decompose(const int num_procs) 
{
    assert(num_procs > 0);

    // vector of subdomains
    std::vector<SubDomain> subdomains;
    subdomains.reserve(num_procs);

    int best_px = 1;
    int best_py = num_procs;
    double best_score = std::numeric_limits<double>::infinity();

    // find factor pair with best aspect ratio
    for (int px = 1; px <= num_procs; ++px)
    {
        if (num_procs % px != 0) continue;
        int py = num_procs / px; // compute py from px

        // compute aspect ratio
        // we want this close to 1
        double block_aspect =
            (static_cast<double>(Nx) * static_cast<double>(py)) /
            (static_cast<double>(Ny) * static_cast<double>(px));

        // compoute 'score': squared distance of curr aspect from ideal unity
        double score = std::pow(block_aspect - 1, 2);

        // update current best px and py if score lower
        if (score < best_score)
        {
            best_score = score;
            best_px = px;
            best_py = py;
        }
    }

    // set px and py with best found grid layout
    const int px = best_px;
    const int py = best_py;

    // compute local portions of Nx/ Ny and remaining elements
    const size_t nx_base = Nx / static_cast<size_t>(px);
    const size_t nx_rem = Nx % static_cast<size_t>(px);
    const size_t ny_base = Ny / static_cast<size_t>(py);
    const size_t ny_rem = Ny % static_cast<size_t>(py);

    // create SubDomain objects
    for (int rank = 0; rank < num_procs; ++rank)
    {
        // compute grid pos
        const size_t grid_i = static_cast<size_t>(rank % px);
        const size_t grid_j = static_cast<size_t>(rank / px);

        // start from the even split, then add one extra cell for the first few subdomains
        size_t extra_x = 0;
        if (grid_i < nx_rem) extra_x = 1;

        size_t extra_y = 0;
        if (grid_j < ny_rem) extra_y = 1;

        const size_t nx_local = nx_base + extra_x;
        const size_t ny_local = ny_base + extra_y;

        // compute the global start index by adding the base blocks
        // also add any earlier remainder cells
        size_t x0_global = grid_i * nx_base;
        if (grid_i < nx_rem) x0_global += grid_i;
        else x0_global += nx_rem;

        size_t y0_global = grid_j * ny_base;
        if (grid_j < ny_rem) y0_global += grid_j;
        else y0_global += ny_rem;

        // add a SubDomain element to the vector
        subdomains.push_back(SubDomain(rank, grid_i, grid_j, 
            static_cast<size_t>(px), static_cast<size_t>(py), 
            nx_local, ny_local, x0_global, y0_global));
    }
    return subdomains;
}

// utility namepsace for encapsulation
namespace utils
{
    // print metadata on computed decomposition
    void print_metadata(std::vector<SubDomain>& metadata) {
        std::cout << "Proc 0 decomposed the domain:" << std::endl;
        std::cout.flush();
        std::cout << "\t- Total: " << metadata.size() << " domains." << std::endl;
        std::cout.flush();
        std::cout << "\t- Shape: " << metadata[0].px << " by " << metadata[0].py << std::endl << std::endl;
        std::cout.flush();
        for (auto sd : metadata) {
            std::cout << "Proc " << sd.rank << " is " << sd.nx_local << " by " << sd.ny_local <<
                        " with grid identifier " << sd.grid_i << ", " << sd.grid_j <<  std::endl;
            std::cout.flush();
            std::cout << "\tNeighbouring procs:";
            // print each neighbour only if it exists (rank != -1)
            if (sd.left_rank != MPI_PROC_NULL) {
                const auto &n = metadata[sd.left_rank];
                std::cout << "\t" << sd.left_rank << " (left: " << n.grid_i << "," << n.grid_j << ")";
            }
            if (sd.right_rank != MPI_PROC_NULL) {
                const auto &n = metadata[sd.right_rank];
                std::cout << "\t" << sd.right_rank << " (right: " << n.grid_i << "," << n.grid_j << ")";
            }
            if (sd.up_rank != MPI_PROC_NULL) {
                const auto &n = metadata[sd.up_rank];
                std::cout << "\t" << sd.up_rank << " (above: " << n.grid_i << "," << n.grid_j << ")";
            }
            if (sd.down_rank != MPI_PROC_NULL) {
                const auto &n = metadata[sd.down_rank];
                std::cout << "\t" << sd.down_rank << " (below: " << n.grid_i << "," << n.grid_j << ")";
            }
            std::cout << std::endl << std::endl;
            std::cout.flush();
        }
        std::cout.flush();
    }

    Domain build_local_domain(const SubDomain& sd, bool is_P=false, double P_max=0.0, const bool test_exchange=false) {
        // allocate local domain with halos

        // fill with rank n for excahnge testing purposes, else init to zero
        double fill_value = 0.0;
        if (test_exchange) fill_value = static_cast<double>(sd.rank);  
        
        // create Domain object
        Domain local(sd.nx_local, sd.ny_local, true, fill_value);

        // if on global inlet (left edge), set inlet pressure on physical boundary
        if (sd.has_bc_left && is_P && !test_exchange) {
            for (size_t j = 1; j <= sd.ny_local; j++) {
                local(1, j) = P_max;  // index 1 is first interior cell
            }
        }
        return local;  // return by value
    }

    void exchange_halos(Domain& field, const SubDomain& sd) {
        using CustomTypes::MPI_Domainrow;
        using CustomTypes::MPI_Domaincol;

        // one request per halo transfer: 4 receives + 4 sends
        MPI_Request requests[8];
        int req_count = 0;

        // comms tags for debugging
        constexpr int TAG_HALO_ROW = 100;
        constexpr int TAG_HALO_COL = 200;

        /*Post receives*/
        // left/ right receive
        MPI_Irecv(&field(0, 1), 1, MPI_Domaincol, sd.left_rank, TAG_HALO_COL, MPI_COMM_WORLD, &requests[req_count++]);
        MPI_Irecv(&field(sd.nx_local + 1, 1), 1, MPI_Domaincol, sd.right_rank, TAG_HALO_COL, MPI_COMM_WORLD, &requests[req_count++]);

        MPI_Irecv(&field(1, 0), 1, MPI_Domainrow, sd.down_rank, TAG_HALO_ROW, MPI_COMM_WORLD, &requests[req_count++]);
        MPI_Irecv(&field(1, sd.ny_local + 1), 1, MPI_Domainrow, sd.up_rank, TAG_HALO_ROW, MPI_COMM_WORLD, &requests[req_count++]);

        // then send the interior faces
        MPI_Isend(&field(1, 1), 1, MPI_Domaincol, sd.left_rank, TAG_HALO_COL, MPI_COMM_WORLD, &requests[req_count++]);
        MPI_Isend(&field(sd.nx_local, 1), 1, MPI_Domaincol, sd.right_rank, TAG_HALO_COL, MPI_COMM_WORLD, &requests[req_count++]);
        MPI_Isend(&field(1, 1), 1, MPI_Domainrow, sd.down_rank, TAG_HALO_ROW, MPI_COMM_WORLD, &requests[req_count++]);
        MPI_Isend(&field(1, sd.ny_local), 1, MPI_Domainrow, sd.up_rank, TAG_HALO_ROW, MPI_COMM_WORLD, &requests[req_count++]);

        MPI_Waitall(req_count, requests, MPI_STATUSES_IGNORE);
    }
    
    void check_is_open(std::fstream& f, const std::stringstream& fname, const SubDomain& sd) {
        // check file opens
        if (!f.is_open()) {
        std::cerr << "Rank " << sd.rank << " failed to open " << fname.str() << ": " << strerror(errno) << "\n";
        std::cout.flush();
        MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }

    // utility to prep write files by writing subdomain metadata
    void prep_files(const SubDomain& sd, const std::string root="./out/extension") {

        // check dirs exist
        bool P_dir, u_dir, v_dir;
        P_dir = std::filesystem::exists(root +"/P");
        u_dir = std::filesystem::exists(root + "/u");
        v_dir = std::filesystem::exists(root + "/v");

        // create if they don't
        if (!(P_dir)) std::filesystem::create_directories(root + "/P");
        if (!(u_dir)) std::filesystem::create_directories(root + "/u");
        if (!(v_dir)) std::filesystem::create_directories(root + "/v");

        std::stringstream fname;
	    std::fstream f1;

        fname << root + "/P/P_rank" << sd.rank << ".dat";
	    f1.open(fname.str().c_str(), std::ios_base::out);
        check_is_open(f1, fname, sd);
        
        f1 << "# Pressure history" << '\n'
            << "# rank: " << sd.rank << '\n'
            << "# grid i: " << sd.grid_i << '\n' 
            << "# grid j: " << sd.grid_j << '\n'
            << "# i0: " << sd.x0_global << '\n'
            << "# j0: " << sd.y0_global << '\n'
            << "# nx: " << sd.nx_local << '\n'
            << "# ny: " << sd.ny_local << '\n' << '\n';
        f1.close();
        fname.str("");

        fname << root + "/u/u_rank" << sd.rank << ".dat";
	    f1.open(fname.str().c_str(), std::ios_base::out);
        check_is_open(f1, fname, sd);

        f1 << "# Velocity (u) history" << '\n'
            << "# rank: " << sd.rank << '\n'
            << "# grid i: " << sd.grid_i << '\n' 
            << "# grid j: " << sd.grid_j << '\n'
            << "# i0: " << sd.x0_global << '\n'
            << "# j0: " << sd.y0_global << '\n'
            << "# nx: " << sd.nx_local << '\n'
            << "# ny: " << sd.ny_local << '\n' << '\n';
        f1.close();
        fname.str("");

        fname << root + "/v/v_rank" << sd.rank << ".dat";
	    f1.open(fname.str().c_str(), std::ios_base::out);
        check_is_open(f1, fname, sd);
        
        f1 << "# Velocity (v) history" << '\n'
            << "# rank: " << sd.rank << '\n'
            << "# grid i: " << sd.grid_i << '\n' 
            << "# grid j: " << sd.grid_j << '\n'
            << "# i0: " << sd.x0_global << '\n'
            << "# j0: " << sd.y0_global << '\n'
            << "# nx: " << sd.nx_local << '\n'
            << "# ny: " << sd.ny_local << '\n' << '\n';
        f1.close();
    }
}

namespace extension {
    Domain build_square_obstacle(const SubDomain& sd, const int half_w, const int nx_global, const int ny_global)
    {
        Domain obstacle(sd.nx_local, sd.ny_local, true, 0.0);

        // compute global centre indices and half width
        const int ic = nx_global / 2;
        const int jc = ny_global / 2;

        // compute bounds of square
        const int i_min = ic - half_w;
        const int i_max = ic + half_w;
        const int j_min = jc - half_w;
        const int j_max = jc + half_w;

        // loop over local domain
        for (int i = 1; i <= sd.nx_local; ++i) {
            const int gi = sd.x0_global + (i - 1); // global i index

            for (int j = 1; j <= sd.ny_local; ++j) {
                const int gj = sd.y0_global + (j - 1); // global j index
                
                // if within the square, set value to 1 as an indicator
                if (gi >= i_min && gi <= i_max && gj >= j_min && gj <= j_max) {
                    obstacle(i, j) = 1.0;
                }
            }
        }
        return obstacle;
    }
}