//
// Created by Raunak Farhaz on 07.09.22.

//// Solving the Hartree-Fock equations Be atom
// note :: psi is a vector consisting of phi which are the spin orbitals

#include <madness/mra/mra.h>
#include <madness/mra/vmra.h>
#include <madness/mra/operator.h>
#include <madness/mra/nonlinsol.h>

using namespace madness;
const double L = 80.0;

// struct 1 :: unknown
template<typename T, std::size_t NDIM>
struct allocator {
    World &world;
    const int n;

    /// @param[in]	world	the world
    /// @param[in]	nn		the number of functions in a given vector
    allocator(World &world, const int nn) :
            world(world), n(nn) {
    }

    /// allocate a vector of n empty functions
    std::vector<Function<T, NDIM> > operator()() {
        return zero_functions<T, NDIM>(world, n);
    }
};

//guess functions :: sto-2g turbomole
static double guess0(const coord_3d &coords) {
    const double x = coords[0], y = coords[1], z = coords[2];
    const double r = sqrt(x * x + y * y + z * z );
    return 0.43 * exp(-11.53 * r * r) + 0.67 * exp(-2.05 * r * r);
}

static double guess1(const coord_3d &coords) {
    const double x = coords[0], y = coords[1], z = coords[2];
    const double r = sqrt(x * x + y * y + z * z );
    return 0.049 * exp(-0.5 * r * r) + 0.96 * exp(-0.13 * r * r);
}

// Electron-Nuclear attraction potential
static double nuclear_potential(const coord_3d &coords) {
    double Z = 4.0; // Atomic Charge of the desired element
    double x = coords[0], y = coords[1], z = coords[2];
    double r = sqrt(x * x + y * y + z * z + 1e-6);
    return -Z / r;
}

// computing electron density
real_function_3d compute_density(World &world,
                                 const std::vector<real_function_3d> &psi) {
    real_function_3d rho = real_factory_3d(world);
    for (int i = 0; i < psi.size(); ++i) {
        rho += 2.0 * square(psi[i]);
    }
    return rho;
}

//computing overlap
double compute_overlap(World &world,
                       const real_function_3d &bra,
                       const real_function_3d &ket) {
    double s = inner(bra, ket);
    return s;
}

//computing overlap matrix
Tensor<double> Smat(World &world,
                    const std::vector<real_function_3d> &bra,
                    const std::vector<real_function_3d> &ket) {
    int n = bra.size();
    Tensor<double> S(n, n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            S(i, j) = compute_overlap(world, bra[i], ket[j]);
        }
    }
    return S;
}

// compute orbital kinetic energy
double compute_orbital_kinetic_energy(World &world,
                                      const real_function_3d &bra,
                                      const real_function_3d &ket) {
    double ke_orb = 0.0;
    for (int axis = 0; axis < 3; axis++) {
        real_derivative_3d D = free_space_derivative<double, 3>(world, axis);
        real_function_3d dbra = D(bra);
        real_function_3d dket = D(ket);
        ke_orb += 0.5 * inner(dbra, dket);
    }
    return ke_orb;
}

// build coulomb operator
real_function_3d build_coulomb(World &world,
                               const std::vector<real_function_3d> &psi,
                               const real_function_3d &target) {
    real_convolution_3d op = CoulombOperator(world, 0.001, 1e-6);
    real_function_3d rho = compute_density(world, psi);
    real_function_3d Jphi = apply(op, rho).truncate() * target;
    return Jphi;
}

//build exchange operator
real_function_3d build_exchange(World &world,
                                const std::vector<real_function_3d> &psi,
                                const real_function_3d &target) {
    real_convolution_3d op = CoulombOperator(world, 0.001, 1e-6);
    real_function_3d Kphi = real_factory_3d(world);
    for (int i = 0; i < psi.size(); i++) {
        Kphi += psi[i] * op(psi[i] * target);
    }
    return Kphi;
}

// build the overall potential :: check(x)
real_function_3d build_Vtot(World &world,
                            const real_function_3d &Vnuc,
                            const std::vector<real_function_3d> &psi,
                            const real_function_3d &phi) {
    real_function_3d Vtot = real_factory_3d(world);
    real_function_3d V_ne_phi = Vnuc * phi;
    V_ne_phi.truncate();
    real_function_3d J_phi = build_coulomb(world, psi, phi);
    real_function_3d K_phi = build_exchange(world, psi, phi);
    Vtot =(V_ne_phi + J_phi - K_phi).truncate();
    return Vtot;
}


// computing fockian element :: check(x)
double compute_fock_element(World &world,
                            const std::vector<real_function_3d> &psi,
                            const real_function_3d &bra,
                            const real_function_3d &ket,
                            const real_function_3d &Vnuc) {
    double kinetic = compute_orbital_kinetic_energy(world, bra, ket);
    real_function_3d Vtot_ket = build_Vtot(world, Vnuc, psi, ket);
    double potential = inner(bra, Vtot_ket);
    return kinetic + potential;
}


// compute orbital coulomb energy
double compute_orbital_coulomb_energy(World &world,
                                      const std::vector<real_function_3d> &psi,
                                      const real_function_3d &bra,
                                      const real_function_3d &ket) {
    real_function_3d Jket = build_coulomb(world, psi, ket);
    double coulomb_en = inner(bra, Jket); // <bra | J | ket>
    return coulomb_en;
}

//compute orbital exchange energy
double compute_orbital_exchange_energy(World &world,
                                       const std::vector<real_function_3d> &psi,
                                       const real_function_3d &bra,
                                       const real_function_3d &ket) {
    real_function_3d Kket = build_exchange(world, psi, ket);
    double exchange_en = inner(bra, Kket); // <bra | K | ket>
    return exchange_en;
}

// compute orbital El-Nuc potential
double compute_orbital_enuc_energy(World &world,
                                   const real_function_3d &Vnuc,
                                   const real_function_3d &bra,
                                   const real_function_3d &ket) {
    real_function_3d Uket = Vnuc * ket;
    Uket.truncate();
    double ENUC_en = inner(bra, Uket);
    return ENUC_en;
}

// compute orbital energy
double compute_orbital_energy(World &world,
                              const std::vector<real_function_3d> &psi,
                              const real_function_3d &bra,
                              const real_function_3d &ket,
                              const real_function_3d &Vnuc) {
    double kinetic = compute_orbital_kinetic_energy(world, bra, ket);
    double Unuc = compute_orbital_enuc_energy(world, Vnuc, bra, ket);
    double coulomb = compute_orbital_coulomb_energy(world, psi, bra, ket);
    double exchange = compute_orbital_exchange_energy(world, psi, bra, ket);
    double total = kinetic + Unuc + coulomb - exchange;
    return total;
}

// compute hamiltonian matrix
Tensor<double> Hmat(World &world,
                    const std::vector<real_function_3d> &psi,
                    const std::vector<real_function_3d> &bra,
                    const std::vector<real_function_3d> &ket,
                    const real_function_3d &Vnuc) {
    int nbra = bra.size();
    int nket = ket.size();
    Tensor<double> H(nbra, nket);
    for (int i = 0; i < nbra; ++i) {
        for (int j = 0; j < nket; ++j) {
            H(i, j) = compute_orbital_energy(world, psi, bra[i], ket[j], Vnuc);
        }
    }
    return H;
}


// diagonalisation
void diagonalise(World &world,
                 Tensor<double> &H,
                 Tensor<double> &S,
                 Tensor<double> &V,
                 Tensor<double> &e) {
    sygv(H, S, 1, V, e);
}

// iterate function
std::vector<double> iterate(World &world,
                            std::vector<real_function_3d> &psi,
                            const real_function_3d &Vnuc,
                            Tensor<double> &eps,
                            XNonlinearSolver<std::vector<real_function_3d>, double, allocator<double, 3>> vectorsolver,
                            const int iter) {
    std::vector<real_function_3d> Vnucpsi = zero_functions_compressed<double, 3>(world, psi.size());
    std::vector<real_function_3d> Jpsi = zero_functions_compressed<double, 3>(world, psi.size());
    std::vector<real_function_3d> Kpsi = zero_functions_compressed<double, 3>(world, psi.size());
    std::vector<real_function_3d> Vpsi(psi.size());
    for (int i = 0; i < psi.size(); ++i) {
        Vnucpsi[i] = Vnuc * psi[i];
        Vnucpsi[i].truncate();
        Jpsi[i] = build_coulomb(world, psi, psi[i]);
        Kpsi[i] = build_exchange(world, psi, psi[i]);
    }

    std::vector<real_function_3d> residual(psi.size());
    std::vector<double> rnorm(psi.size());
    std::vector<double> eps_new(psi.size());
    std::vector<double> norm(psi.size());
    for (int i = 0; i < psi.size(); ++i) {
        real_convolution_3d op3d = BSHOperator<3>(world,
                                                  sqrt(-2 * std::min(-0.05,eps[i])),
                                                  0.001,
                                                  1.e-6);
        Vpsi[i] = Vnucpsi[i] + Jpsi[i] - Kpsi[i];
        Vpsi[i].truncate();
        residual[i] = (psi[i] + 2.0 * op3d(Vpsi[i])).truncate();
        rnorm[i] = residual[i].norm2();
    }
    if (iter > 5) {
        psi = truncate(vectorsolver.update(psi, residual));
    } else {
        for (int j = 0; j < psi.size(); ++j) {
            psi[j] = psi[j] - residual[j];
        }
    }

    for (int i = 0; i < psi.size(); i++) {
        norm[i] = psi[i].norm2();
        psi[i].scale(1.0 / norm[i]);
        eps_new[i] = compute_orbital_energy(world, psi, psi[i], psi[i], Vnuc);
        print("eps_new for orbital: ", i+1, " is ", eps_new[i]);
        if (world.rank() == 0) {
            print(" norm= ", norm[i], " Energy= ", eps_new[i], " err(psi)=", rnorm[i],
                  " err(eps)= ", eps_new[i] - eps[i]);
        }

        eps[i] = eps_new[i];
    }
    return rnorm;
}

// compute total energy
double compute_total_energy(World& world,
                            const std::vector<real_function_3d>& psi,
                            const std::vector<real_function_3d>& bra,
                            const std::vector<real_function_3d>& ket,
                            const real_function_3d& Vnuc){
    double total_energy,kinetic=0.0,exchange=0.0,coulomb=0.0,Unuc=0.0;
    for(int i = 0; i < psi.size(); ++i){
        kinetic += 2.0 * compute_orbital_kinetic_energy(world,bra[i],ket[i]);
        Unuc += 2.0 * compute_orbital_enuc_energy(world,Vnuc,bra[i],ket[i]);
        coulomb +=  compute_orbital_coulomb_energy(world,psi,bra[i],ket[i]);
        exchange += compute_orbital_exchange_energy(world,psi,bra[i],ket[i]);
    }
    total_energy = kinetic + Unuc + coulomb - exchange;
    if(world.rank() == 0){
        print(" Total kinetic energy                        :",kinetic);
        print(" Total Electron-Nuclear potential energy     :",Unuc);
        print(" Total coulomb potential energy              :",coulomb);
        print(" Total exchange potential energy             :",exchange);
        print(" Total energy                                :",total_energy);
    }

    return total_energy;
}

int main(int argc, char **argv) {
    initialize(argc, argv);
    World world(SafeMPI::COMM_WORLD);
    startup(world, argc, argv);

    if (world.rank() == 0) printf("Starting at time %.1f\n", wall_time());

    std::cout.precision(6);
    FunctionDefaults<3>::set_k(8);
    FunctionDefaults<3>::set_thresh(1e-6);
    FunctionDefaults<3>::set_refine(true);
    FunctionDefaults<3>::set_initial_level(5);
    FunctionDefaults<3>::set_truncate_mode(1);
    FunctionDefaults<3>::set_cubic_cell(-L / 2, L / 2);

    std::vector<real_function_3d> phi_old(2);
    phi_old[0] = real_factory_3d(world).f(guess0);
    phi_old[1] = real_factory_3d(world).f(guess1);
    std::vector<double> norm(phi_old.size());
    for (int i = 0; i < norm.size(); ++i) { // normalisation of the funcs
        norm[i] = phi_old[i].norm2();
        phi_old[i].scale(1.0 / norm[i]);
    }
    real_function_3d Vnuc = real_factory_3d(world).f(nuclear_potential);
    Vnuc.truncate();
    Tensor<double> evec, eval;

    //// KAIN Stuff :: must declared inside MPI
    typedef allocator<double, 3> allocT;
    XNonlinearSolver<std::vector<real_function_3d>, double, allocT> vectorsolver(
            allocT(world, phi_old.size())); // can pass a vector of functions
    vectorsolver.set_maxsub(6);

    int maxiter = 50;
    // parent loop
    for (int iter = 0; iter < maxiter; ++iter) {
        print("for iteration : ", iter + 1);

        Tensor<double> H = Hmat(world, phi_old, phi_old, phi_old, Vnuc);
        Tensor<double> S = Smat(world, phi_old, phi_old);

        diagonalise(world, H, S, evec, eval);


        std::vector<real_function_3d> phi_new = zero_functions_compressed<double, 3>(world, phi_old.size());
        for (int i = 0; i < phi_old.size(); ++i) {
            for (int j = 0; j < phi_old.size(); ++j) {
                phi_new[i] += evec(j, i) * phi_old[j];
            }
        }


        std::vector<double> error(phi_new.size());

        error = iterate(world, phi_new, Vnuc, eval, vectorsolver, iter);

        for (int i = 0; i < phi_new.size(); ++i) {
            phi_old[i] = phi_new[i];
        }
        double errorsum =0.0;
        for (int i = 0 ; i < error.size(); ++i){
            errorsum += error[i];
        }
        if (errorsum < 1e-5) break;
    }
    print("\n       Summing up      \n");
    compute_total_energy(world,phi_old,phi_old,phi_old,Vnuc);
    if (world.rank() == 0) printf("Finishing at time %.1f\n", wall_time());

    finalize();

    return 0;
}