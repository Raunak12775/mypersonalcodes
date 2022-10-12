//
// Created by Raunak Farhaz on 27.06.22.
//


#include <madness/mra/mra.h>
#include <madness/mra/vmra.h>
#include <madness/mra/operator.h>
#include <madness/mra/nonlinsol.h>

using namespace madness;
const double L = 4.0;
static double shift = 17.0;


// The default constructor for functions does not initialize
// them to any value, but the solver needs functions initialized
// to zero for which we also need the world object.
template<typename T, std::size_t NDIM>
struct allocator {
    World& world;
    const int n;

    /// @param[in]	world	the world
    /// @param[in]	nn		the number of functions in a given vector
    allocator(World& world, const int nn) :
            world(world), n(nn) {
    }

    /// allocate a vector of n empty functions
    std::vector<Function<T, NDIM> > operator()() {
        return zero_functions<T, NDIM>(world, n);
    }
};


// struct for g = guess functions
struct g {

    int i;

    g(const int i) : i(i) {}

    double operator()(const coord_1d &r) const {
        return pow((r[0]), i) * exp(-(r[0] * r[0]));
    }
};


// potential
static double Harmonic_Potential(const coord_1d &r) {

    double x = r[0];
    // double well potential is of the form of V(x) = 1/2*x^2 + \alpha*Exp(-\beta*x^2) 
    // where I chose \alpha = 1.0 and \beta = 3.0
    // refer : /Users/raunakfarhaz/Desktop/madness-codes/HarmonicGaussianDoubleWellPotential-author.nb
    // the results are correct upto 2 decimal places	
    return 0.5 * (x * x)  + 1.0*exp(-3.0 * x * x) - shift; // this is for Double Well potential
}


// Convenience routine for plotting
void plot(const char *filename, const real_function_1d &f) {

    coord_1d lo(0.0), hi(0.0);
    lo[0] = -L;
    hi[0] = L;
    plot_line(filename, 401, lo, hi, f);
}


// computing energy function
double compute_energy(World &world, const real_function_1d &psi, const real_function_1d &V) {

    double kinetic_energy = 0.0;
    real_derivative_1d D = free_space_derivative<double, 1>(world, 0);
    real_function_1d dpsi = D(psi);
    kinetic_energy += 0.5 * inner(dpsi, dpsi);
    double potential_energy = inner(psi, V * psi);
    return kinetic_energy + potential_energy;
}

// computing overlap matrix element
double compute_S_matrixelement(World &world, const real_function_1d &bra, const real_function_1d &ket) {
    double s = inner(bra, ket);
    return s;
}

// computing overlap matrix
Tensor<double> Smat(World &world, const std::vector<real_function_1d> &bra, const std::vector<real_function_1d> &ket) {
    int n = bra.size();
    Tensor<double> S(n, n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            S(i, j) = compute_S_matrixelement(world, bra[i], ket[j]);
        }
    }
    return S;
}


// computing Hamiltonian matrix element
double compute_H_matrixelement(World &world, const real_function_1d &bra, const real_function_1d &ket,
                               const real_function_1d &V) {
    double h = inner(bra, V * ket);
    real_derivative_1d D = free_space_derivative<double, 1>(world, 0);
    real_function_1d dbra = D(bra);
    real_function_1d dket = D(ket);
    h += 0.5 * inner(dbra, dket);
    return h;
}

// computing Hamiltonian matrix
Tensor<double> Hmat(World &world, const std::vector<real_function_1d> &bra, const std::vector<real_function_1d> &ket,
                    const real_function_1d &V) {
    int nbra = bra.size();
    int nket = ket.size();
    Tensor<double> H(nbra, nket);
    for (int i = 0; i < nbra; ++i) {
        for (int j = 0; j < nket; ++j) {
            H(i, j) = compute_H_matrixelement(world, bra[i], ket[j], V);
        }
    }
    return H;
}

// diagonalising the H matrix
void diagolisation(World &world, Tensor<double> &H, Tensor<double> &S, Tensor<double> &V, Tensor<double> &e) {

    //// H = Hamiltonian matrix
    //// S = Overlap Matrix
    //// itype = 1 which means it solves the system like => A*x = (lambda)*B*x
    //// V = the eigenvector
    //// e = all the eigenvalues
    sygv(H, S, 1, V, e);
}


std::vector<double> apply_BSH(World &world,
                              const real_function_1d &V,
                              std::vector<real_function_1d> &psi,
                              Tensor<double> &eps,
                              XNonlinearSolver<std::vector<real_function_1d>,double,allocator<double,1>> vectorsolver/*NonlinearSolverND<1> &solver*/,
                              const int iter) {

    std::vector<real_function_1d> Vpsi(psi.size());
    std::vector<real_function_1d> r(psi.size());
    std::vector<double> rnorm(psi.size());
    std::vector<double> eps_new(psi.size());
    std::vector<double> norm(psi.size());
    for (int i = 0; i < psi.size(); ++i) {
        real_convolution_1d op1d = BSHOperator<1>(world, sqrt(-2 * eps[i]), 0.0001, 1.e-6);
//    real_convolution_1d op1d = CoulombOperator<1>(world, 0.0001,1.e-6);

//    real_function_1d
        Vpsi[i] = (V * psi[i]);
        Vpsi[i].truncate();
//    real_function_1d
        r[i] = psi[i] + 2.0 * op1d(Vpsi[i]);
//    double
        rnorm[i] = r[i].norm2();
    }

    if (iter > 3) {
        psi = vectorsolver.update(psi, r); //// tip : r can also be passed as a vector

    } else {
        for (int j = 0; j < psi.size(); ++j)
            psi[j] = r[j] - psi[j];
    }
    //    double
    for (int i = 0; i < psi.size(); ++i) {
        norm[i] = psi[i].norm2();
        psi[i].scale(1.0 / norm[i]);
//    double
        eps_new[i] = compute_energy(world, psi[i], V);
    if (world.rank() == 0) {
        print("for state : ", i, "norm= ", norm[i], " eps(no shift)= ", eps[i] - shift, " err(psi)=", rnorm[i],
              " err(eps)=", eps_new[i] - eps[i],
              " Energy=", eps_new[i] + shift);
    }
    eps[i] = eps_new[i];
}

return rnorm;

}

int main(int argc, char **argv) {

    //initialising the startup things from the runtime layer (copy paste from old code : ref :: ~/devel/madness/src/example/h2.cc)
    initialize(argc, argv);
    World world(SafeMPI::COMM_WORLD);

    startup(world, argc, argv);
    if (world.rank() == 0) printf("starting at time %.1f\n", wall_time());
    //Tensor<double> &V, &e;

    std::cout.precision(6);
    // ##############################################################################

    // within init and finalise you will write what you want to do : marked the boundaries with ###

    FunctionDefaults<1>::set_k(16);
    FunctionDefaults<1>::set_thresh(1e-6);
    FunctionDefaults<1>::set_refine(true);
    FunctionDefaults<1>::set_initial_level(5);
    FunctionDefaults<1>::set_truncate_mode(1);
    FunctionDefaults<1>::set_cubic_cell(-L, L);

    //    MRA representation of the guess functions in a vector
    long n = 5; //// n = number of guess functions
    std::vector<real_function_1d> guess_funs(n);
//    Tensor<real_function_1d> guess_funs;
    for (int i = 0; i < n; ++i) {
        guess_funs[i] = real_factory_1d(world).functor(g(i));
    }


    real_function_1d v = real_factory_1d(world).f(Harmonic_Potential); // potential
    Tensor<double> evec, eval;
//    std::vector<NonlinearSolverND<1>> solver(n);
//    for (int i=0; i<n; ++i) solver[i].set_maxsub(3);

//    allocT(world, nemo.size()
    typedef allocator<double,1> allocT;
    XNonlinearSolver<std::vector<real_function_1d>, double, allocT> vectorsolver(allocT(world, guess_funs.size())); // can pass a vector of functions
    vectorsolver.set_maxsub(15);
    //vectorsolver.update(guess_funs,guess_funs);

//    NonlinearSolverND<1> solver;
//    std::vector<double> error(guess_funs.size());
//    std::vector<NonlinearSolverND<1>> solver(n);// solver1, solver2;
    // zero_functions_compressed() is used for declaring vector of real_1D_functions
    // creates an empty psi_new vector of real_functions
    // std::vector<real_function_1d> psi_new = zero_functions_compressed<double, 1>(world, guess_funs.size());


    int maxiter = 50;
//    int count = 0;
    //// Parent loop construction
    for (int iter = 0; iter < maxiter; ++iter) {

        Tensor<double> H = Hmat(world, guess_funs, guess_funs, v); // constructed H (1)
        Tensor<double> S = Smat(world, guess_funs, guess_funs); // constructed S (2)
//        print("H matrix: ");
//        print(H);
//        print("S matrix: ");
//        print(S);

        // diagonalise H and get evec(i,j) and eval(i) from H.evec=eval.S.evec (3)
        diagolisation(world, H, S, evec, eval);
        // Call a tensor by indices => D(i,j)
//        print("old eigenvalues are : ");
//        print(eval);
        //      print("eigen vector: ");
        //    print(evec);
////        print("eigen value: ");
        //  print(eval);

        std::vector<real_function_1d> psi_new = zero_functions_compressed<double, 1>(world, guess_funs.size());
        // new guess vector is constructed (4)
        for (int i = 0; i < guess_funs.size(); i++) {
            for (int j = 0; j < guess_funs.size(); ++j) {
                psi_new[i] += evec(j, i) * guess_funs[j];
            }
        }

        for (int i = 0; i < psi_new.size(); ++i) {
            double n = psi_new[i].norm2();
            psi_new[i].scale(1 / n);
        }


        // error vector for 6 guess functions created (6)
        std::vector<double> error(guess_funs.size());
        // BSH applied on the vectors (7)
        print("For iteration : ", iter + 1);

//    for (int i = 0; i < guess_funs.size(); ++i) {
        error = apply_BSH(world, v, psi_new, eval, vectorsolver, iter);
//        }
        print("================================================================================");

        for (int i = 0; i < guess_funs.size(); ++i) {
            guess_funs[i] = psi_new[i]; // old guess function vector is modified to new (5)
        }

        if (error[n-1] < 5e-4) break;

    }

    print("New eigenvalues : ");
    for (int i = 0; i < eval.size(); ++i) {
        print("Energy of state : ", i, " ", eval[i] + shift);
    }

    if (world.rank() == 0) printf("finishing at time %.1f\n", wall_time());

    finalize();

    return 0;
}
