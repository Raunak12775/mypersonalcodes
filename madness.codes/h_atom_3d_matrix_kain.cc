//
// Created by Raunak Farhaz on 27.06.22.
//


#include <madness/mra/mra.h>
#include <madness/mra/vmra.h>
#include <madness/mra/operator.h>
#include <madness/mra/nonlinsol.h>

using namespace madness;
const double L = 30.0;

//// we would not require shift because the energies are bound
//static double shift = 17.0;


// The default constructor for functions does not initialize
// them to any value, but the solver needs functions initialized
// to zero for which we also need the world object.
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


// struct for g = guess functions
struct g {

    int i;

    g(const int i) : i(i) {}

    double operator()(const coord_3d &coord) const {
        double x = coord[0];
        double y = coord[1];
        double z = coord[2];
        double rr = sqrt(x*x + y*y + z*z);
        return rr*exp(-(fabs(rr)/i)+1e-8);
        //return pow((r[0]), i) * exp(-(r[0] * r[0]));
//        return (pow(x, i) * exp(-fabs(x))) + (pow(y, i) * exp(-fabs(y))) + (pow(z, i) * exp(-fabs(z)));
    }
};

//r_squared = pow((atom1.coords[0] - atom2.coords[0]), 2.0) + pow((atom1.coords[1] - atom2.coords[1]), 2.0) +
//            pow((atom1.coords[2] - atom2.coords[2]), 2.0);
//length = pow(r_squared, 0.5);

//static double guess(const coord_3d &r) {
//    double x = r[0];
//    double y = r[1];
//    double z = r[2];
//    return exp(-sqrt(x*x + y*y + z*z + 1e-6));
////    return (x + y + z) * exp(-fabs(x+y+z));// + (y * exp(-y)) + (z * exp(-z));
//}

// potential
static double Potential(const coord_3d &r) {

    double x = r[0];
    double y = r[1];
    double z = r[2];
//    double rr =
    return -1.0/sqrt(x*x + y*y + z*z+ 1e-8);
//    return -(1.0 / (fabs(x)+fabs(y)+fabs(z) + 1e-8)); //- (1.0 / (fabs(y) + 1e-6)) - (1.0 / (fabs(z) + 1e-6));
}


// Convenience routine for plotting
void plot(const char *filename, const real_function_3d &f) {

    coord_3d lo(0.0), hi(0.0);
    lo[0] = -L;
    hi[0] = L;
    plot_line(filename, 401, lo, hi, f);
//    plot();
}


// computing energy function
//double compute_energy1d(World &world, const real_function_1d &psi, const real_function_1d &V) {
//
//    double kinetic_energy = 0.0;
////    for (int axis = 0; axis < 3; axis++) {
//        real_derivative_1d D = free_space_derivative<double, 1>(world, 0);
//        real_function_1d dpsi = D(psi);
//        kinetic_energy += 0.5 * inner(dpsi, dpsi);
////    }
//    double potential_energy = inner(psi, V * psi);
//    return kinetic_energy + potential_energy;
//}

double compute_energy(World &world, const real_function_3d &psi, const real_function_3d &V) {

    double kinetic_energy = 0.0;
    for (int axis = 0; axis < 3; axis++) {
    real_derivative_3d D = free_space_derivative<double, 3>(world, axis);
    real_function_3d dpsi = D(psi);
    kinetic_energy += 0.5 * inner(dpsi, dpsi);
    }
    double potential_energy = inner(psi, V * psi);
    return kinetic_energy + potential_energy;
}

// computing overlap matrix element
double compute_S_matrixelement(World &world, const real_function_3d &bra, const real_function_3d &ket) {
    double s = inner(bra, ket);
    return s;
}

// computing overlap matrix
Tensor<double> Smat(World &world, const std::vector<real_function_3d> &bra, const std::vector<real_function_3d> &ket) {
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
double compute_H_matrixelement(World &world, const real_function_3d &bra, const real_function_3d &ket,
                               const real_function_3d &V) {
    double h = inner(bra, V * ket);
    for (int axis = 0; axis < 3; axis++) {
        real_derivative_3d D = free_space_derivative<double, 3>(world, axis);
        real_function_3d dbra = D(bra);
        real_function_3d dket = D(ket);
        h += 0.5 * inner(dbra, dket);
    }
    return h;
}

// computing Hamiltonian matrix
Tensor<double> Hmat(World &world, const std::vector<real_function_3d> &bra, const std::vector<real_function_3d> &ket,
                    const real_function_3d &V) {
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

//double iterate(World &world,
//               const real_function_3d &V,
//               real_function_3d &psi,
//               double &eps, NonlinearSolverND<3> &solver, int iter) {
//
//    real_convolution_3d op1d = BSHOperator<3>(world, sqrt(-2 * eps), 0.0001, 1.e-6);
//    real_function_3d Vpsi = (V * psi);
//    Vpsi.truncate();
//    real_function_3d r = psi + 2.0 * op1d(Vpsi);
//    double rnorm = r.norm2();
//    if (iter > 5){
//    psi = solver.update(psi, r);}
//    else{
//    psi = r -psi;}
//    double norm = psi.norm2();
//    psi.scale(1.0 / norm);
//    double eps_new = compute_energy(world, psi, V);
//    if (world.rank() == 0) {
//        print("norm=", norm, " eps=", eps, " err(psi)=", rnorm, " err(eps)=", eps_new - eps);
//    }
//    eps = eps_new;
//    return rnorm;
//}


std::vector<double> apply_BSH(World &world,
                              const real_function_3d &V,
                              std::vector<real_function_3d> &psi,
                              Tensor<double> &eps,
                              XNonlinearSolver<std::vector<real_function_3d>, double, allocator<double, 3>> vectorsolver/*NonlinearSolverND<1> &solver*/,
                              const int iter) {

    std::vector<real_function_3d> Vpsi(psi.size());
    std::vector<real_function_3d> r(psi.size());
    std::vector<double> rnorm(psi.size());
    std::vector<double> eps_new(psi.size());
    std::vector<double> norm(psi.size());
    for (int i = 0; i < psi.size(); ++i) {
        real_convolution_3d op3d = BSHOperator<3>(world, sqrt(-2 * eps[i]), 0.0001, 1.e-6);
//    real_convolution_1d op1d = CoulombOperator<1>(world, 0.0001,1.e-6);

        Vpsi[i] = (V * psi[i].truncate());
        Vpsi[i].truncate();
        r[i] = psi[i] + 2.0 * op3d(Vpsi[i]);
        rnorm[i] = r[i].norm2();
    }

    if (iter > 3) {
        psi = vectorsolver.update(psi, r); //// tip : r can also be passed as a vector
    } else {
        for (int j = 0; j < psi.size(); ++j)
            psi[j] = r[j] - psi[j];
    }
    for (int i = 0; i < psi.size(); ++i) {
        norm[i] = psi[i].norm2();
        psi[i].scale(1.0 / norm[i]);
        eps_new[i] = compute_energy(world, psi[i], V);
        if (world.rank() == 0) {
            print("for state : ", i, "norm= ", norm[i], " Energy= ", eps_new[i], " err(psi)=", rnorm[i],
                  " err(eps)=", eps_new[i] - eps[i]);
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
    std::cout.precision(6);

    FunctionDefaults<3>::set_k(6);
    FunctionDefaults<3>::set_thresh(1e-6);
    FunctionDefaults<3>::set_refine(true);
    FunctionDefaults<3>::set_initial_level(5);
    FunctionDefaults<3>::set_truncate_mode(1);
    FunctionDefaults<3>::set_cubic_cell(-L, L);

//    real_function_3d psi = real_factory_3d(world).f(guess); // guess
//    real_function_3d Vnuc = real_factory_3d(world).f(Potential); // potential
//    psi.truncate();
//    double norm = psi.norm2();
//    psi.scale(1.0/norm);
//    double eps = -0.3;
//    NonlinearSolverND<3> solver;
//    double err = 0.0;
//    for (int iter = 0; iter < 50; iter++){
//        err = iterate(world,Vnuc,psi,eps,solver,iter);
//        if (err < 5e-6) break;
//    }
//
//    print(" energy is : ", eps);



    //    MRA representation of the guess functions in a vector
    long n = 3; //// n = number of guess functions
    std::vector<real_function_3d> guess_funs(n);
//    for (int i = 0; i < n; ++i) {
//        guess_funs[i] = real_factory_3d(world).functor(g(i+1));
//    }
    guess_funs[0] = real_factory_3d(world).functor(g(1));
    guess_funs[1] = real_factory_3d(world).functor(g(2));
    guess_funs[2] = real_factory_3d(world).functor(g(3));
//    guess_funs[3] = real_factory_3d(world).functor(g(4));
//    guess_funs[4] = real_factory_3d(world).functor(g(5));
//    for (int i = 1; i <= n; ++i){
//        guess_funs.emplace_back(real_factory_3d(world).functor(g(i)));
//    }
    char file[256];
    real_function_3d v = real_factory_3d(world).f(Potential); // potential
    Tensor<double> evec, eval;
    typedef allocator<double,3> allocT;
    XNonlinearSolver<std::vector<real_function_3d>, double, allocT> vectorsolver(allocT(world, guess_funs.size())); // can pass a vector of functions
    vectorsolver.set_maxsub(6);
    //vectorsolver.update(guess_funs,guess_funs);

    // zero_functions_compressed() is used for declaring vector of real_1D_functions
    // creates an empty psi_new vector of real_functions
    // std::vector<real_function_1d> psi_new = zero_functions_compressed<double, 1>(world, guess_funs.size());


    int maxiter = 50;
    //// Parent loop construction
    for (int iter = 0; iter < maxiter; ++iter) {

        Tensor<double> H = Hmat(world, guess_funs, guess_funs, v); // constructed H (1)
        Tensor<double> S = Smat(world, guess_funs, guess_funs); // constructed S (2)

        // diagonalise H and get evec(i,j) and eval(i) from H.evec=eval.S.evec (3)
        diagolisation(world, H, S, evec, eval);
        // Call a tensor by indices => D(i,j)
        //print("eigenvalues : \n",eval);

        //print("eigenvector : \n",evec);

        std::vector<real_function_3d> psi_new = zero_functions_compressed<double, 3>(world, guess_funs.size());
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

        error = apply_BSH(world, v, psi_new, eval, vectorsolver, iter);
        print("================================================================================");

        for (int i = 0; i < guess_funs.size(); ++i) {
            guess_funs[i] = psi_new[i]; // old guess function vector is modified to new (5)
        }

        if (error[0] < 5e-4) break;

    }

    print("New eigenvalues : ");
    for (int i = 0; i < eval.size(); ++i) {
        print("Energy of state : ", i, " ", eval[i]);
    }

//    for (int i = 0; i < n ; ++i){
//        sprintf(file,"psi-%3.3d.dat", i);
//        plot(file,guess_funs[i]);
//    }

    if (world.rank() == 0) printf("finishing at time %.1f\n", wall_time());

    finalize();

    return 0;
}
