#include <iostream>
#include <complex>
#include <cmath>
#include <fstream>
#include <string>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <vector>

using cd = std::complex<double>;
using Eigen::MatrixXcd;
using Eigen::VectorXcd;
using T = Eigen::Triplet<cd>;

typedef Eigen::SparseMatrix<std::complex<double>> SpComplexMatrix;

cd pi = cd(3.14159265359, 0);

struct sim{
    int nt;
    int nx;
    cd dt;
    cd dx;
    cd sx;
    cd V0;
    cd hbar;
    cd m;
    cd t;
    cd E;       // Particle's energy
    cd r;       // auxiliary for iteration
};

struct wp{
    cd sigma;
    cd x0;
    cd y0;
    cd kx;
    cd ky;
};

// ======================
// Potentials
// ======================
void potentialBarrier(VectorXcd& A, const sim& p) {
    // Potential Barrier at x = 0
    for (int i = 0; i < p.nx; i++) {
        for (int j = p.nx/2; j < p.nx; j++) {
            A(i*p.nx + j) = p.V0;
        }
    }
}

void potentialSlit(VectorXcd& A, const sim& p) {
    // Potential Barrier at x = 0
    for (int i = 0; i < p.nx; i++) {
        if (i > p.nx/2 + p.nx/16 || i < p.nx/2 - p.nx/16) {
            A(i*p.nx + p.nx/2) = cd(1000, 0);
        }
    }
}

void potentialDoubleSlit(VectorXcd& A, const sim& p) {
    // Potential Barrier at x = 0
    for (int i = 0; i < p.nx; i++) {
        if (i > p.nx/2 + p.nx/20 || i < p.nx/2 - p.nx/20) {
            A(i*p.nx + p.nx/2) = cd(-1000, 0);
        }
    }

    for (int i = 0; i < p.nx; i++) {
        if (i < p.nx/2 + p.nx/(10*3) && i > p.nx/2 - p.nx/(10*3)) {
            A(i*p.nx + p.nx/2) = cd(-1000, 0);
        }
    }
}

void absorbingFrontiers(VectorXcd& A, const sim& p) {
    double k = -1000;
    for (int i = 0; i < p.nx; i++) {
        A(i*p.nx) = cd(0, k);
        A(i*p.nx + p.nx-1) = cd(0, k);
        A(i) = cd(0, k);
        A((p.nx-1)*p.nx + i) = cd(0, k);
    }
}

// ======================
// Sources
// ======================
void wavePacket(VectorXcd& psi, VectorXcd& X, VectorXcd& Y, const sim& p, const wp& w) {
    // Set a normalized Gaussian wavepacket
    int k;
    cd x;
    cd y;
    cd N = cd(std::pow(pi.real()*w.sigma.real()*w.sigma.real(), -1.0/4.0), 0);
    for (int i = 0; i < p.nx; i++) {
        for (int j = 0; j < p.nx; j ++) {
            k = i*p.nx + j;
            x = X(k);
            y = Y(k);
            psi(k) =  N * std::exp(-((x-w.x0)*(x-w.x0) + (y-w.y0)*(y-w.y0)) / (cd(2, 0)*w.sigma*w.sigma)) * 
                          std::exp(cd(0, -1) * (w.kx*(x-w.x0) + w.ky*(y-w.y0)));
        }
     } 
}


// ======================
// Coordinate system
// ======================
// For making the meshgrid
// of the computational domain

void meshGrid(VectorXcd& X, VectorXcd& Y, const sim& p) {
    for (int i = 0; i < p.nx; i++) {
        for (int j = 0; j < p.nx; j++) {
            X(i*p.nx + j) = -p.sx/2.0 + (p.sx/(cd(p.nx, 0) - cd(1, 0)))*cd(j, 0);
            Y(i*p.nx + j) = -p.sx/2.0 + (p.sx/(cd(p.nx, 0) - cd(1, 0)))*cd(i, 0);
        }
    }
}

// ======================
// Auxiliary functions
// ======================
void saveParams2File(std::ofstream& file, const sim& p) {
    file << "Time steps: " << p.nt << std::endl;
    file << "Spatial Resolution: " << p.nx << std::endl;
    file << "dt: " << p.dt.real() << std::endl;
    file << "dx: " << p.dx.real() << std::endl;
    file << "Cell size sx: " << p.sx.real() << std::endl;
    file << "Barrier at x = 0. V0: " << p.V0.real() << std::endl;
    file << "hbar: " << p.hbar.real() << std::endl;
    file << "Mass: " << p.m.real() << std::endl;
    file << "Energy: " << p.E.real() << std::endl;
}


void saveWavefunction2file(std::ofstream& file, VectorXcd& psi, const sim& p) {
    for (int i = 0; i < p.nx; i++) {
        for (int j = 0; j < p.nx; j++) {
            // file << psi1[i*p.nx + j] << "\t";
            file << psi(i*p.nx + j).real() << "\t"<< psi(i*p.nx + j).imag() << "\t";
        }
        file << std::endl;
    }
}


void saveMatrix2file(std::ofstream& file, SpComplexMatrix& psi, const sim& p) {

    MatrixXcd A = MatrixXcd(psi);

    for (int i = 0; i < p.nx * p.nx; i++) {
        for (int j = 0; j < p.nx * p.nx; j++) {
            // file << psi1[i*p.nx + j] << "\t";
            file << A(i, j).real() << "\t"<< A(i, j).imag() << "\t";
        }
        file << std::endl;
    }
}

void spySparse(std::ofstream& file, const SpComplexMatrix& mat, int displaySize) {
    for (int i = 0; i < displaySize; ++i) {
        for (int j = 0; j < displaySize; ++j) {
            // Check if there is a non-zero value at (i, j)
            if (std::abs(mat.coeff(i, j)) > 1e-9) {
                file << "# "; // Non-zero entry
            } else {
                file << ". "; // Zero
            }
        }
        file << std::endl;
    }
}


// ======================
// Simulation's functions
// ======================

void AB_matrices(SpComplexMatrix& A, SpComplexMatrix& B, VectorXcd& V, const sim& p) {
    std::vector<T> A_triplets;
    std::vector<T> B_triplets;
    cd a;
    cd b;
    int k;

    for (int i = 0; i < p.nx; i++) {
        for (int j = 0; j < p.nx; j++) {
            k = i*p.nx + j;

            a = cd(cd(1, 0) + cd(4, 0)*p.r + ((cd(0, 1)*p.dt*V(k))/(cd(2, 0)*p.hbar)));
            b = cd(cd(1, 0) - cd(4, 0)*p.r - ((cd(0, 1)*p.dt*V(k))/(cd(2, 0)*p.hbar)));

            A_triplets.emplace_back(k, k, a);
            B_triplets.emplace_back(k, k, b);

            if (k + p.nx < p.nx * p.nx && i+1 < p.nx) {
                A_triplets.emplace_back(k, k + p.nx, -p.r);
                B_triplets.emplace_back(k, k + p.nx, p.r);
            }

            if (k - p.nx >= 0 && i-1 >= 0) {
                A_triplets.emplace_back(k, k - p.nx, -p.r);
                B_triplets.emplace_back(k, k - p.nx, p.r);
            }

            if (k + 1 < p.nx * p.nx && j + 1 < p.nx) {
                A_triplets.emplace_back(k, k + 1, -p.r);
                B_triplets.emplace_back(k, k + 1, p.r);
            }

            if (k - 1 >= 0 && j - 1 >= 0) {
                A_triplets.emplace_back(k, k - 1, -p.r);
                B_triplets.emplace_back(k, k - 1, p.r);
            }
        }
    }

    A.setFromTriplets(A_triplets.begin(), A_triplets.end());
    B.setFromTriplets(B_triplets.begin(), B_triplets.end());
}


void AB_matrices2(SpComplexMatrix& A, SpComplexMatrix& B, VectorXcd& V, const sim& p) {
    std::vector<T> A_triplets;
    std::vector<T> B_triplets;
    cd a;
    cd b;
    int i;
    int j;

    for (int k = 0; k < p.nx * p.nx; k++) {
        i = k / p.nx;
        j = k % p.nx;

        a = cd(1, 0) + cd(4, 0)*p.r + ((cd(0, 1) * p.dt * V(k)) / (cd(2, 0) * p.hbar));
        b = cd(1, 0) - cd(4, 0)*p.r - ((cd(0, 1) * p.dt * V(k)) / (cd(2, 0) * p.hbar));

        A_triplets.emplace_back(k, k, a);
        B_triplets.emplace_back(k, k, b);

        // y + 1
        if (i != p.nx - 1) {
            A_triplets.emplace_back(k, (i+1)*p.nx + j, -p.r);
            B_triplets.emplace_back(k, (i+1)*p.nx + j, p.r);
        
        }

        // y - 1
        if (i != 0) {
            A_triplets.emplace_back(k, (i-1)*p.nx + j, -p.r);
            B_triplets.emplace_back(k, (i-1)*p.nx + j, p.r);
        
        }

        // x + 1
        if (j != p.nx - 1) {
            A_triplets.emplace_back(k, i*p.nx + j + 1, -p.r);
            B_triplets.emplace_back(k, i*p.nx + j + 1, p.r);
        
        }

        // x - 1
        if (j != 0) {
            A_triplets.emplace_back(k, i*p.nx + j - 1, -p.r);
            B_triplets.emplace_back(k, i*p.nx + j - 1, p.r);
        
        }

    }
    

    A.setFromTriplets(A_triplets.begin(), A_triplets.end());
    B.setFromTriplets(B_triplets.begin(), B_triplets.end());
}

// ======================
// Main function
// ======================
int main() {
    std::string fileP = "Params.txt";
    std::string fileE = "TimeEvolution.txt";
    std::string fileM = "Matrix2.txt";


    // ======================
    // Simulation parameters
    // ======================
    sim s;
    s.nt = 1000;
    s.nx = 160;
    s.dt = cd(0.00125, 0);
    s.sx = cd(10, 0);
    s.V0 = cd(-20, 0);
    s.hbar = cd(1, 0);
    s.m = cd(1, 0);
    s.dx = s.sx/cd(s.nx - 1, 0);
    s.E = 2;
    s.r = (cd(0, 1)*s.hbar*s.dt)/(cd(4, 0)*s.m*s.dx*s.dx);

    wp w;
    w.sigma = cd(0.2, 0);
    w.x0 = cd(0, 0);
    w.y0 = cd(0, 0);
    w.kx = cd(0, 0)*pi;
    w.ky = cd(0, 0)*pi;


    // ======================
    // Simulation Objects
    // ======================
    VectorXcd V(s.nx * s.nx);                   // Potential matrix [row*nx + column]
    SpComplexMatrix A(s.nx*s.nx, s.nx*s.nx);    // A matrix for time iteration
    SpComplexMatrix M(s.nx*s.nx, s.nx*s.nx);    // M matrix for time iteration
    VectorXcd psi0(s.nx * s.nx);                // Wavefunction before
    VectorXcd psi1(s.nx * s.nx);                // Wavefunction after

    VectorXcd X(s.nx * s.nx);                   // x-coordinates
    VectorXcd Y(s.nx * s.nx);                   // y-corrdinates

    // Initializing simulation objects
    V.setZero();
    A.setZero();
    M.setZero();
    psi0.setZero();
    psi1.setZero();
    X.setZero();
    Y.setZero();

    // Select Potential
    // potentialBarrier(V, s);
    // potentialDoubleSlit(V, s);
    // absorbingFrontiers(V, s);
    meshGrid(X, Y, s);

    // Computing A and B
    // AB_matrices(A, M, V, s);
    AB_matrices2(A, M, V, s);
    
    // Files for saving data
    std::ofstream paramsFile(fileP);
    std::ofstream evolutionFile(fileE);
    std::ofstream matrixFile(fileM);

    saveParams2File(paramsFile, s);
    // spySparse(matrixFile, A, s.nx*s.nx);

    // Time evolution loop
    Eigen::SparseLU<SpComplexMatrix> solver;
    solver.compute(A);

    if (solver.info() != Eigen::Success) {
        std::cout << "Decomposition Failed";
    } else {

        s.t = cd(0, 0);                 // initial time
        wavePacket(psi0, X, Y, s, w);

        for (int i = 0; i < s.nt; i++) {
            std::cout << 1.0*i / (1.0*s.nt) << std::endl;

            saveWavefunction2file(evolutionFile, psi0, s);

            psi1 = solver.solve(M * psi0);

            psi0 = psi1;
            s.t = s.t + cd(s.dt.real(), 0);
        }
    }

    paramsFile.close();
    evolutionFile.close();
    return 0;
}
