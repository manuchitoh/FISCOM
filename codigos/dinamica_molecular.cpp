#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <fstream>

// Parámetros de la simulación
const int Nc = 30; // Número de celdas
const int N = 900; // Número de partículas
const double L = 55.0; // Longitud de la caja
const double h = 0.005; // Paso de tiempo
const int steps = 1; // Número de pasos

// Variables de cada partícula
struct Atom {
    double x, y; // Posiciones
    double vx, vy; // Velocidades
    double fx, fy; // Fuerzas
};

// Potencial de Lennard-Jones
// Devuelve el valor del potencial y actualiza la fuerza
double potential(double r2, double& force) {
    double rm2 = 1.0 / r2;
    double rm6 = rm2*rm2*rm2;
    double rm12 = rm6*rm6;
    force = 24*(2*rm12 - rm6)*rm2;
    return 4*(rm12 - rm6);
}

// Inicializo las partículas:
// * posiciones iniciales (na,ma), con a=L/(Nc+1) y n,m=1,...,Nc
// * velocidades iniciales (+1.1, 0) o (-1.1, 0) elegidos al azar
void initializeAtoms(std::vector<Atom>& atoms) {
    double a = L/(Nc+1);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(0, 1);
    for (int i = 0; i < N; i++) {
        atoms[i].x = a*(1 + (i%Nc));
        atoms[i].y = a*(1 + (i/Nc));
        atoms[i].vx = dis(gen) < 0.5 ? 1.1 : -1.1;
        atoms[i].vy = 0;
    }
}
// Condición inicial de red triangular que ocupe toda la caja
void initializeAtoms_triangular(std::vector<Atom>& atoms) {
    double a = 1.122462048;
    double h = a*sqrt(3)/2;
    double sen = sqrt(3)/2;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(0, 1);
    for (int i = 0; i < N; i++) {
        atoms[i].x = a*(1 + (i%Nc)) + (i/Nc)%2 * a/2 + 10.38;
        atoms[i].y = h*(1 + (i/Nc)) + 12.92;
        atoms[i].vx = dis(gen) < 0.5 ? 0.5*0.5 : -0.5*0.5;
        atoms[i].vy = dis(gen) < 0.5 ? 0.5*sen : -0.5*sen;
    }
}

// Obtener las fuerzas y la energía potencial
double getForces(std::vector<Atom>& atoms) {
    double potentialEnergy = 0.0;
    for (auto& at : atoms) at.fx = at.fy = 0.0;

    for (size_t i = 0; i < atoms.size(); ++i) {
        for (size_t j = i + 1; j < atoms.size(); ++j) {
            double dx = atoms[i].x - atoms[j].x;
            double dy = atoms[i].y - atoms[j].y;

            // Imagen mínima para condiciones de contorno periódicas
            dx -= L * round(dx/L);
            dy -= L * round(dy/L);

            // PREGUNTAR SOBRE S(r) para contnuidad
            double r2 = dx*dx + dy*dy;
            if (r2 < 6.25) {    // 2.5^2 = 6.25
                double force;
                potentialEnergy += potential(r2, force) + 0.016363168911; // Corrección de corte
                double fx = force*dx;
                double fy = force*dy;
                atoms[i].fx += fx;
                atoms[i].fy += fy;
                atoms[j].fx -= fx;
                atoms[j].fy -= fy;
            }
        }
    }
    return potentialEnergy;
}

// Método Verlet de velocidades
double velocityVerlet(std::vector<Atom>& atoms) {
    for (auto& at : atoms) {
        // Posición
        at.x += at.vx * h + 0.5 * at.fx * h*h;
        at.y += at.vy * h + 0.5 * at.fy * h*h;
        // Condiciones periódicas
        at.x -= L * floor(at.x/L);
        at.y -= L * floor(at.y/L);
    }

    // Guardo las fuerzas anteriores y calculo las nuevas
    std::vector<double> fx_old(N), fy_old(N);
    for (size_t i = 0; i < atoms.size(); ++i) {
        fx_old[i] = atoms[i].fx;
        fy_old[i] = atoms[i].fy;
    }
    
    double potentialEnergy = getForces(atoms);

    for (size_t i = 0; i < atoms.size(); ++i) {
        // Velocidad
        atoms[i].vx += 0.5 * (atoms[i].fx + fx_old[i]) * h;
        atoms[i].vy += 0.5 * (atoms[i].fy + fy_old[i]) * h;
    }
    return potentialEnergy;
}

// Obtener la energía cinética
double getKineticEnergy(const std::vector<Atom>& atoms) {
    double kineticEnergy = 0.0;
    for (const auto& at : atoms) {
        kineticEnergy += 0.5 * (at.vx * at.vx + at.vy * at.vy);
    }
    return kineticEnergy;
}


int main() {
    std::vector<Atom> atoms(N);
    initializeAtoms_triangular(atoms);
    getForces(atoms);

    // Abro los archivos
    std::ofstream outfileE("energia_dinamica_molecular_triang.dat");
    std::ofstream outfileV("velocidades_triang.dat");
    std::ofstream outfileR("posiciones_triang.dat");

    if (!outfileE || !outfileV || !outfileR) {
        std::cerr << "Error al abrir los archivos." << std::endl;
        return 1;
    }
    outfileE << "Paso\tEnergía Total\tEnergía Cinética\tEnergía Potencial\n";

    for (int i = 0; i < steps; i++) {
        std::cout << "Iteración " << i << "..." << std::endl;

        double potentialEnergy = velocityVerlet(atoms);
        double kineticEnergy = getKineticEnergy(atoms);
        double totalEnergy = kineticEnergy + potentialEnergy;

        outfileE << i << "\t" << totalEnergy << "\t" << kineticEnergy << "\t" << potentialEnergy << "\n";

        if (i == 0 || i == 200 || i == steps-1) {
            outfileV << "Paso " << i << ":\n";
            for (auto& at : atoms) {
                outfileV << at.vx << "\t" << at.vy << "\n";
            }
        }    

        outfileR << "Paso " << i << ":\n";
        for (auto& at : atoms) {
            outfileR << at.x << "\t" << at.y << "\n";
        }
    }

    // Cierro los archivos
    outfileE.close();
    outfileV.close();
    outfileR.close();

    return 0;
}
