#include <iostream>
#include <cmath>
#include <vector>

constexpr int size = 10; // Size of the flow data
constexpr double dt = 1.0e-5;

class Flow{
public:
    double gamma = 1.4; // Specific heat ratio for air
    double R = 287.05; // Specific gas constant for air in J/(kg·K)
    double c_v = 717.0; // Specific heat at constant volume for air in J/(kg·K)
    std::vector<double> rho;
    std::vector<double> rhoU;
    std::vector<double> energyDensity;
    Flow(int);
    double pressure(double rho, double rhoU, double energyDensity);
    double velocity(double rhoU, double rho);
    double soundSpeed(double rho, double rhoU, double energyDensity);
    double temperature(double energyDensity, double rhoU, double rho);
    double totalEnergy(double rho, double rhoU, double energyDensity);
    double totalPressure(double rho, double rhoU, double energyDensity);
    double calculateEnergyDensity(double pressure, double rho, double rhoU);
};

double Flow::pressure(double rho, double rhoU, double energyDensity) {
    return (gamma - 1.0) * (energyDensity - 0.5 * rhoU * rhoU / rho);
}

double Flow::velocity(double rhoU, double rho) {
    return rhoU / rho;
}

double Flow::soundSpeed(double rho, double rhoU, double energyDensity) {
    return std::sqrt(gamma * pressure(rho, rhoU, energyDensity) / rho);
}

double Flow::temperature(double energyDensity, double rhoU, double rho) {
    return (energyDensity - 0.5 * rhoU * rhoU / rho) / (c_v * rho);
}

double Flow::totalEnergy(double rho, double rhoU, double energyDensity) {
    return energyDensity + 0.5 * rhoU * rhoU / rho;
}

double Flow::totalPressure(double rho, double rhoU, double energyDensity) {
    return pressure(rho, rhoU, energyDensity) + 0.5 * rhoU * rhoU / rho;
}

double Flow::calculateEnergyDensity(double pressure, double rho, double rhoU) {
    return pressure / (gamma - 1.0) + 0.5 * rhoU * rhoU / rho;
}

Flow::Flow(int size) : rho(size), rhoU(size), energyDensity(size) {
    for (int i = 0; i < size; ++i) {
        rho[i] = 1.29;
        rhoU[i] = 0.0;
        energyDensity[i] = calculateEnergyDensity(101325.0, rho[i], 0.0); // 101325 Pa
        if (i > size / 2) {
            rho[i] = 0.7 * 1.29;
            energyDensity[i] = calculateEnergyDensity(0.7 * 101325.0, rho[i], 0.0);
        }
    }
}

int main() {

    Flow flow(size);

    for(int i = 0; i < size; ++i) {
        std::cout << "Cell " << i << ": "
                  << "Density = " << flow.rho[i] << ", "
                  << "EnergyDensity = " << flow.energyDensity[i] << ", "
                  << "Pressure = " << flow.pressure(flow.rho[i], flow.rhoU[i], flow.energyDensity[i]) << ", "
                  << "Velocity = " << flow.velocity(flow.rhoU[i], flow.rho[i]) << ", "
                  << "Sound Speed = " << flow.soundSpeed(flow.rho[i], flow.rhoU[i], flow.energyDensity[i]) << ", "
                  << "Temperature = " << flow.temperature(flow.energyDensity[i], flow.rhoU[i], flow.rho[i]) << ", "
                  << "Total Energy = " << flow.totalEnergy(flow.rho[i], flow.rhoU[i], flow.energyDensity[i]) << ", "
                  << "Total Pressure = " << flow.totalPressure(flow.rho[i], flow.rhoU[i], flow.energyDensity[i]) 
                  << std::endl;
    }
    
    std::cout << "Hello, World!" << std::endl;
    return 0;
}