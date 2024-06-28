#ifndef PROTON_H
#define PROTON_H
/**
 * @brief A class for a proton that is detected in the Roman Pots,
 * 
 */
class Proton
{
public:
    const double p = 6500; // Momentum of the photon. We assume the change is basically nonexistant and stays the same after the collition
    double Thx; // Scattering angle on the xz-plane
    double Thy; // Scattering angle on the yz-plane
    double px; // X component of momentum
    double py; // Y component of momentum

    /**
     * @brief Construct a new Proton object
     * 
     * @param Thx ThxL/R
     * @param Thy ThyL/R
     */
    Proton(double Thx, double Thy);

private:
    /**
     * @brief Calculate the x and y components of the momentum
     * 
     */
    void calculate_momenta();
};
#endif