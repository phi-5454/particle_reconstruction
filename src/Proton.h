#ifndef PROTON_H
#define PROTON_H
/**
 * @brief A class for a proton that is detected in the Roman Pots,
 * 
 */
class Proton
{
public:
    const float p = 6500; // Momentum of the photon. We assume the change is basically nonexistant and stays the same after the collition
    float Thx; // Scattering angle on the xz-plane
    float Thy; // Scattering angle on the yz-plane
    float px; // X component of momentum
    float py; // Y component of momentum

    /**
     * @brief Construct a new Proton object
     * 
     * @param Thx ThxL/R
     * @param Thy ThyL/R
     */
    Proton(float Thx, float Thy);

    /**
     * @brief Creates a proton and calculates its momenta
     * 
     * @param Thx ThxL/R
     * @param Thy ThyL/R
     */
    void initialize_proton(float Thx, float Thy);

private:
    /**
     * @brief Calculate the x and y components of the momentum
     * 
     */
    void calculate_momenta();
};
#endif