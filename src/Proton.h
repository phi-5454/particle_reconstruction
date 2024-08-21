#ifndef PROTON_H
#define PROTON_H
/**
 * @brief A class for a proton that is detected in the Roman Pots.
 * 
 */
class Proton
{
public:
    // Old ntuple parameters
    double Thx; // Scattering angle on the xz-plane
    double Thy; // Scattering angle on the yz-plane
    double pr_px; // X component of momentum
    double pr_py; // Y component of momentum
    double pr_pz = 6500; // Z component of momentum. We assume the momentum change is basically nonexistant compared to the initial momentum.

    // New ntuple parameters
    double pr_ptx; // X component of momentum. Redundant
    double pr_pty; // Y component of momentum. Redundant
    double pr_ptx_sigma; // Error of X component of momentum
    double pr_pty_sigma; // Error of Y component of momentum
    double pr_posx; // Proton vertex X position
    double pr_posy; // Proton vertex Y position
    double pr_posx_sigma; // Error of the vertex X position
    double pr_posy_sigma; // Error of the vertex Y position


    /**
     * @brief Construct a new Proton object with just the scattering angles.
     *
     * @param Thx ThxL/R
     * @param Thy ThyL/R
     */
    Proton(double Thx, double Thy);

    /**
     * @brief Construct a new Proton object.
     * 
     * @param Thx ThxL/R
     * @param Thy ThyL/R
     * @param pr_px pr_px_a/b
     * @param pr_py pr_py_a/b
     * @param pr_pz pr_pz_a/b
     * @param pr_ptx pr_ptx_a/b
     * @param pr_pty pr_pt_a/b
     * @param pr_ptx_sigma pr_ptx_sigma_a/b
     * @param pr_pty_sigma pr_pty_sigma_a/b
     * @param pr_posx pr_posx_a/b
     * @param pr_posy pr_posy_a/b
     * @param pr_posx_sigma pr_posx_sigma_a/b
     * @param pr_posy_sigma pr_posy_sigma_a/b
     */
    Proton(double Thx, double Thy, double pr_px, double pr_py, double pr_pz,
           double pr_ptx, double pr_pty, double pr_ptx_sigma, double pr_pty_sigma,
           double pr_posx, double pr_posy, double pr_posx_sigma, double pr_posy_sigma);

private:
    /**
     * @brief Calculate the x and y components of the momentum based on the Z component and the scattering angles.
     *
     */
    void calculate_momenta();
};
#endif