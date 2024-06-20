/**
 * @brief Class for an individual particle
 * 
 */
#ifndef PARTICLE_H
#define PARTICLE_H
class Particle
{
public:
    float p; // Momentum
    float pt; // Transverse (= perpendicular to the beam) momentum
    float px; // X component of momentum, towards the center of the ring
    float py; // Y component of momentum, "upwards"
    float pz; // Z component of momentum, parallel to the beam
    float eta; // Pseudorapidity
    float theta; // Scattering angle, from z-axis.
    float phi; // Azimuthal angle, from positive x-axis
    int q; // Charge
    float dxy; // Particle's track's smallest deviation in xy-plane from the primary vertex 
    float dz; // Particle's track's smallest deviation along z-axis from the primary vertex
    float mass; // Expected mass of the particle

    /**
     * @brief Construct a new Particle object
     * 
     */
    Particle();

    /**
     * @brief Calculate the different components of momentum
     * 
     */
    void calculate_momenta();

    /**
     * @brief Calculate the scattering angle of the particle
     * 
     */
    void calculate_theta();

};
#endif