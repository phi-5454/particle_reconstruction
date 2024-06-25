#ifndef PARTICLE_H
#define PARTICLE_H
/**
 * @brief A class for a particle that's detected in the CMS
 *
 */
class Particle {
public:
  float p;     // Momentum
  float pt;    // Transverse (= perpendicular to the beam) momentum
  float px;    // X component of momentum, towards the center of the ring
  float py;    // Y component of momentum, "upwards"
  float pz;    // Z component of momentum, parallel to the beam
  float eta;   // Pseudorapidity
  float theta; // Scattering angle, from z-axis.
  float phi;   // Azimuthal angle, from positive x-axis
  int q;       // Charge
  float dxy;   // Particle's track's smallest deviation in xy-plane from the
               // primary vertex
  float dz;    // Particle's track's smallest deviation along z-axis from the
               // primary vertex
  float mass;  // Expected invariant mass of the particle
  float E;     // Expected energy of the particle

  /**
   * @brief Construct a new Particle object with mass.
   *
   * @param p trk_p
   * @param pt trk_pt
   * @param eta trk_eta
   * @param phi trk_phi
   * @param q trk_q
   * @param dxy trk_dxy
   * @param dz trk_dz
   * @param mass expected mass of the particle
   */
  Particle(float p, float pt, float eta, float phi, int q, float dxy, float dz,
           float mass);

  /**
   * @brief Construct a new Particle object without mass
   *
   * @param p trk_p
   * @param pt trk_pt
   * @param eta trk_eta
   * @param phi trk_phi
   * @param q trk_q
   * @param dxy trk_dxy
   * @param dz trk_dz
   */
  Particle(float p, float pt, float eta, float phi, int q, float dxy, float dz);

  /**
   * @brief Calculate the energy of the particle
   *
   */
  void calculate_energy();

  /**
   * @brief Calculate the mass of the particle
   *
   */
  void calculate_mass();

  /**
   * @brief Initiates a particle object and calculates its attributes
   *
   * @param p trk_p
   * @param pt trk_pt
   * @param eta trk_eta
   * @param phi trk_phi
   * @param q trk_q
   * @param dxy trk_dxy
   * @param dz trk_dz
   * @param mass expected mass of the particle
   */
  void init_particle(float p, float pt, float eta, float phi, int q, float dxy,
                     float dz, float mass);

private:
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
