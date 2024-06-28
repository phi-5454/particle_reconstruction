#ifndef PARTICLE_H
#define PARTICLE_H
/**
 * @brief A class for a particle that's detected in the CMS
 *
 */
class Particle {
public:
  double p;     // Momentum
  double pt;    // Transverse (= perpendicular to the beam) momentum
  double ptErr; // Transverse momentum error
  double px;    // X component of momentum, towards the center of the ring
  double py;    // Y component of momentum, "upwards"
  double pz;    // Z component of momentum, parallel to the beam
  double eta;   // Pseudorapidity
  double theta; // Scattering angle, from z-axis.
  double phi;   // Azimuthal angle, from positive x-axis
  int q;       // Charge
  double dxy;   // Particle's track's smallest deviation in xy-plane from the
               // primary vertex
  double dxyErr;// xy-deviation error
  double dz;    // Particle's track's smallest deviation along z-axis from the
               // primary vertex
  double dzErr; // z-deviation error
  double mass;  // Expected invariant mass of the particle
  double E;     // Expected energy of the particle

  /**
   * @brief Construct a new Particle object with known momenta along axis' 
   * 
   * @param p Momentum
   * @param pt Transverse momentun
   * @param px X component of momentun
   * @param py Y component of momentum
   * @param pz Z conponent of momentun
   * @param eta Pseudorapidity
   * @param phi Azimuthal angle
   * @param q Charge
   * @param dxy Assumed deviation from the primary vertex in xy-plane
   * @param dz Assumed deviation from the primary vertex in z-axis
   * @param mass Mass
   * @param E Energy
   * @param ptErr Transerse momentim error
   * @param dxyErr xy-deviation error
   * @param dzErr z-deviation error
   */

  Particle(double p, double pt, double px, double py, double pz, double eta, double phi,
           int q, double dxy, double dz, double mass, double E, double ptErr, double dxyErr, double dzErr);

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
   * @param ptErr Transerse momentim error
   * @param dxyErr xy-deviation error
   * @param dzErr z-deviation error
   */
  Particle(double p, double pt, double eta, double phi, int q, double dxy, double dz,
           double mass, double ptErr, double dxyErr, double dzErr);

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
   * @param ptErr Transerse momentim error
   * @param dxyErr xy-deviation error
   * @param dzErr z-deviation error
   */
  Particle(double p, double pt, double eta, double phi, int q, double dxy, double dz,
           double ptErr, double dxyErr, double dzErr);

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
   * @brief Calculate the total momentum
   * 
   */
  void calculate_momentum();

  /**
   * @brief Prints the particle information
   * 
   */
  void print();

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
