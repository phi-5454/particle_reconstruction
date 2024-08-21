#ifndef EVENT_H
#define EVENT_H

#include <algorithm>
#include <vector>
#include "TMath.h"
#include "Particle.h"
#include "Proton.h"
/**
 * @brief A class for an event registered both in the CMS and in the Roman Pots
 * 
 */
class Event
{
public:
    // Old ntuple / MC data parameters
    int ntracks; // Amount of particles in the event, not including protons
    std::vector<std::vector<std::vector<Particle *>>> particles; // The particles, initial and reconstructed, associated with an event
    /// particles[a][b][c]:
    /// a: The reconstruction layer
    /// b: The alternate combination of lower-level particles on the given reconstruction layer
    /// c: The particle, reconstructed or track-based.
    std::vector<Proton*> protons; // The protons associated with an event
    long EventNum; // Event number
    double zPV; // Z coordinate of the primary vertex of the event

    // New ntuple parameters
    double xPV; // X coordinate of the primary vertex of the event
    double yPV; // Y coordinate of the primary vertex of the event

    /**
     * @brief Construct a new Event object from new ntuples.
     * 
     * @param ntracks ntrk
     * @param zPV zPV
     * @param xPV xPV
     * @param yPV yPV
     * @param eventNum EventNum
     */
    Event(int ntracks, double zPV, double xPV, double yPV, long eventNum);

    /**
     * @brief Construct a new Event object from old ntuples.
     * 
     * @param ntracks ntrk
     * @param zPV zPV
     * @param eventNum EventNum
     */
    Event(int ntracks, double zPV, long eventNum);

    /**
     * @brief Adds a particle with mass to the event.
     * 
     * @param p trk_p
     * @param pt trk_pt
     * @param eta trk_eta
     * @param phi trk_phi
     * @param q trk_q
     * @param dxy trk_dxy
     * @param dz trk_dz
     * @param mass Expected mass of the particle
     * @param ptErr trk_pterr
     * @param dxyErr trk_dxyerr
     * @param dzErr trk_dzerr
     * @param i Particle's iteration level
     * @param j Particles' permutation
     */
    void add_particle(double p, double pt, double eta, double phi, int q, double dxy, double dz, double mass,
                      double ptErr, double dxyErr, double dzErr, int i, int j);

    /**
     * @brief Adds a particle without mass to the event.
     * 
     * @param p trk_p
     * @param pt trk_pt
     * @param eta trk_eta
     * @param phi trk_phi
     * @param q trk_q
     * @param dxy trk_dxy
     * @param dz trk_dz
     * @param ptErr trk_pterr
     * @param dxyErr trk_dxyerr
     * @param dzErr trk_dzerr
     * @param i Particle's iteration level
     * @param j Particles' permutation
     */
    void add_particle(double p, double pt, double eta, double phi, int q, double dxy, double dz, double ptErr,
                      double dxyErr, double dzErr, int i, int j);

    /**
     * @brief Adds a proton to the event.
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
    void add_proton(double Thx, double Thy, double pr_px, double pr_py, double pr_pz,
                    double pr_ptx, double pr_pty, double pr_ptx_sigma, double pr_pty_sigma,
                    double pr_posx, double pr_posy, double pr_posx_sigma, double pr_posy_sigma);

    /**
     * @brief Adds a proton with scattering angles to the event.
     * 
     * @param Thx ThxL/R
     * @param Thy ThyL/R
     */
    void add_proton(double Thx, double Thy);

    /**
     * @brief Get the particle object at location [i][j][k].
     * 
     * @param i Particle's iteration level
     * @param j Particles' permutation
     * @param k Particle's position in the vector
     * @return Particle 
     */
    Particle* get_particle(int i, int j, int k);

    /**
     * @brief Get the proton object.
     * 
     * @param i Proton's position in the vector
     * @return Proton 
     */
    Proton* get_proton(int i);

    /**
     * @brief Set the masses and energies of the particles.
     * 
     * @param mass Mass to be given
     */
    void set_masses_and_energies(double mass);

    /**
     * @brief Filters the individual particles in an event based on the lambda function
     * 
     * @tparam F A lambda function
     * @param lambda Function used to filter
     */
    template <typename F> void filter_tracks(F &&lambda) {
        auto helper = std::vector<Particle *>(particles[0][0].size());
        auto it = std::copy_if(particles[0][0].begin(), particles[0][0].end(), helper.begin(), lambda);
        helper.resize(it - helper.begin());
        particles[0][0] = helper;
        ntracks = helper.size();
    }

    /**
     * @brief Filters the reconstructed particle pairs based on the lambda function.
     * 
     * @tparam F A lambda function
     * @param lambda Function used to filter
     */
    template <typename F> void filter_reco(F &&lambda) {
        auto helper = std::vector<std::vector<Particle *>>(particles[1].size());
        auto it = std::copy_if(particles[1].begin(), particles[1].end(), helper.begin(), lambda);
        helper.resize(it - helper.begin());
        particles[1] = helper;
    }

    /**
     * @brief Filters the twice reconstructed particle (supposed glueball) based on the lambda function.
     *
     * @tparam F A lambda function
     * @param lambda Function used to filter
     */
    template <typename F> void filter_orig(F &&lambda) {
        auto helper = std::vector<std::vector<Particle *>>(particles[2].size());
        auto it = std::copy_if(particles[2].begin(), particles[2].end(), helper.begin(), lambda);
        helper.resize(it - helper.begin());
        particles[2] = helper;
    }

    /**
     * @brief Relativistic p-wave Breit-Wigner function for fitting (from https://cds.cern.ch/record/1156140). 
     * 
     * @param x X
     * @param par Parameters: Half-width, location, mass of particle produced in decay (pion) and scale
     * @return double Value of the function with given parameters
     */
    static double CauchyRelaDist(double *x, double *par);

    /**
     * @brief The coefficient of the term f_i(m) in the Söding model (from https://cds.cern.ch/record/1156140).
     * 
     * @param x X
     * @param par Parameters: Half-width from signal (fixed), location from signal (fixed),
     * mass of particle produced in decay (pion) and free parameter C.
     * @return double Value of the coefficient with given parameters
     */
    static double SodingFit(double *x, double *par);

    /**
     * @brief Additional fit function from Söding model (from https://cds.cern.ch/record/1156140).
     * 
     * @param x X
     * @param par Parameters: BW half-width, location and scale, Söding half-width, location, decay mass and free parameter C
     * @return double Value of the function with given parameters
     */
    static double CauchyRelaSodingFit(double *x, double *par);

    /**
     * @brief Reconstructs a new particle from two given particles.
     * 
     * @param p1 First particle
     * @param p2 Second particle
     * @return Particle* Reconstructed particle
     */
    Particle* reconstruct_particle(Particle* p1, Particle* p2);

    /**
     * @brief Reconstructs all possible particles from different particle pairs where q_tot = 0 and adds them to the particles vector.
     * 
     */
    void reconstruct();

    /**
     * @brief Prints info of the event.
     * 
     */
    void print();
};
#endif