#ifndef EVENT_H
#define EVENT_H

#include <algorithm>
#include <vector>
#include "Particle.h"
#include "Proton.h"
/**
 * @brief A class for an event registered both in the CMS and in the Roman Pots
 * 
 */
class Event
{
public:
    int ntracks; // Amount of particles in the event, not including protons
    // The particles, initial and reconstructed, associated with an event
    /// particles[a][b][c]:
    /// a: The reconstruction layer
    /// b: The alternate combination of lower-level particles on the given reconstruction layer
    /// c: The particle, reconstructed or track-based.
    std::vector<std::vector<std::vector<Particle *>>> particles;
    std::vector<Proton*> protons; // The protons associated with an event
    double zPV; // Z coordinate of the primary vertex of the event
    double xPV; // X coordinate of the primary vertex of the event
    double yPV; // Y coordinate of the primary vertex of the event
    long EventNum; // Event number

    /**
     * @brief Construct a new Event object
     * 
     * @param ntracks ntrk
     * @param zPV zPV
     * @param xPV xPV
     * @param yPV yPV
     * @param eventNum EventNum
     */
    Event(int ntracks, double zPV, double xPV, double yPV, long eventNum);

    /**
     * @brief Construct a new Event object
     * 
     * @param ntracks ntrk
     * @param zPV zPV
     * @param eventNum EventNum
     */
    Event(int ntracks, double zPV, long eventNum);

    /**
     * @brief Adds a particle with mass to the event
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
     * @brief Adds a particle without mass to the event
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

    
    void add_proton(double Thx, double Thy, double px, double py);


    /**
     * @brief Adds a proton to the event
     * 
     * @param Thx ThxL/R
     * @param Thy ThyL/R
     */
    void add_proton(double Thx, double Thy) {
        add_proton(Thx, Thy, 0, 0);
    }


    /**
     * @brief Get the particle object at location [i][j][k]
     * 
     * @param i Particle's iteration level
     * @param j Particles' permutation
     * @param k Particle's position in the vector
     * @return Particle 
     */
    Particle* get_particle(int i, int j, int k);

    /**
     * @brief Get the proton object
     * 
     * @param i Proton's position in the vector
     * @return Proton 
     */
    Proton* get_proton(int i);

    /**
     * @brief Set the masses and energies of the particles
     * 
     * @param mass Mass to be given
     */
    void set_masses_and_energies(double mass);

    /**
     * @brief Filters the particles in an event based on the lambda
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
     * @brief Filters the reconstructed particle pairs based on the lambda
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
     * @brief Filters the original, reconstructed particle based on the lambda
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
     * @brief Reconstructs a new particle from two given particles
     * 
     * @param p1 First particle
     * @param p2 Second particle
     * @return Particle* Reconstructed particle
     */
    Particle* reconstruct_particle(Particle* p1, Particle* p2);

    /**
     * @brief Reconstructs one new particle per two exiting ones and adds them to the particles vector 
     * 
     */
    void reconstruct();

    /**
     * @brief Prints info of the event
     * 
     */
    void print();

    void
    add_proton(double Thx, double Thy, double old_px, double old_py, double pr_px, double pr_py, double pr_pz, double pr_ptx,
               double pr_pty, double pr_ptx_sigma, double pr_pty_sigma, double pr_posx, double pr_posy,
               double pr_posx_sigma, double pr_posy_sigma);
};
#endif