#ifndef EVENT_H
#define EVENT_H

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
    std::vector<std::vector<Particle*>> particles; // The particles, initial and reconstructed, associated with an event
    std::vector<Proton*> protons; // The protons associated with an event
    float zPV; // Z coordinate of the primary vertex of the event

    /**
     * @brief Construct a new Event object
     * 
     * @param ntracks ntrk
     * @param zPV zPV
     */
    Event(int ntracks, float zPV);

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
     * @param i Particle's iteration level 
     */
    void add_particle(float p, float pt, float eta, float phi, int q, float dxy, float dz, float mass, int i);

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
     * @param i Particle's iteration level
     */
    void add_particle(float p, float pt, float eta, float phi, int q, float dxy, float dz, int i);

    /**
     * @brief Adds a proton to the event
     * 
     * @param Thx ThxL/R
     * @param Thy ThyL/R
     */
    void add_proton(float Thx, float Thy);

    /**
     * @brief Removes a particle from the event 
     * 
     * @param i Particle's iteration level
     * @param j Particle's position in the vector
     */
    void remove_particle(int i, int j);

    /**
     * @brief Get the particle object at location [i][j]
     * 
     * @param i Particle's iteration level
     * @param j Particle's position in the vector
     * @return Particle 
     */
    Particle* get_particle(int i, int j);

    /**
     * @brief Get the proton object
     * 
     * @param i Proton's position in the vector
     * @return Proton 
     */
    Proton* get_proton(int i);
};
#endif