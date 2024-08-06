#include <vector>
#include <iostream>
#include "TMath.h"
#include "Event.h"
#include "Particle.h"
#include "Proton.h"

Event::Event(int ntracks, double zPV, double xPV, double yPV, long eventNum) {
    this->ntracks = ntracks;
    this->zPV = zPV;
    this->xPV = xPV;
    this->yPV = yPV;
    this->EventNum = eventNum;
    this->particles = std::vector<std::vector<std::vector<Particle *>>>{};
    this->protons = std::vector<Proton *>{};
}

Event::Event(int ntracks, double zPV, long eventNum) : Event(ntracks, zPV, 0, 0, eventNum)
{
    // empty
}

void Event::add_particle(double p, double pt, double eta, double phi, int q, double dxy, double dz, double mass, double ptErr,
                      double dxyErr, double dzErr, int i, int j)
{
    auto* part = new Particle(p, pt, eta, phi, q, dxy, dz, mass, ptErr, dxyErr, dzErr);
    particles[i][j].push_back(part);
}

void Event::add_particle(double p, double pt, double eta, double phi, int q, double dxy, double dz, double ptErr,
                      double dxyErr, double dzErr, int i, int j)
{
    add_particle(p, pt, eta, phi, q, dxy, dz, 0, ptErr, dxyErr, dzErr, i, j);
}

void Event::add_proton(double Thx, double Thy, double px, double py)
{
    auto* proton = new Proton(Thx, Thy, px, py);
    protons.push_back(proton);
}

void Event::add_proton(double Thx, double Thy, double px, double py,
                       double pr_px,
                       double pr_py,
                       double pr_pz,
                       double pr_ptx,
                       double pr_pty,
                       double pr_ptx_sigma,
                       double pr_pty_sigma,
                       double pr_posx,
                       double pr_posy,
                       double pr_posx_sigma,
                       double pr_posy_sigma)
{
    auto* proton = new Proton(Thx, Thy, px, py, pr_px, pr_py, pr_pz, pr_ptx, pr_pty, pr_ptx_sigma, pr_pty_sigma, pr_posx, pr_posy, pr_posx_sigma, pr_posy_sigma);
    protons.push_back(proton);
}


Particle* Event::get_particle(int i, int j, int k)
{
    return particles[i][j][k];
}

Proton* Event::get_proton(int i)
{
    return protons.at(i);
}

void Event::set_masses_and_energies(double mass)
{
    for (Particle* &part : particles[0][0])
    {
        part->mass = mass;
        part->calculate_energy();
    }
}

Particle* Event::reconstruct_particle(Particle* p1, Particle* p2)
{
    double E = p1->E + p2->E;
    double px = p1->px + p2->px;
    double py = p1->py + p2->py;
    double pz = p1->pz + p2->pz;
    double p = sqrt(pow(px, 2) + pow(py, 2) + pow(pz, 2));
    double pt = sqrt(pow(px, 2) + pow(py, 2));
    double mass = sqrt(pow(E, 2) - pow(p, 2));
    double eta = TMath::ATanH(pz / p);
    double phi = atan(py / px);
    int q = p1->q + p2->q;
    
    double dxy = sqrt(pow(p1->dxy, 2) + pow(p2->dxy, 2));
    double dz = sqrt(pow(p1->dz, 2) + pow(p2->dz, 2));

    return new Particle(p, pt, px, py, pz, eta, phi, q, dxy, dz, mass, E, 0, 0, 0);
}

void Event::reconstruct()
{
    std::vector<std::vector<Particle *>> init_particles = particles[particles.size() - 1];
    particles.push_back(std::vector<std::vector<Particle*>>{});

    for (int j = 0; j < init_particles.size(); ++j)
    { 
        if (init_particles[0].size() == 2) {
            Particle* part = reconstruct_particle(init_particles[j][0], init_particles[j][1]);
            particles[particles.size() - 1].push_back(std::vector<Particle *>{part});
        }
        else {
            Particle* part1 = reconstruct_particle(init_particles[j][0], init_particles[j][1]);
            Particle* part2 = reconstruct_particle(init_particles[j][2], init_particles[j][3]);
            if (part1->q == 0 && part2->q == 0) { 
                particles[particles.size() - 1].push_back(std::vector<Particle *>{part1, part2});
            }

            part1 = reconstruct_particle(init_particles[j][0], init_particles[j][2]);
            part2 = reconstruct_particle(init_particles[j][1], init_particles[j][3]);
            if (part1->q == 0 && part2->q == 0) {
                particles[particles.size() - 1].push_back(std::vector<Particle *>{part1, part2});
            }

            part1 = reconstruct_particle(init_particles[j][0], init_particles[j][3]);
            part2 = reconstruct_particle(init_particles[j][1], init_particles[j][2]);
            if (part1->q == 0 && part2->q == 0) {
                particles[particles.size() - 1].push_back(std::vector<Particle *>{part1, part2});
            }
        }
    }
}

void Event::print() {
    std::cout << "A new event" << " " << this->EventNum << std::endl;
    for (int i = 0; i < 4; ++i)
        std::cout << this->particles[0][0][i]->p << "\t" << this->particles[0][0][i]->E << std::endl;
}