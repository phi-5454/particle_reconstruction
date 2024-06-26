#include <vector>
#include <iostream>
#include "TMath.h"
#include "Event.h"
#include "Particle.h"
#include "Proton.h"

Event::Event(int ntracks, float zPV)
{
    this->ntracks = ntracks;
    this->zPV = zPV;
    this->particles = std::vector<std::vector<std::vector<Particle *>>>{};
    this->protons = std::vector<Proton *>{};
}

void Event::add_particle(float p, float pt, float eta, float phi, int q, float dxy, float dz, float mass, int i, int j)
{
    Particle* part = new Particle(p, pt, eta, phi, q, dxy, dz, mass);
    particles[i][j].push_back(part);
}

void Event::add_particle(float p, float pt, float eta, float phi, int q, float dxy, float dz, int i, int j)
{
    add_particle(p, pt, eta, phi, q, dxy, dz, 0, i, j);
}

void Event::add_proton(float Thx, float Thy)
{
    Proton* proton = new Proton(Thx, Thy);
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

void Event::set_masses_and_energies(float mass)
{
    for (Particle* &part : particles[0][0])
    {
        part->mass = mass;
        part->E = sqrt(pow(part->mass, 2) + pow(part->p, 2));
    }
}

Particle* Event::reconstruct_particle(Particle* p1, Particle* p2)
{
    float E = p1->E + p2->E;
    float px = p1->px + p2->px;
    float py = p1->py + p2->py;
    float pz = p1->pz + p2->pz;
    float p = sqrt(pow(px, 2) + pow(py, 2) + pow(pz, 2));
    float pt = sqrt(pow(px, 2) + pow(py, 2));
    float mass = sqrt(pow(E, 2) - pow(p, 2));
    float eta = TMath::ATanH(pz / p);
    float phi = atan(py / px);
    int q = p1->q + p2->q;
    
    float dxy = sqrt(pow(p1->dxy, 2) + pow(p2->dxy, 2));
    float dz = sqrt(pow(p1->dz, 2) + pow(p2->dz, 2));

    return new Particle(p, pt, px, py, pz, eta, phi, q, dxy, dz, mass, E);
}

void Event::reconstruct()
{
    std::vector<std::vector<Particle *>> init_particles = particles[particles.size() - 1];
    particles.push_back(std::vector<std::vector<Particle*>>{});

    for (int j = 0; j < init_particles.size(); ++j)
    { 
        for (int i = 0; i < init_particles[0].size() - 1; ++i)
        {
            if (init_particles[0].size() == 2)
            {
                Particle* part = reconstruct_particle(init_particles[j][i], init_particles[j][i + 1]);
                particles[particles.size() - 1].push_back(std::vector<Particle *>{part});
            }
            else
            {
                Particle* part1 = reconstruct_particle(init_particles[j][i], init_particles[j][i + 1]); /*Works only with <=4 particles*/
                Particle* part2 = reconstruct_particle(init_particles[j][(i + 2) % 4], init_particles[j][(i + 3) % 4]);
                if (part1->q ==0 && part2->q == 0) {
                    particles[particles.size() - 1].push_back(std::vector<Particle *>{part1, part2});
                }
            }
        }
    }
}