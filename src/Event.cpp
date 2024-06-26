#include <vector>
#include <iostream>
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