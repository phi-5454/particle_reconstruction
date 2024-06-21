#include <vector>
#include <iostream>
#include "Event.h"
#include "Particle.h"
#include "Proton.h"

Event::Event(int ntracks, float zPV)
{
    this->ntracks = ntracks;
    this->zPV = zPV;
    this->particles = std::vector<std::vector<Particle*>>{};
    this->protons = std::vector<Proton*>{};
}

void Event::add_particle(float p, float pt, float eta, float phi, int q, float dxy, float dz, float mass, int i)
{
    Particle* part = new Particle(p, pt, eta, phi, q, dxy, dz, mass);
    this->particles[i].push_back(part);
}

void Event::add_particle(float p, float pt, float eta, float phi, int q, float dxy, float dz, int i)
{
    Particle* part = new Particle(p, pt, eta, phi, q, dxy, dz, 0);
    auto j = this->particles[i];
    this->particles[i].push_back(part);
}

void Event::add_proton(float Thx, float Thy)
{
    Proton* proton = new Proton(Thx, Thy);
    this->protons.push_back(proton);
}

void Event::remove_particle(int i, int j)
{
    std::vector<Particle*>::iterator it = this->particles[i].begin() + j;
    this->particles[i].erase(it);
}

Particle* Event::get_particle(int i, int j)
{
    return this->particles[i][j];
}

Proton* Event::get_proton(int i)
{
    return this->protons.at(i);
}