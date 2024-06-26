#include "Particle.h"
#include "math.h"
#include <iostream>

Particle::Particle(float p, float pt, float px, float py, float pz, float eta,
                   float phi, int q, float dxy, float dz, float mass, float E)
                   {
  this->p = p;
  this->pt = pt;
  this->px = px;
  this->py = py;
  this->pz = pz;
  this->eta = eta;
  this->phi = phi;
  this->q = q;
  this->dxy = dxy;
  this->dz = dz;
  this->mass = mass;
  this->E = E;

  this->calculate_theta();
}

Particle::Particle(float p, float pt, float eta, float phi, int q, float dxy, float dz,
                   float mass) : Particle(p, pt, 0, 0, 0, eta, phi, q, dxy, dz, mass, 0) {
  this->calculate_momenta();
}

Particle::Particle(float p, float pt, float eta, float phi, int q, float dxy,
                   float dz) : Particle(p, pt, 0, 0, 0, eta, phi, q, dxy, dz, 0, 0) {
  this->calculate_momenta();
}

void Particle::calculate_momenta() {
  px = cos(phi) * pt;
  py = sin(phi) * pt;
  pz = sqrt(pow(p, 2) - pow(pt, 2));
}

void Particle::calculate_momentum() { p = sqrt(pow(px, 2) + pow(py, 2) + pow(pz, 2)); }

void Particle::calculate_theta() { theta = asin(pt / p); }

void Particle::calculate_mass() { mass = sqrt(pow(E, 2) - pow(p, 2)); }

void Particle::calculate_energy() { E = sqrt(pow(mass, 2) + pow(p, 2)); }

void Particle::print()
{
  std::cout << this->dxy << "\t" << this->p << std::endl;
}
