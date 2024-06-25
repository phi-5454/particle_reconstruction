#include "Particle.h"
#include "math.h"

Particle::Particle(float p, float pt, float eta, float phi, int q, float dxy,
                   float dz, float mass) {
  this->p = p;
  this->pt = pt;
  this->phi = phi;
  this->q = q;
  this->dxy = dxy;
  this->dz = dz;
  this->mass = mass;
}

Particle::Particle(float p, float pt, float eta, float phi, int q, float dxy,
                   float dz) {
  Particle(p, pt, eta, phi, q, dxy, dz, 0);
}

void Particle::calculate_momenta() {
  px = cos(phi) * pt;
  py = sin(phi) * pt;
  pz = sqrt(pow(p, 2) - pow(pt, 2));
}

void Particle::calculate_theta() { theta = asin(pt / p); }

void Particle::calculate_mass() { mass = sqrt(pow(E, 2) - pow(p, 2)); }

void Particle::calculate_energy() { E = sqrt(pow(mass, 2) + pow(p, 2)); }

void Particle::init_particle(float p, float pt, float eta, float phi, int q,
                             float dxy, float dz, float mass) {
  Particle(p, pt, eta, phi, q, dxy, dz, mass);
  calculate_momenta();
  calculate_theta();
}
