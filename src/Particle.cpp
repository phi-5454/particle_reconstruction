#include "Particle.h"
#include "math.h"
#include <iostream>

Particle::Particle(double p, double pt, double px, double py, double pz, double eta,
                   double phi, int q, double dxy, double dz, double mass, double E,
                   double ptErr, double dxyErr, double dzErr)
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
  this->ptErr = ptErr;
  this->dxyErr = dxyErr;
  this->dzErr = dzErr;

  this->calculate_theta();
}

Particle::Particle(double p, double pt, double eta, double phi, int q, double dxy, double dz,
                   double mass, double ptErr, double dxyErr, double dzErr) : Particle(
                    p, pt, 0, 0, 0, eta, phi, q, dxy, dz, mass, 0, ptErr, dxyErr, dzErr) {
  this->calculate_momenta();
}

Particle::Particle(double p, double pt, double eta, double phi, int q, double dxy,
                   double dz, double ptErr, double dxyErr, double dzErr) : Particle(
                    p, pt, 0, 0, 0, eta, phi, q, dxy, dz, 0, 0, ptErr, dxyErr, dzErr) {
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
