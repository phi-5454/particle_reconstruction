#include "Proton.h"


void Proton::calculate_momenta()
{
    px = Thx * p;
    py = -Thy * p; //Minus to make the axis directions the same as in CMS tracker
}

Proton::Proton(double Thx, double Thy, double px, double py, double pr_px, double pr_py, double pr_pz, double pr_ptx,
               double pr_pty, double pr_ptx_sigma, double pr_pty_sigma, double pr_posx, double pr_posy,
               double pr_posx_sigma, double pr_posy_sigma) {

    this->Thx = Thx;
    this->Thy = Thy;
    this->px = px;
    this->py = py;

    this->pr_px = pr_px;
    this->pr_py = pr_py;
    this->pr_pz = pr_pz;
    this->pr_ptx = pr_ptx;
    this->pr_pty = pr_pty;
    this->pr_ptx_sigma = pr_ptx_sigma;
    this->pr_pty_sigma = pr_pty_sigma;
    this->pr_posx = pr_posx;
    this->pr_posy = pr_posy;
    this->pr_posx_sigma = pr_posx_sigma;
    this->pr_posy_sigma = pr_posy_sigma;
}

Proton::Proton(double Thx, double Thy, double px, double py) :Proton(Thx, Thy, px, py, 0,0,0,0,0,0,0,0,0,0,0){
//calculate_momenta();
}