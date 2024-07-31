#include "Proton.h"

Proton::Proton(double Thx, double Thy, double px, double py)
{
    this->Thx = Thx;
    this->Thy = Thy;
    this->px = px;
    this->py = py;

    //calculate_momenta();
}

void Proton::calculate_momenta()
{
    px = Thx * p;
    py = -Thy * p; //Minus to make the axis directions the same as in CMS tracker
}