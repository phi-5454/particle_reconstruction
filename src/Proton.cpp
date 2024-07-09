#include "Proton.h"

Proton::Proton(double Thx, double Thy)
{
    this->Thx = Thx;
    this->Thy = Thy;

    calculate_momenta();
}

void Proton::calculate_momenta()
{
    px = Thx * p;
    py = -Thy * p; //Minus to make the axis directions the same as in CMS tracker
}