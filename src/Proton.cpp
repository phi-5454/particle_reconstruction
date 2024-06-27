#include "Proton.h"

Proton::Proton(float Thx, float Thy)
{
    this->Thx = Thx;
    this->Thy = Thy;

    calculate_momenta();
}

void Proton::calculate_momenta()
{
    px = Thx * p;
    py = Thy * p;
}