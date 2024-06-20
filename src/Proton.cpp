#include "Proton.h"


Proton::Proton(float Thx, float Thy)
{
    Thx = Thx;
    Thy = Thy;
}

void Proton::calculate_momenta()
{
    px = Thx * p;
    py = Thy * p;
}

void Proton::initialize_proton(float Thx, float Thy)
{
    Proton(Thx, Thy);
    calculate_momenta();
}