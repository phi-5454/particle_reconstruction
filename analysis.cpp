#include "EventCollector.h"

int main(int argc, char** argv)
{
    EventCollector evc;

    evc.initialize_events();
    evc.filter_initial_events();
    evc.analyze("hist0.pdf");
    evc.filter();
    evc.analyze("hist1.pdf");

    return 0;
}