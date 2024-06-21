#ifndef EVENTCOLLECTOR_H
#define EVENTCOLLECTOR_H

#include <vector>
#include <string>
#include "Event.h"

/**
 * @brief Handles all the inspected events.
 * 
 */
class EventCollector
{
public:
    std::string filepath = "/eos/cms/store/group/phys_diffraction/CMSTotemLowPU2018/ntuples/data/";
    std::string results = "/afs/cern.ch/user/p/ptuomola/private/particle_reconstruction_results";
    std::vector<Event*> events;

    /**
     * @brief Construct a new Event Collector object
     * 
     */
    EventCollector();

    /**
     * @brief Creates and assings all the events (and thus particles)
     * 
     */
    void initialize_events();

    /**
     * @brief Filters the inial events to fit the provided criteria
     * 
     */
    void filter_initial_events();

    /**
     * @brief Analyzes the data
     * 
     */
    void analyze();
};
#endif