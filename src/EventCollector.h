#ifndef EVENTCOLLECTOR_H
#define EVENTCOLLECTOR_H

#include <vector>
#include <string>
#include "Event.h"
#include "TFitResultPtr.h"
#include "TH1.h"

/**
 * @brief Handles all the inspected events.
 * 
 */
class EventCollector
{
public:
    std::string filepath = "/eos/cms/store/group/phys_diffraction/CMSTotemLowPU2018/ntuples/data/";
                            /*/eos/user/y/yelberke/TOTEM_2018_ADDEDVARS_OUT*/
    std::string results = "/afs/cern.ch/user/p/ptuomola/private/particle_reconstruction_results.root";
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
     * @brief Create a 1D histogram
     * 
     * @tparam F A lambda function
     * @param lambda Property to be fitted as a lambda for event. Needs to return a float.
     * @param bins Amount of bins
     * @param low Lower limit of bins
     * @param high Upper limit of bins
     * @param title Title of the histogram
     * @param draw Whether or not to draw the histogram
     * @return TH1F* 
     */
    template <typename F>
    TH1F* create_1Dhistogram(F&& lambda, int bins, float low, float high, std::string title, bool draw);

    /**
     * @brief Create a 1D histogram for fitting
     * 
     * @tparam F A lambda function
     * @param lambda Property to be fitted as a lambda for event. Needs to return a float.
     * @param draw Whether or not to draw the histogram
     * @return TH1F* 
     */
    template <typename F>
    TH1F* create_1Dhistogram(F&& lambda, bool draw);

    /**
     * @brief Create a histogram and fit it
     * 
     * @tparam F A lambda function
     * @param lambda Property to be fitted as a lambda
     * @param distr The distribution used in the fit
     * @return TFitResultPtr 
     */
    template <typename F>
    TFitResultPtr create_histogram_fit(F&& lambda, std::string distr);

    /**
     * @brief Filters the events based on the given lambda function
     * 
     * @tparam F A lambda function
     * @param lambda The filter function
     */
    template <typename F>
    void filter_events(F&& lambda);

    /**
     * @brief Filters the events based on the distribution
     * 
     * @tparam F A Lambda function
     * @param lambda The property to be fitted in the histogram that is then fitted
     * @param distr The distribution function
     * @param sigma The selected area around mean value, +- sigmas. Eg. 3 = -+ 3 sigmas.
     */
    template <typename F>
    void filter_events_distribution(F&& lambda, std::string distr, float sigma);

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