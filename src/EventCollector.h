#ifndef EVENTCOLLECTOR_H
#define EVENTCOLLECTOR_H

#include "Event.h"
#include "TH1.h"
#include "TH2.h"
#include <string>
#include <vector>


/**
 * @brief Handles all the inspected events.
 *
 */
class EventCollector {
public:
//    std::string filepath =
//            "/eos/cms/store/group/phys_diffraction/CMSTotemLowPU2018/ntuples/data/";
//    std::string filepath =
//            "/eos/user/y/yelberke/TOTEM_2018_ADDEDVARS_OUT/";
     std::string results =
     "/afs/cern.ch/user/p/ptuomola/private/particle_reconstruction_results.root";
    //std::string results = "../res.root";
    std::vector<Event *> events;

    /**
     * @brief Construct a new Event Collector object
     *
     */
    EventCollector();

    /**
     * @brief Creates and assings all the events (and thus particles)
     *
     * @param isNew Whether the used data is new or old
     */
    void initialize_events(bool isNew);

    /**
     * @brief Finds and returns the lowest and highest values of given lambda
     * function
     *
     * @tparam F A lambda function
     * @param lambda Property to be measured
     * @return float The lowest and highest values of the property
     */
    template <typename F> std::tuple<float, float> find_min_max(F &&lambda);

    /**
     * @brief Create a 1D histogram
     *
     * @tparam F A lambda function
     * @param lambda Property to be fitted as a lambda for event. Needs to return
     * a vector of floats.
     * @param bins Amount of bins
     * @param low Lower limit of bins
     * @param high Upper limit of bins
     * @param title Title of the histogram
     * @param draw Whether or not to draw the histogram
     * @return TH1F*
     */
    template <typename F>
    TH1F *create_1Dhistogram(F &&lambda, int bins, float low, float high,
                             std::string title, bool draw);

    /**
     * @brief Create a 1D histogram for fitting
     *
     * @tparam F A lambda function
     * @param lambda Property to be fitted as a lambda for event. Needs to return
     * a vector of floats.
     * @param draw Whether or not to draw the histogram
     * @return TH1F*
     */
    template <typename F> TH1F *create_1Dhistogram(F &&lambda, bool draw);

    /**
     * @brief Create a histogram with a fit
     *
     * @tparam F A lambda function
     * @param lambda Property to be fitted as a lambda
     * @param bins Amount of bins
     * @param low Lower limit of bins
     * @param high Upper limit of bins
     * @param title Title of the histogram
     * @param distr The distribution used in the fit
     * @return TF1* The fitted function
     */
    template <typename F>
    TF1 *create_1Dhistogram_fit(F &&lambda, int bins, float low, float high,
                                std::string title, std::string distr);

    /**
     * @brief Create a histogram with a fit
     *
     * @tparam F A lambda function
     * @param lambda Property to be fitted as a lambda
     * @param distr The distribution used in the fit
     * @return TF1* The fitted function
     */
    template <typename F>
    TF1 *create_1Dhistogram_fit(F &&lambda, std::string distr);


    /**
     * @brief Create a 2D histogram. MAKE SURE TO PLOT COMPATIBLE VARIABLES
     *
     * @tparam F lambda function type
     * @param lambda_x Property to be fitted as a lambda for event on x axis. Needs to return
     * a vector of floats.
     * @param lambda_x Property to be fitted as a lambda for event on y axis. Needs to return a vector of floats
     * @param bins Amount of bins
     * @param low Lower limit of bins
     * @param high Upper limit of bins
     * @param title Title of the histogram
     * @param draw Whether or not to draw the histogram
     * @return TH2F*
     */
    template <typename F1, typename F2>
    TH2F *create_2Dhistogram(F1 &&lambda_x, F2 &&lambda_y, int bins_x, float low_x, float high_x,
                             int bins_y, float low_y, float high_y,
                             const std::string& title, bool draw);


    /**
     * @brief Create a 1D histogram for fitting
     *
     * @tparam F A lambda function
     * @param lambda Property to be fitted as a lambda for event. Needs to return
     * a vector of floats.
     * @param draw Whether or not to draw the histogram
     * @return TH1F*
     */
    template <typename F1, typename F2> TH2F *create_2Dhistogram(F1 &&lambda_x, F2 &&lambda_y, const std::string& title, bool draw);

    /**
     * @brief Filters the events based on the given lambda function
     *
     * @tparam F A lambda function
     * @param lambda The filter function
     */
    template <typename F> void filter_events(F &&lambda);

    /**
     * @brief Filters the events based on the distribution
     *
     * @tparam F A Lambda function
     * @param lambda The property to be fitted
     * @param distr The distribution function
     * @param sigma The selected area around mean value, +- sigmas. Eg. 3 = -+ 3
     * sigmas.
     * @param bins Amount of bins
     * @param low Lower limit of bins
     * @param high Upper limit of bins
     * @param title Title of the histogram
     */
    template <typename F>
    void filter_events_distribution(F &&lambda, std::string distr, float sigma,
                                    int bins, float low, float high,
                                    std::string title);

    /**
     * @brief Filters the events based on the distribution
     *
     * @tparam F A Lambda function
     * @param lambda The property to be fitted in the histogram that is then
     * fitted
     * @param distr The distribution function
     * @param sigma The selected area around mean value, +- sigmas. Eg. 3 = -+ 3
     * sigmas.
     */
    template <typename F>
    void filter_events_distribution(F &&lambda, std::string distr, float sigma);

    /**
     * @brief Filters the inial events to fit the provided criteria
     *
     */
    void filter_initial_events();

    /**
     * @brief Analyzes the data
     *
     * @param filename Name of the file (pdf) where the histograms are drawn
     */
    void analyze(std::string filename);

    /**
     * @brief Analyzes the data through the new variables 
     *
     * @param filename Name of the file (pdf) where the histograms are drawn
     */
    void analyze_new(std::string filename);

    /**
     * @brief Analyzes the data of recreated particles
     *
     * @param filename Name of the file (pdf) where the histograms are drawn
     */
    void analyze_reco(std::string filename);

    /**
     * @brief Filters the data
     *
     */
    void filter();

    /**
     * @brief Initializes the masses of the particles to given value and energy
     * based on the mass
     *
     * @param mass Particle's mass
     */
    void init_masses_and_energy(float mass);

    /**
     * @brief Reconstructs one particle per two particles in an event.
     *
     */
    void reconstruct_particles();
};
#endif
