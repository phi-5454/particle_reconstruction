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
    std::string filepath;
    std::string results;
    std::vector<Event *> events;

    /**
     * @brief Construct a new Event Collector object
     *
     * @param input The path to the data file/s with trees
     * @param output The path to the created .root file and histograms
     */
    EventCollector(std::string input, std::string output);

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
     * @return double The lowest and highest values of the property
     */
    template <typename F> std::tuple<double, double> find_min_max(F &&lambda);

    /**
     * @brief Create a 1D histogram
     *
     * @tparam F A lambda function
     * @param lambda Property to be fitted as a lambda for event. Needs to return
     * a vector of doubles.
     * @param bins Amount of bins
     * @param low Lower limit of bins
     * @param high Upper limit of bins
     * @param title Title of the histogram
     * @param draw Whether or not to draw the histogram
     * @return TH1F*
     */
    template <typename F>
    TH1F *create_1Dhistogram(F &&lambda, int bins, double low, double high,
                             std::string title, bool draw);

    /**
     * @brief Create a 1D histogram for fitting
     *
     * @tparam F A lambda function
     * @param lambda Property to be fitted as a lambda for event. Needs to return
     * a vector of doubles.
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
    TF1 *create_1Dhistogram_fit(F &&lambda, int bins, double low, double high,
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
     * a vector of doubles.
     * @param lambda_x Property to be fitted as a lambda for event on y axis. Needs to return a vector of doubles
     * @param bins Amount of bins
     * @param low Lower limit of bins
     * @param high Upper limit of bins
     * @param title Title of the histogram
     * @param draw Whether or not to draw the histogram
     * @return TH2F*
     */
    template <typename F1, typename F2>
    TH2F *create_2Dhistogram(F1 &&lambda_x, F2 &&lambda_y, int bins_x, double low_x, double high_x,
                             int bins_y, double low_y, double high_y,
                             const std::string& title, bool draw);


    /**
     * @brief Create a 1D histogram for fitting
     *
     * @tparam F A lambda function
     * @param lambda Property to be fitted as a lambda for event. Needs to return
     * a vector of doubles.
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
    void filter_events_distribution(F &&lambda, std::string distr, double sigma,
                                    int bins, double low, double high,
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
    void filter_events_distribution(F &&lambda, std::string distr, double sigma);

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
     * @brief Initializes the masses of the particles to given value and energy
     * based on the mass
     *
     * @param mass Particle's mass
     */
    void init_masses_and_energy(double mass);

    /**
     * @brief Reconstructs one particle per two particles in an event.
     *
     */
    void reconstruct_particles();
};
#endif
