#ifndef EVENTCOLLECTOR_H
#define EVENTCOLLECTOR_H

#include "Event.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
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
    template <typename F> std::tuple<double, double> find_min_max(F &&lambda) {
        double min = HUGE_VALF;
        double max = -HUGE_VALF;
        for (Event *&event : events) {
            std::vector<double> values = lambda(event);
            for (double value : values) {
            if (value < min)
                min = value;
            if (value > max)
                max = value;
            }
        }
        return std::make_tuple(min, max);
    };

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
                             std::string title, bool draw) {
          TH1F *hist = new TH1F("hist", title.c_str(), bins, low, high);
        for (Event *&event : events) {
            std::vector<double> values = lambda(event);
            for (double value : values)
                hist->Fill(value);
        }
        if (draw)
            hist->Draw("E");
        return hist;
    };

    /**
     * @brief Create a 1D histogram for fitting
     *
     * @tparam F A lambda function
     * @param lambda Property to be fitted as a lambda for event. Needs to return
     * a vector of doubles.
     * @param draw Whether or not to draw the histogram
     * @return TH1F*
     */
    template <typename F> TH1F *create_1Dhistogram(F &&lambda, bool draw) {
        auto [min, max] = find_min_max(lambda);
        return create_1Dhistogram(lambda, round(sqrt(2 * events.size())), min, max,
                                  "A histogram for a fit", draw);
    };

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
                                std::string title, std::string distr) {
        TH1F *h1 = create_1Dhistogram(lambda, bins, low, high, title, true);
        h1->Fit(distr.c_str());
        return h1->GetFunction(distr.c_str());
    };

    /**
     * @brief Create a histogram with a fit
     *
     * @tparam F A lambda function
     * @param lambda Property to be fitted as a lambda
     * @param distr The distribution used in the fit
     * @return TF1* The fitted function
     */
    template <typename F>
    TF1 *create_1Dhistogram_fit(F &&lambda, std::string distr) {
        auto [min, max] = find_min_max(lambda);
        return create_1Dhistogram_fit(lambda, round(sqrt(events.size() / 2)), min,
                                      max, "A histogram for a fit", distr);
    };


    /**
     * @brief Create a 2D histogram. MAKE SURE TO PLOT COMPATIBLE VARIABLES
     *
     * @tparam F lambda function type
     * @param lambda_x Property to be fitted as a lambda for event on x axis. Needs to return
     * a vector of doubles.
     * @param lambda_y Property to be fitted as a lambda for event on y axis. Needs to return a vector of doubles
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
                             const std::string& title, bool draw) {
        TH2F *hist = new TH2F("hist", title.c_str(), bins_x, low_x, high_x, bins_y, low_y, high_y);
        for (Event *&event : events) {
            std::vector<double> values_x = lambda_x(event);
            std::vector<double> values_y = lambda_y(event);
            for (int i = 0; i < values_x.size(); ++i) {
                hist->Fill(values_x[i], values_y[i]);
            }
        }
        if (draw)
            hist->Draw("Colz");
        return hist;
    }


    /**
     * @brief Create a 2D histogram
     *
     * @tparam F A lambda function
     * @param lambda_x Property to be fitted as a lambda for event on x axis. Needs to return
     * a vector of doubles.
     * @param lambda_y Property to be fitted as a lambda for event on y axis. Needs to return a vector of doubles
     * a vector of doubles.
     * @param title Title of the histogram.
     * @param draw Whether or not to draw the histogram
     * @return TH2F*
     */
    template <typename F1, typename F2> TH2F *create_2Dhistogram(F1 &&lambda_x, F2 &&lambda_y, const std::string& title, bool draw) {
        auto [min_x, max_x] = find_min_max(lambda_x);
        auto [min_y, max_y] = find_min_max(lambda_y);
        return create_2Dhistogram(lambda_x, lambda_y, round(pow(2 * events.size(), 0.5)), min_x, max_x,
                                round(pow(2 * events.size(), 0.5)), min_y, max_y,
                                title, draw);
    }

    /**
     * @brief Filters the events based on the given lambda function
     *
     * @tparam F A lambda function
     * @param lambda The filter function
     */
    template <typename F> void filter_events(F &&lambda) {
        auto helper = std::vector<Event *>(events.size());
        auto it = std::copy_if(events.begin(), events.end(), helper.begin(), lambda);
        helper.resize(it - helper.begin());
        events = helper;
    };

    /**
     * @brief Filters the events based on the distribution
     *
     * @tparam F A Lambda function
     * @param lambda The property to be fitted
     * @param distr The distribution function
     * @param sigmaMulti The selected area around mean value, +- sigmas. Eg. 3 = -+ 3
     * sigmas.
     * @param bins Amount of bins
     * @param low Lower limit of bins
     * @param high Upper limit of bins
     * @param title Title of the histogram
     */
    template <typename F>
    void filter_events_distribution(F &&lambda, std::string distr, double sigmaMulti,
                                    int bins, double low, double high,
                                    std::string title) {
        TF1 *fit = create_1Dhistogram_fit(lambda, bins, low, high, title, distr);
        double mean = fit->GetParameter("Mean");
        double sigma = fit->GetParameter("Sigma");
        filter_events([&](Event *event) {
            std::vector<double> values = lambda(event);
            for (double value : values) {
            if (value < mean - sigmaMulti * sigma ||
                value > mean + sigmaMulti * sigma)
                return false;
            }
            return true;
        });
    };

    /**
     * @brief Filters the events based on the distribution
     *
     * @tparam F A Lambda function
     * @param lambda The property to be fitted in the histogram that is then
     * fitted
     * @param distr The distribution function
     * @param sigmaMulti The selected area around mean value, +- sigmas. Eg. 3 = -+ 3
     * sigmas.
     */
    template <typename F>
    void filter_events_distribution(F &&lambda, std::string distr, double sigmaMulti) {
        auto [min, max]= find_min_max(lambda);
        filter_events_distribution(lambda, distr, sigmaMulti,
                                   round(sqrt((double)events.size() / 2)), min, max,
                                   "A histogram for a fit");
    };

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
