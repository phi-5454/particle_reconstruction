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
    std::string filepath; // Path to the files with the data in TTrees
    std::string results; // Path to the resulting .root file where the histograms are drawn
    std::vector<Event *> events; // Vector of all the events that are analyzed

    /**
     * @brief Construct a new Event Collector object.
     *
     * @param input The path to the data file/s with trees
     * @param output The path to the created .root file and histograms
     */
    EventCollector(std::string input, std::string output);

    /**
     * @brief Creates and assings all the events (and thus particles).
     *
     * @param new_ntuples Whether the used ntuples are new or old / MC
     * @param new_protons Whether the used protons are new or old
     */
    void initialize_events(bool new_ntuples, bool new_protons);

    /**
     * @brief Finds and returns the lowest and highest values of given lambda function.
     *
     * @tparam F A lambda function
     * @param lambda Property to be measured
     * @return double The lowest and highest values of the property in the event
     */
    template <typename F> 
    std::tuple<double, double> find_min_max(F &&lambda) {
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
     * @brief Create a one-dimensional histogram.
     * 
     * @tparam F A lambda function
     * @param lambda Property to be fitted as a lambda for event. Needs to return a vector of doubles.
     * @param bins Amount of bins
     * @param low Lower limit of bins
     * @param high Upper limit of bins
     * @param title Title of the histogram
     * @param draw Whether or not to draw the histogram
     * @param drawOpt Drawing options for the histogram
     * @param scaleFac The value to scale the second histogram values drawn on top of the first one.
     * @param xtitle Title of the x-axis
     * @param ytitle Title of the y-axis
     * @return TH1F* Drawn histogram
     */
    template <typename F>
    TH1F *create_1Dhistogram(F &&lambda, int bins, double low, double high, std::string title, bool draw,
                             std::string drawOpt, double scaleFac, std::string xtitle, std::string ytitle) {
        TH1F *hist = new TH1F("hist", title.c_str(), bins, low, high);
        for (Event *&event : events) {
            std::vector<double> values = lambda(event);
            for (double value : values)
                hist->Fill(value);
        }
        hist->GetXaxis()->SetTitle(xtitle.c_str());
        hist->GetYaxis()->SetTitle(ytitle.c_str());
        
        // This has to do with the Monte Carlo drawings. Basically it scales the histograms 
        // to have the same max value for presentation purposes.
        if (drawOptions.find("SAME") != -1 ) {
            double max = 1.07*hist->GetMaximum();
            double scale = scaleFac / max;
            hist->Scale(scale);
        }
        if (draw)
            hist->Draw(drawOpt.c_str());
        return hist;
    }

    /**
     * @brief Create a one-dimensional histogram.
     *
     * @tparam F A lambda function
     * @param lambda Property to be fitted as a lambda for event. Needs to return a vector of doubles.
     * @param bins Amount of bins
     * @param low Lower limit of bins
     * @param high Upper limit of bins
     * @param title Title of the histogram
     * @param draw Whether or not to draw the histogram
     * @param drawOpt Drawing options for the histogram
     * @param xtitle Title of the x-axis
     * @param ytitle Title of the y-axis
     * @return TH1F* Drawn histogram
     */
    template <typename F>
    TH1F *create_1Dhistogram(F &&lambda, int bins, double low, double high, std::string title, bool draw,
                             std::string drawOpt, std::string xtitle, std::string ytitle) {
        return create_1Dhistogram(lambda, bins, low, high, title, draw, drawOpt, 0, xtitle, ytitle);
    };

    /**
     * @brief Create a one-dimensional histogram.
     *
     * @tparam F A lambda function
     * @param lambda Property to be fitted as a lambda for event. Needs to return a vector of doubles.
     * @param bins Amount of bins
     * @param low Lower limit of bins
     * @param high Upper limit of bins
     * @param title Title of the histogram
     * @param draw Whether or not to draw the histogram
     * @param xtitle Title of the x-axis
     * @param ytitle Title of the y-axis
     * @return TH1F* Drawn histogram
     */
    template <typename F>
    TH1F *create_1Dhistogram(F &&lambda, int bins, double low, double high, std::string title, bool draw,
                             std::string xtitle, std::string ytitle) {
        return create_1Dhistogram(lambda, bins, low, high, title, draw, "E", 0, xtitle, ytitle);
    };

    /**
     * @brief Create a one-dimensional histogram.
     *
     * @tparam F A lambda function
     * @param lambda Property to be fitted as a lambda for event. Needs to return
     * a vector of doubles.
     * @param bins Amount of bins
     * @param low Lower limit of bins
     * @param high Upper limit of bins
     * @param title Title of the histogram
     * @param draw Whether or not to draw the histogram
     * @return TH1F* Drawn histogram
     */
    template <typename F>
    TH1F *create_1Dhistogram(F &&lambda, int bins, double low, double high,
                             std::string title, bool draw) {
        return create_1Dhistogram(lambda, bins, low, high, title, draw, "E", 0, "GeV", "Events");
    };

    /**
     * @brief Create a one-dimensional histogram.
     *
     * @tparam F A lambda function
     * @param lambda Property to be fitted as a lambda for event. Needs to return
     * a vector of doubles.
     * @param draw Whether or not to draw the histogram
     * @return TH1F* Drawn histogram
     */
    template <typename F> 
    TH1F *create_1Dhistogram(F &&lambda, bool draw) {
        auto [min, max] = find_min_max(lambda);
        return create_1Dhistogram(lambda, round(sqrt(2 * events.size())), min, max,
                                  "A histogram", draw, "E", 0);
    };

    /**
     * @brief Create a histogram with a fit.
     *
     * @tparam F A lambda function
     * @param lambda Property to be fitted as a lambda
     * @param bins Amount of bins
     * @param low Lower limit of bins
     * @param high Upper limit of bins
     * @param title Title of the histogram
     * @param distr The distribution used in the fit. Has to be a predefined function in ROOT (such as "gaus")
     * @return TF1* The fitted function
     */
    template <typename F>
    TF1 *create_1Dhistogram_fit(F &&lambda, int bins, double low, double high,
                                std::string title, std::string distr) {
        TH1F *h1 = create_1Dhistogram(lambda, bins, low, high, title, false);
        h1->Fit(distr.c_str());
        return h1->GetFunction(distr.c_str());
    };

    /**
     * @brief Create a histogram with a fit.
     *
     * @tparam F A lambda function
     * @param lambda Property to be fitted as a lambda
     * @param distr The distribution used in the fit. Has to be a predefined function in ROOT (such as "gaus")
     * @return TF1* The fitted function
     */
    template <typename F>
    TF1 *create_1Dhistogram_fit(F &&lambda, std::string distr) {
        auto [min, max] = find_min_max(lambda);
        return create_1Dhistogram_fit(lambda, round(sqrt(events.size() / 2)), min,
                                      max, "A histogram for a fit", distr);
    };

    /**
     * @brief Create a 2D histogram. MAKE SURE TO PLOT COMPATIBLE VARIABLES.
     *
     * @tparam F A lambda function
     * @param lambda_x Property to be fitted as a lambda for event on x axis. Needs to return a vector of doubles
     * @param lambda_y Property to be fitted as a lambda for event on y axis. Needs to return a vector of doubles
     * @param bins_x Amount of bins on x axis
     * @param low_x Lower limit of bins on x axis
     * @param high_x Upper limit of bins
     * @param bins_y Amount of bins on y axis
     * @param low_y Lower limit of bins on y axis
     * @param high_y Upper limit of bins on y axis
     * @param title Title of the histogram
     * @param draw Whether or not to draw the histogram
     * @param xtitle Title of the x-axis
     * @param ytitle Title of the y-axis
     * @return TH2F* Drawn histogram
     */
    template <typename F1, typename F2>
    TH2F *create_2Dhistogram(F1 &&lambda_x, F2 &&lambda_y, int bins_x, double low_x, double high_x,
                             int bins_y, double low_y, double high_y, const std::string title,
                             bool draw, std::string xtitle, std::string ytitle) {
        TH2F *hist = new TH2F("hist", title.c_str(), bins_x, low_x, high_x, bins_y, low_y, high_y);
        for (Event *&event : events) {
            std::vector<double> values_x = lambda_x(event);
            std::vector<double> values_y = lambda_y(event);
            for (int i = 0; i < values_x.size(); ++i) {
                hist->Fill(values_x[i], values_y[i]);
            }
        }
        hist->GetXaxis()->SetTitle(xtitle.c_str());
        hist->GetYaxis()->SetTitle(ytitle.c_str());
        if (draw)
            hist->Draw("Colz");
        return hist;
    }

    /**
     * @brief Create a 2D histogram. MAKE SURE TO PLOT COMPATIBLE VARIABLES
     *
     * @tparam F A lambda function
     * @param lambda_x Property to be fitted as a lambda for event on x axis. Needs to return a vector of doubles
     * @param lambda_y Property to be fitted as a lambda for event on y axis. Needs to return a vector of doubles
     * @param bins_x Amount of bins on x axis
     * @param low_x Lower limit of bins on x axis
     * @param high_x Upper limit of bins
     * @param bins_y Amount of bins on y axis
     * @param low_y Lower limit of bins on y axis
     * @param high_y Upper limit of bins on y axis
     * @param title Title of the histogram
     * @param draw Whether or not to draw the histogram
     * @return TH2F* Drawn histogram
     */
    template <typename F1, typename F2>
    TH2F *create_2Dhistogram(F1 &&lambda_x, F2 &&lambda_y, int bins_x, double low_x, double high_x,
                             int bins_y, double low_y, double high_y,
                             const std::string& title, bool draw) {
        return create_2Dhistogram(lambda_x, lambda_y, bins_x, low_x, high_x, bins_y, low_y, high_y, title, draw,
                                  "X-axis", "Y-axis");
    }

    /**
     * @brief Create a 2D histogram
     *
     * @tparam F A lambda function
     * @param lambda_x Property to be fitted as a lambda for event on x axis. Needs to return a vector of doubles
     * @param lambda_y Property to be fitted as a lambda for event on y axis. Needs to return a vector of doubles
     * @param title Title of the histogram
     * @param draw Whether or not to draw the histogram
     * @return TH2F* Drawn histogram
     */
    template <typename F1, typename F2> 
    TH2F *create_2Dhistogram(F1 &&lambda_x, F2 &&lambda_y, const std::string& title, bool draw) {
        auto [min_x, max_x] = find_min_max(lambda_x);
        auto [min_y, max_y] = find_min_max(lambda_y);
        return create_2Dhistogram(lambda_x, lambda_y, round(pow(2 * events.size(), 0.5)), min_x, max_x,
                                round(pow(2 * events.size(), 0.5)), min_y, max_y,
                                title, draw);
    }

    /**
     * @brief Filters whole events based on the given lambda function
     *
     * @tparam F A lambda function
     * @param lambda The filter function
     */
    template <typename F> 
    void filter_events(F &&lambda) {
        auto helper = std::vector<Event *>(events.size());
        auto it = std::copy_if(events.begin(), events.end(), helper.begin(), lambda);
        helper.resize(it - helper.begin());
        events = helper;
    };

    /**
     * @brief Filters individual tracks in an event based on the given lambda function
     * 
     * @tparam F A lambda function
     * @param lambda The filter function
     */
    template <typename F>
    void filter_tracks(F &&lambda) {
        for (Event* &event: events)
            event->filter_tracks(lambda);
    }

    /**
     * @brief Filters the events based on the distribution.
     *
     * @tparam F A lambda function
     * @param lambda The property to be filtered with
     * @param distr The distribution function. Has to be a predefined function in ROOT (such as "gaus")
     * @param sigmaMulti The selected area around mean value, +- sigmas. Eg. 3 = -+ 3 sigmas.
     * @param bins Amount of bins
     * @param low Lower limit of bins
     * @param high Upper limit of bins
     * @param title Title of the histogram
     */
    template <typename F>
    void filter_events_distribution(F &&lambda, std::string distr, double sigmaMulti,
                                    int bins, double low, double high, std::string title) {
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
     * @brief Filters the events based on the distribution.
     *
     * @tparam F A lambda function
     * @param lambda The property to be filtered with
     * @param distr The distribution function. Has to be a self-made function
     * @param sigmaMulti The selected area around mean value, +- sigmas. Eg. 3 = -+ 3 sigmas.
     * @param bins Amount of bins
     * @param low Lower limit of bins
     * @param high Upper limit of bins
     * @param title Title of the histogram
     */
    template <typename F>
    void filter_events_distribution(F &&lambda, TF1* func, double sigmaMulti,int bins,
                                    double low, double high, std::string title) {
        TH1 *hist = create_1Dhistogram(lambda, bins, low, high, title, false);
        hist->Fit(func, "0");
        TF1 *fit = hist->GetFunction("fit");
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
        delete hist;
    };

    /**
     * @brief Filters the events based on the distribution
     *
     * @tparam F A Lambda function
     * @param lambda The property to be filtered with
     * @param distr The distribution function. Has to be a predefined function in ROOT (such as "gaus")
     * @param sigmaMulti The selected area around mean value, +- sigmas. Eg. 3 = -+ 3 sigmas.
     */
    template <typename F>
    void filter_events_distribution(F &&lambda, std::string distr, double sigmaMulti) {
        auto [min, max]= find_min_max(lambda);
        filter_events_distribution(lambda, distr, sigmaMulti,
                                   round(sqrt((double)events.size() / 2)), min, max,
                                   "A histogram for a fit");
    };

    /**
     * @brief Filters the reconstructed particle pairs
     * 
     * @tparam F A lambda function
     * @param lambda Function used to filter
     */
    template <typename F>
    void filter_reconstruction(F &&lambda) {
        for (Event* &event : events)
            event->filter_reco(lambda);
    }

    /**
     * @brief Filters the twice reconstructed particles (supposed glueballs) based on the lambda function.
     *
     * @tparam F A lambda function
     * @param lambda Function used to filter
     */
    template <typename F>
    void filter_original(F &&lambda) {
        for (Event* &event : events)
            event->filter_orig(lambda);
    }

    /**
     * @brief Initializes the masses of the particles to given value and calculates 
     * energy based on the mass
     *
     * @param mass Particle's mass
     */
    void init_masses_and_energy(double mass);

    /**
     * @brief Reconstructs one particle from two particles with all permutations in events.
     *
     */
    void reconstruct_particles();
};
#endif
