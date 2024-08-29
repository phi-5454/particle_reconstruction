#include "TCanvas.h"
#include "TFile.h"
#include "TMath.h"
#include "TLine.h"
#include "TEllipse.h"
#include "TStyle.h"
#include "TROOT.h"
#include <iostream>
#include <math.h>
#include <fstream>
#include "EventCollector.h"
#include "TApplication.h"

const float RHO_MASS = 0.760;
const float RHO_WIDTH = 0.070;
const float PHI_MASS = 1.021;
const float PHI_WIDTH = 0.034;
const float PION_MASS = 0.13957;
const float KAON_MASS = 0.49367;

/**
 * @brief Create an object for each event and fill it with particles. Also give all particles a certain mass.
 * 
 * @param evc EventCollector
 * @param part Particle type (pion or kaon)
 * @param new_ntuples Whether the used ntuples are new or old / MC
 * @param new_protons Whether the used protons are new or old
 */
void initialize_particles(EventCollector& evc, std::string part, bool new_ntuples, bool new_protons) {
    std::cout << "Initializing events." << std::endl;

    evc.initialize_events(new_ntuples, new_protons);

    if (part == "pion") evc.init_masses_and_energy(PION_MASS);
    else evc.init_masses_and_energy(KAON_MASS);

    std::cout << "Finished initializing events." << std::endl;
}

/**
 * @brief Cauchy or Breit-Wigner distribution function.
 * 
 * @param x X
 * @param par Parameters: Width, location and scale 
 * @return Double_t Value of the function with given parameters
 */
Double_t CauchyDist(Double_t *x, Double_t *par) {
    return par[2] * par[0] / (2 * (pow(par[0], 2) / 4 + pow(x[0] - par[1], 2)));
}

/**
 * @brief Relativistic p-wave Breit-Wigner function for fitting (from https://cds.cern.ch/record/1156140). 
 * 
 * @param x X
 * @param par Parameters: Half-width, location, mass of particle produced in decay (pion) and scale
 * @return double Value of the function with given parameters
 */
Double_t CauchyRelaDist(double *x, double *par) {
    double q = sqrt(pow(x[0], 2) / 4 - pow(par[2], 2));
    double q0 = sqrt(pow(par[1], 2) / 4 - pow(par[2], 2));
    double sigma = par[0] * pow(q / q0, 3) * 2 * pow(q0, 2) / (pow(q0, 2) + pow(q, 2));
    return par[3] * x[0] * par[1] * sigma / (pow(pow(x[0], 2) - pow(par[1], 2), 2) + pow(par[1], 2) * pow(sigma, 2));
}

/**
 * @brief Landau distribution function.
 * 
 * @param x X
 * @param par Parameters: Width, location and scale
 * @return Double_t Value of the function with given parameters
 */
Double_t LandauDist(Double_t *x, Double_t *par) {
    return par[2] * TMath::Landau(x[0], par[0], par[1]);
}

/**
 * @brief Combined Cauchy and Landau distributions.
 * 
 * @param x X
 * @param par Parameters: Cauchy width, location and scale, Landau width, location and scale
 * @return Double_t Value of the function with given parameters
 */
Double_t CauchyLandauDist(Double_t *x, Double_t *par) {
    Double_t p1[3] = {par[0],par[1],par[2]};
    Double_t p2[3] = {par[3],par[4],par[5]};
    return CauchyDist(x, p1) + LandauDist(x, p2) + par[6];
}

/**
 * @brief The coefficient of the term f_i(m) in the Söding model (from https://cds.cern.ch/record/1156140).
 * 
 * @param x X
 * @param par Parameters: Half-width from signal (fixed), location from signal (fixed),
 * mass of particle produced in decay (pion) and free parameter C.
 * @return Double_t Value of the coefficient with given parameters
 */
Double_t SodingFit(double *x, double *par) {
    double q = sqrt(pow(x[0], 2) / 4 - pow(par[2], 2));
    double q0 = sqrt(pow(par[1], 2) / 4 - pow(par[2], 2));
    double sigma = par[0] * pow(q / q0, 3) * 2 * pow(q0, 2) / (pow(q0, 2) + pow(q, 2));
    return par[3] * (pow(par[1], 2) - pow(x[0], 2)) / (x[0] * sigma);
}

/**
 * @brief Additional fit function from Söding model (from https://cds.cern.ch/record/1156140).
 * 
 * @param x X
 * @param par Parameters: BW half-width, location and scale, Söding half-width, location, decay mass and free parameter C
 * @return double Value of the function with given parameters
 */
Double_t CauchyRelaSodingFit(double *x, double *par) {
    double p1[4] = {par[0],par[1],par[5],par[2]};
    double p2[4] = {par[3],par[4],par[5],par[6]};
    return CauchyRelaDist(x, p1) * SodingFit(x, p2);
}

/**
 * @brief Filters the events based on given criteria.
 * 
 * @param evc EventCollector
 */
void filter(EventCollector& evc, bool new_ntuples, bool new_protons) {
    std::cout << "Filtering events." << std::endl;
    std::cout << "Current events: " << evc.events.size() << std::endl;

    // Filter out events with two or less tracks
    evc.filter_events([](Event *event) { return event->ntracks > 2; });
    std::cout << "Without two-tracks: " << evc.events.size() << std::endl;

    // Filter out looper tracks (multiple detections for a single particle)
    evc.filter_events(
        [](Event *event) {
            for (int i = 0; i < event->ntracks - 1; ++i) {
                Particle* part1 = event->get_particle(0, 0, i);
                for (int j = i + 1; j < event->ntracks; ++j) {
                    Particle* part2 = event->get_particle(0, 0, j);
                    if(sqrt(pow(part1->px + part2->px, 2) + pow(part1->py + part2->py, 2) + pow(part1->pz + part2->pz, 2)) < 0.05 /*GeV*/)
                        return false;
                }
            }
            return true;
        }
    );
    std::cout << "Without loopers: " << evc.events.size() << std::endl;

    // Filter our badly reconstructed tracks based on primary vertex Z position
    evc.filter_events_distribution(
        [](Event *event) {
            std::vector<double> values = {event->zPV};
            return values;
        },
        "gaus", 3);
    std::cout << "After ZPV filter: " << evc.events.size() << std::endl;

    // Filter out background particles and bad tracks based on particle's 
    // smallest distance from the primary vertex in xy-plane
    evc.filter_tracks(
        [](Particle* part) {
            return abs(part->dxy) < 0.0435328; // Two sigmas
        }
    );

    // Filter out background particles and bad tracks based on particle's
    // smallest distance from the primary vertex in z-axis
    evc.filter_tracks(
        [](Particle* part) {
            return abs(part->dz) < 0.03 + abs(0.01*part->eta); // Two-dimensional with pseudorapidity. Three sigmas 0.0773883
        }
    );

    if (new_ntuples) { 
    // Filter our badly reconstructed tracks or long-lived particles based on primary collision vertex XY position
    evc.filter_events_distribution(
        [](Event *event) {
            std::vector<double> values = {sqrt(pow(event->xPV, 2) + pow(event->yPV, 2))};
            return values;
        }, "gaus", 3);
    std::cout << "After XYPV filter: " << evc.events.size() << std::endl;
    }

    if (new_protons) {
    /*
    // Filter badly reconstructed proton tracks based on proton vertex x position
    evc.filter_events([](Event * event) {
        const double threshold = 0.05;
        auto prtx_a = event->get_proton(0)->pr_posx;
        auto prtx_b = event->get_proton(1)->pr_posx;
        return abs(prtx_a - prtx_b) < threshold;
    });
    */
/*
    // Filter out elastic proton collisions
    evc.filter_events([](Event * event) {
        double px = 0;
        double py = 0;
        for (int i = 0; i < 2; ++i) {
            Proton *prot = event->get_proton(i);
            px += prot->pr_px;
            py += prot->pr_py;
        }
        return sqrt(px*px / 0.0176815 + py*py / 0.003880805) > 1; // An ellipse, axle length is half a sigma
    });

    std::cout << "Without elastic protons: " << evc.events.size() << std::endl;
     */
    // Filter out badly reconstructed tracks based on proton-track momentum matching.
    // This is currently done when creating the ntuples, this just describes the process and filter values
    evc.filter_events([](Event * event) {
        double trk_px = 0;
        double trk_py = 0;
        double pr_px = 0;
        double pr_py = 0;
        for (int i = 0; i < 2; ++i) {
            Proton *prot = event->get_proton(i);
            pr_px += prot->pr_px;
            pr_py += prot->pr_py;
        }
        for (int i = 0; i < event->particles[0][0].size(); ++i) {
            Particle *part = event->get_particle(0,0,i);
            trk_px += part->px;
            trk_py += part->py;
        }
        auto sumx = trk_px + pr_px;
        auto sumy = trk_py + pr_py;
        const double threshold_x = 0.13;
        const double threshold_y = 0.06;
        bool a = -threshold_x < sumx && sumx < threshold_x;
        bool b = -threshold_y < sumy && sumy < threshold_y;
        return a && b;
    });
    std::cout << "After momentum matcing filter: " << evc.events.size() << std::endl;
    }

    // Filter out everything but four-track events
    evc.filter_events([](Event *event) { return event->ntracks == 4; });
    std::cout << "After four-track + dxy + dz: " << evc.events.size() << std::endl;

    // Filter out everything but net zero charge events. This gives two positive and two negative particles.
    evc.filter_events([](Event *event) {
        int j = 0;
        for (int i = 0; i < event->ntracks; ++i) {
        j += event->get_particle(0, 0, i)->q;
        }
        return j == 0;
    });
    std::cout << "After zero charge: " << evc.events.size() << std::endl;

    std::cout << "Finished filtering events." << std::endl;
}

/**
 * @brief Reconstructs the particles.
 * 
 * @param evc EventCollector
 * @param useMCCoupling Whether to use Monte Carlo coupling and not all possible combinations
 */
void reconstruct(EventCollector& evc, bool useMCCoupling) {
    std::cout << "Reconstructing particles." << std::endl;
    evc.reconstruct_particles(useMCCoupling);
    std::cout << "Finished reconstructing." << std::endl;
}

/**
 * @brief Draw visual lines indicating std deviations onto a histogram fitted with a gaussian distribution
 * 
 * @param hist Histogram to be modified
 * @param sigmas How many sigmas away are the lines drawn
 */
void draw_limits(TH1* hist, double sigmas) {
    hist->Fit("gaus");
    TF1* fit = hist->GetFunction("gaus");
    double max = hist->GetMaximum();
    double down = fit->GetParameter(1) - sigmas * fit->GetParameter(2);
    double up = fit->GetParameter(1) + sigmas * fit->GetParameter(2);
    TLine* l1 = new TLine(down, 0, down, max);
    TLine* l2 = new TLine(up, 0, up, max);
    l1->SetLineStyle(9);
    l1->SetLineColor(7);
    l2->SetLineStyle(9);
    l2->SetLineColor(7);
    l1->Draw();
    l2->Draw();
}

/**
 * @brief Draw histograms of dxy and dz
 * 
 * @param evc EventCollector
 * @param filename Result pdf filename
 * @param drawOpt Drawing options for the histograms
 * @param c1 The canvas the histograms are drawn onto. Has to be divided into four parts
 * @param scales The scales used to scale the histograms drawn on top of each other
 */
void analyze_impact(EventCollector& evc, std::string filename, std::string drawOpt, TCanvas* c1, std::vector<double> scales) {
    c1->cd(1);
    // Particle smallest distance from the primary vertex in xy-plane
    TH1* h11 = evc.create_1Dhistogram(
        [](Event *event) {
            std::vector<double> values(4);
            for (int i = 0; i < 4; ++i)
                values[i] = event->get_particle(0, 0, i)->dxy;
            return values;
        },
        200, -0.3, 0.3, "Distance from primary vertex in xy-plane", true, drawOpt, scales[0], "Distance (mm)", "Events/3 um");

    c1->cd(2);
    // Particle smallest distance from the primary vertex in z-axis
    TH1F* h12 = evc.create_1Dhistogram(
        [](Event *event) {
            std::vector<double> values(4);
            for (int i = 0; i < 4; ++i)
                values[i] = event->get_particle(0, 0, i)->dz;
            return values;
        },
        200, -0.3, 0.3, "Distance from primary vertex in z-axis", true, drawOpt, scales[1], "Distance (mm)", "Events/3 um");

    c1->cd(3);
    // dxy standard deviation
    TH1* h13 = evc.create_1Dhistogram(
        [](Event *event) {
            double avg = 0;
            double value = 0;
            for (int i = 0; i < 4; ++i)
                avg += event->get_particle(0, 0, i)->dxy / 4;
            for (int i = 0; i < 4; ++i)
                value += pow(avg - event->get_particle(0, 0, i)->dxy, 2);
            return std::vector<double>{sqrt(value)};
        }, 100, 0, 0.1, "Std for dxy of events", true, drawOpt, scales[2], "Distance (mm)", "Events");

    c1->cd(4);
    // dz standard deviation
    TH1* h14 = evc.create_1Dhistogram(
        [](Event *event) {
            double avg = 0;
            double value = 0;
            for (int i = 0; i < 4; ++i)
                avg += event->get_particle(0, 0, i)->dz / 4;
            for (int i = 0; i < 4; ++i)
                value += pow(avg - event->get_particle(0, 0, i)->dz, 2);
            return std::vector<double>{sqrt(value)};
        }, 100, 0, 0.1, "Std for dz of events", true, drawOpt, scales[3], "Distance (mm)", "Events");

    c1->SaveAs((filename + "_data_impact.pdf").c_str());
}

/**
 * @brief Draw histograms of phi and eta
 *
 * @param evc EventCollector
 * @param filename Result pdf filename
 * @param drawOpt Drawing options for the histograms
 * @param c1 The canvas the histograms are drawn onto. Has to be divided into four parts
 * @param scales The scales used to scale the histograms drawn on top of each other
 */
void analyze_angles(EventCollector& evc, std::string filename, std::string drawOpt, TCanvas* c1, std::vector<double> scales) {
    c1->cd(1);
    // Particle dxy vs. azimuthal angle
        TH2* h22 = evc.create_2Dhistogram(
        [](Event *event) {
            std::vector<double> values(4);
            for (int i = 0; i < 4; ++i)
                values[i] = event->get_particle(0, 0, i)->dxy;
            return values;
        }, [](Event *event) {
            std::vector<double> values(4);
            for (int i = 0; i < 4; ++i)
                values[i] = event->get_particle(0, 0, i)->phi;
            return values;
        }, 100, -0.2, 0.2, 100, -3.2, 3.2, "Particle dxy vs. azimuthal angle", true, "dxy (mm)", "Azimuthal angle (rad)");
    h22->SetMinimum(5);
    
    c1->cd(2);
    // Particle dz versus pseudorapidity
    TH2* h21 = evc.create_2Dhistogram(
        [](Event *event) {
            std::vector<double> values(4);
            for (int i = 0; i < 4; ++i)
                values[i] = event->get_particle(0, 0, i)->dz;
            return values;
        }, [](Event *event) {
            std::vector<double> values(4);
            for (int i = 0; i < 4; ++i)
                values[i] = event->get_particle(0, 0, i)->eta;
            return values;
        }, 100, -0.1, 0.1, 100, -2.8, 2.8, "Particle dz vs. pseudorapidity", true, "dz (mm)", "Pseudorapidity");
    h21->SetMinimum(5);

    c1->cd(3);
    // Particle azimuthal angle
    TH1* h23 = evc.create_1Dhistogram(
        [](Event *event) {
            std::vector<double> values(4);
            for (int i = 0; i < 4; ++i)
                values[i] = event->get_particle(0, 0, i)->phi;
            return values;
        }, 110, -3.3, 3.3, "Particle track azimuthal angle", true, drawOpt, scales[2], "Azimuthal angle (rad)", "Events/60 mrad");

    c1->cd(4);
    // Particle pseudorapidity
    TH1* h24 = evc.create_1Dhistogram(
        [](Event *event) {
            std::vector<double> values(4);
            for (int i = 0; i < 4; ++i)
                values[i] = event->get_particle(0, 0, i)->eta;
            return values;
        }, 100, -3, 3, "Particle track pseudorapidity", true, drawOpt, scales[3], "Pseudorapidity", "Events/0.06");

    c1->SaveAs((filename + "_data_angles.pdf").c_str());
}

/**
 * @brief Draw histograms of particle and proton matching
 * 
 * @param evc EventCollector
 * @param filename Result pdf filename
 */
void analyze_trackmatch(EventCollector& evc, std::string filename) {
    TCanvas *c12 = new TCanvas("c12", "c12");
    c12->DivideSquare(8);
    c12->Draw();

    c12->cd(1);
    // Proton/track momentum matching
    TH2* h22 = evc.create_2Dhistogram(
            [](Event *event) {
                double px_trk = 0;
                for (int i = 0; i < event->ntracks; ++i)
                    px_trk += event->get_particle(0, 0, i)->px;
                return std::vector{px_trk};
            }, [](Event *event) {
                double px_prt = 0;
                for (int i = 0; i < 2; ++i)
                    px_prt += event->get_proton(i)->pr_px;
                return std::vector{px_prt};
            }, 100, -1, 1, 100, -1, 1, "Track x momentum sum vs. proton x momentum sum ", true, "Momentum (GeV)", "Momentum (GeV)");

    c12->cd(2);
    // Proton/track momentum matching
    TH2* h23 = evc.create_2Dhistogram(
            [](Event *event) {
                std::vector<double> values(4);
                double py_trk = 0;
                for (int i = 0; i < event->ntracks; ++i)
                    py_trk += event->get_particle(0, 0, i)->py;
                return std::vector{py_trk};
            }, [](Event *event) {
                double py_prt = 0;
                for (int i = 0; i < 2; ++i)
                    py_prt += event->get_proton(i)->pr_py;
                return std::vector{py_prt};
            }, 100, -1, 1, 100, -1, 1, "Track y momentum sum vs. proton y momentum sum", true, "Momentum (GeV)", "Momentum (GeV)");

    c12->cd(3);
    // Sum of proton and track x momentum
    TH1* h24 = evc.create_1Dhistogram(
            [](Event *event) {
                double sigma_x_trk = 0;
                for (int i = 0; i < event->ntracks; ++i)
                    sigma_x_trk += event->get_particle(0,0,i)->px;
                double sigma_x_prt = 0;
                for (int i = 0; i < 2; ++i)
                    sigma_x_prt += event->get_proton(i)->pr_px;
                return std::vector{sigma_x_prt+sigma_x_trk};
            },
            200, -1.5, 1.5, "Sum between track x momenta and proton x momenta", true, "Momentum (GeV)", "Momentum (GeV)");

    c12->cd(4);
    TH1* h25 = evc.create_1Dhistogram(
            [](Event *event) {
                double sigma_x_trk = 0;
                for (int i = 0; i < event->ntracks; ++i)
                    sigma_x_trk += event->get_particle(0,0,i)->px;
                double sigma_x_prt = 0;
                for (int i = 0; i < 2; ++i)
                    sigma_x_prt += event->get_proton(i)->pr_px;
                return std::vector{sigma_x_prt+sigma_x_trk};
            },
            200, -1.5, 1.5, "Sum between track y momenta and proton y momenta", true, "Momentum (GeV)", "Momentum (GeV)");

    c12->cd(5);
    // Proton/track momentum matching
    TH2* h26 = evc.create_2Dhistogram(
            [](Event *event) {
                double sum_y_trk = 0;
                for (int i = 0; i < event->ntracks; ++i)
                    sum_y_trk += event->get_particle(0, 0, i)->py;
                return std::vector{sum_y_trk};
            }, [](Event *event) {
                double sum_y_prt = 0;
                for (int i = 0; i < 2; ++i)
                    sum_y_prt += event->get_proton(i)->pr_py;
                return std::vector{sum_y_prt};
            }, 100, -1, 1, 100, -1,1, "track y momentum sum vs. proton y momentum sum ", true, "pr_py", "trk_py");
    //h22->SetMinimum(5);
    h26->Draw();


    c12->cd(6);
    // Proton/track momentum matching
    TH2* h27 = evc.create_2Dhistogram(
            [](Event *event) {
                double sum_y_trk = 0;
                for (int i = 0; i < event->ntracks; ++i)
                    sum_y_trk += event->get_particle(0, 0, i)->py;
                return std::vector{sum_y_trk};
            }, [](Event *event) {
                double sum_y_prt = 0;
                for (int i = 0; i < 2; ++i)
                    sum_y_prt += event->get_proton(i)->pr_py;
                return std::vector{sum_y_prt};
            }, 100, -1, 1, 100, -1,1, "track y momentum sum vs. new proton y momentum sum ", true, "pr_py", "trk_py");
    //h22->SetMinimum(5);
    h27->Draw();


    c12->cd(7);
    TH1* h28 = evc.create_1Dhistogram(
            [](Event *event) {
                double sigma_y_trk = 0;
                for (int i = 0; i < event->ntracks; ++i)
                    sigma_y_trk += event->get_particle(0,0,i)->py;
                double sigma_y_prt = 0;
                for (int i = 0; i < 2; ++i)
                    sigma_y_prt += event->get_proton(i)->pr_py;
                return std::vector{sigma_y_prt+sigma_y_trk};
            },
            200, -1.5, 1.5, "Sum between track y momenta and proton y momenta", true, "", "");
    h28->Draw();

    c12->cd(8);
    TH1* h29 = evc.create_1Dhistogram(
            [](Event *event) {
                double sigma_y_trk = 0;
                for (int i = 0; i < event->ntracks; ++i)
                    sigma_y_trk += event->get_particle(0,0,i)->py;
                double sigma_y_prt = 0;
                for (int i = 0; i < 2; ++i)
                    sigma_y_prt += event->get_proton(i)->pr_py;
                return std::vector{sigma_y_prt+sigma_y_trk};
            },
            200, -1.5, 1.5, "Sum between track y momenta and proton y momenta", true, "", "");
    h28->Draw();

    c12->SaveAs((filename + "_data_prot_trac.pdf").c_str());
}

/**
 * @brief Draw histograms of proton vertex matching
 * 
 * @param evc 
 * @param filename 
 */
void analyze_vertmatch(EventCollector& evc, std::string filename) {
    TCanvas *c20 = new TCanvas("c20", "c20");
    c20->DivideSquare(8);
    c20->Draw();

    c20->cd(1);
    // Proton/track momentum matching
    TH2* h22 = evc.create_2Dhistogram(
            [](Event *event) {
                auto prt_x = event->get_proton(0)->pr_posx;
                return std::vector{prt_x};
            }, [](Event *event) {
                auto prt_x = event->get_proton(1)->pr_posx;
                return std::vector{prt_x};
            }, 100, -1, 1, 100, -1, 1, "Reconstructed vertex x, proton b vs. a", true, "pr_px", "trk_px");
    //h22->SetMinimum(5);
    h22->Draw();

    c20->cd(2);
    // Proton/track momentum matching
    TH2* h23 = evc.create_2Dhistogram(
            [](Event *event) {
                auto prt_y = event->get_proton(0)->pr_posy;
                return std::vector{prt_y};
            }, [](Event *event) {
                auto prt_y = event->get_proton(1)->pr_posy;
                return std::vector{prt_y};
            }, 100, -1, 1, 100, -1, 1, "Reconstructed vertex y, proton b vs. a", true, "pr_py", "pr_py");
    //h22->SetMinimum(5);
    h23->Draw();

    c20->SaveAs((filename + "_data_vertex.pdf").c_str());
}

/**
 * @brief Draw histograms of protons
 *
 * @param evc EventCollector
 * @param filename Result pdf filename
 */
void analyze_protons(EventCollector& evc, std::string filename) {
    TCanvas *c13 = new TCanvas("c13", "c13");
    c13->DivideSquare(4);
    c13->Draw();

    c13->cd(1);
    // X momentums of the two protons
    TH2* h31 = evc.create_2Dhistogram(
        [](Event* event) {
            return std::vector<double>{event->get_proton(0)->pr_px};
        }, [](Event *event) {
            return std::vector<double>{event->get_proton(1)->pr_px};
        }, 300, -1.5, 1.5, 300, -1.5, 1.5, "Proton 1 X momentum vs. proton 2 X momentum", true,
        "Px of proton 1 (GeV)", "Px of proton 2 (GeV)");

    c13->cd(2);
    // Y momentums of the two protons
    TH2* h32 = evc.create_2Dhistogram(
        [](Event* event) {
            return std::vector<double>{event->get_proton(0)->pr_py};
        }, [](Event *event) {
            return std::vector<double>{event->get_proton(1)->pr_py};
        }, 400, -1, 1, 400, -1, 1, "Proton 1 Y momentum vs. proton 2 Y momentum", true,
        "Py of proton 1 (GeV)", "Py of proton 2 (GeV)");

    c13->cd(3);
    // Transverse momentums of the two protons
    TH2* h33 = evc.create_2Dhistogram(
        [](Event* event) {
            return std::vector<double>{sqrt(pow(event->get_proton(0)->pr_py, 2) + pow(event->get_proton(0)->pr_px, 2))};
        }, [](Event *event) {
            return std::vector<double>{sqrt(pow(event->get_proton(1)->pr_py, 2) + pow(event->get_proton(1)->pr_px, 2))};
        }, 300, 0, 1.5, 300, 0, 1.5, "Proton 1 pt vs. proton 2 pt", true,
        "Pt of proton 1 (GeV)", "Pt of proton 2 (GeV)");

    c13->cd(4);
    // Proton X momentum sum vs. Y momentum sum
    TH2* h34 = evc.create_2Dhistogram(
        [](Event* event) {
            return std::vector<double>{event->get_proton(0)->pr_px + event->get_proton(1)->pr_px};
        }, [](Event *event) {
            return std::vector<double>{event->get_proton(0)->pr_py + event->get_proton(1)->pr_py};
        }, 300, -1.5, 1.5, 300, -1.5, 1.5, "Proton X momentum sum vs. Y momentum sum", true,
        "Px of protons (GeV)", "Py of protons (GeV)");

    TEllipse* elli = new TEllipse(0, 0, 0.132972, 0.0623);
    elli->SetFillStyle(0);
    elli->SetLineColor(2);
    elli->Draw();
    
    c13->SaveAs((filename + "_data_protons.pdf").c_str());
}

/**
 * @brief Draw histograms used in filtering
 *
 * @param evc EventCollector
 * @param filename Result pdf filename
 * @param drawOpt Drawing options for the histograms
 * @param c1 The canvas the histograms are drawn onto. Has to be divided into four parts
 * @param scales The scales used to scale the histograms drawn on top of each other
 */
void analyze_filters(EventCollector& evc, std::string filename, std::string drawOpt, TCanvas* c1, std::vector<double> scales) {
    c1->cd(1);
    // Check for loopers in an event
    TH1* h41 = evc.create_1Dhistogram(
        [](Event *event) {
            std::vector<double> values;
            for (int i = 0; i < event->ntracks - 1; ++i) {
                Particle* part1 = event->get_particle(0, 0, i);
                for (int j = i + 1; j < event->ntracks; ++j) {
                    Particle* part2 = event->get_particle(0, 0, j);
                    values.push_back(sqrt(pow(part1->px + part2->px, 2) + pow(part1->py + part2->py, 2) + pow(part1->pz + part2->pz, 2)));
                }
            }
            return values;
        }, 100, 0, 0.5, "Sum of two tracks' total momentum", true, drawOpt, scales[0], "Sum of total momentum (GeV)", "Events/5 MeV");

    TLine* l1 = new TLine(0.05, 0, 0.05, h41->GetMaximum());
    l1->SetLineStyle(9);
    l1->SetLineColor(7);
    l1->Draw();

    c1->cd(2);
    // Z Location of the primary vertex
    TH1* h42 = evc.create_1Dhistogram(
        [](Event *event) {
            return std::vector<double>{event->zPV};
        }, 100, -15, 15, "Z position of event primary vertex", true, drawOpt, scales[1], "Location of the vertex (mm)", "Events/0.3 mm");

    c1->SaveAs((filename + "_data_filters.pdf").c_str());
}

/**
 * @brief Draw histograms of min and max dxy and dz
 *
 * @param evc EventCollector
 * @param filename Result pdf filename
 * @param drawOpt Drawing options for the histograms
 * @param c1 The canvas the histograms are drawn onto. Has to be divided into four parts
 * @param scales The scales used to scale the histograms drawn on top of each other
 */
void analyze_impact_minmax(EventCollector& evc, std::string filename, std::string drawOpt, TCanvas* c1, std::vector<double> scales) {
    c1->cd(1);
    // Min dxy of event
    TH1* h51 = evc.create_1Dhistogram(
        [](Event* event) {
            std::vector<double> values;
            double min = abs(event->get_particle(0, 0, 0)->dxy);
            for (int i = 1; i < 4; ++i) {
                double dxy = abs(event->get_particle(0, 0, i)->dxy);
                if (dxy < min) min = dxy;
            }
            return std::vector<double>{min};
        }, 150, 0, 0.03, "Min dxy of event", true, drawOpt, scales[0], "Distance (mm)", "Events/0,003 mm");

    c1->cd(2);
    // Min dz of event
    TH1* h52 = evc.create_1Dhistogram(
        [](Event* event) {
            std::vector<double> values;
            double min = abs(event->get_particle(0, 0, 0)->dz);
            for (int i = 1; i < 4; ++i) {
                double dz = abs(event->get_particle(0, 0, i)->dz);
                if (dz < min) min = dz;
            }
            return std::vector<double>{min};
        }, 140, 0, 0.03, "Min dz of event", true, drawOpt, scales[1], "Distance (mm)" ,"Events/0,005 mm");

    c1->cd(3);
    // Max dxy of event
    TH1* h53 = evc.create_1Dhistogram(
        [](Event* event) {
            std::vector<double> values;
            double max = abs(event->get_particle(0, 0, 0)->dxy);
            for (int i = 1; i < 4; ++i) {
                double dxy = abs(event->get_particle(0, 0, i)->dxy);
                if (dxy > max) max = dxy;
            }
            return std::vector<double>{max};;
        }, 90, 0, 0.045, "Max dxy of event", true, drawOpt, scales[2], "Distance (mm)", "Events/0,05 mm");

    c1->cd(4);
    // Max dz of event
    TH1* h54 = evc.create_1Dhistogram(
        [](Event* event) {
            std::vector<double> values;
            double max = abs(event->get_particle(0, 0, 0)->dz);
            for (int i = 1; i < 4; ++i) {
                double dz = abs(event->get_particle(0, 0, i)->dz);
                if (dz > max) max = dz;
            }
            return std::vector<double>{max};
        }, 100, 0, 0.06, "Max dz of events", true, drawOpt, scales[3], "Distance (mm)", "Events/0,05 mm");

    c1->SaveAs((filename + "_data_impact_minmax.pdf").c_str());
}

/**
 * @brief Analyzes the initial data that hasn't been reconstructed
 * 
 * @param evc EventCollector
 * @param filename Result pdf filename
 * @param drawOpt Drawing options for the histograms
 * @param c1 The canvas the histograms are drawn onto. Has to be divided into four parts
 */
void analyze_data(EventCollector& evc, std::string filename, std::string drawOpt, std::vector<TCanvas*> c) {
    std::cout << "Analyzing data." << std::endl;
    TFile *results = TFile::Open(evc.results.c_str(), "");

    std::vector<std::vector<Double_t>> scales(c.size(), std::vector<Double_t>(4, 0));
    for (int j = 0; j < c.size(); ++j) {
        for (int i = 1; i < 5; ++i)
            scales[j][i-1] = c[j]->GetPad(i)->GetUymax();
    }

    analyze_impact(evc, filename, drawOpt, c[0], scales[0]);
    analyze_impact_minmax(evc, filename, drawOpt, c[1], scales[1]);
    analyze_angles(evc, filename, drawOpt, c[2], scales[2]);
    analyze_filters(evc, filename, drawOpt, c[3], scales[3]);
    //analyze_protons(evc, filename);

    std::cout << "Finished analyzing data." << std::endl;
}

/**
 * @brief Analyzes the first reconstruction data
 * 
 * @param evc EventCollector
 * @param filename Result pdf filename
 * @param type Type of initial particle
 * @param drawOpt Drawing options for the histograms
 * @param c1 The canvas the histograms are drawn onto. Has to be divided into four parts
 */
void analyze_reco1(EventCollector& evc, std::string filename, std::string type, std::string drawOpt, std::vector<TCanvas*> c) {
    std::cout << "Analyzing the first iteration of recreated particles." << std::endl;
    TFile *results = TFile::Open(evc.results.c_str(), "UPDATE");

    std::vector<std::vector<Double_t>> scales(c.size(), std::vector<Double_t>(4, 0));
    for (int j = 0; j < 1; ++j) {
        for (int i = 1; i < 5; ++i) {
            scales[j][i-1] = c[j]->GetPad(i)->GetUymax();
        }
    }
    scales[1][0] = c[1]->GetUymax();

    int bins = 120;

    float min = 0.9;
    float max = 1.2;
    float partMass = PHI_MASS;
    float partWidth = PHI_WIDTH;

    if (type == "pion" ) {
        min = 0.2;
        max = 1.4;
        partMass = RHO_MASS;
        partWidth = RHO_WIDTH;
    }

    c[0]->cd(1);
    //Azimuthal angle between particles
    /*TH1* h23 = evc.create_1Dhistogram(
        [](Event *event) {
            std::vector<double> values(event->particles[1].size());
            for (int i = 0; i < event->particles[1].size(); ++i) {
                values[i] = abs(event->get_particle(1, i, 0)->phi - event->get_particle(1, i, 1)->phi);
                if (values[i] > TMath::PiOver2()) values[i] = TMath::Pi() - values[i];
            }
            return values;
        }, 160, 0, 1.6, "Azimuthal angle between two particles", true, drawOpt, scales[0][0], "Azimuthal angle (rad)", "Events/10 mrad");
    */

    // Dxy for the particles
    TH1* h23 = evc.create_1Dhistogram(
        [](Event *event) {
            std::vector<double> values(event->particles[1].size() * 2);
            for (int i = 0; i < event->particles[1].size(); ++i)
                for (int j = 0; j < 2; ++j)
                    values[2 * i + j] = event->get_particle(1, i, j)->dxy;
            return values;
        }, 300, 0, 0.3, "Dxy of recreated particles", true, drawOpt, scales[0][0], "Distance (mm)", "Events/1 um");

    c[0]->cd(2);
    // Longitunal angle between particles
    /*TH1* h27 = evc.create_1Dhistogram(
        [](Event *event) {
            std::vector<double> values(event->particles[1].size());
            for (int i = 0; i < event->particles[1].size(); ++i) {
                double ang1 = event->get_particle(1, i, 0)->theta;
                double ang2 = event->get_particle(1, i, 1)->theta;
                if (event->get_particle(1, i, 0)->py < 0) ang1 = -ang1;
                if (event->get_particle(1, i, 1)->py < 0) ang2 = -ang2;
                values[i] = abs(event->get_particle(1, i, 0)->theta - event->get_particle(1, i, 1)->theta);
                if (values[i] > TMath::PiOver2()) values[i] = TMath::Pi() - values[i];
            }
            return values;
        }, 160, 0, 1.6, "Longitunal angle between two particles", true, drawOpt, scales[0][1], "Longitunal angle (rad)", "Events/10 mrad");
    */

    // Dz for the particles
    TH1* h27 = evc.create_1Dhistogram(
        [](Event *event) {
            std::vector<double> values(event->particles[1].size() * 2);
            for (int i = 0; i < event->particles[1].size(); ++i)
                for (int j = 0; j < 2; ++j)
                    values[2 * i + j] = event->get_particle(1, i, j)->dz;
            return values;
        }, 250, 0, 0.5, "Dz of recreated particles", true, drawOpt, scales[0][1], "Distance (mm)", "Events/2 um");

    c[0]->cd(3);
    //  Pseudorapidity of particles
    TH1* h26 = evc.create_1Dhistogram(
        [](Event* event) {
            std::vector<double> values(event->particles[1].size() * 2);
            for (int i = 0; i < event->particles[1].size(); ++i)
                for (int j = 0; j < 2; ++j)
                    values[2*i+j] = event->get_particle(1, i, j)->eta;
            return values;
        }, 100, -5, 5, "Pseudorapidity of recreated particles", false, drawOpt, scales[0][2], "Pseudorapidity", "Events/0,1");

    // Transverse momentum of reconstructed particles
    c[0]->cd(4);
    TH1* h25 = evc.create_1Dhistogram(
        [](Event* event) {
            std::vector<double> values(event->particles[1].size() * 2);
            for (int i = 0; i < event->particles[1].size(); ++i)
                for (int j = 0; j < 2; ++j)
                    values[2*i+j] = event->get_particle(1, i, j)->pt;
            return values;
        }, 200, 0, 2, "Transverse momentum of recreated particles", false, drawOpt, scales[0][3], "Momentum (GeV)", "Events/10 MeV");
    
    c[0]->SaveAs((filename + "_reco1A.pdf").c_str());

    c[1]->cd();
    // Mass distribution of reconstructed particles
    TH1* h24 = evc.create_1Dhistogram(
        [](Event* event) {
            std::vector<double> values(event->particles[1].size() * 2);
            for (int i = 0; i < event->particles[1].size(); ++i)
                for (int j = 0; j < 2; ++j)
                    values[2*i+j] = event->get_particle(1, i, j)->mass;
            return values;
        }, 120, min, max, "Mass of rho meson from diagonal events", true, drawOpt, scales[1][0], "Mass (GeV)", "Events");

    h24->SetName("Mass_rho_diag");
    h24->Write();

    //c[1]->cd();
    // Mass distribution of each reconstructed particle pair
    TH2 *h21 = evc.create_2Dhistogram(
        [](Event *event) {
            std::vector<double> values(event->particles[1].size());
            for (int i = 0; i < event->particles[1].size(); ++i)
                values[i] = event->get_particle(1, i, 0)->mass;
            return values;
        }, [](Event *event) {
            std::vector<double> values(event->particles[1].size());
            for (int i = 0; i < event->particles[1].size(); ++i)
                values[i] = event->get_particle(1, i, 1)->mass;
            return values;
        }, bins, min, max, bins, min, max, "Mass of recreated particles, assumed " + type + "s",
        false, "m_p1 (GeV)", "m_p2 (GeV)");

    // Mass distribution of second particle if first is assumed rho/phi
    int lowbin = bins * (partMass - partWidth - min) / (max - min);
    int highbin = bins * (partMass + partWidth - min) / (max - min);

    TH1* h22 = h21->ProjectionX("px", lowbin, highbin);
    h22->SetTitle("Mass of second particle when first has mass 0.760 +- 0.070 GeV from diagonal events");
    h22->SetName("Mass_rho_diag_second");

    if (drawOpt.find("SAME") != -1 ) {
        double max = 1.07*h22->GetMaximum();
        double scale = scales[1][2] / max;
        h22->Scale(scale);
    }

    //c[1]->cd();
    //h22->Draw(drawOpt.c_str());
    h22->Write();

    //c[1]->cd(3);
    c[1]->SaveAs((filename + "_reco1B.pdf").c_str());

    // Different kinds of fit functions that can be used for fitting
    TF1* f1 = new TF1("CauchyFit", CauchyDist, -15, 15, 3);
    f1->SetParameters(0.15, 0.77, 340);
    //f1->SetParameters(0.01, 1.02, 100);
    f1->SetParNames("Sigma", "Mean", "Scale");
    //h22->Fit("CauchyFit", "", "", 0.7, 0.8);

    TF1* f2 = new TF1("CauchyLandau", CauchyLandauDist, -15, 15, 7);
    f2->SetParameters(0.15, 0.77, 1400, 0.6, 0.12, 90000, 0);
    //f2->SetParameters(0.05, 1.02, 1400, 1.2, 0.05, 150000, 0);
    f2->SetParNames("SigmaC", "MeanC", "ScaleC", "MeanL", "SigmaL", "ScaleL", "Const");

    TF1* f3 = new TF1("CauchySodingFit", CauchyRelaSodingFit, -15, 15, 4);
    f3->SetParNames("Sigma", "Mean", "Scale", "Const");
    f3->FixParameter(0, 0.180);
    f3->FixParameter(1, 0.7551);
    f3->SetParameter("Scale", 1000);
    f3->SetParameter("Const", 0.3);
    f3->SetLineColor(kAzure);
    
    //h24->Fit("CauchySodingFit", "+", "", 0.72, 0.82); // C = -0.316786
    //h24->Fit("CauchyFit", "", "", 0.7, 0.85);
    //h23->Fit("CauchyFit", "", "", 1.014, 1.026);
    //h23->Fit("CauchyLandau", "", "", 0, 2);

    results->Close();
    std::cout << "Finished analyzing the first iteration of recreated particles." << std::endl;
}

/**
 * @brief Filters particles reconstructed from two tracks.
 *
 * @param evc EventCollector
 */
void filter_reco1(EventCollector& evc, std::string part) {
    evc.filter_reconstruction(
        [&](std::vector<Particle*> parts) {
            if (parts.size() == 0) return false;

            double partMass = PHI_MASS;
            double partWidth = PHI_WIDTH;
            double pt_cut = 0.3; // 0.3 seems to be the best for diagonal phis. Doesn't work for parallel,
            // use filter_reco2 eta constraint instead.
            if (part == "pion") {
                partMass = RHO_MASS;
                partWidth = RHO_WIDTH;
                pt_cut = 0.8; // 0.8 from the 2015 paper
            }

            // Mass constraint on the particles
            for (int i = 0; i < parts.size(); ++i) {
                double mass = parts[i]->mass;
                if (mass < partMass - partWidth || mass > partMass + partWidth)
                    return false;
            }

            // Constraint on total pt
            double pt_sum = 0;
            double pt_x = 0;
            double pt_y = 0;
            for (auto & part : parts) {
                pt_sum += part->pt; // For some reason this works for the rho signal
                pt_x += part->px;
                pt_y += part->py;
            }

            //pt_sum = sqrt(pow(pt_x, 2) + pow(pt_y, 2));
            if(pt_sum > pt_cut) return false;

            // Constraint on the azimuthal angle between the particles
            double angle = abs(parts[0]->phi - parts[1]->phi);
            if (angle > TMath::PiOver2()) angle = TMath::Pi() - angle;
            //if (angle < 0.6) return false;

            return true;
        }
    );
}

/**
 * @brief Filters the twice reconstructed particles (supposed glueballs).
 *
 * @param evc EventCollector
 */
void filter_reco2(EventCollector& evc) {
    std::cout << "Filtering the second iteration of recreated particles." << std::endl;
    evc.filter_original(
            [](std::vector<Particle*> part) {
                if(abs(part[0]->eta) > 1.25) return false;
                return true;
            }
    );
    std::cout << "Finished filtering the second iteration of recreated particles." << std::endl;
}

/**
 * @brief Analyzes the twice reconstructed data
 * 
 * @param evc EventCollector
 * @param filename Name of the created histogram pdf file
 * @param drawOpt Drawing options for the histograms
 * @param c1 The canvases the histograms are drawn onto. Have to be divided into four parts
 */
void analyze_reco2(EventCollector& evc, std::string filename, std::string drawOpt, std::vector<TCanvas*> c1) {
    std::cout << "Analyzing the second iteration of recreated particles." << std::endl;
    TFile *results = TFile::Open(evc.results.c_str(), "UPDATE");
    
    std::vector<std::vector<Double_t>> scales(c1.size(), std::vector<Double_t>(4, 0));
    for (int j = 0; j < 3; j = j + 2) {
        for (int i = 1; i < 5; ++i) {
            scales[j][i-1] = c1[j]->GetPad(i)->GetUymax();
        }
    }

    // Reconstructed particle momentum phase space
    TH2* h30 = evc.create_2Dhistogram(
        [](Event *event) {
            std::vector<double> values(event->particles[1].size() * 2);
            for (int i = 0; i < event->particles[1].size(); ++i) {
                for (int j = 0; j < 2; ++j)
                    values[2*i + j] = event->get_particle(1, i, j)->pz;
            }
            return values;
        }, [](Event *event) {
            std::vector<double> values(event->particles[1].size() * 2);
            for (int i = 0; i < event->particles[1].size(); ++i) {
                for (int j = 0; j < 2; ++j)
                    values[2*i + j] = event->get_particle(1, i, j)->pt;
            }
            return values;
        }, 200, -10, 10, 100, 0, 1.5, "Particle pz vs pt of rho mesons", false, "Longitudinal momentum (GeV)", "Transverse momentum (GeV)");

    c1[3]->cd();
    // Reconstructed particle mass
    TH1 *h31 = evc.create_1Dhistogram(
        [](Event *event) {
            std::vector<double> values(event->particles[2].size());
            for (int i = 0; i < event->particles[2].size(); ++i) {
                values[i] = event->get_particle(2, i, 0)->mass;
            }
            return values;
        },
        100, 1, 3, "Mass of the recreated particle from diagonal rho events", true, drawOpt, c1[3]->GetUymax(), "Mass (GeV)", "Events");
    h31->SetName("Mass_glue_from_rho_diag");
    h31->Write();

    TF1* f1 = new TF1("CauchyFit", CauchyDist, 2.1, 2.3, 3);
    f1->SetParameters(0.02, 2.22, 6);
    f1->SetParNames("Sigma", "Mean", "Scale");
    //h31->Fit("CauchyFit", "", "", 2.17, 2.25);

    TF1* f2 = new TF1("CauchyLandau", CauchyLandauDist, -15, 15, 7);
    f2->SetParameters(0.02, 2.22, 10, 0.6, 2.4, 10);
    f2->SetParNames("SigmaC", "MeanC", "ScaleC", "MeanL", "SigmaL", "ScaleL", "Const");
    //h31->Fit("CauchyLandau", "", "", 2.1, 3);

    c1[2]->cd(2);
    // Reconstructed particle pseudorapidity
    TH1* h32 = evc.create_1Dhistogram(
        [](Event *event) {
            std::vector<double> values(event->particles[2].size());
            for (int i = 0; i < event->particles[2].size(); ++i) {
                values[i] = event->get_particle(2, i, 0)->eta;
            }
            return values;
        }, 50, -5, 5, "Pseudorapidity of the recreated particle", true, drawOpt, scales[2][1], "Pseudorapidity", "Events");
    
    c1[2]->cd(3);
    // Reconstructed particle azimuthal angle
    TH1* h33 = evc.create_1Dhistogram(
        [](Event *event) {
            std::vector<double> values(event->particles[2].size());
            for (int i = 0; i < event->particles[2].size(); ++i) {
                values[i] = event->get_particle(2, i, 0)->phi;
            }
            return values;
        }, 64, -3.2, 3.2, "Azimuthal angle of the recreated particle", true, drawOpt, scales[2][2], "Angle (rad)", "Events");

    c1[2]->cd(4);
    // Reconstructed particle transverse momentum
    TH1* h34 = evc.create_1Dhistogram(
        [](Event *event) {
            std::vector<double> values(event->particles[2].size());
            for (int i = 0; i < event->particles[2].size(); ++i) {
                values[i] = event->get_particle(2, i, 0)->pt;
            }
            return values;
        }, 100, 0, 5, "Transverse momentum of the recreated particle", true, drawOpt, scales[2][3], "Momentum (GeV)", "Events");

    c1[2]->SaveAs((filename + "_reco2A.pdf").c_str());

    c1[3]->cd();
    // Particle mass vs. eta
    TH2* h35 = evc.create_2Dhistogram(
        [](Event* event) {
            std::vector<double> values(event->particles[2].size());
            for (int i = 0; i < event->particles[2].size(); ++i)
                values[i] = event->get_particle(2, i, 0)->mass;
            return values;
        }, [](Event* event) {
            std::vector<double> values(event->particles[2].size());
            for (int i = 0; i < event->particles[2].size(); ++i)
                values[i] = event->get_particle(2, i, 0)->eta;
            return values;
        }, 100, 1, 3, 100, -5, 5, "Particle mass vs. eta", false, "Mass (GeV)", "Pseudorapidity");

    c1[3]->SaveAs((filename + "_reco2B.pdf").c_str());

    results->Close();

    std::cout << "Finished analyzing the second iteration of recreated particles." << std::endl;
}

int main()
{
    TApplication app("app", nullptr, nullptr);

    // Set the assumed initial particle here
    const std::string part_type = "pion";

    // EventCollector for the actual data
    EventCollector evc_data(
                //"/eos/cms/store/group/phys_diffraction/CMSTotemLowPU2018/ntuples/data/TOTEM*.root?#tree"
                "/eos/cms/store/group/phys_diffraction/CMSTotemLowPU2018/YounesNtuples/TOTEM2*.root?#tree"
                ,"/afs/cern.ch/user/p/ptuomola/private/particle_reconstruction_results.root"
                //,"particle_reconstruction_results.root"
    );
    // EventCollector for the Monte Carlo data
    EventCollector evc_mc(
                "/eos/cms/store/group/phys_diffraction/CMSTotemLowPU2018/ntuples/mc/phi.root?#tree"
                ,"/afs/cern.ch/user/p/ptuomola/private/particle_reconstruction_mc.root"
    );
    // Canvases used to draw the histograms onto
    TCanvas* c100 = new TCanvas("c100", "c100");
    TCanvas* c101 = new TCanvas("c101", "c101");
    TCanvas* c102 = new TCanvas("c102", "c102");
    TCanvas* c103 = new TCanvas("c103", "c103");
    //TCanvas* c104 = new TCanvas("c104", "c104");

    c100->DivideSquare(4);
    //c101->DivideSquare(4);
    c102->DivideSquare(4);
    //c103->DivideSquare(4);
    //c104->DivideSquare(4);

    std::vector<TCanvas *> c = {c100, c101, c102, c103};

    // Histogram style for the data
    TStyle* data = new TStyle("Data", "Data style");
    data->SetHistFillColor(kBlue);
    data->SetMarkerColor(kBlue);
    data->SetHistLineColor(kBlue);

    // Histogram style for the Monte Carlo
    TStyle* MC = new TStyle("MC", "Monte Carlo style");
    MC->SetHistFillColor(2);
    MC->SetMarkerColor(2);
    MC->SetHistLineColor(2);

    // First initialize the particles for data and MC
    initialize_particles(evc_data, part_type, true, true);
    //initialize_particles(evc_mc, part_type, false, false);

    // Do the initial filtering
    filter(evc_data, true, true);
    //filter(evc_mc, false, false);

    // Draw histograms analyzing the data if wanted
    //data->cd();
    //analyze_data(evc_data, "data1", "E", c);
    //MC->cd();
    //analyze_data(evc_mc, "MC", "E SAME", c);

    // Do the first reconstruction to phi/rho
    reconstruct(evc_data, false);
    //reconstruct(evc_mc, true);


    // Draw histograms analyzing the reconstructed particles if wanted
    //data->cd();
    analyze_reco1(evc_data, "data1", part_type, "E", c);
    //MC->cd();
    //analyze_reco1(evc_mc, "MC", part_type, "E SAME", c);

    // Filter the reconstructed data
    filter_reco1(evc_data, part_type);
    //filter_reco1(evc_mc, part_type);

    // Do the second reconstruction to (hopefully) glueballs
    reconstruct(evc_data, false);
    //reconstruct(evc_mc, false);

    // Filter the final data
    //filter_reco2(evc_data);

    // Analyze the final data
    //data->cd();
    analyze_reco2(evc_data, "data1", "E", c);
    //MC->cd();
    //analyze_reco2(evc_mc, "MC", "E SAME", c);

    // Write the data to a csv file if wanted
    //write_to_csv("testcsv.csv", evc_data);

    app.Run();

    return 0;
}