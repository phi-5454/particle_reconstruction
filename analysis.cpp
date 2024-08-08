#include "TCanvas.h"
#include "TFile.h"
#include "TMath.h"
#include "TLine.h"
#include "TEllipse.h"
#include <iostream>
#include <math.h>
#include <fstream>
#include "EventCollector.h"
#include "TApplication.h"

const float RHO_MASS = 0.760;
//const float RHO_WIDTH = 0.236;
//const float RHO_MASS = 0.770;
const float RHO_WIDTH = 0.070;
const float PHI_MASS = 1.021;
const float PHI_WIDTH = 0.034;

/**
 * @brief Create object for each event and fill it with particles. Also give all particles a certain mass.
 * 
 * @param evc EventCollector
 * @param part Particle type (pion or kaon)
 */
void initialize_particles(EventCollector& evc, std::string part) {
    std::cout << "Initializing events." << std::endl;

    evc.initialize_events(true);

    if (part == "pion") evc.init_masses_and_energy(0.13957039);
    else evc.init_masses_and_energy(0.493667);

    std::cout << "Finished initializing events." << std::endl;
}

/**
 * @brief Cauchy or Breit-Wigner distribution function
 * 
 * @param x X
 * @param par Parameters: Width, location and scale 
 * @return Double_t Value of the function with given parameters
 */
Double_t CauchyDist(Double_t *x, Double_t *par) {
    return par[2] * par[0] / (2 * (pow(par[0], 2) / 4 + pow(x[0] - par[1], 2)));
}

/**
 * @brief Landau distribution function
 * 
 * @param x X
 * @param par Parameters: Width, location and scale
 * @return Double_t Value of the function with given parameters
 */
Double_t LandauDist(Double_t *x, Double_t *par) {
    return par[2] * TMath::Landau(x[0], par[0], par[1]);
}

/**
 * @brief Combined Cauchy and Landau distributions
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
 * @brief Filters the events based on given criteria
 * 
 * @param evc EventCollector
 */
void filter(EventCollector& evc) {
    std::cout << "Filtering events." << std::endl;

    // Non-two-track events
    evc.filter_events([](Event *event) { return event->ntracks > 2; });
    //evc.filter_events([](Event *event) { return event->ntracks == 4; });

    // No loopers
    evc.filter_events(
        [](Event *event) {
            for (int i = 0; i < event->ntracks - 1; ++i) {
                Particle* part1 = event->get_particle(0, 0, i);
                for (int j = i + 1; j < event->ntracks; ++j) {
                    Particle* part2 = event->get_particle(0, 0, j);
                    if(sqrt(pow(part1->px + part2->px, 2) + pow(part1->py + part2->py, 2) + pow(part1->pz + part2->pz, 2)) < 0.05)
                        return false;
                }
            }
            return true;
        }
    );

    // Primary vertex XY position
    evc.filter_events_distribution(
        [](Event *event) {
            std::vector<double> values = {sqrt(pow(event->xPV, 2) + pow(event->yPV, 2))};
            return values;
        },
        "gaus", 3);

    // Primary vertex Z position
    evc.filter_events_distribution(
        [](Event *event) {
            std::vector<double> values = {event->zPV};
            return values;
        },
        "gaus", 3);

    // Particle smallest distance from the primary vertex in xy-plane
    evc.filter_tracks(
        [](Particle* part) {
            return abs(part->dxy) < 0.0435328; // Two sigmas
        }
    );

    // Particle smallest distance from the primary vertex in z-axis
    evc.filter_tracks(
        [](Particle* part) {
            return abs(part->dz) < 0.03 + abs(0.01*part->eta); // Two-dimensional with pseudorapidity. Three sigmas 0.0773883
        }
    );

    // No elastic protons
    /*
    evc.filter_events([](Event * event) {
        double px = 0;
        double py = 0;
        for (int i = 0; i < 2; ++i) {
            Proton *prot = event->get_proton(i);
            px += prot->old_px;
            py += prot->old_py;
        }
        return sqrt(px*px / 0.0176815 + py*py / 0.003880805) > 1; // Axle length is half a sigma
    });
     */


    // Proton-track momentum matching
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
        for (int i = 0; i < 4; ++i) {
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

    /*
// Proton vertex x position
    evc.filter_events([](Event * event) {
        const double threshold = 0.05;
        auto prtx_a = event->get_proton(0)->pr_posx;
        auto prtx_b = event->get_proton(1)->pr_posx;
        return abs(prtx_a - prtx_b) < threshold;
    });
     */

    // Four track events
    evc.filter_events([](Event *event) { return event->ntracks == 4; });

    // Net zero charge events
    evc.filter_events([](Event *event) {
        int j = 0;
        for (int i = 0; i < event->ntracks; ++i) {
        j += event->get_particle(0, 0, i)->q;
        }
        return j == 0;
    });

    std::cout << "Finished filtering events." << std::endl;
}

/**
 * @brief Reconstructs the particles from pairs
 * 
 * @param evc EventCollector
 */
void reconstruct(EventCollector& evc) {
    std::cout << "Reconstructing particles." << std::endl;
    evc.reconstruct_particles();
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
 */
void analyze_impact(EventCollector& evc, std::string filename) {
    TCanvas *c11 = new TCanvas("c11", "c11");
    c11->DivideSquare(4);
    c11->Draw();

    c11->cd(1);
    // Particle smallest distance from the primary vertex in xy-plane
    TH1* h11 = evc.create_1Dhistogram(
        [](Event *event) {
            std::vector<double> values(4);
            for (int i = 0; i < 4; ++i)
                values[i] = event->get_particle(0, 0, i)->dxy;
            return values;
        },
        200, -0.3, 0.3, "Distance from primary vertex in xy-plane", true, "Distance (mm)", "Events/3 um");

    draw_limits(h11, 3);

    c11->cd(2);
    // Particle smallest distance from the primary vertex in z-axis
    TH1F* h12 = evc.create_1Dhistogram(
        [](Event *event) {
            std::vector<double> values(4);
            for (int i = 0; i < 4; ++i)
                values[i] = event->get_particle(0, 0, i)->dz;
            return values;
        },
        200, -0.3, 0.3, "Distance from primary vertex in z-axis", true, "Distance (mm)", "Events/3 um");

    draw_limits(h12, 3);

    c11->cd(3);
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
        }, 100, 0, 0.1, "Std for dxy of events", true, "Distance (mm)", "Events");

    c11->cd(4);
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
        }, 100, 0, 0.1, "Std for dz of events", true, "Distance (mm)", "Events");

    c11->SaveAs((filename + "_data_impact.pdf").c_str());
}

/**
 * @brief Draw histograms of phi and eta
 *
 * @param evc EventCollector
 * @param filename Result pdf filename
 */
void analyze_angles(EventCollector& evc, std::string filename) {
    TCanvas *c21 = new TCanvas("c21", "c21");
    //c12->DivideSquare(4);
    c21->Draw();

    c21->cd(1);
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
    
    c21->cd(2);
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

    c21->cd(3);
    // Particle azimuthal angle
    TH1* h23 = evc.create_1Dhistogram(
        [](Event *event) {
            std::vector<double> values(4);
            for (int i = 0; i < 4; ++i)
                values[i] = event->get_particle(0, 0, i)->phi;
            return values;
        }, 110, -3.3, 3.3, "Particle track azimuthal angle", true, "Azimuthal angle (rad)", "Events/20 mrad");

    c21->cd(4);
    // Particle pseudorapidity
    TH1* h24 = h21->ProjectionY();
    h24->SetTitle("Particle pseudorapidity");
    h24->Draw("E");

    c21->SaveAs((filename + "_data_angles.pdf").c_str());
}

void analyze_trackmatch(EventCollector& evc, std::string filename) {
    TCanvas *c12 = new TCanvas("c12", "c12");
    c12->DivideSquare(8);
    c12->Draw();

    c12->cd(1);
    // Proton/track momentum matching
    TH2* h22 = evc.create_2Dhistogram(
            [](Event *event) {
                double sigma_x_trk = 0;
                for (int i = 0; i < event->ntracks; ++i)
                    sigma_x_trk += event->get_particle(0,0,i)->px;
                return std::vector{sigma_x_trk};
            }, [](Event *event) {
                double sigma_x_prt = 0;
                sigma_x_prt += event->get_proton(0)->old_px;
                sigma_x_prt += event->get_proton(1)->old_px;
                return std::vector{sigma_x_prt};
            }, 100, -1, 1, 100, -1, 1, "track x momentum sum vs. proton x momentum sum ", true, "pr_px", "trk_px");
    //h22->SetMinimum(5);
    h22->Draw();
    c12->cd(2);
    // Proton/track momentum matching
    TH2* h23 = evc.create_2Dhistogram(
            [](Event *event) {
                std::vector<double> values(4);
                double sum_x_trk = 0;
                for (int i = 0; i < event->ntracks; ++i)
                    sum_x_trk += event->get_particle(0, 0, i)->px;
                return std::vector{sum_x_trk};
            }, [](Event *event) {
                double sum_x_prt = 0;
                for (int i = 0; i < 2; ++i)
                    sum_x_prt += event->get_proton(i)->pr_px;
                return std::vector{sum_x_prt};
            }, 100, -1, 1, 100, -1,1, "track x momentum sum vs. new proton x momentum sum", true, "new pr_py", "trk_py");
    //h22->SetMinimum(5);
    h23->Draw();

    c12->cd(3);
    TH1* h24 = evc.create_1Dhistogram(
            [](Event *event) {
                double sigma_x_trk = 0;
                for (int i = 0; i < event->ntracks; ++i)
                    sigma_x_trk += event->get_particle(0,0,i)->px;
                double sigma_x_prt = 0;
                for (int i = 0; i < 2; ++i)
                    sigma_x_prt += event->get_proton(i)->old_px;
                return std::vector{sigma_x_prt+sigma_x_trk};
            },
            200, -1.5, 1.5, "Sum between track x momenta and proton x momenta", true, "", "");
    h24->Draw();

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
            200, -1.5, 1.5, "Sum between track x momenta and new proton x momenta", true, "", "");
    h25->Draw();

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
                    sum_y_prt += event->get_proton(i)->old_py;
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
                    sigma_y_prt += event->get_proton(i)->old_py;
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

    c12->SaveAs((filename + "_data_angles.pdf").c_str());
}

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

    c20->SaveAs((filename + "_data_angles.pdf").c_str());
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
            return std::vector<double>{event->get_proton(0)->old_px};
        }, [](Event *event) {
            return std::vector<double>{event->get_proton(1)->old_px};
        }, 300, -1.5, 1.5, 300, -1.5, 1.5, "Proton 1 X momentum vs. proton 2 X momentum", true,
        "Px of proton 1 (GeV)", "Px of proton 2 (GeV)");

    c13->cd(2);
    // Y momentums of the two protons
    TH2* h32 = evc.create_2Dhistogram(
        [](Event* event) {
            return std::vector<double>{event->get_proton(0)->old_py};
        }, [](Event *event) {
            return std::vector<double>{event->get_proton(1)->old_py};
        }, 400, -1, 1, 400, -1, 1, "Proton 1 Y momentum vs. proton 2 Y momentum", true,
        "Py of proton 1 (GeV)", "Py of proton 2 (GeV)");

    c13->cd(3);
    // Transverse momentums of the two protons
    TH2* h33 = evc.create_2Dhistogram(
        [](Event* event) {
            return std::vector<double>{sqrt(pow(event->get_proton(0)->old_py, 2) + pow(event->get_proton(0)->old_px, 2))};
        }, [](Event *event) {
            return std::vector<double>{sqrt(pow(event->get_proton(1)->old_py, 2) + pow(event->get_proton(1)->old_px, 2))};
        }, 300, 0, 1.5, 300, 0, 1.5, "Proton 1 pt vs. proton 2 pt", true,
        "Pt of proton 1 (GeV)", "Pt of proton 2 (GeV)");

    c13->cd(4);
    // Proton X momentum sum vs. Y momentum sum
    TH2* h34 = evc.create_2Dhistogram(
        [](Event* event) {
            return std::vector<double>{event->get_proton(0)->old_px + event->get_proton(1)->old_px};
        }, [](Event *event) {
            return std::vector<double>{event->get_proton(0)->old_py + event->get_proton(1)->old_py};
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
 */
void analyze_filters(EventCollector& evc, std::string filename) {
    TCanvas *c14 = new TCanvas("c14", "c14");
    c14->DivideSquare(4);
    c14->Draw();

    c14->cd(1);
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
        }, 100, 0, 0.5, "Sum of two tracks' total momentum", true, "Sum of total momentum (GeV)", "Events/5 MeV");

    TLine* l1 = new TLine(0.05, 0, 0.05, h41->GetMaximum());
    l1->SetLineStyle(9);
    l1->SetLineColor(7);
    l1->Draw();

    c14->SaveAs((filename + "_data_filters.pdf").c_str());
}

/**
 * @brief Draw histograms of min and max dxy and dz
 *
 * @param evc EventCollector
 * @param filename Result pdf filename
 */
void analyze_impact_minmax(EventCollector& evc, std::string filename) {
    TCanvas* c15 = new TCanvas("c15", "c15");
    c15->DivideSquare(4);
    c15->Draw();

    c15->cd(1);
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
        }, 150, 0, 0.03, "Min dxy of event", true, "Distance (mm)", "Events/0,003 mm");

    c15->cd(2);
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
        }, 140, 0, 0.03, "Min dz of event", true, "Distance (mm)" ,"Events/0,005 mm");

    c15->cd(3);
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
        }, 90, 0, 0.045, "Max dxy of event", true, "Distance (mm)", "Events/0,05 mm");

    c15->cd(4);
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
        }, 100, 0, 0.06, "Max dz of events", true, "Distance (mm)", "Events/0,05 mm");

    c15->SaveAs((filename + "_data_impact_minmax.pdf").c_str());
}

/**
 * @brief Analyzes the initial data that hasn't been reconstructed
 * 
 * @param evc EventCollector
 * @param filename Name of the created histogram pdf file
 */
void analyze_data(EventCollector& evc, std::string filename) {
    std::cout << "Analyzing data." << std::endl;
    TFile *results = TFile::Open(evc.results.c_str(), "");

    TF1* f1 = new TF1("fit", CauchyDist, -15, 15, 3);
    f1->SetParNames("Sigma", "Mean", "Scale");

    analyze_impact(evc, filename);
    analyze_impact_minmax(evc, filename);
    analyze_angles(evc, filename);
    analyze_protons(evc, filename);
    analyze_filters(evc, filename);

    std::cout << "Finished analyzing data." << std::endl;
}

/**
 * @brief Analyzes the first reconstruction data
 * 
 * @param evc EventCollector
 * @param filename Name of the created histogram pdf file
 * @param type Type of initial particle
 */
void analyze_reco1(EventCollector& evc, std::string filename, std::string type) {
    std::cout << "Analyzing the first iteration of recreated particles." << std::endl;
    TFile *results = TFile::Open(evc.results.c_str(), "");

    float min = 0.9;
    float max = 1.2;
    if (type == "pion" ) {
        min = 0.2;
        max = 1.4;
    }

    TCanvas *c21 = new TCanvas("c21", "c21");
    c21->Draw();

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
        }, 120, min, max, 120, min, max, "Mass of recreated particles, assumed " + type + "s",
        true, "m_p1 (GeV)", "m_p2 (GeV)");
    
    c21->SaveAs((filename + "_reco1A.pdf").c_str());

    TCanvas* c22 = new TCanvas("c22", "c22");
    c22->Draw();

    // Mass distribution of second particle if first is assumed rho/phi
    TH1* h22 = evc.create_1Dhistogram(
        [](Event* event) {
            std::vector<std::vector<Particle*>> parts = event->particles[1];
            std::vector<double> values(parts.size());
            for (int i = 0; i < parts.size(); ++i) {
                for (int j = 0; j < 2; ++j) {
                    if (event->get_particle(1, i, j)->mass > RHO_MASS - RHO_WIDTH && event->get_particle(1, i, j)->mass < RHO_MASS + RHO_WIDTH) {
                        values[i] = event->get_particle(1, i, (j+1)%2)->mass;
                        break;
                    }
                }
            }
            return values;
        }, 120, min, max, "Mass of another particle when first is assumed rho meson", true, "Mass (GeV)", "Events");

    // Mass distribution of reconstructed particles
    TH1* h23 = evc.create_1Dhistogram(
        [](Event* event) {
            std::vector<double> values(event->particles[1].size() * 2);
            for (int i = 0; i < event->particles[1].size(); ++i)
                for (int j = 0; j < 2; ++j)
                    values[2*i+j] = event->get_particle(1, i, j)->mass;
            return values;
        }, 120, min, max, "Mass of recreated particles", false, "Mass (GeV)", "Events");

    TF1* f1 = new TF1("CauchyFit", CauchyDist, -15, 15, 3);
    f1->SetParameters(0.15, 0.77, 340);
    //f1->SetParameters(0.01, 1.02, 100);
    f1->SetParNames("Sigma", "Mean", "Scale");

    TF1* f2 = new TF1("CauchyLandau", CauchyLandauDist, -15, 15, 7);
    f2->SetParameters(0.15, 0.77, 1400, 0.6, 0.12, 90000, 0);
    //f2->SetParameters(0.05, 1.02, 1400, 1.2, 0.05, 150000, 0);
    f2->SetParNames("SigmaC", "MeanC", "ScaleC", "MeanL", "SigmaL", "ScaleL", "Const");

    //h23->Fit("CauchyFit", "", "", 0.69, 0.85);
    //h23->Fit("CauchyFit", "", "", 1.014, 1.026);
    //h23->Fit("CauchyLandau", "", "", 0, 2);

    c22->SaveAs((filename + "_reco1B.pdf").c_str());

    std::cout << "Finished analyzing the first iteration of recreated particles." << std::endl;
}

/**
 * @brief Filters particles reconstructed from two tracks
 *
 * @param evc EventCollector
 */
void filter_reco1(EventCollector& evc) {
    evc.filter_reconstruction(
        [](std::vector<Particle*> parts) {
            if (parts.size() == 0) return false;
            double mass_sum = 0;
            for (int i = 0; i < parts.size(); ++i) {
                double mass = parts[i]->mass;
                double pt = parts[i]->pt;
                mass_sum += mass;
                //if (mass < 1.021 - 0.034 || mass > 1.021 + 0.034)
                //if (mass > PHI_MASS - PHI_WIDTH && mass < PHI_MASS + PHI_WIDTH)
                //if (mass < PHI_MASS - PHI_WIDTH || mass > PHI_MASS + PHI_WIDTH)
                if (mass < RHO_MASS - RHO_WIDTH || mass > RHO_MASS + RHO_WIDTH)
                //if (mass > RHO_MASS - RHO_WIDTH && mass < RHO_MASS + RHO_WIDTH)
                //if (mass < 0.754 - 0.062 || mass > 0.754 + 0.064)
                    return false;
            }
            //if(mass_sum <= 1.867 || mass_sum >= 2.243) return false;

            // Constraint on total pt
            // TODO: Study as to what constraints must apply for pions that have come from a rho.
            //  Any place in phase space which would accentuate Fj(2220)?
            double pt_sum = 0;
            for (auto & part : parts) {
                pt_sum += part->pt;
            }

            // From the 2015 paper
             if(pt_sum > 0.800) return false;

            return true;
        }
    );
}


/**
 * @brief Filters on the original, reconstructed particle
 *
 * @param evc EventCollector
 */
void filter_reco2(EventCollector& evc) {
    evc.filter_original(
            [](std::vector<Particle*> part) {
                if (part.size() == 0) return false;
                //std::cerr << part[0]->mass;
                //if(part[0]->pt > 0.200) return false;
                //if(abs(part[0]->eta) >= 1) return false;
                return true;
            }
    );
}

/**
 * @brief Analyzes the second reconstructed data
 * 
 * @param evc EventCollector
 * @param filename Name of the created histogram pdf file
 */
void analyze_reco2(EventCollector& evc, std::string filename) {
    std::cout << "Analyzing the second iteration of recreated particles." << std::endl;
    TFile *results = TFile::Open(evc.results.c_str(), "");
/*
    evc.filter_events(
        [](Event *event) {
            for (int i = 0; i < event->particles[2].size(); ++i) {
                if (event->get_particle(2, i, 0)->eta > 1)
                    return false;
            }
            return true;
        }
    );
*/
    TCanvas *c31 = new TCanvas("c31", "c31");
    c31->Draw();

    c31->SetLogz();
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
        }, 200, -10, 10, 100, 0, 1.5, "Particle pz vs pt of rho mesons", true, "Longitudinal momentum (GeV)", "Transverse momentum (GeV)");
    h30->SetMinimum(1);
    h30->Draw("Colz");

    TH1 *h31 = evc.create_1Dhistogram(
        [](Event *event) {
            std::vector<double> values(4);
            for (int i = 0; i < event->particles[2].size(); ++i) {
                values[i] = event->get_particle(2, i, 0)->mass;
            }
            return values;
        },
        100, 1, 3, "Mass of the recreated particle",
        true);

    TF1* f1 = new TF1("CauchyFit", CauchyDist, 2.1, 2.3, 3);
    f1->SetParameters(0.02, 2.22, 6);
    f1->SetParNames("Sigma", "Mean", "Scale");
    //h31->Fit("CauchyFit", "", "", 2.2, 2.24);
    
    c31->SaveAs((filename + "_reco2A.pdf").c_str());

    TCanvas *c32 = new TCanvas("c32", "c32");
    c32->Draw();

    TH2* h32 = evc.create_2Dhistogram(
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
        }, 100, 1, 3, 100, -5, 5, "Particle mass vs. eta", true, "Mass (GeV)", "Pseudorapidity");

    c32->SaveAs((filename + "_reco2B.pdf").c_str());

    std::cout << "Finished analyzing the second iteration of recreated particles." << std::endl;
}

int write_to_csv(const std::string& filename, const EventCollector& ec){
    std::ofstream file(filename);

    //TODO: error checking on file open/close

    // Write the header
    file
            << "pa1_p" << ","
            << "pa1_pt" << ","
            << "pa1_ptErr" << ","
            << "pa1_px" << ","
            << "pa1_py" << ","
            << "pa1_pz" << ","
            << "pa1_eta" << ","
            << "pa1_theta" << ","
            << "pa1_phi" << ","
            << "pa1_q" << ","
            << "pa1_dxy" << ","
            << "pa1_dxyErr" << ","
            << "pa1_dz" << ","
            << "pa1_dzErr" << ","
            << "pa1_mass" << ","
            << "pa1_E" << ",";
    file
            << "pa2_p" << ","
            << "pa2_pt" << ","
            << "pa2_ptErr" << ","
            << "pa2_px" << ","
            << "pa2_py" << ","
            << "pa2_pz" << ","
            << "pa2_eta" << ","
            << "pa2_theta" << ","
            << "pa2_phi" << ","
            << "pa2_q" << ","
            << "pa2_dxy" << ","
            << "pa2_dxyErr" << ","
            << "pa2_dz" << ","
            << "pa2_dzErr" << ","
            << "pa2_mass" << ","
            << "pa2_E" << ",";
    file
            << "paf_p" << ","
            << "paf_pt" << ","
            << "paf_ptErr" << ","
            << "paf_px" << ","
            << "paf_py" << ","
            << "paf_pz" << ","
            << "paf_eta" << ","
            << "paf_theta" << ","
            << "paf_phi" << ","
            << "paf_q" << ","
            << "paf_dxy" << ","
            << "paf_dxyErr" << ","
            << "paf_dz" << ","
            << "paf_dzErr" << ","
            << "paf_mass" << ","
            << "paf_E" << ",";
    file
            << "pr1_p" << ","
            << "pr1_Thx" << ","
            << "pr1_Thy" << ","
            << "pr1_px" << ","
            << "pr1_py" << ","
            << "pr1_pr_px" << ","
            << "pr1_pr_py" << ","
            << "pr1_pr_pz" << ","
            << "pr1_pr_ptx" << ","
            << "pr1_pr_pty" << ","
            << "pr1_pr_ptx_sigma" << ","
            << "pr1_pr_pty_sigma" << ","
            << "pr1_pr_posx" << ","
            << "pr1_pr_posy" << ","
            << "pr1_pr_posx_sigma" << ","
            << "pr1_pr_posy_sigma" << "," ;
    file
            << "pr2_p" << ","
            << "pr2_Thx" << ","
            << "pr2_Thy" << ","
            << "pr2_px" << ","
            << "pr2_py" << ","
            << "pr2_pr_px" << ","
            << "pr2_pr_py" << ","
            << "pr2_pr_pz" << ","
            << "pr2_pr_ptx" << ","
            << "pr2_pr_pty" << ","
            << "pr2_pr_ptx_sigma" << ","
            << "pr2_pr_pty_sigma" << ","
            << "pr2_pr_posx" << ","
            << "pr2_pr_posy" << ","
            << "pr2_pr_posx_sigma" << ","
            << "pr2_pr_posy_sigma" << "," ;
    file << "\n";

    // Write each struct as a CSV row
    for (const auto& event : ec.events) {
        for (int i = 0; i < event->particles[1].size(); ++i) {
            auto pa1 = event->get_particle(1,i,0);
            auto pa2 = event->get_particle(1,i,1);
            auto paf = event->get_particle(2,i,0);
            auto pr1 = event->get_proton(0);
            auto pr2 = event->get_proton(1);
            // CMS particles
            file
                    << pa1->p << ","
                    << pa1->pt << ","
                    << pa1->ptErr << ","
                    << pa1->px << ","
                    << pa1->py << ","
                    << pa1->pz << ","
                    << pa1->eta << ","
                    << pa1->theta << ","
                    << pa1->phi << ","
                    << pa1->q << ","
                    << pa1->dxy << ","
                    << pa1->dxyErr << ","
                    << pa1->dz << ","
                    << pa1->dzErr << ","
                    << pa1->mass << ","
                    << pa1->E << ",";

            file
                    << pa2->p << ","
                    << pa2->pt << ","
                    << pa2->ptErr << ","
                    << pa2->px << ","
                    << pa2->py << ","
                    << pa2->pz << ","
                    << pa2->eta << ","
                    << pa2->theta << ","
                    << pa2->phi << ","
                    << pa2->q << ","
                    << pa2->dxy << ","
                    << pa2->dxyErr << ","
                    << pa2->dz << ","
                    << pa2->dzErr << ","
                    << pa2->mass << ","
                    << pa2->E << ",";

            file
                    << paf->p << ","
                    << paf->pt << ","
                    << paf->ptErr << ","
                    << paf->px << ","
                    << paf->py << ","
                    << paf->pz << ","
                    << paf->eta << ","
                    << paf->theta << ","
                    << paf->phi << ","
                    << paf->q << ","
                    << paf->dxy << ","
                    << paf->dxyErr << ","
                    << paf->dz << ","
                    << paf->dzErr << ","
                    << paf->mass << ","
                    << paf->E << ",";

            // Protons
            file
                    << pr1->p << ","
                    << pr1->Thx << ","
                    << pr1->Thy << ","
                    << pr1->old_px << ","
                    << pr1->old_py << ","
                    << pr1->pr_px << ","
                    << pr1->pr_py << ","
                    << pr1->pr_pz << ","
                    << pr1->pr_ptx << ","
                    << pr1->pr_pty << ","
                    << pr1->pr_ptx_sigma << ","
                    << pr1->pr_pty_sigma << ","
                    << pr1->pr_posx << ","
                    << pr1->pr_posy << ","
                    << pr1->pr_posx_sigma << ","
                    << pr1->pr_posy_sigma << "," ;

            file
                    << pr2->p << ","
                    << pr2->Thx << ","
                    << pr2->Thy << ","
                    << pr2->old_px << ","
                    << pr2->old_py << ","
                    << pr2->pr_px << ","
                    << pr2->pr_py << ","
                    << pr2->pr_pz << ","
                    << pr2->pr_ptx << ","
                    << pr2->pr_pty << ","
                    << pr2->pr_ptx_sigma << ","
                    << pr2->pr_pty_sigma << ","
                    << pr2->pr_posx << ","
                    << pr2->pr_posy << ","
                    << pr2->pr_posx_sigma << ","
                    << pr2->pr_posy_sigma << "," ;

            file << "\n";
        }

    }

    file.close();
    std::cout << "Wrote data to " << filename << std::endl;
    return 0;
}

int main()
{
    TApplication app("app", nullptr, nullptr);

    const std::string part_type = "pion";
    EventCollector evc(
//             "/eos/cms/store/group/phys_diffraction/CMSTotemLowPU2018/ntuples/data/TOTEM*.root?#tree"
                //"/eos/user/y/yelberke/TOTEM_2018_ADDEDVARS_OUT/combined/TOTEM*.root?#tree"
               //,"/afs/cern.ch/user/p/ptuomola/private/particle_reconstruction_results.root"
            "/home/younes/totemdata/combined/TOTEM2*.root?#tree"
            ,"particle_reconstruction_results.root"
               );

    initialize_particles(evc, part_type);
    filter(evc);
    //analyze_data(evc, "histogram1");
    reconstruct(evc);
    //analyze_reco1(evc, "histogram1", part_type);
    filter_reco1(evc);
    //analyze_trackmatch(evc, "histogram1");
    //analyze_vertmatch(evc, "histogram22");
    reconstruct(evc);
    //filter_reco2(evc);
    analyze_reco2(evc, "histogram1");

    //write_to_csv("testcsv.csv", evc);

    app.Run();

    return 0;
}