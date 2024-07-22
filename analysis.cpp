#include "TCanvas.h"
#include "TFile.h"
#include <iostream>
#include <math.h>
#include "EventCollector.h"

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

void initialize_protons(EventCollector& evc) {
    evc.initialize_protons();
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
 * @brief Filters the events based on given criteria
 * 
 * @param evc EventCollector
 */
void filter(EventCollector& evc) {
    std::cout << "Filtering events." << std::endl;

    // Cauchy distribution for filtering
    TF1* f1 = new TF1("fit", CauchyDist, -15, 15, 3);
    f1->SetParNames("Sigma", "Mean", "Scale");

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
            return abs(part->dxy) < 0.3;
        }
    );

    // Same as above, just with event-based filtering instead of particle-based
/*
    f1->SetParameters(0, 0.1, 50);
    evc.filter_events_distribution(
        [](Event *event) {
            std::vector<double> values(4);
            for (int i = 0; i < event->ntracks; ++i)
            values[i] = event->get_particle(0, 0, i)->dxy;
            return values;
        },
        f1, 3, 100, -0.5, 0.5, "Title");
*/
    // Particle smallest distance from the primary vertex in z-axis
    evc.filter_tracks(
        [](Particle* part) {
            return abs(part->dz) < 0.3;
        }
    );

    // Same as above, just with event-based filtering instead of particle-based
/*
    evc.filter_events_distribution(
        [](Event *event) {
            std::vector<double> values(4);
            for (int i = 0; i < event->ntracks; ++i)
            values[i] = event->get_particle(0, 0, i)->dz;
            return values;
        },
        "gaus", 3);
*/
    // No elastic protons
    evc.filter_events([](Event * event) {
        double px = 0;
        double py = 0;
        for (int i = 0; i < 2; ++i) {
            Proton *prot = event->get_proton(i);
            px += prot->px;
            py += prot->py;
        }
        return (abs(px) > 0.05 && abs(py) > 0.05);
    });

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

    TCanvas *c11 = new TCanvas("c11", "c11");
    c11->DivideSquare(4);
    c11->Draw();

    c11->cd(1);
    // Primary vertex Z position
    evc.create_1Dhistogram(
        [](Event *event) {
            std::vector<double> values = {event->zPV};
            return values;
        },
        60, -15, 15, "Primary vertex Z position", true, "Distance (mm)", "Events/0,5 mm");

    c11->cd(2);

    // Primary vertex XY position
    TH1* h12 = evc.create_1Dhistogram(
        [](Event *event) {
            std::vector<double> values = {sqrt(pow(event->xPV, 2) + pow(event->yPV, 2))};
            return values;
        },
        70, 0.12, 0.19, "Primary vertex radial position", true, "Distance (mm)", "Events/1 μm");
/*
    // Proton X momentum sum vs. Y momentum sum
    TH2* h12 = evc.create_2Dhistogram(
        [](Event* event) {
            return std::vector<double>{event->get_proton(0)->px + event->get_proton(1)->px};
        }, [](Event *event) {
            return std::vector<double>{event->get_proton(0)->py + event->get_proton(1)->py};
        }, 100, -2, 2, 75, -1.5, 1.5, "Proton X momentum sum vs. Y momentum sum", true,
        "Px of protons (GeV)", "Py of protons (GeV)");
*/
    c11->cd(3);
    // Particle smallest distance from the primary vertex in xy-plane
    TH1* h13 = evc.create_1Dhistogram(
        [](Event *event) {
            std::vector<double> values(4);
            for (int i = 0; i < 4; ++i)
                values[i] = event->get_particle(0, 0, i)->dxy;
            return values;
        },
        200, -0.3, 0.3, "Distance from primary vertex in xy-plane", true, "Distance (mm)", "Events/3 μm");

    c11->cd(4);
    // Particle smallest distance from the primary vertex in z-axis
    evc.create_1Dhistogram(
        [](Event *event) {
            std::vector<double> values(4);
            for (int i = 0; i < 4; ++i)
                values[i] = event->get_particle(0, 0, i)->dz;
            return values;
        },
        200, -0.3, 0.3, "Distance from primary vertex in z-axis", true, "Distance (mm)", "Events/3 μm");

    c11->SaveAs((filename + "_dataA.pdf").c_str());

    TCanvas *c12 = new TCanvas("c12", "c12");
    c12->DivideSquare(4);
    c12->Draw();

    c12->cd(1);
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
    
    c12->cd(2);
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
        }, 100, -0.3, 0.3, 100, -2.8, 2.8, "Particle dz vs. pseudorapidity", true, "dz (mm)", "Pseudorapidity");
    h21->SetMinimum(5);

    c12->cd(3);
    // Particle azimuthal angle
    TH1* h23 = evc.create_1Dhistogram(
        [](Event *event) {
            std::vector<double> values(4);
            for (int i = 0; i < 4; ++i)
                values[i] = event->get_particle(0, 0, i)->phi;
            return values;
        }, 330, -3.3, 3.3, "Particle track azimuthal angle", true, "Azimuthal angle (rad)", "Events/20 mrad");

    c12->cd(4);
    // Event transverse momentum
    evc.create_1Dhistogram(
        [](Event *event) {
            double px = 0;
            double py = 0;
            for (int i = 0; i < 4; ++i) {
                Particle* part = event->get_particle(0, 0, i);
                px += part->px;
                py += part->py;
            }
            return std::vector<double>{sqrt(pow(px, 2) + pow(py, 2))};
        },
        100, 0, 2, "Event transverse momentum", true, "Transverse momentum (GeV)", "Events/20 MeV");

    c12->SaveAs((filename + "_dataB.pdf").c_str());

    TCanvas *c13 = new TCanvas("c13", "c13");
    c13->DivideSquare(4);
    c13->Draw();

    c13->cd(1);
    // Particle dxy / dxy error
    evc.create_1Dhistogram(
        [](Event *event) {
            std::vector<double> values(4);
            for (int i = 0; i < 4; ++i)
                values[i] = event->get_particle(0, 0, i)->dxy / event->get_particle(0, 0, i)->dxyErr;
            return values;
        },
        100, -5, 5, "Particle dxy / dxyErr", true, "dxy/dxy error", "Events/0,1");

    c13->cd(2);
    // Particle dz / dz error
    evc.create_1Dhistogram(
        [](Event *event) {
            std::vector<double> values(4);
            for (int i = 0; i < 4; ++i)
                values[i] = event->get_particle(0, 0, i)->dz / event->get_particle(0, 0, i)->dzErr;
            return values;
        },
        200, -5, 5, "Particle dz / dzErr", true, "dz/dz error", "Events/0,05");
    
    c13->cd(3);
    // Particle pt error / pt
    evc.create_1Dhistogram(
        [](Event *event) {
            std::vector<double> values(4);
            for (int i = 0; i < 4; ++i)
                values[i] = event->get_particle(0, 0, i)->ptErr / event->get_particle(0, 0, i)->pt;
            return values;
        },
        100, 0, 0.1, "Particle ptErr / pt", true, "pt error / pt", "Events/0,1%");

    c13->cd(4);
    // Check for loopers in an event
    evc.create_1Dhistogram(
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

    c13->SaveAs((filename + "_dataC.pdf").c_str());

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

    TCanvas *c21 = new TCanvas("c21", "c21");
    c21->Draw();

    float min = 0.2;
    float max = 2.6;
    if (type == "pion" ) {
        min = 0.2;
        max = 2.6;
    }

    TH2 *h21 = evc.create_2Dhistogram(
        [](Event *event) {
            std::vector<double> values(2);
            for (int i = 0; i < event->particles[1].size(); ++i)
                values[i] = event->get_particle(1, i, 0)->mass;
            return values;
        }, [](Event *event) {
            std::vector<double> values(2);
            for (int i = 0; i < event->particles[1].size(); ++i)
                values[i] = event->get_particle(1, i, 1)->mass;
            return values;
        }, 240, min, max, 240, min, max, "Mass of recreated particles, assumed " + type + "s",
        true, "m_p1 (GeV)", "m_p2 (GeV)");

    TF1* f1 = new TF1("CauchyFit", CauchyDist, -15, 15, 3);
    f1->SetParameters(0.1, 1.02, 30);
    f1->SetParNames("Sigma", "Mean", "Scale");

    //h21->Fit("CauchyFit", "", "", 1.01, 1.03);
    
    c21->SaveAs((filename + "_reco1A.pdf").c_str());

    TCanvas* c22 = new TCanvas("c2", "c2");
    c22->Draw();

    std::cout << "Finished analyzing the first iteration of recreated particles." << std::endl;
    
    TH1 *h22 = h21->ProjectionX();
    h22->Add(h21->ProjectionY());
    h22->GetXaxis()->SetTitle("Mass (GeV)");
    h22->GetYaxis()->SetTitle("Events/10 MeV");
    h22->Draw("E");

    evc.filter_reconstruction(
        [](std::vector<Particle*> parts) {
            if (parts.size() == 0) return false;
            for (int i = 0; i < parts.size(); ++i) {
                double mass = parts[i]->mass;
                if (mass < 0.77 - 0.15 || mass > 0.77 + 0.15)
                    return false;
            }
            return true;
        }
    );

    c22->SaveAs((filename + "_reco1B.pdf").c_str());
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
            for (int i = 0; i < 2; ++i) {
                if (event->get_particle(2, i, 0)->eta > 0.65)
                    return false;
            }
            return true;
        }
    );
*/
    TCanvas *c31 = new TCanvas("c31", "c31");
    c31->Draw();

    c31->cd(1);
    TH1 *h31 = evc.create_1Dhistogram(
        [](Event *event) {
            std::vector<double> values(4);
            for (int i = 0; i < event->particles[2].size(); ++i) {
                values[i] = event->get_particle(2, i, 0)->mass;
            }
            return values;
        },
        200, 2, 3, "Mass of the recreated particle",
        true);
    
    c31->SaveAs((filename + "_reco2.pdf").c_str());

    std::cout << "Finished analyzing the second iteration of recreated particles." << std::endl;
}

int main()
{
    const std::string part_type = "pion";
    EventCollector evc(
//             "/eos/cms/store/group/phys_diffraction/CMSTotemLowPU2018/ntuples/data/TOTEM*.root?#tree"
//               "/eos/user/y/yelberke/TOTEM_2018_ADDEDVARS_OUT/minimal/TOTEM*.root?#tree"
                "/eos/user/y/yelberke/TOTEM_2018_ADDEDVARS_OUT/indiv_partial2/TOTEM*.root?#tree"
               ,"/afs/cern.ch/user/p/ptuomola/private/particle_reconstruction_results.root");

    initialize_particles(evc, part_type);
    filter(evc);
//    analyze_data(evc, "histogram1");
    reconstruct(evc);
    analyze_reco1(evc, "histogram1", part_type);
//    analyze_reco1(evc, "histogram2", part_type);
    reconstruct(evc);
    analyze_reco2(evc, "histogram1");

    return 0;
}