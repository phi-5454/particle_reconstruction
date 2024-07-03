#include "EventCollector.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include <string>

EventCollector::EventCollector(std::string in, std::string out) {
  this->filepath = in;
  this->results = out;
}

void EventCollector::initialize_events(bool isNew) {
  TChain *chain = new TChain("hugetree");

  chain->Add(this->filepath.c_str());

  TTreeReader myReader(chain);
  TTreeReaderValue<float> zPV(myReader, "zPV");
  TTreeReaderValue<Int_t> ntrk(myReader, "ntrk");
  TTreeReaderArray<float> p(myReader, "trk_p");
  TTreeReaderArray<float> pt(myReader, "trk_pt");
  TTreeReaderArray<float> eta(myReader, "trk_eta");
  TTreeReaderArray<float> phi(myReader, "trk_phi");
  TTreeReaderArray<Int_t> q(myReader, "trk_q");
  TTreeReaderArray<float> dxy(myReader, "trk_dxy");
  TTreeReaderArray<float> dz(myReader, "trk_dz");
  TTreeReaderValue<float> ThxR(myReader, "ThxR");
  TTreeReaderValue<float> ThxL(myReader, "ThxL");
  TTreeReaderValue<float> ThyR(myReader, "ThyR");
  TTreeReaderValue<float> ThyL(myReader, "ThyL");

/*  TTreeReaderValue<float> xPV(myReader, "xPV");
  TTreeReaderValue<float> yPV(myReader, "yPV");
  TTreeReaderArray<float> dxyErr(myReader, "trk_dxyerr");
  TTreeReaderArray<float> dzErr(myReader, "trk_dzerr");
  TTreeReaderArray<float> ptErr(myReader, "trk_pterr");
*/
  std::cout << "Initializing events" << std::endl;

  while (myReader.Next()) {
    Event *ev;
    if (isNew) {
//      ev = new Event(*ntrk, *zPV, *xPV, *yPV);
    }
    else {
      ev = new Event(*ntrk, *zPV);
    }
    events.push_back(ev);
    ev->particles.push_back(std::vector<std::vector<Particle *>>{});
    ev->particles[0].push_back(std::vector<Particle *>{});
    if (isNew) {
      for (int i = 0; i < *ntrk; ++i) {
//        ev->add_particle(p[i], pt[i], eta[i], phi[i], q[i], dxy[i], dz[i], ptErr[i], dxyErr[i], dzErr[i], 0, 0);
      }
    }
     else {
      for (int i = 0; i < *ntrk; ++i) {
        ev->add_particle(p[i], pt[i], eta[i], phi[i], q[i], dxy[i], dz[i], 0, 0, 0, 0, 0);
      }
    }
    ev->add_proton(*ThxR, *ThyR);
    ev->add_proton(*ThxL, *ThyL);
  }
  std::cout << "Finished initializing events" << std::endl;

  // h->Close();
}

/// 2D Histograms

template <typename F1, typename F2>
TH2F *EventCollector::create_2Dhistogram(F1 &&lambda_x, F2 &&lambda_y, int bins_x, double low_x, double high_x,
                             int bins_y, double low_y, double high_y,
                             const std::string& title, bool draw) {
    TH2F *hist = new TH2F("hist", title.c_str(), bins_x, low_x, high_x, bins_y, low_y, high_y);
    for (Event *&event : events) {
        std::vector<double> values_x = lambda_x(event);
        std::vector<double> values_y = lambda_y(event);
        // Ass. len(values_x) == len(values_y)
        for (int i = 0; i < values_x.size(); ++i) {
            hist->Fill(values_x[i], values_y[i]);
        }
    }
    if (draw)
        hist->Draw("Colz");
    return hist;
}


template <typename F1, typename F2>
TH2F *EventCollector::create_2Dhistogram(F1 &&lambda_x, F2 &&lambda_y, const std::string& title,  bool draw) {
    auto [min_x, max_x] = find_min_max(lambda_x);
    auto [min_y, max_y] = find_min_max(lambda_y);
    return create_2Dhistogram(lambda_x, lambda_y, round(pow(2 * events.size(), 0.5)), min_x, max_x,
                              round(pow(2 * events.size(), 0.5)), min_y, max_y,
                              title, draw);
}

void EventCollector::analyze(std::string filename) {
  std::cout << "Analyzing events" << std::endl;
  TFile *results = TFile::Open(this->results.c_str(), "RECREATE");

  TCanvas *c1 = new TCanvas("c1", "c1");
  c1->DivideSquare(4);
  c1->Draw();

  c1->cd(1);
  TH1 *h1 = create_1Dhistogram(
      [](Event *event) {
        std::vector<double> values = {event->zPV};
        return values;
      },
      200, -15, 15, "Primary vertex Z position", true);

  c1->cd(3);
  TH1 *h2 = create_1Dhistogram(
      [](Event *event) {
        std::vector<double> values(4);
        for (int i = 0; i < 4; ++i) {
          values[i] = event->get_particle(0, 0, i)->dxy;
        }
        return values;
      },
      300, -1, 1, "Particle distance from primary vertex in xy-plane",
      true);

  c1->cd(4);
  TH1 *h3 = create_1Dhistogram(
      [](Event *event) {
        std::vector<double> values(4);
        for (int i = 0; i < 4; ++i) {
          values[i] = event->get_particle(0, 0, i)->dz;
        }
        return values;
      },
      300, -1, 1, "Particle distance from primary vertex along z-axis",
      true);

  TCanvas *c2 = new TCanvas("c2", "c2");
  c2->DivideSquare(4);
  c2->Draw();

  c2->cd(1);
  auto l1 = [](Event *event) {
      double px = 0;
      for (int i = 0; i < 4; ++i) {
          px += event->get_particle(0, 0, i)->px;
      }
      return std::vector<double>{px};
  };
  auto l2 = [](Event *event) {
      double px = 0;
      for (int i = 0; i < 2; ++i) {
          px += event->get_proton(i)->px;
      }
      return std::vector<double>{px};
  };
  TH2 *h4 = create_2Dhistogram(
          l1, l2, 100, -1, 1, 100, -1, 1,
          "Track px vs proton px",
          true);

  c2->cd(2);
  auto l3 = [](Event *event) {
      double py = 0;
      for (int i = 0; i < 4; ++i) {
          py += event->get_particle(0, 0, i)->py;
      }
      return std::vector<double>{py};
  };
  auto l4 = [](Event *event) {
      double py = 0;
      for (int i = 0; i < 2; ++i) {
          py += event->get_proton(i)->py;
      }
      return std::vector<double>{py};
  };
  TH2 *h5 = create_2Dhistogram(
          l3, l4, 100, -1, 1, 100, -1, 1,
          "Track py vs proton py",
          true);

  c2->cd(4);
  auto l5 = [](Event *event) {
      std::vector<double> values(4);
      for (int i = 0; i < 4; ++i) {
          values[i] = event->get_particle(0, 0, i)->dz;
      }
      return values;
  };
  auto l6 = [](Event *event) {
      std::vector<double> values(4);
      for (int i = 0; i < 4; ++i) {
          values[i] = event->get_particle(0, 0, i)->eta;
      }
      return values;
  };
  TH2 *h6 = create_2Dhistogram(
          l5, l6, 200, -0.4, 1.1, 200, -3, 3,
          "Particle dz vs. eta",
          true);

  c2->cd(3);
  auto l7 = [](Event *event) {
      double px = 0;
      for (int i = 0; i < 2; ++i) {
          px += event->get_proton(i)->px;
      }
      return std::vector<double>{px};
  };
  auto l8 = [](Event *event) {
      double py = 0;
      for (int i = 0; i < 2; ++i) {
          py += event->get_proton(i)->py;
      }
      return std::vector<double>{py};
  };
  TH2 *h7 = create_2Dhistogram(
          l7, l8, 200, -1.5, 1.5, 100, -0.6, 0.6,
          "Proton px vs. py",
          true);

  TCanvas* c3 = new TCanvas("c3", "c3");
  c3->Draw();

  auto l9 = [](Event *event) {
    std::vector<double> values(4);
    for (int i = 0; i < 4; ++i) {
      values[i] = event->get_particle(0, 0, i)->dxy;
    }
    return values;
  };

  auto l10 = [](Event *event) {
    std::vector<double> values(4);
    for (int i = 0; i < 4; ++i) {
      values[i] = event->get_particle(0, 0, i)->phi;
    }
    return values;
  };

  TH2 *h8 = create_2Dhistogram(
            l9, l10, 200, -1, 1, 200, -3.2, 3.2,
            "Particle dxy vs. phi", true);
  
  h1->Write();
  h2->Write();
  h3->Write();
  h4->Write();
  h5->Write();
  h6->Write();
  h7->Write();
  c1->SaveAs((filename + "A.pdf").c_str());
  c2->SaveAs((filename + "B.pdf").c_str());
  c3->SaveAs((filename + "C.pdf").c_str());
  results->Close();

  std::cout << "Finished analyzing" << std::endl;
}

void EventCollector::analyze_new(std::string filename) {
  TFile *results = TFile::Open(this->results.c_str(), "RECREATE");

  TCanvas *c1 = new TCanvas("c1", "c1");
  c1->DivideSquare(4);
  c1->Draw();

  c1->cd(1);
  TH1 *h1 = create_1Dhistogram(
      [](Event *event) {
        std::vector<double> values = {event->xPV};
        return values;
      },
      20, 0.06, 0.14, "Primary vertex X position", true);

  c1->cd(2);
  TH1 *h2 = create_1Dhistogram(
      [](Event *event) {
        std::vector<double> values = {event->yPV};
        return values;
      },
      20, 0.08, 0.16, "Primary vertex Y position", true);

  c1->cd(3);
  TH1 *h3 = create_1Dhistogram(
    [](Event *event) {
      std::vector<double> values(4);
      for (int i = 0; i < 4; ++i) {
        values[i] = event->get_particle(0, 0, i)->dz / event->get_particle(0, 0, i)->dzErr;
      }
      return values;
    },
    50, -5, 5, "dz divided by dzErr", true);

  c1->cd(4);
  TH1 *h4 = create_1Dhistogram(
    [](Event *event) {
      std::vector<double> values(4);
      for (int i = 0; i < 4; ++i) {
        values[i] = event->get_particle(0, 0, i)->dxy / event->get_particle(0, 0, i)->dxyErr;
      }
      return values;
    },
    50, -5, 5, "dxy divided by dxyErr", true);

  TCanvas *c2 = new TCanvas("c2", "c2");
  c2->DivideSquare(4);
  c2->Draw();
  c2->cd(1);

  TH1 *h5 = create_1Dhistogram(
    [](Event *event) {
      std::vector<double> values(4);
      for (int i = 0; i < 4; ++i) {
        values[i] = event->get_particle(0, 0, i)->ptErr / event->get_particle(0, 0, i)->pt;
      }
      return values;
    },
    100, 0, 0.2, "ptErr divided by pt", true);

  c2->cd(2);
  TH1 *h6 = create_1Dhistogram(
      [](Event *event) {
        std::vector<double> values = {(double)sqrt(pow(event->xPV, 2) + pow(event->yPV, 2))};
        return values;
      },
      30, 0.12, 0.18, "Primary vertex R position", true);

  c1->SaveAs((filename + "A.pdf").c_str());
  c2->SaveAs((filename + "B.pdf").c_str());
  results->Close();
}

void EventCollector::init_masses_and_energy(double mass) {
  for (Event *&event : events)
    for (int i = 0; i < 4; ++i) {
      event->set_masses_and_energies(mass);
    }
}

void EventCollector::reconstruct_particles() {
  for (Event *&event : events) {
    event->reconstruct();
  }
}

void EventCollector::analyze_reco(std::string filename) {
  std::cout << "Analyzing recreated events" << std::endl;
  TFile *results = TFile::Open(this->results.c_str(), "RECREATE");

  TCanvas *c1 = new TCanvas("c1", "c1");

  c1->Draw();
/*  TH1 *h1 = create_1Dhistogram(
      [](Event *event) {
        std::vector<double> values(4);
        for (int i = 0; i < 2; ++i) {
          for (int j = 0; j < 2; ++j) {
            values[2 * i + j] = event->get_particle(1, i, j)->mass;
          }
        }
        return values;
      },
      200, 0.25, 1, "Mass of recreated particles, assumed pions",
      true);*/

  TH2 *h2 = create_2Dhistogram([](Event* event) {return std::vector<double>{event->get_particle(1, 0, 0)->mass, event->get_particle(1, 1, 1)->mass};},
                               [](Event* event) {return std::vector<double>{event->get_particle(1, 1, 0)->mass, event->get_particle(1, 0, 1)->mass};},
                               100, 1, 1.5, 100, 1, 1.5, "Masses of two recreated particles", true);
/*    TH1 *h2 = create_1Dhistogram([](Event* event) {return std::vector<double>{event->get_particle(2, 0, 0)->mass,
                                 event->get_particle(2, 1, 0)->mass};}, 50, 1.5, 2.5, "Glueball", true );*/
  h2->Write();

  TH1 *h3 = h2->ProjectionX();
  TH1 *h4 = h2->ProjectionY();
  h3->Add(h4);
  h3->Draw("E");
  c1->SaveAs((filename + ".pdf").c_str());
  results->Close();

  std::cout << "Finished analyzing" << std::endl;
}
