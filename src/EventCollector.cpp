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
  TTreeReaderValue<unsigned long long> eventN(myReader, "EventNum");
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

  while (myReader.Next()) {
    Event *ev;
    if (isNew) {
//      ev = new Event(*ntrk, *zPV, *xPV, *yPV, *eventN);
    }
    else {
      ev = new Event(*ntrk, *zPV, *eventN);
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

void EventCollector::init_masses_and_energy(double mass) {
  for (Event *&event : events)
    for (int i = 0; i < event->ntracks; ++i) {
      event->set_masses_and_energies(mass);
    }
}

void EventCollector::reconstruct_particles() {
  for (Event *&event : events)
    event->reconstruct();
}
