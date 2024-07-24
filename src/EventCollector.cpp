#include "EventCollector.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include <string>

EventCollector::EventCollector(std::string in, std::string out) {
  this->filepath = in;
  this->results = out;
}

void EventCollector::initialize_events(bool isNew) {
  TChain *chain_part = new TChain("tree_part");
  TChain *chain_prot = new TChain("tree_prot");

  chain_part->Add(this->filepath.c_str());
  chain_prot->Add("/eos/user/y/yelberke/TOTEM_2018_ADDEDVARS_OUT/improved_protons_small/TOTEM20.root?#tree");
  chain_prot->BuildIndex("EventNum");
  chain_part->AddFriend(chain_prot);

  TTreeReader myReader(chain_part);
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
/*
  TTreeReaderValue<double> PtxR(myReader, "pr_ptx_a");
  TTreeReaderValue<double> PtxL(myReader, "pr_ptx_b");
  TTreeReaderValue<double> PtyR(myReader, "pr_pty_a");
  TTreeReaderValue<double> PtyL(myReader, "pr_pty_b");
*/
  TTreeReaderValue<float> xPV(myReader, "xPV");
  TTreeReaderValue<float> yPV(myReader, "yPV");
  TTreeReaderArray<float> dxyErr(myReader, "trk_dxyerr");
  TTreeReaderArray<float> dzErr(myReader, "trk_dzerr");
  TTreeReaderArray<float> ptErr(myReader, "trk_pterr");

  while (myReader.Next()) {
    Event *ev;
    if (isNew) {
      ev = new Event(*ntrk, *zPV, *xPV, *yPV, *eventN);
    }
    else {
      ev = new Event(*ntrk, *zPV, *eventN);
    }
    events.push_back(ev);
    ev->particles.push_back(std::vector<std::vector<Particle *>>{});
    ev->particles[0].push_back(std::vector<Particle *>{});
    if (isNew) {
      for (int i = 0; i < *ntrk; ++i) {
        ev->add_particle(p[i], pt[i], eta[i], phi[i], q[i], dxy[i], dz[i], ptErr[i], dxyErr[i], dzErr[i], 0, 0);
      }
    }
     else {
      for (int i = 0; i < *ntrk; ++i) {
        ev->add_particle(p[i], pt[i], eta[i], phi[i], q[i], dxy[i], dz[i], 0, 0, 0, 0, 0);
      }
    }
//    ev->add_proton(*ThxR, *ThyR, *PtxR, *PtyR);
//    ev->add_proton(*ThxL, *ThyL, *PtxL, *PtyL);
    ev->add_proton(*ThxR, *ThyR);
    ev->add_proton(*ThxL, *ThyL);
  }
}

void EventCollector::initialize_protons() {
  TChain *chain = new TChain("hugetree");
  chain->Add("/eos/user/y/yelberke/TOTEM_2018_ADDEDVARS_OUT/improved_protons/TOTEM*.root?#tree");

  TTreeReader myReader(chain);
  TTreeReaderValue<unsigned long long> EV(myReader, "EventNum");
  TTreeReaderValue<double> ThxR(myReader, "pr_ptx_a");
  TTreeReaderValue<double> ThxL(myReader, "pr_ptx_b");
  TTreeReaderValue<double> ThyR(myReader, "pr_pty_a");
  TTreeReaderValue<double> ThyL(myReader, "pr_pty_b");

  for (Event* &event : events) {
    while (myReader.Next()) {
      if (event->EventNum == *EV) {
        event->get_proton(0)->px = *ThxR;
        event->get_proton(0)->py = *ThyR;
        event->get_proton(1)->px = *ThxL;
        event->get_proton(1)->py = *ThyL;
        break;
      }
    }
  }
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
