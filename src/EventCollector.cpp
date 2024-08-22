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

void EventCollector::initialize_events(bool new_ntuples, bool new_protons) {
    // Add all the TTrees from the files to the TChain.
    TChain *chain_part = new TChain("tree_part");
    chain_part->Add(this->filepath.c_str());

    // Initialize the values from the TTrees to be read.
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

    TTreeReaderValue<float> *xPV{nullptr};
    TTreeReaderValue<float> *yPV{nullptr};
    TTreeReaderArray<float> *dxyErr{nullptr};
    TTreeReaderArray<float> *dzErr{nullptr};
    TTreeReaderArray<float> *ptErr{nullptr};

    if(new_ntuples){
        xPV   = new TTreeReaderValue<float>{myReader, "xPV"};
        yPV   = new TTreeReaderValue<float>{myReader, "yPV"};
        dxyErr= new TTreeReaderArray<float>{myReader, "trk_dxyerr"};
        dzErr = new TTreeReaderArray<float>{myReader, "trk_dzerr"};
        ptErr = new TTreeReaderArray<float>{myReader, "trk_pterr"};
    }

    TTreeReaderValue<double> *pr_px_a{nullptr};
    TTreeReaderValue<double> *pr_px_b{nullptr};
    TTreeReaderValue<double> *pr_py_a{nullptr};
    TTreeReaderValue<double> *pr_py_b{nullptr};
    TTreeReaderValue<double> *pr_pz_a{nullptr};
    TTreeReaderValue<double> *pr_pz_b{nullptr};
    TTreeReaderValue<double> *pr_ptx_a{nullptr};
    TTreeReaderValue<double> *pr_ptx_b{nullptr};
    TTreeReaderValue<double> *pr_pty_a{nullptr};
    TTreeReaderValue<double> *pr_pty_b{nullptr};
    TTreeReaderValue<double> *pr_ptx_sigma_a{nullptr};
    TTreeReaderValue<double> *pr_ptx_sigma_b{nullptr};
    TTreeReaderValue<double> *pr_pty_sigma_a{nullptr};
    TTreeReaderValue<double> *pr_pty_sigma_b{nullptr};
    TTreeReaderValue<double> *pr_posx_a{nullptr};
    TTreeReaderValue<double> *pr_posx_b{nullptr};
    TTreeReaderValue<double> *pr_posy_a{nullptr};
    TTreeReaderValue<double> *pr_posy_b{nullptr};
    TTreeReaderValue<double> *pr_posx_sigma_a{nullptr};
    TTreeReaderValue<double> *pr_posx_sigma_b{nullptr};
    TTreeReaderValue<double> *pr_posy_sigma_a{nullptr};
    TTreeReaderValue<double> *pr_posy_sigma_b{nullptr};

    if(new_protons){
        pr_px_a = new TTreeReaderValue<double>(myReader, "pr_px_a");
        pr_px_b = new TTreeReaderValue<double>(myReader, "pr_px_b");
        pr_py_a = new TTreeReaderValue<double>(myReader, "pr_py_a");
        pr_py_b = new TTreeReaderValue<double>(myReader, "pr_py_b");
        pr_pz_a = new TTreeReaderValue<double>(myReader, "pr_pz_a");
        pr_pz_b = new TTreeReaderValue<double>(myReader, "pr_pz_b");
        pr_ptx_a = new TTreeReaderValue<double>(myReader, "pr_ptx_a");
        pr_ptx_b = new TTreeReaderValue<double>(myReader, "pr_ptx_b");
        pr_pty_a = new TTreeReaderValue<double>(myReader, "pr_pty_a");
        pr_pty_b = new TTreeReaderValue<double>(myReader, "pr_pty_b");
        pr_ptx_sigma_a = new TTreeReaderValue<double>(myReader, "pr_ptx_sigma_a");
        pr_ptx_sigma_b = new TTreeReaderValue<double>(myReader, "pr_ptx_sigma_b");
        pr_pty_sigma_a = new TTreeReaderValue<double>(myReader, "pr_pty_sigma_a");
        pr_pty_sigma_b = new TTreeReaderValue<double>(myReader, "pr_pty_sigma_b");
        pr_posx_a = new TTreeReaderValue<double>(myReader, "pr_posx_a");
        pr_posx_b = new TTreeReaderValue<double>(myReader, "pr_posx_b");
        pr_posy_a = new TTreeReaderValue<double>(myReader, "pr_posy_a");
        pr_posy_b = new TTreeReaderValue<double>(myReader, "pr_posy_b");
        pr_posx_sigma_a = new TTreeReaderValue<double>(myReader, "pr_posx_sigma_a");
        pr_posx_sigma_b = new TTreeReaderValue<double>(myReader, "pr_posx_sigma_b");
        pr_posy_sigma_a = new TTreeReaderValue<double>(myReader, "pr_posy_sigma_a");
        pr_posy_sigma_b = new TTreeReaderValue<double>(myReader, "pr_posy_sigma_b");
    }

    // Read one event at a time and create required objects for it.
    while (myReader.Next()) {
        // Create the event
        Event *ev;
        if (new_ntuples) {
            ev = new Event(*ntrk, *zPV, **xPV, **yPV, *eventN);
        }
        else {
            ev = new Event(*ntrk, *zPV, *eventN);
        }
        events.push_back(ev);

        // Create the particles
        ev->particles.push_back(std::vector<std::vector<Particle *>>{});
        ev->particles[0].push_back(std::vector<Particle *>{});
        if (new_ntuples) {
            for (int i = 0; i < *ntrk; ++i) {
                ev->add_particle(p[i], pt[i], eta[i], phi[i], q[i], dxy[i], dz[i], (*ptErr)[i], (*dxyErr)[i], (*dzErr)[i], 0, 0);
            }
        }
        else {
            for (int i = 0; i < *ntrk; ++i) {
                ev->add_particle(p[i], pt[i], eta[i], phi[i], q[i], dxy[i], dz[i], 0, 0, 0, 0, 0);
            }
        }

        // Create the protons
        if(new_protons) {
            ev->add_proton(*ThxR, *ThyR, **pr_px_b, **pr_py_b, **pr_pz_b, **pr_ptx_b, **pr_pty_b, **pr_ptx_sigma_b,
                           **pr_pty_sigma_b, **pr_posx_b, **pr_posy_b, **pr_posx_sigma_b, **pr_posy_sigma_b);
            ev->add_proton(*ThxL, *ThyL, **pr_px_a, **pr_py_a, **pr_pz_a, **pr_ptx_a, **pr_pty_a, **pr_ptx_sigma_a,
                           **pr_pty_sigma_a, **pr_posx_a, **pr_posy_a, **pr_posx_sigma_a, **pr_posy_sigma_a);
        }
        else{
            ev->add_proton(*ThxR, *ThyR);
            ev->add_proton(*ThxL, *ThyL);
        }
    }
}

void EventCollector::init_masses_and_energy(double mass) {
    for (Event *&event : events)
        for (int i = 0; i < event->ntracks; ++i) {
            event->set_masses_and_energies(mass);
        }
}

void EventCollector::reconstruct_particles(bool useMCCoupling) {
    for (Event *&event : events)
        event->reconstruct(useMCCoupling);
}
