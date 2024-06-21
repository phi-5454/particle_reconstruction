#include "EventCollector.h"
#include "TFile.h"
#include "TChain.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include <string>

EventCollector::EventCollector()
{
    // empty
}

void EventCollector::initialize_events()
{
    std::string files = this->filepath + "TOTEM2*.root?#tree";

    TFile* h = TFile::Open(this->filepath.c_str());
    TChain* chain = new TChain("hugetree");
    chain->Add(files.c_str());

    TTreeReader myReader(chain);
    TTreeReaderValue<Float_t> zPV(myReader, "zPV");
    TTreeReaderValue<Int_t> ntrk(myReader, "ntrk");
    TTreeReaderArray<Float_t> p(myReader, "trk_p");
    TTreeReaderArray<Float_t> pt(myReader, "trk_pt");
    TTreeReaderArray<Float_t> eta(myReader, "trk_eta");
    TTreeReaderArray<Float_t> phi(myReader, "trk_phi");
    TTreeReaderArray<Int_t> q(myReader, "trk_q");
    TTreeReaderArray<Float_t> dxy(myReader, "trk_dxy");
    TTreeReaderArray<Float_t> dz(myReader, "trk_dz");
    TTreeReaderArray<Float_t> ThxR(myReader, "ThxR");
    TTreeReaderArray<Float_t> ThxL(myReader, "ThxL");
    TTreeReaderArray<Float_t> ThyR(myReader, "ThyR");
    TTreeReaderArray<Float_t> ThyL(myReader, "ThyL");

    std::cout << "Starting to initialize events" << std::endl;

    while(myReader.Next())
    {
        Event* ev = new Event(*ntrk, *zPV);
        this->events.push_back(ev);

        for (int i = 0; i < *ntrk; ++i)
        {
            ev->add_particle(p[i], pt[i], eta[i], phi[i], q[i], dxy[i], dz[i]);
        }
    }

    std::cout << "Finished initializing events" << std::endl;
}

void EventCollector::filter_initial_events()
{
    for (Event* &event : events)
    {
        if (event->ntracks != 4) events.erase(std::remove(events.begin(), events.end(), event));
    }
}