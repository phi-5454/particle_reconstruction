#include "EventCollector.h"
#include "TFile.h"
#include "TChain.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TH1.h"
#include "TCanvas.h"
#include <string>

EventCollector::EventCollector()
{
    // empty
}

void EventCollector::initialize_events()
{
    
    std::string infile = this->filepath + "TOTEM20.root";
//    std::string files = this->filepath + "TOTEM2*.root?#tree";
    std::string files = this->filepath + "TOTEM20.root?#tree";

    TFile* h = TFile::Open(infile.c_str());

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
        ev->particles.push_back(std::vector<Particle*>{});
        for (int i = 0; i < *ntrk; ++i)
        {
            ev->add_particle(p[i], pt[i], eta[i], phi[i], q[i], dxy[i], dz[i], 0);
        }
    }
 
    std::cout << "Finished initializing events" << std::endl;

    h->Close();
}

void EventCollector::filter_initial_events()
{
    std::cout << "Starting filtering" << std::endl;
    
    auto helper = std::vector<Event*>(events.size());
    auto it = std::copy_if(events.begin(), events.end(), helper.begin(), [](Event* e) {return e->ntracks == 4;});
    size_t len = it - helper.begin();
    helper.resize(len);

    int rem = 0;
/*    for (int i = len - 1; i > -1; --i)
    {
        int charge = 0;
        for (int j = 0; i < 4; ++i)
        {
            charge += events[i]->get_particle(0, j)->q;
        }
        if (charge != 0)
        {
            events.erase(events.begin() + i);
            ++rem;
        }
    }
*/
    helper.resize(len - rem);
    events = helper;
    std::cout << "Finished filtering" << std::endl;
}

void EventCollector::analyze()
{
    std::cout << "Starting analyzing" << std::endl;
    TFile* results = TFile::Open(this->results.c_str(), "RECREATE");
    TH1F* h1 = new TH1F("h1", "h1", 400, -20, 20);
    TCanvas* c1 = new TCanvas("c1", "c1");

    for (Event* &event : events)
    {
        h1->Fill(event->zPV);
    }

    h1->Write();
    h1->Draw();
    results->Close();

    std::cout << "Finished analyzing" << std::endl;
}