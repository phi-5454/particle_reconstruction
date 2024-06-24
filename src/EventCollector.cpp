#include "EventCollector.h"
#include "TFile.h"
#include "TChain.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TF1.h"
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

    std::cout << "Initializing events" << std::endl;


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

template <typename F>
TH1F* EventCollector::create_1Dhistogram(F&& lambda, int bins, float low, float high, std::string title, bool draw)
{
    TH1F* hist = new TH1F("h1", title.c_str(), bins, low, high);
    for (Event* &event : events)
    {
        hist->Fill(lambda(event));
    }
    if (draw) hist->Draw();
    return hist;
}

template <typename F>
TH1F* EventCollector::create_1Dhistogram(F&& lambda, bool draw)
{
    float min = 0;
    float max = 1;
    for (Event* &event : events)
    {
        if (lambda(event) < min) min = lambda(event);
        if (lambda(event) > max) max = lambda(event);
    }
    return create_1Dhistogram(lambda, round(sqrt(events.size() / 2)), min, max, "A histogram for a fit", draw);
}

template <typename F>
TF1* EventCollector::create_1Dhistogram_fit(F&& lambda, int bins, float low, float high, std::string title, std::string distr)
{
    TH1F* h1 = create_1Dhistogram(lambda, bins, low, high, title, false);
    h1->Fit(distr.c_str());
    return h1->GetFunction("gaus");
}

template <typename F>
TF1* EventCollector::create_1Dhistogram_fit(F&& lambda, std::string distr)
{
    TH1F* h1 = create_1Dhistogram(lambda, false);
    h1->Fit(distr.c_str());
    return h1->GetFunction("gaus");
}

template <typename F>
void EventCollector::filter_events(F&& lambda)
{
    auto helper = std::vector<Event*>(events.size());
    auto it = std::copy_if(events.begin(), events.end(), helper.begin(), lambda);
    helper.resize(it - helper.begin());
    events = helper;
}

template <typename F>
void EventCollector::filter_events_distribution(F&& lambda, std::string distr, float sigmaMulti, int bins, float low, float high, std::string title, bool isPart)
{
    TF1* fit = create_1Dhistogram_fit(lambda, bins, low, high, title, distr);
    float mean = fit->GetParameter("Mean");
    float sigma = fit->GetParameter("Sigma");
    filter_events([&](Event* e) {
        if (isPart)
        {
            for (int i = 0; i < 4; ++i) {if (lambda(e) < mean - sigmaMulti * sigma || lambda(e) > mean + sigmaMulti * sigma) return false;} return true;
        }
        return (lambda(e) > mean - sigmaMulti * sigma && lambda(e) < mean + sigmaMulti * sigma);});
}

template <typename F>
void EventCollector::filter_events_distribution(F&& lambda, std::string distr, float sigmaMulti, bool isPart)
{
    TF1* fit = create_1Dhistogram_fit(lambda, distr);
    float mean = fit->GetParameter("Mean");
    float sigma = fit->GetParameter("Sigma");
    filter_events([&](Event* e) {
        if (isPart)
        {
            for (int i = 0; i < 4; ++i) {if (lambda(e) < mean - sigmaMulti * sigma || lambda(e) > mean + sigmaMulti * sigma) return false;} return true;
        }
        return (lambda(e) > mean - sigmaMulti * sigma && lambda(e) < mean + sigmaMulti * sigma);});
}

void EventCollector::filter_initial_events()
{
    std::cout << "Filtering events" << std::endl;
    filter_events([](Event* e) {return e->ntracks == 4;});
    filter_events([](Event* e) {int j = 0; for (int i = 0; i < 4; ++i) { j+= e->get_particle(0, i)->q;} return j == 0;});

    std::cout << "Finished filtering" << std::endl;
}

void EventCollector::analyze()
{
    std::cout << "Analyzing events" << std::endl;
    TCanvas* c1 = new TCanvas("c1", "c1"); 
    c1->Draw();
    TFile* results = TFile::Open(this->results.c_str(), "RECREATE");
    filter_events_distribution([](Event* event) {return event->zPV;}, "gaus", 1, false);
    TH1* h1 = create_1Dhistogram([](Event* event) {return event->zPV;}, true);
//    h1->Write();
    c1->SaveAs("hist1.pdf");
    results->Close();

    std::cout << "Finished analyzing" << std::endl;
}