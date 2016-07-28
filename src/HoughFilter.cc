#include <cstdlib>
#include <fstream>
#include <vector>
#include <set>
#include <map>
#include <iostream>
#include <sstream>      
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TText.h>
#include <TLine.h>
#include <TLegend.h>
#include <TH2F.h>
#include "TROOT.h"
#include "TRint.h"
#include "../interface/Hough.h"
#include "../interface/Stub.h"
#include "../interface/Settings.h"
#include <vector>
#ifdef __MAKECINT__
#pragma link C++ class vector< std::vector<int> >+;
#endif

using namespace std;

bool insideEta(double r, double z){
    double chosenRofZ_ = 45.;
    double etaMin_ = -0.37;
    double etaMax_ = 0.37;
    double zOuterMin_ = chosenRofZ_ / tan( 2. * atan(exp(-etaMin_)) );
    double zOuterMax_ = chosenRofZ_ / tan( 2. * atan(exp(-etaMax_)) );
    double beamWindowZ_ = 15.;

    double zMin = ( zOuterMin_ * r - beamWindowZ_ * fabs(r - chosenRofZ_) ) / chosenRofZ_;
    // Calculate z coordinate of upper edge of this eta region, evaluated at radius of stub.
    double zMax = ( zOuterMax_ * r + beamWindowZ_ * fabs(r - chosenRofZ_) ) / chosenRofZ_;

    // zMin = ( zRangeMin * stub->r() - beamWindowZ_ * fabs(stub->r() - rOuterMin_) ) / rOuterMin_;
    // zMax = ( zRangeMax * stub->r() + beamWindowZ_ * fabs(stub->r() - rOuterMax_) ) / rOuterMax_;

    bool inside = (z > zMin && z < zMax);
    return inside;
}

bool insidePhi( double phi, double r ){

  // N.B. The logic here for preventing a stub being assigned to > 2 sectors seems overly agressive.
  // But attempts at improving it have failed ...

  bool okPhi    = true;

  double phiCentre_ = 2.*M_PI * (0.5 + 26) / 32 - M_PI;
  double invPtToDphi = 3.8112*(3.0E8/2.0E11);
  double sectorHalfWidth_ = M_PI / 32; // Sector half width excluding overlaps.

    float delPhi = phi - phiCentre_; // Phi difference between stub & sector in range -PI to +PI.
    float tolerancePhi = fabs(r)*invPtToDphi/3.;

    float outsidePhi = fabs(delPhi) - sectorHalfWidth_ - tolerancePhi; // If > 0, then stub is not compatible with being inside this sector. 
    if (outsidePhi > 0) okPhi = false;
  
    return okPhi;
}

double GetYmax(TH1F* h1, TH1F *h2){
    double y1 = h1->GetBinContent(h1->GetMaximumBin());
    double y2 = h2->GetBinContent(h2->GetMaximumBin());
    return max(y1,y2);
}

double GetYmin(TH1F* h1, TH1F *h2){
    double y1 = h1->GetBinContent(h1->GetMinimumBin());
    double y2 = h2->GetBinContent(h2->GetMinimumBin());
    return min(y1,y2);
}

void RatioPlot(TH1F* h1, TH1F* h2, string title, string name, int nevents){
    
    TH1F *h1copy = (TH1F*)h1->Clone("h1copy");
    TH1F *h2copy = (TH1F*)h2->Clone("h2copy");

    // Define the Canvas
    TCanvas *cratio = new TCanvas("cratio", "canvas", 1024, 800);
    // gStyle->SetTitleAlign(11); //Align to left-top

    // Upper plot will be in pad1
    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.305, 1, 1.0);
    pad1->SetBottomMargin(0.02); // Upper and lower plot are joined
    pad1->SetGridx();         // Vertical grid
    
    // if(logy)
    //  pad1->SetLogy();

    pad1->Draw();             // Draw the upper pad: pad1


    pad1->cd();               // pad1 becomes the current pad

    double ymax = GetYmax(h1,h2)*1.1;
    double ymin = GetYmin(h1,h2)*0.9;

    h1copy->SetStats(0);          // No statistics on upper plot
    h1copy->Draw("P");            // Draw h1copy
    h1copy->SetTitle("");
    h2copy->Draw("hist same");         // Draw h2copy on top of h1copy
    h1copy->SetMarkerStyle(20);

    h1copy->SetTitleSize(0.1);
    h1copy->GetXaxis()->SetLabelSize(0.);

    h1copy->GetYaxis()->SetLabelSize(15);
    h1copy->GetYaxis()->SetLabelFont(43);
    h1copy->GetYaxis()->SetTitle(title.c_str());
    h1copy->GetYaxis()->SetRangeUser(ymin,ymax);

    TLegend *legratio = new TLegend(0.7,0.905,0.9,1);
    legratio->SetNColumns(2);
    legratio->SetFillColor(0);
    legratio->SetBorderSize(0);
    // legratio->SetTextAlign(21);
    legratio->SetTextSize(0.035);
    legratio->SetTextFont(42);
    legratio->AddEntry(h1copy, "pre-filtered", "p");
    legratio->AddEntry(h2copy, "filtered", "p");

    legratio->Draw();

    ostringstream titolo;
    titolo.str("");
    titolo << "AM HoughFilter TTbar+PU140 " <<nevents << " Events"; 
    
    TText *t = new TText(.1,.92,titolo.str().c_str() );
    t->SetNDC();
    t->SetTextFont(103);
    t->SetTextSize(24);
    t->Draw();

    // lower plot will be in pad
    cratio->cd();          // Go back to the main canvas before defining pad2
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.3);
    pad2->SetGridx(); // vertical grid
    pad2->Draw();
    pad2->cd();       // pad2 becomes the current pad

    // Define the ratio plot
    TH1F *h3 = (TH1F*)h1copy->Clone("h3");
    h3->SetLineColor(kBlack);
    h3->SetMinimum(0.8);  // Define Y ..
    h3->SetMaximum(1.25); // .. range
    h3->Sumw2();
    h3->SetStats(0);      // No statistics on lower plot
    h3->Divide(h2copy);
    h3->SetMarkerStyle(21);
    h2->SetLineColor(kBlue);
    h3->Draw("p");       // Draw the ratio plot

    TLine* unity = new TLine(h3->GetXaxis()->GetBinLowEdge(1),1,h3->GetXaxis()->GetBinUpEdge(h3->GetNbinsX()),1);
    cout << h3->GetMinimum() << ", " << h3->GetMaximum() << endl;
    unity->SetLineColor(kRed);
    unity->SetLineWidth(2);
    unity->Draw();
    h3->Draw("p same");

    // h1copy settings
    h2copy->SetLineColor(kBlue+1);
    h2copy->SetLineWidth(2);

    // Y axis h1copy plot settings
    h1copy->GetYaxis()->SetTitleSize(20);
    h1copy->GetYaxis()->SetTitleFont(43);
    h1copy->GetYaxis()->SetTitleOffset(1.45);

    // h2copy settings
    h1copy->SetLineColor(kBlack);
    // h2copy->SetLineWidth(2);

    // Ratio plot (h3) settings
    h3->SetTitle(""); // Remove the ratio title

    // Y axis ratio plot settings
    h3->GetYaxis()->SetTitle("ratio filtered/not-filtered");

    h3->GetYaxis()->SetNdivisions(505);
    h3->GetYaxis()->SetTitleSize(20);
    h3->GetYaxis()->SetTitleFont(43);
    h3->GetYaxis()->SetTitleOffset(1.45);
    h3->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    h3->GetYaxis()->SetLabelSize(15);
    h3->GetYaxis()->CenterTitle();
    // X axis ratio plot settings
    h3->GetXaxis()->SetTitleSize(20);
    h3->GetXaxis()->SetTitleFont(43);
    h3->GetXaxis()->SetTitleOffset(4.);
    h3->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    h3->GetXaxis()->SetLabelSize(15);
    h3->SetLineColor(kBlue);

    cratio->Print(name.c_str());
    cout << "Saving file: "<< name.c_str() << endl;

    cratio->Close();
}



int main(int argc, char* argv[]){
  gROOT->ProcessLine("#include <vector>");

  if ( argc<2 ){
    cout << "\n*** No arguments provided! HoughFilter will run with the default configuration.\n" << endl;
    }

    string input_file = "extracted.root";
    string output_file = "output.root";
    string res_folder = "resolutions";
    long long int n_events = -1;
    unsigned htType = 0;
    unsigned xbins = 10;
    unsigned ybins = 10;
    unsigned averageRes = 0;
    unsigned debug = 0;

    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i],"-maxEvents")) {
            n_events = atoi(argv[i+1]);
        }
        if (!strcmp(argv[i],"-inputFile")) {
            input_file = argv[i+1];
        }
        if (!strcmp(argv[i],"-outputFile")) {
            output_file = argv[i+1];
        }
        if (!strcmp(argv[i],"-HTtype")) {
            htType = atoi(argv[i+1]);
            if(htType > 1){
                cout << "Value of HTtype out of range (max 1)" << endl;
                return EXIT_SUCCESS;
            }
        }
        if (!strcmp(argv[i],"-Xbins")) {
            xbins = atoi(argv[i+1]);
        }
        if (!strcmp(argv[i],"-Ybins")) {
            ybins = atoi(argv[i+1]);
        }
        if (!strcmp(argv[i],"-averageRes")) {
            averageRes = atoi(argv[i+1]);
            if(averageRes > 1){
                cout << "Value of averageRes out of range (max 1)" << endl;
            return EXIT_SUCCESS;
            }
        }
        if (!strcmp(argv[i],"-resFolder")) {
            res_folder = argv[i+1];
        }
        if (!strcmp(argv[i],"-debug")) {
            debug = atoi(argv[i+1]);
        }
        if (!strcmp(argv[i],"-h")) {
            cout << "-h \t print this help and exit" << endl 
            << "-maxEvents # \t set number of events (default = -1)" << endl 
            << "-inputFile <filename> \t input root file containing extracted information (default = extracted.root)" << endl
            << "-outputFile <filename> \t output root file containing the resulting histograms" << endl
            << "-HTtype # \t type of Hough Transform to be processed (0: rphi, 1: rz) (default: 0)"<< endl
            << "-Xbins # \t number of bins of HT array along x axis (default 10)" << endl
            << "-Ybins # \t number of bins of HT array along y axis (default 10)" << endl
            << "-averageRes # \t if set to 1 use the average pattern resolution instead of the individual road one (default 0)" << endl
            << "-resFolder <path> \t path to the folder containing the parameter resolution values (default: resolution)" << endl
            << "-debug # \t verbosity level (default = 0):" << endl
            << "         0: no printouts" << endl
            << "         1: print some infos" << endl;
            return EXIT_SUCCESS;
        }
    }

    Settings* settings = new Settings(htType, xbins, ybins, averageRes, debug);

	// =============================================================================================
	//  Get the trees from the root files
	// =============================================================================================

    TChain *m_PATT = new TChain("L1tracks");  // Tree containing the L1 track reco info
	TChain *m_STUBS = new TChain("TkStubs");   // Tree containing the stub info 
	TChain *m_MC   = new TChain("MC");        // Tree containing the true MC info 

    // Case 1: it's a root file
    m_PATT->Add(input_file.c_str());
    m_STUBS->Add(input_file.c_str());
    m_MC->Add(input_file.c_str());

    if(n_events == -1)
        n_events = m_STUBS->GetEntries();

    n_events = std::min(m_STUBS->GetEntries(), n_events);

    int eventId = -1;
	int n_patterns = 0;  // Number of patterns in event


    // --- Book the histogram arrays:
    TH1F *hPatterns = new TH1F("hPatterns","Number of Patterns per event",150,0,150);
    TH1F *hStubsPerPattern = new TH1F("hStubsPerPattern","Number of Stubs per Pattern",64,0,64);
    TH1F *hFilteredPatterns = new TH1F("hFilteredPatterns","Number of Filtered Patterns per event",150,0,150);
    TH1F *hStubsPerFiltPattern = new TH1F("hStubsPerFiltPattern","Number of Stubs per Filtered Pattern",64,0,64);
    TH1F *hRecoTracks = new TH1F("hRecoTracks", "Number of reconstructed tracks per event", 50, 0, 50);
    TH1F *hFiltRecoTracks = new TH1F("hFiltRecoTracks", "Number of filtered tracks per event", 50, 0, 50);
    TH1F *hComb = new TH1F("hComb", "Number of Combinations per pattern", 100, 0, 1000);
    TH1F *hFiltComb = new TH1F("hFiltComb", "Number of Combinations per filtered pattern", 100, 0, 1000);
    TH1F *hSpecComb = new TH1F("hSpecComb", "Number of Special Combinations per pattern", 100, 0, 1000);
    TH1F *hFiltSpecComb = new TH1F("hFiltSpecComb", "Number of Special Combinations per filtered pattern", 100, 0, 1000);
    TGraph *gComb = new TGraph();
    gComb->SetTitle("Number of Combinations vs. number of filtered combinations per pattern");
    gComb->SetName("gComb");
    TGraph *gCombRatio = new TGraph();
    gCombRatio->SetTitle("Number of Combinations per pattern vs. ratio filtered combinations/number of combinations"); 

    // m_PATT->SetBranchAddress("L1TRK_n", &n_tracks);
    // m_PATT->SetBranchAddress("L1TRK_links", &track_links);

    // m_PATT->SetBranchAddress("L1TRK_secid", &track_secId);
    // m_PATT->SetBranchAddress("L1TRK_pt", &track_pT);
    // m_PATT->SetBranchAddress("L1TRK_phi", &track_phi);
    // m_PATT->SetBranchAddress("L1TRK_z", &track_z);
    // m_PATT->SetBranchAddress("L1TRK_eta", &track_eta);

    int n_pu = 0; 
    std::vector<int>     *pu_pdg  = NULL;
    std::vector<float>   *pu_px   = NULL;
    std::vector<float>   *pu_py   = NULL;
    std::vector<float>   *pu_pz   = NULL;
    std::vector<float>   *pu_pt   = NULL;

    std::vector<float>   *pu_eta  = NULL;
    std::vector<float>   *pu_theta  = NULL;

    std::vector<float>   *pu_phi  = NULL;
    std::vector<float>   *pu_vx   = NULL;
    std::vector<float>   *pu_vy   = NULL;
    std::vector<float>   *pu_vz   = NULL;
    std::vector<int>   *pu_charge   = NULL;

    int n_gen = 0; 
    std::vector<int>     *gen_pdg = NULL;
    std::vector<float>   *gen_px  = NULL;
    std::vector<float>   *gen_py  = NULL;
    std::vector<float>   *gen_pz  = NULL;
    std::vector<float>   *gen_vx  = NULL;
    std::vector<float>   *gen_vy  = NULL;
    std::vector<float>   *gen_vz  = NULL;

    m_MC->SetBranchAddress("gen_n",   &n_gen);
    m_MC->SetBranchAddress("gen_pdg", &gen_pdg);
    m_MC->SetBranchAddress("gen_px",  &gen_px);
    m_MC->SetBranchAddress("gen_py",  &gen_py);
    m_MC->SetBranchAddress("gen_pz",  &gen_pz);
    m_MC->SetBranchAddress("gen_x",   &gen_vx);
    m_MC->SetBranchAddress("gen_y",   &gen_vy);
    m_MC->SetBranchAddress("gen_z",   &gen_vz);  
    m_MC->SetBranchAddress("subpart_n",     &n_pu);
    m_MC->SetBranchAddress("subpart_pdgId", &pu_pdg);
    m_MC->SetBranchAddress("subpart_px",    &pu_px);
    m_MC->SetBranchAddress("subpart_py",    &pu_py);
    m_MC->SetBranchAddress("subpart_pt",    &pu_pt);

    m_MC->SetBranchAddress("subpart_pz",    &pu_pz);
    m_MC->SetBranchAddress("subpart_eta",   &pu_eta);
    m_MC->SetBranchAddress("subpart_theta",   &pu_theta);

    m_MC->SetBranchAddress("subpart_phi",   &pu_phi);
    m_MC->SetBranchAddress("subpart_x",     &pu_vx);
    m_MC->SetBranchAddress("subpart_y",     &pu_vy);
    m_MC->SetBranchAddress("subpart_z",     &pu_vz);
    m_MC->SetBranchAddress("subpart_charge",     &pu_charge);

    int n_stubs = 0;
    std::vector<int> *stub_layer = NULL;
    std::vector<int> *stub_tp = NULL;
    std::vector<double> *stub_phi = NULL;
    std::vector<double> *stub_r = NULL;
    std::vector<double> *stub_z = NULL;

    m_STUBS->SetBranchAddress("L1TkSTUB_n",         &n_stubs);
    m_STUBS->SetBranchAddress("L1TkSTUB_layer",     &stub_layer);
    m_STUBS->SetBranchAddress("L1TkSTUB_tp", &stub_tp);
    m_STUBS->SetBranchAddress("L1TkSTUB_phi", &stub_phi);
    m_STUBS->SetBranchAddress("L1TkSTUB_r", &stub_r);
    m_STUBS->SetBranchAddress("L1TkSTUB_z", &stub_z);


    std::vector< std::vector<int> > *patt_links = NULL; // Links to the stub ids of the pattern
    std::vector<int>                *patt_secId = NULL; // Link to the sector ids of the pattern
    std::vector<std::string>    *patt_id = NULL; // Link to the unique pattern id
    std::vector<int>                *patt_nmiss = NULL;

    // std::vector< std::vector<int> > *track_links = NULL;
    // std::vector<int>                 *track_secId  = NULL; 
    // std::vector<double>          *track_pT  = NULL; 
    // std::vector<double>          *track_phi  = NULL; 
    // std::vector<double>          *track_z  = NULL; 
    // std::vector<double>          *track_eta  = NULL; 

    m_PATT->SetBranchAddress("L1evt",            &eventId); 
    m_PATT->SetBranchAddress("L1PATT_n",     &n_patterns);
    m_PATT->SetBranchAddress("L1PATT_links", &patt_links);
    m_PATT->SetBranchAddress("L1PATT_secid", &patt_secId);
    m_PATT->SetBranchAddress("L1PATT_id", &patt_id);
    m_PATT->SetBranchAddress("L1PATT_nmiss", &patt_nmiss);

    // Maps containing the pattern windows 
    map<std::string, std::pair<double, double> > phiLimits;
    map<std::string, std::pair<double, double> > phi65Limits;
    map<std::string, std::pair<double, double> > etaLimits;
    map<std::string, std::pair<double, double> > z0Limits;
    map<std::string, std::pair<double, double> > ptLimits;
    map<std::string, std::pair<double, double> > thetaLimits;
    map<std::string, double> roadCharge;
    // Open files containing pattern resolutions
    ostringstream filename;
    filename.str("");
    filename.clear();
    filename << res_folder << "/z0resolution.txt";

    std::ifstream ifs_z0(filename.str().c_str(), std::ifstream::in);
    if(ifs_z0.good()==false){
        cout << "ERROR: file "<< filename.str() << " not found!" <<endl;
        return EXIT_SUCCESS;
    }
    filename.str("");
    filename.clear();
    filename << res_folder << "/phiresolution.txt";
    std::ifstream ifs_phi(filename.str().c_str(), std::ifstream::in);
    if(ifs_phi.good()==false){
        cout << "ERROR: file "<< filename.str() << " not found!" <<endl;
        return EXIT_SUCCESS;
    }

    filename.str("");
    filename.clear();
    filename << res_folder << "/phi65resolution.txt";
    std::ifstream ifs_phi65(filename.str().c_str(), std::ifstream::in);
    if(ifs_phi65.good()==false){
        cout << "ERROR: file "<< filename.str() << " not found!" <<endl;
        return EXIT_SUCCESS;
    }

    filename.str("");
    filename.clear();
    filename << res_folder << "/thetaresolution.txt";
    std::ifstream ifs_theta(filename.str().c_str(), std::ifstream::in);
    if(ifs_theta.good()==false){
        cout << "ERROR: file "<< filename.str() << " not found!" <<endl;
        return EXIT_SUCCESS;
    }

    filename.str("");
    filename.clear();
    filename << res_folder << "/pTresolution.txt";
    std::ifstream ifs_pT(filename.str().c_str(), std::ifstream::in);
    if(ifs_pT.good()==false){
        cout << "ERROR: file "<< filename.str() << " not found!" <<endl;
        return EXIT_SUCCESS;
    }

    filename.str("");
    filename.clear();
    filename << res_folder << "/Charge.txt";
    std::ifstream ifs_charge(filename.str().c_str(), std::ifstream::in);
    if(ifs_charge.good()==false){
        cout << "ERROR: file "<< filename.str() << " not found!" <<endl;
        return EXIT_SUCCESS;
    }

    std::string key, header;
    double min, max;
    cout << "here" << endl;
    getline( ifs_z0, header);
    while(ifs_z0 >> key >> min >> max){
        z0Limits[key] = std::make_pair(min,max);
    }
    getline( ifs_phi, header);
    while(ifs_phi >> key >> min >> max){
        phiLimits[key] = std::make_pair(min,max);
    }
    getline( ifs_phi65, header);
    while(ifs_phi65 >> key >> min >> max){
        phi65Limits[key] = std::make_pair(min,max);
    }
    getline( ifs_pT, header);
    while(ifs_pT >> key >> min >> max){
        ptLimits[key] = std::make_pair(min,max);
    }
    getline(ifs_theta, header);
    while(ifs_theta >> key >> min >> max){
        thetaLimits[key] = std::make_pair(min,max);
    }
    getline(ifs_charge, header);
    while(ifs_charge >> key >> min){
        roadCharge[key] = min;
    }
    unsigned int TracksOutOfRange = 0;
    unsigned int nPatternsInAMsector = 0;
    unsigned int RecoTPinAMsector = 0;

    unsigned int nFiltPatternsInAMsector = 0;
    unsigned int FiltRecoTPinAMsector = 0;
    unsigned int MissedTrack = 0;
    unsigned int point = 0;
    // Events loop
    for (long long int ientry=0; ientry<n_events; ientry++){
        m_PATT->GetEntry(ientry);
        m_STUBS->GetEntry(ientry);
        m_MC->GetEntry(ientry);
        unsigned n_prepatterns = 0;

        
	    // --- Skip events with no roads:
        if ( n_patterns == 0 ) continue;

        // hPatterns->Fill(n_patterns); // Fill histogram of number of patterns per event before filtering
        int nFilteredPatterns=0;
        cout << "event "<< ientry << "/"<<n_events << endl;
        set<int> RecoTracks;
        set<int> RecoFiltTracks;
        set<int> RecoTracksInAMsector;
        set<int> RecoFiltTracksInAMsector;

        for (int iroad=0; iroad<n_patterns; ++iroad){
            bool patternInAMsector = false;
            if(patt_nmiss->at(iroad)!=0) continue;

		    hStubsPerPattern->Fill(patt_links->at(iroad).size()); // Fill No. Stubs per pattern histogram
            Hough* ht = new Hough();
            n_prepatterns++;
            // cout << "Pattern loop "<< endl;
            double xmin, xmax, ymin, ymax, charge;
            if(htType==0){
                xmin = ptLimits[patt_id->at(iroad)].first;
                xmax = ptLimits[patt_id->at(iroad)].second;
                
                if(settings->chosenRofPhi()==0){
                    ymin = phiLimits[patt_id->at(iroad)].first ;
                    ymax = phiLimits[patt_id->at(iroad)].second ;
                } else{
                    ymin = phi65Limits[patt_id->at(iroad)].first ;
                    ymax = phi65Limits[patt_id->at(iroad)].second ;
                }

                // cout << "limits assigned"<< endl;
                // std::map<std::string, std::pair<double, double> >::iterator it;
                // it = ptLimits.find(patt_id->at(iroad));
                if(xmin == 0 && xmax == 0 && ymin == 0 && ymax == 0){
                    xmin = -3;
                    xmax = +3;
                    ymin = 1.1;
                    ymax = 2.9;
                } else if(ymax - ymin < 0.001){
                    double ymean = (ymax+ymin)/2.;
                    ymin = ymean - .035;
                    ymax = ymean + .035;
                } else{
                    ymin = ymin - 0.005;
                    ymax = ymax + 0.005;
                }

                if(xmax - xmin < 0.05){
                    double xmean = (xmax + xmin)/2.;
                    xmin = xmean - 1.;
                    xmax = xmean + 1.;
                }
            } else{
                xmin = thetaLimits[patt_id->at(iroad)].first;
                xmax = thetaLimits[patt_id->at(iroad)].second;
                ymin = z0Limits[patt_id->at(iroad)].first;
                ymax = z0Limits[patt_id->at(iroad)].second;
                
                if(xmin == 0 && xmax == 0 && ymin == 0 && ymax == 0){
                    ymin = -15.;
                    ymax = 15.;
                    xmin = 1;
                    xmax = 2;
                } else if(ymax - ymin < 0.5){
                    double ymean = (ymax+ymin)/2.;
                    ymin = ymean - 5.;
                    ymax = ymean + 5.;
                } else{
                    ymin = ymin -5.;
                    ymax = ymax + 5;
                }
                if(xmax - xmin < 0.01){
                    double xmean = (xmax + xmin)/2.;
                    xmin = xmean - 0.03;
                    xmax = xmean + 0.03;
                } else{
                    xmin = xmin - 0.03;
                    xmax = xmax + 0.03;
                }
            }

            charge = roadCharge[patt_id->at(iroad)];
            
            ht->init(settings, xmin, xmax, ymin, ymax, charge);
            // cout << "HT initialised "<< endl;
            map< int, vector<Stub*> > StubsPerTP;
            map< int, vector<Stub*> > filtStubsPerTP;

            std::vector<unsigned int> vLayers(30,0);
            std::vector<bool> vPhiLayers(30,0); 

            for (unsigned int istub=0; istub<patt_links->at(iroad).size(); ++istub ){
                int stubRef = patt_links->at(iroad).at(istub);                
                Stub* stub = new Stub(stub_layer->at(stubRef), stub_phi->at(stubRef), stub_r->at(stubRef), stub_z->at(stubRef), stub_tp->at(stubRef));
                if(stub_tp->at(stubRef)>-1){
                    StubsPerTP[stub_tp->at(stubRef)].push_back(stub);
                    ht->store(stub);
                    vLayers[stub_layer->at(stubRef)-5]++;
                    double phi = stub_phi->at(stubRef);
                    double z = stub_z->at(stubRef);
                    double r = stub_r->at(stubRef);

                    if(insideEta(r,z) && insidePhi(phi, r)){
                        vPhiLayers[stub_layer->at(stubRef)-5] = true;
                    }
                }
            }

            unsigned int nPhiLayer = 0;
            for(bool activeLayer : vPhiLayers){
                if(activeLayer)
                    nPhiLayer++;
            }

            if(nPhiLayer > 4){
                nPatternsInAMsector++;
                patternInAMsector = true;
            }
            // cout <<"layers loop ended"<< endl;

            unsigned int comb = 1;
            unsigned int comb3 = 1;
            unsigned int sum = 0;
            for(unsigned i = 0; i < vLayers.size(); ++i) {
                if(vLayers[i]!=0){
                    comb = comb*vLayers[i];
                    if(i<3)
                        comb3 = comb3*vLayers[i];
                    if(i>2)
                        sum = sum + vLayers[i];
                }       
            }
            // cout <<"comb loop "<< endl;
            comb3 = comb3*sum;
            hSpecComb->Fill(comb3);
            hComb->Fill(comb);


            for(map<int, vector<Stub*> >::iterator it = StubsPerTP.begin(); it != StubsPerTP.end(); it++){
                const int maxLayerID(30);
                vector<bool> foundLayers(maxLayerID, false);
                for(Stub* stub: it->second){
                    foundLayers[stub->layerId()] = true;
                }
                // cout << "loop 0 "<< endl;
                unsigned int ncount = 0;
                for (const bool& found: foundLayers) {
                    if (found) ncount++;
                }
                // cout << "loop 1"<< endl;
                // cout << "ncount "<< ncount << " tp "<< it->first << endl;
                if(ncount > 4 && pu_pt->at(it->first) > 3.){
                    if(patternInAMsector) RecoTracksInAMsector.insert(it->first);
                    RecoTracks.insert(it->first);
                    // cout << "Track "<< it->first << endl;
                    if(ht->trackCandFound()){
                       RecoFiltTracks.insert(it->first);
                       if(patternInAMsector)RecoFiltTracksInAMsector.insert(it->first);
                    } else{
                        // cout <<"Stubs surviving HT "<<ht->maxNumStubs() <<  endl;
                        if(htType==0){
                            double phi65 = pu_phi->at(it->first) - 3.8112*(3.0E8/2.0E11)*65*pu_charge->at(it->first)/pu_pt->at(it->first);
                            if( double(pu_charge->at(it->first))/pu_pt->at(it->first) < -1/3. || double(pu_charge->at(it->first))/pu_pt->at(it->first) > 1/.3 || (settings->chosenRofPhi() == 0. && (pu_phi->at(it->first) < ht->yMin() || pu_phi->at(it->first) > ht->yMax()) ) ||  (settings->chosenRofPhi() > 0. && ( phi65 < ht->yMin() || phi65 > ht->yMax() ) ) )
                                TracksOutOfRange++;
                            else
                                MissedTrack++;
                            if(debug==2)
                                cout << "Track q/pt "<< double(pu_charge->at(it->first))/pu_pt->at(it->first) << " phi "<< pu_phi->at(it->first) << " array y: ["<< ht->yMin() << " , "<< ht->yMax()<<"]"<< " array x : ["<< ht->xMin() << ", "<< ht->xMax() << "]"<< " xmin "<<xmin << " xmax "<< xmax << " charge " << charge << " pattId "<<patt_id->at(iroad) << endl;
                        }

                        if(htType==1){
                            if(averageRes==0){
                                if( 1/tan(pu_theta->at(it->first)) < ht->xMin() || 1/tan(pu_theta->at(it->first)) > ht->xMax() || pu_vz->at(it->first) < ht->yMin() || pu_vz->at(it->first) > ht->yMax() )
                                TracksOutOfRange++;
                            else
                                MissedTrack++;
                        } else{
                            if( 1/tan(pu_theta->at(it->first)) < 1/tan(2.) || 1/tan(pu_theta->at(it->first)) > 1/tan(0.8) || pu_vz->at(it->first) < ht->yMin() || pu_vz->at(it->first) > ht->yMax() ){
                                TracksOutOfRange++;
                                cout << "Track cotan(theta) " <<1/tan(pu_theta->at(it->first))<< " array x: ["<<1/tan(2.) << " , "<< 1/tan(0.8)<< " z0: "<< pu_vz->at(it->first) << " y : ["<< ht->yMin() << " , " << ht->yMax() << endl;
                            }
                            else
                                MissedTrack++;
                        }

                            // cout << "Track q/pt "<< double(pu_charge->at(it->first))/pu_pt->at(it->first) << " phi "<< pu_phi->at(it->first) <<endl;
                        }

                        // else
                        //     cout << "Track z0 "<< pu_vz->at(it->first) << " cotan(theta) "<< 1/tan(pu_theta->at(it->first) ) <<endl;

                        // cout << "Hough range x: "<<ht->xMin()<< ", "<< ht->xMax() <<" y: "<<   ht->yMin()<< ", "<< ht->yMax() << endl;
                        // cout << "input xmin " << xmin <<" xmax "<< xmax << endl;
                        //             cout << "number of stubs in ht "<< ht->maxNumStubs() << endl;

                    }
                }
                // cout << "if 0 "<< endl;
   
            }

            // cout << "Reconstructed tracks loop"<< endl;

            unsigned int fcomb = 0;

            if(ht->trackCandFound()){
                nFilteredPatterns++;
                if(patternInAMsector) nFiltPatternsInAMsector++;
                hStubsPerFiltPattern->Fill(ht->maxNumStubs());
                std::vector<unsigned int> vFiltLayers(30,0);
                for(Stub* stub: ht->FilteredStubs()){
                    vFiltLayers[stub->layerId()-5]++;   
                }
                fcomb = 1;
                unsigned int fcomb3 = 1;

                unsigned int fsum = 0;
                for(unsigned i = 0; i < vFiltLayers.size(); ++i) {
                    if(vFiltLayers[i]!=0){
                        fcomb = fcomb*vFiltLayers[i];
                        if(i<3)
                            fcomb3 = fcomb3*vFiltLayers[i];
                        if(i>2)
                        fsum = fsum + vFiltLayers[i];
                    }       
                }

                fcomb3 = fcomb3*fsum;
                hFiltComb->Fill(fcomb);
                hFiltSpecComb->Fill(fcomb3);
            }
            gComb->SetPoint(point, comb, fcomb);
            // cout << "fcomb/comb "<< double(fcomb)/double(comb) <<" ("<<fcomb<<"/"<<comb<<")"<< endl;
            gCombRatio->SetPoint(point, comb, double(fcomb)/double(comb));
            ++point;
        }

        hRecoTracks->Fill(RecoTracks.size());
        hFiltRecoTracks->Fill(RecoFiltTracks.size());
        hFilteredPatterns->Fill(nFilteredPatterns);
        hPatterns->Fill(n_prepatterns);

        RecoTPinAMsector = RecoTPinAMsector + RecoTracksInAMsector.size();
        FiltRecoTPinAMsector = FiltRecoTPinAMsector + RecoFiltTracksInAMsector.size();

    }

    // Save the histograms
    TFile f(output_file.c_str(),"recreate");

    hPatterns->Write(); 
    hStubsPerPattern->Write();  
    hFilteredPatterns->Write(); 
    hStubsPerFiltPattern->Write();
    hRecoTracks->Write();
    hFiltRecoTracks->Write();
    hComb->Write();
    hFiltComb->Write();
    hSpecComb->Write();
    hFiltSpecComb->Write();

    gComb->GetYaxis()->SetTitle("n. filtered combinations");
    gComb->GetXaxis()->SetTitle("n. combinations");
    gComb->GetXaxis()->SetRangeUser(0,100);
    gComb->GetYaxis()->SetRangeUser(0,70);
    gComb->SetMarkerStyle(20);
    gComb->SetMarkerColor(kRed);
    gComb->Write();

    gCombRatio->GetYaxis()->SetTitle("n. filtered combinations/n.combinations");
    gCombRatio->GetXaxis()->SetTitle("n. combinations");
    gCombRatio->GetXaxis()->SetRangeUser(0,100);
    gCombRatio->GetYaxis()->SetRangeUser(0,1.1);
    gCombRatio->SetMarkerStyle(20);
    gCombRatio->SetMarkerColor(kRed);
    gCombRatio->Write();

    string name = "plots/npatterns.pdf";
    string title = "no. patterns";

    RatioPlot(hPatterns, hFilteredPatterns, title, name, n_events);

    cout << "========= SUMMARY =======" << endl;
    cout << "Average number of Patterns per event "<< hPatterns->GetMean() << endl;
    cout << "Average number of filtered Patterns per event "<< hFilteredPatterns->GetMean()<< endl;

    cout << "Average number of Stubs per Pattern "<< hStubsPerPattern->GetMean() << endl;
    cout << "Average number of Stubs per filtered Pattern "<< hStubsPerFiltPattern->GetMean()<< endl;

    cout << "Average number of Reconstructed Tracks per event "<< hRecoTracks->GetMean() << endl;
    cout << "Average number of Reconstructed Tracks per event after filtering "<< hFiltRecoTracks->GetMean() <<"( "<<hFiltRecoTracks->GetMean()*100./hRecoTracks->GetMean() <<"%)"<< endl;

    cout << "Average number of Combinations per Pattern "<< hComb->GetMean() << endl;
    cout << "Average number of Combinations per filtered Pattern "<< hFiltComb->GetMean()<< endl;
    cout << "Average number of Special Combinations per Pattern "<< hSpecComb->GetMean() << endl;
    cout << "Average number of Special Combinations per filtered Pattern "<< hFiltSpecComb->GetMean()<< endl;

    cout << "Tracks out of array ranges: "<<TracksOutOfRange<< ", missed tracks: "<< MissedTrack << endl;

    cout << "==== COMPARISON in TMT sector (4,26) =====" << endl;
    cout << "Number of Patterns "<< nPatternsInAMsector << endl;
    cout << "Number of Filtered Patterns "<< nFiltPatternsInAMsector << endl;
    cout << "Number of Reco TPs "<< RecoTPinAMsector << endl;
    cout << "Number of Filtered Reco TPs "<< FiltRecoTPinAMsector << endl;


    f.Close();
}