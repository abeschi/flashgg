/* To compile
g++ -Wall -o ZMuMuGammaPlotter `root-config --cflags --glibs` -L $ROOTSYS/lib -lFoam -lMinuit -lMathMore CMS_lumi.C tdrstyle.C  rochcor2016.cc RoccoR.cc ZMuMuGammaPlotter.cpp
*/


#ifndef CMS_LUMI_H
#include "CMS_lumi.h"
#endif

#ifndef CMS_STYLE
#include "tdrstyle.h"
#endif

#ifndef ElectroWeakAnalysis_RoccoR
#include "RoccoR.h"
#endif

#ifndef ROCHCOR2016_H
#include "rochcor2016.h"
#endif

#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TPaveStats.h"
#include "TLegend.h"
#include "TChain.h"
#include "TVirtualFitter.h"
#include "TLatex.h"
#include "TObjArray.h"
#include "TMinuit.h"
#include "TRandom.h"
#include "TLorentzVector.h"

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <ctime>
#include <map>
#include <algorithm>
#include <cmath>
#include <vector>
#include <cstdlib>


using namespace std;

float compute_DeltaPhi(float phi_1, float phi_2)
{
	float delta_phi = fabs(phi_1 - phi_2);
	if(delta_phi > 3.14159265359)
		delta_phi = 2*3.14159265359 - delta_phi;
	if(phi_1 - phi_2 >= 0)	
		return delta_phi;
	else
		return -delta_phi;
}



void HistoPlotter(TH1F* MC, TH1F* Data, TString name)
{
	TCanvas* c = new TCanvas();

	Data -> SetMarkerStyle(20);
	Data -> SetMarkerSize(1);
	MC -> SetFillStyle(1001);
	MC->SetFillColor(kYellow);

	MC -> Scale((Data->GetSumOfWeights()/MC->GetSumOfWeights()));

	double max = Data->GetBinContent(Data->GetMaximumBin());
	if(MC->GetBinContent(MC->GetMaximumBin()) > max )
		max = MC->GetBinContent(MC->GetMaximumBin());

	MC-> GetYaxis() -> SetRangeUser(0, 1.2*max);
	MC-> GetYaxis() -> SetTitleSize(0.05);
	MC-> GetYaxis() -> SetTitleFont(42);
	MC-> GetYaxis() -> SetLabelSize(0.045);
	MC-> GetYaxis() -> SetLabelFont(42);
	MC-> GetXaxis() -> SetLabelSize(0);

	TPad *pad1 = new TPad("pad1", "pad1", 0, 0.35, 1, 0.97);
	pad1->SetBottomMargin(0.035);
	pad1->Draw();
	pad1->cd();

	MC -> Draw("histo");
	Data -> Draw("EP SAME");


	TLegend* leg = new TLegend(0.65, 0.70, 0.8, 0.85);
	leg -> AddEntry(Data, "Data", "p");
	leg -> AddEntry(MC, "MC", "f");
	leg -> Draw("SAME");


	c->cd();
	TPad *pad2 = new TPad("pad2", "pad2", 0, 0.10, 1, 0.35);
	pad2->SetTopMargin(0);
	pad2->SetBottomMargin(0);
	pad2->Draw();
	pad2->cd();

	TH1F *h = (TH1F*)Data->Clone("h");
	h->SetLineColor(kBlack);
	h->SetMinimum(0.5);  // Define Y ..
	h->SetMaximum(1.5); // .. range
	h->Sumw2();
	h->SetStats(0);      // No statistics on lower plot
	h->Divide(MC);
	h->SetMarkerStyle(21);
	h -> SetTitle("");
	h-> GetYaxis() -> SetTitle("Data/MC");


	// Y axis ratio plot settings
	h->GetYaxis()->SetNdivisions(-6);
	h->GetYaxis()->SetTitleSize(0.13);
	h->GetYaxis()->SetTitleFont(42);
	h->GetYaxis()->SetTitleOffset(0.5);
	h->GetYaxis()->SetLabelFont(42); // Absolute font size in pixel (precision 3)
	h->GetYaxis()->SetLabelSize(0.12);

	// X axis ratio plot settings
	h->GetXaxis()->SetTitleSize(0.15);
	h->GetXaxis()->SetTitleFont(42);
	h->GetXaxis()->SetLabelFont(42); // Absolute font size in pixel (precision 3)
	h->GetXaxis()->SetLabelSize(0.12);

	h->Draw("EP");

	CMS_lumi(c, 0, 10);


	c -> SaveAs(name + ".png");
	c -> SaveAs(name + ".pdf");

	delete c;
	delete leg;

	return;
}







int main(int argc, char *argv[])
{
	writeExtraText = true;       // if extra text
	extraText  = "Preliminary";  // default extra text is "Preliminary"
	lumi_sqrtS = "36.73 fb^{-1} (13 TeV)";       // used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)

	setTDRStyle();
	gStyle -> SetOptFit(0);
	gStyle -> SetOptStat(0);


	//Final plots

	TH1F* MCMass = new TH1F("MCMass", "; M_{Z} (GeV); Counts", 40 , 70, 110);
	TH1F* MCMuMuMass = new TH1F("MCMuMuMass", "; M_{#mu#mu} (GeV); Counts", 45 , 0, 90);
	TH1F* MCMuMuPt = new TH1F("MCMuMuPt", "; P_{T #mu#mu} (GeV); Counts", 50 , 0, 100);
	TH1F* MCLeadingMuP = new TH1F("MCLeadingMuP", "; P (GeV); Counts", 70 , 10, 150);
	TH1F* MCLeadingMuPt = new TH1F("MCLeadingMuPt", "; P_{T} (GeV); Counts",50 , 0, 100);
	TH1F* MCLeadingMuEta = new TH1F("MCLeadingMuEta", "; #eta; Counts", 50 , -3, 3);
	TH1F* MCLeadingMuPhi = new TH1F("MCLeadingMuPhi", "; #varphi; Counts", 50 , -3.15, 3.15);
	TH1F* MCLeadingMuCharge = new TH1F("MCLeadingCharge", "; Charge; Counts", 10 , -2, 2);
	TH1F* MCLeadingMuNtrk = new TH1F("MCLeadingMuNtrk", "; N_{tkr}; Counts", 31 , -0.5, 30.5);        // For Rochester corrections
	TH1F* MCSubleadingMuP = new TH1F("MCSubleadingMuP", "; P (GeV); Counts", 50 , 0, 100);
	TH1F* MCSubleadingMuPt = new TH1F("MCSubleadingMuPt", "; P_{T} (GeV); Counts", 75 , 0, 150);
	TH1F* MCSubleadingMuEta = new TH1F("MCSubleadingMuEta", "; #eta; Counts", 50 , -3, 3);
	TH1F* MCSubleadingMuPhi = new TH1F("MCSubleadingMuPhi", "; #varphi; Counts", 50 , -3.15, 3.15);
	TH1F* MCSunleadingMuCharge = new TH1F("MCSubleadingCharge", "; Charge; Counts", 10 , -2, 2);
	TH1F* MCSubleadingMuNtrk = new TH1F("MCSubleadingMuNtrk", "; N_{tkr}; Counts",  31 , -0.5, 30.5);        // For Rochester corrections
	TH1F* MCPhotonE = new TH1F("MCPhotonE", "; Energy (GeV); Counts", 45 , 10, 100);
	TH1F* MCPhotonPt = new TH1F("MCPhotonPt", "; P_{T} (GeV); Counts", 45 , 10, 100);
	TH1F* MCPhotonEta = new TH1F("MCPhotonEta", "; #eta; Counts", 50 , -2.5, 2.5);
	TH1F* MCPhotonPhi = new TH1F("MCPhotonPhi", "; #varphi; Counts", 50 , -3.15, 3.15);
	TH1F* MCPhotonR9 = new TH1F("MCPhotonR9", "; R_{9}; Counts", 25 , 0.5, 1);
	TH1F* MCPhotonS4 = new TH1F("MCPhotonS4", "; S_{4}; Counts", 25 , 0, 1);
	TH1F* MCPhotonEtaWidth = new TH1F("MCPhotonEtaWidth", "; EtaWidth; Counts", 25 , 0, 0.05);

	TH1F* DataMass = new TH1F("DataMass", "; M_{Z} (GeV); Counts", 40 , 70, 110);
	TH1F* DataMuMuMass = new TH1F("DataMuMuMass", "; M_{#mu#mu} (GeV); Counts", 45 , 0, 90);
	TH1F* DataMuMuPt = new TH1F("DataMuMuPt", "; P_{T #mu#mu} (GeV); Counts", 50 , 0, 100);
	TH1F* DataLeadingMuP = new TH1F("DataLeadingMuP", "; P (GeV); Counts",140 , 10, 150);
	TH1F* DataLeadingMuPt = new TH1F("DataLeadingMuPt", "; P_{T} (GeV); Counts", 50 , 0, 100);
	TH1F* DataLeadingMuEta = new TH1F("DataLeadingMuEta", "; #eta; Counts", 50 , -3, 3);
	TH1F* DataLeadingMuPhi = new TH1F("DataLeadingMuPhi", "; #varphi; Counts", 50 , -3.15, 3.15);
	TH1F* DataLeadingMuCharge = new TH1F("DataLeadingCharge", "; Charge; Counts", 10 , -2, 2);
	TH1F* DataLeadingMuNtrk = new TH1F("DataLeadingMuNtrk", "; N_{tkr}; Counts", 31 , -0.5, 30.5);        // For Rochester corrections
	TH1F* DataSubleadingMuP = new TH1F("DataSubleadingMuP", "; P (GeV); Counts", 75 , 0, 150);
	TH1F* DataSubleadingMuPt = new TH1F("DataSubleadingMuPt", "; P_{T} (GeV); Counts", 50 , 0, 100);
	TH1F* DataSubleadingMuEta = new TH1F("DataSubleadingMuEta", "; #eta; Counts", 50 , -3, 3);
	TH1F* DataSubleadingMuPhi = new TH1F("DataSubleadingMuPhi", "; #varphi; Counts", 50 , -3.15, 3.15);
	TH1F* DataSunleadingMuCharge = new TH1F("DataSubleadingCharge", "; Charge; Counts", 10 , -2, 2);
	TH1F* DataSubleadingMuNtrk = new TH1F("DataSubleadingMuNtrk", "; N_{tkr}; Counts", 31 , -0.5, 30.5);        // For Rochester corrections
	TH1F* DataPhotonE = new TH1F("DataPhotonE", "; Energy (GeV); Counts", 45 , 10, 100);
	TH1F* DataPhotonPt = new TH1F("DataPhotonPt", "; P_{T} (GeV); Counts", 45 , 10, 100);
	TH1F* DataPhotonEta = new TH1F("DataPhotonEta", "; #eta; Counts", 50 , -2.5, 2.5);
	TH1F* DataPhotonPhi = new TH1F("DataPhotonPhi", "; #varphi; Counts", 50 , -3.15, 3.15);
	TH1F* DataPhotonR9 = new TH1F("DataPhotonR9", "; R_{9}; Counts", 25 , 0.5, 1);
	TH1F* DataPhotonS4 = new TH1F("DataPhotonS4", "; S_{4}; Counts", 25 , 0, 1);
	TH1F* DataPhotonEtaWidth = new TH1F("DataPhotonEtaWidth", "; EtaWidth; Counts", 25 , 0, 0.05);


	TString MCFolder = "/afs/cern.ch/work/a/abeschi/MuMuGammaMC/";
	TString DataFolder = "/afs/cern.ch/work/a/abeschi/MuMuGamma/";

	int nentries;

	float weight = 0;
	float mass_mmg = 0;
	float mass_mumu = 0;
	float pt_mumu = 0;
	float leadMuonP = 0;
	float leadMuonPt = 0;
	float leadMuonEta = 0;
	float leadMuonPhi = 0;
	float leadMuonCharge = 0;
	float leadMuonNtrk = 0;
	float subleadMuonP = 0;
	float subleadMuonPt = 0;
	float subleadMuonEta = 0;
	float subleadMuonPhi = 0;
	float subleadMuonCharge = 0;
	float subleadMuonNtrk = 0;
	float photonE = 0;
	float photonPt = 0;
	float photonEta = 0;
	float photonPhi = 0;
	float photonR9 = 0;
	float photonS4 = 0;
	float photonEtaWidth = 0;

	float muMass = 0.1057;
	int evData = 0;
	int evMC = 0;
	float sumW = 0;

	TChain* DataTree;
	TChain* MCTree;

	float lumiFactor = 36.46;

	DataTree = new TChain("mumugammaDumper/trees/tree");
	DataTree -> Add(DataFolder + "MuMuGammaData_*.root");		
	MCTree = new TChain("mumugammaDumper/trees/tree");
	MCTree -> Add(MCFolder + "MuMuGammaMC_*.root");

	MCTree -> SetBranchAddress("weight", &weight);
       	MCTree -> SetBranchAddress("mass_mmg", &mass_mmg);
	MCTree -> SetBranchAddress("mass_mumu", &mass_mumu);
	MCTree -> SetBranchAddress("pt_mumu", &pt_mumu);
	MCTree -> SetBranchAddress("leadMuonP", &leadMuonP);
	MCTree -> SetBranchAddress("leadMuonPt", &leadMuonPt);
	MCTree -> SetBranchAddress("leadMuonEta", &leadMuonEta);
	MCTree -> SetBranchAddress("leadMuonPhi", &leadMuonPhi);
	MCTree -> SetBranchAddress("leadMuonCharge", &leadMuonCharge);
	MCTree -> SetBranchAddress("leadMuonNtrk", &leadMuonNtrk);
	MCTree -> SetBranchAddress("subleadMuonP", &subleadMuonP);
	MCTree -> SetBranchAddress("subleadMuonPt", &subleadMuonPt);
	MCTree -> SetBranchAddress("subleadMuonEta", &subleadMuonEta);
	MCTree -> SetBranchAddress("subleadMuonPhi", &subleadMuonPhi);
	MCTree -> SetBranchAddress("subleadMuonCharge", &subleadMuonCharge);
	MCTree -> SetBranchAddress("subleadMuonNtrk", &subleadMuonNtrk);
	MCTree -> SetBranchAddress("photonE", &photonE);
	MCTree -> SetBranchAddress("photonPt", &photonPt);
	MCTree -> SetBranchAddress("photonEta", &photonEta);
	MCTree -> SetBranchAddress("photonPhi", &photonPhi);
	MCTree -> SetBranchAddress("photonR9", &photonR9);
	MCTree -> SetBranchAddress("photonS4", &photonS4);
	MCTree -> SetBranchAddress("photonEtaWidth", &photonEtaWidth);

	DataTree -> SetBranchAddress("mass_mmg", &mass_mmg);
	DataTree -> SetBranchAddress("mass_mumu", &mass_mumu);
	DataTree -> SetBranchAddress("pt_mumu", &pt_mumu);
	DataTree -> SetBranchAddress("leadMuonP", &leadMuonP);
	DataTree -> SetBranchAddress("leadMuonPt", &leadMuonPt);
	DataTree -> SetBranchAddress("leadMuonEta", &leadMuonEta);
	DataTree -> SetBranchAddress("leadMuonPhi", &leadMuonPhi);
	DataTree -> SetBranchAddress("leadMuonCharge", &leadMuonCharge);
	DataTree -> SetBranchAddress("leadMuonNtrk", &leadMuonNtrk);
	DataTree -> SetBranchAddress("subleadMuonP", &subleadMuonP);
	DataTree -> SetBranchAddress("subleadMuonPt", &subleadMuonPt);
	DataTree -> SetBranchAddress("subleadMuonEta", &subleadMuonEta);
	DataTree -> SetBranchAddress("subleadMuonPhi", &subleadMuonPhi);
	DataTree -> SetBranchAddress("subleadMuonCharge", &subleadMuonCharge);
	DataTree -> SetBranchAddress("subleadMuonNtrk", &subleadMuonNtrk);
	DataTree -> SetBranchAddress("photonE", &photonE);
	DataTree -> SetBranchAddress("photonPt", &photonPt);
	DataTree -> SetBranchAddress("photonEta", &photonEta);
	DataTree -> SetBranchAddress("photonPhi", &photonPhi);
	DataTree -> SetBranchAddress("photonR9", &photonR9);
	DataTree -> SetBranchAddress("photonS4", &photonS4);
	DataTree -> SetBranchAddress("photonEtaWidth", &photonEtaWidth);

	nentries = MCTree -> GetEntries();
	rochcor2016 *rmcor = new rochcor2016();

	for(int i=0; i<nentries; i++)
	{
		MCTree -> GetEntry(i);
		weight = weight*lumiFactor;

		sumW += weight;
		evMC++;

		if(i%1000==0)
			cout << "Processed " << i << " events out of " << nentries << endl; 

		TLorentzVector mu1, mu2;
		mu1.SetPtEtaPhiM(leadMuonPt, leadMuonEta, leadMuonPhi, muMass);
    		mu2.SetPtEtaPhiM(subleadMuonPt, subleadMuonEta, subleadMuonPhi, muMass);
		float qter = 1.0;

		rmcor->momcor_mc(mu1, leadMuonCharge, leadMuonNtrk, qter); 
		rmcor->momcor_mc(mu2, subleadMuonCharge, subleadMuonNtrk, qter); 

		float mass_mumu_corr = mass_mumu * sqrt ( mu1.P()*mu2.P() / leadMuonP/subleadMuonP);

		MCMass -> Fill(mass_mmg, weight);
		MCMuMuMass -> Fill(mass_mumu, weight);
		MCMuMuPt -> Fill(pt_mumu, weight);
		MCLeadingMuP -> Fill(leadMuonP, weight);
		MCLeadingMuPt -> Fill(leadMuonPt, weight);
		MCLeadingMuEta -> Fill(leadMuonEta, weight);
		MCLeadingMuPhi -> Fill(leadMuonPhi, weight);
		MCLeadingMuCharge -> Fill(leadMuonCharge, weight);
		MCLeadingMuNtrk -> Fill(leadMuonNtrk, weight);
		MCSubleadingMuP -> Fill(subleadMuonP, weight);
		MCSubleadingMuPt -> Fill(subleadMuonPt, weight);
		MCSubleadingMuEta -> Fill(subleadMuonEta, weight);
		MCSubleadingMuPhi -> Fill(subleadMuonPhi, weight);
		MCSunleadingMuCharge -> Fill(subleadMuonCharge, weight);
		MCSubleadingMuNtrk -> Fill(subleadMuonNtrk, weight);
		MCPhotonE -> Fill(photonE, weight);
		MCPhotonPt -> Fill(photonPt, weight);
		MCPhotonEta -> Fill(photonEta, weight);
		MCPhotonPhi -> Fill(photonPhi, weight);
		MCPhotonR9 -> Fill(photonR9, weight);
		MCPhotonS4 -> Fill(photonS4, weight);
		MCPhotonEtaWidth -> Fill(photonEtaWidth, weight);
	}	


	nentries = DataTree -> GetEntries();

	for(int i=0; i<nentries; i++)
	{
		DataTree -> GetEntry(i);
		weight = 1;

		evData++;
		if(i%1000==0)
			cout << "Processed " << i << " events out of " << nentries << endl;

		TLorentzVector mu1, mu2;
		mu1.SetPtEtaPhiM(leadMuonPt, leadMuonEta, leadMuonPhi, muMass);
    		mu2.SetPtEtaPhiM(subleadMuonPt, subleadMuonEta, subleadMuonPhi, muMass);
		float qter = 1.0;

		rmcor->momcor_data(mu1, leadMuonCharge, 0, qter); 
		rmcor->momcor_data(mu2, subleadMuonCharge, 0, qter); 

		float mass_mumu_corr = mass_mumu * sqrt ( mu1.P()*mu2.P() / leadMuonP/subleadMuonP);

		DataMass -> Fill(mass_mmg, weight);
		DataMuMuMass -> Fill(mass_mumu, weight);
		DataMuMuPt -> Fill(pt_mumu, weight);
		DataLeadingMuP -> Fill(leadMuonP, weight);
		DataLeadingMuPt -> Fill(leadMuonPt, weight);
		DataLeadingMuEta -> Fill(leadMuonEta, weight);
		DataLeadingMuPhi -> Fill(leadMuonPhi, weight);
		DataLeadingMuCharge -> Fill(leadMuonCharge, weight);
		DataLeadingMuNtrk -> Fill(leadMuonNtrk, weight);
		DataSubleadingMuP -> Fill(subleadMuonP, weight);
		DataSubleadingMuPt -> Fill(subleadMuonPt, weight);
		DataSubleadingMuEta -> Fill(subleadMuonEta, weight);
		DataSubleadingMuPhi -> Fill(subleadMuonPhi, weight);
		DataSunleadingMuCharge -> Fill(subleadMuonCharge, weight);
		DataSubleadingMuNtrk -> Fill(subleadMuonNtrk, weight);
		DataPhotonE -> Fill(photonE, weight);
		DataPhotonPt -> Fill(photonPt, weight);
		DataPhotonEta -> Fill(photonEta, weight);
		DataPhotonPhi -> Fill(photonPhi, weight);
		DataPhotonR9 -> Fill(photonR9, weight);
		DataPhotonS4 -> Fill(photonS4, weight);
		DataPhotonEtaWidth -> Fill(photonEtaWidth, weight);
	}	

	HistoPlotter(MCMass, DataMass, "Mass");
	HistoPlotter(MCMuMuMass, DataMuMuMass, "MuMuMass");
	HistoPlotter(MCMuMuPt, DataMuMuPt, "MuMuPt");
	HistoPlotter(MCLeadingMuP, DataLeadingMuP, "LeadingMuP");
	HistoPlotter(MCLeadingMuPt, DataLeadingMuPt, "LeadingMuPt");
	HistoPlotter(MCLeadingMuEta, DataLeadingMuEta, "LeadingMuEta");
	HistoPlotter(MCLeadingMuPhi, DataLeadingMuPhi, "LeadingMuPhi");
	HistoPlotter(MCLeadingMuCharge, DataLeadingMuCharge, "LeadingMuCharge");
	HistoPlotter(MCLeadingMuNtrk, DataLeadingMuNtrk, "LeadingMuNtrk");
	HistoPlotter(MCSubleadingMuP, DataSubleadingMuP, "SubleadingMuP");
	HistoPlotter(MCSubleadingMuPt, DataSubleadingMuPt, "SubleadingMuPt");
	HistoPlotter(MCSubleadingMuEta, DataSubleadingMuEta, "SubleadingMuEta");
	HistoPlotter(MCSubleadingMuPhi, DataSubleadingMuPhi, "SubleadingMuPhi");
	HistoPlotter(MCSunleadingMuCharge, DataSunleadingMuCharge, "SunleadingMuCharge");
	HistoPlotter(MCSubleadingMuNtrk, DataSubleadingMuNtrk, "SubleadingMuNtrk");
	HistoPlotter(MCPhotonE, DataPhotonE, "PhotonE");
	HistoPlotter(MCPhotonPt, DataPhotonPt, "PhotonPt");
	HistoPlotter(MCPhotonEta, DataPhotonEta, "PhotonEta");
	HistoPlotter(MCPhotonPhi, DataPhotonPhi, "PhotonPhi");
	HistoPlotter(MCPhotonR9, DataPhotonR9, "PhotonR9");
	HistoPlotter(MCPhotonS4, DataPhotonS4, "PhotonS4");
	HistoPlotter(MCPhotonEtaWidth, DataPhotonEtaWidth, "PhotonEtaWidth");	


	system("mv *.pdf ~/www/MuMuGamma");
	system("mv *.png ~/www/MuMuGamma");



	cout << "Data events: " << evData << endl;
	cout << "MC events: " << sumW << " (" << evMC << " unweighted)" << endl; 

}





















