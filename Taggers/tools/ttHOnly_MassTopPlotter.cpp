/*
g++ -Wall -o ttHOnly_MassTopPlotter `root-config --cflags --glibs` -L $ROOTSYS/lib -lRooFitCore -lFoam -lMinuit -lMathMore CMS_lumi.C tdrstyle.C ttHOnly_MassTopPlotter.cpp
*/

#ifndef CMS_LUMI_H
#include "CMS_lumi.h"
#endif

#ifndef CMS_STYLE
#include "tdrstyle.h"
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
#include "TAxis.h"
#include "TMath.h"
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


using namespace std;


float compute_R(float eta_1, float eta_2, float phi_1, float phi_2)
{
	float delta_eta_squared = (eta_1 - eta_2)*(eta_1 - eta_2);
	float delta_phi = fabs(phi_1 - phi_2);
		if(delta_phi > 3.14159265359)
			delta_phi =2*3.14159265359 - delta_phi;

	return sqrt(delta_eta_squared + delta_phi*delta_phi);
}


std::vector<int> genMatching(float* top1_eta, float* top1_phi, float* top1_pdgId, float* top2_eta, float* top2_phi, float* top2_pdgId, float* jet_eta, float* jet_phi, float* jet_pt, TH1F* h)
{
	float tops_eta[6];
	float tops_phi[6];
	float tops_pdgId[6];

	for(int i=0; i<3; i++)
	{
		tops_eta[i] = top1_eta[i];
		tops_phi[i] = top1_phi[i];
		tops_pdgId[i] = top1_pdgId[i];
	}

	for(int i=0; i<3; i++)
	{
		tops_eta[i+3] = top2_eta[i];
		tops_phi[i+3] = top2_phi[i];
		tops_pdgId[i+3] = top2_pdgId[i];
	}

	std::vector<int> Indexes;
	std::vector<int> MatchedIndexes;
	for(int i=0; i<9; i++)
	{
		if(jet_pt[i]>0.)
			Indexes.push_back(i);
	}
	
	for(unsigned int i=0; i<6 ; i++)
	{
		float minDr = 100.;
		int minDrIdx = -1;
		float maxDrThreshold = 0.15;

		//if(Indexes.size()<6) break;
		for(unsigned int j=0; j<Indexes.size(); j++)
		{
			float deltaR = compute_R(tops_eta[i], jet_eta[Indexes[j]], tops_phi[i], jet_phi[Indexes[j]]);

			if(deltaR<minDr)
			{
				minDr = deltaR;
				minDrIdx = Indexes[j];
				h -> Fill(minDr);
			}
			
		}

		if(minDr>maxDrThreshold) break;
		MatchedIndexes.push_back(minDrIdx);

		std::vector<int> tmp;
		for(unsigned int idx=0; idx<Indexes.size(); idx++)
		{
			if(Indexes[idx]!=minDrIdx)
				tmp.push_back(Indexes[idx]);
		}

		Indexes.clear();
		Indexes = tmp;
		tmp.clear();
	}



/*	if(MatchedIndexes.size() == 6)
	{
		for(unsigned int i=0; i<MatchedIndexes.size(); i++)
		{
			cout << i << " " << MatchedIndexes[i] << endl;
			cout << tops_eta[i] << " " << tops_phi[i] << " " << tops_pdgId[i] << endl; 
			cout << jet_eta[MatchedIndexes[i]] << " " << jet_phi[MatchedIndexes[i]] << " " << jet_pt[MatchedIndexes[i]] << endl; 
		}
	}
*/
	return MatchedIndexes;
}




float NormalizationFactor(TH1F* h)
{

	float n = 1;

	TF1* f = new TF1("f", "expo", 100, 180);

	int i=0;
	int counter = 0;

	while(counter<10 && i!=0)
	{
		TFitResultPtr p = h -> Fit("f", "RQN");
		i = p;
		counter++;
	}

	
	if(i==0)
		n = 1. - f->Integral(115, 135)/f->Integral(100, 180);

	return n;
}




void MakePlot(TH1F* histo, TString title)
{
	TCanvas* c = new TCanvas;
	c -> cd();

	histo -> SetLineWidth(2);
	histo -> SetLineColor(kBlue);
	histo -> GetYaxis() -> SetRangeUser(0, 1.2*histo -> GetMaximum());
	histo -> Draw("histo");

	CMS_lumi(c, 0, 0);

	c -> SaveAs(title + ".pdf");
	c -> SaveAs(title + ".png");

	return;
}

void MakePlot2D(TH2F* histo, TString title)
{
	TCanvas* c = new TCanvas;
	c -> cd();

	histo -> Draw("COLZ");

	CMS_lumi(c, 0, 0);

	c -> SaveAs(title + ".pdf");
	c -> SaveAs(title + ".png");

	return;
}

void Make2Plot(TH1F* histo1, TH1F* histo2, TString title, TString l1, TString l2, bool Normalize=0)
{
	TCanvas* c = new TCanvas;
	c -> cd();

	histo1 -> SetLineWidth(3);
	histo1 -> SetLineColor(kRed+1);

	histo2 -> SetLineWidth(3);
	histo2 -> SetLineColor(kGreen+2);

	float Max = 0;
	if(!Normalize)
		Max = max(histo1->GetMaximum(), histo2->GetMaximum());
	else
		Max = max(histo1->GetMaximum()/histo1->Integral(), histo2->GetMaximum()/histo2->Integral());

	histo1 -> GetYaxis() -> SetRangeUser(0, 1.2*Max);
	if(Normalize)
	{
		histo1 -> DrawNormalized("histo");
		histo2 -> DrawNormalized("histo SAME");	
	}	

	else
	{
		histo1 -> Draw("histo");
		histo2 -> Draw("histo SAME");	
	}	

	TLegend* leg = new TLegend(0.7, 0.7, 0.9, 0.9);
	leg -> AddEntry(histo1, l1, "l");
	leg -> AddEntry(histo2, l2, "l");
	leg -> Draw("SAME");

	CMS_lumi(c, 0, 0);

	c -> SaveAs(title + ".pdf");
	c -> SaveAs(title + ".png");

	return;
}


int main(int argc, char *argv[])
{
	writeExtraText = true;       // if extra text
	extraText  = "Preliminary Simulation";  // default extra text is "Preliminary"
	lumi_sqrtS = "(13 TeV)";       // used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)

	setTDRStyle();
	gStyle -> SetOptFit(0);
	gStyle -> SetOptStat(0);

	TString SignalFold = "/afs/cern.ch/work/a/abeschi/ttH_MassTop_v3/";

	TChain* ttH;
	
	float lumiFactor = 35.9;
	float bDiscriminantThreshold = 0.6;
	float topPtThreshold = 100;
	float topMassLowThreshold = 140;
	float topMassHighThreshold = 200;
	float topDeltaEtaThreshold = 2.;
	float topDeltaPhiThreshold = 2.5;
	float WDeltaEtaThreshold = 1.5;
	float WDeltaPhiThreshold = 2.;
	float WMassLowThreshold = 60;
	float WMassHighThreshold = 110;
	float WPtThreshold = 80;

	bool isLept = (argc!=1 ? 0 : 1);

	if(isLept)
	{	cout << "Processing leptonic tag" << endl;

		ttH = new TChain("TTHLeptonicDumper/trees/tth_13TeV_all");
		ttH -> Add(SignalFold + "ttHJetToGG_M125_13TeV_*.root");		
	}

	else
	{	cout << "Processing hadronic tag" << endl;

		ttH = new TChain("TTHHadronicDumper/trees/tth_13TeV_all");
		ttH -> Add(SignalFold + "ttHJetToGG_M125_13TeV_*.root");		
	}

	TH1F* Dr = new TH1F("Dr", "; #DeltaR^{Min}; Counts", 50, 0, 1);
	TH1F* WMass1 = new TH1F("WMass1", "; mass^{W} (GeV); Counts", 160, 40, 200);
	TH1F* WMass2 = new TH1F("WMass2", "; mass^{W} (GeV); Counts", 160, 40, 200);
	TH1F* TopMass1 = new TH1F("TopMass1", "; mass^{top} (GeV); Counts", 120, 110, 230);
	TH1F* TopMass2 = new TH1F("TopMass2", "; mass^{top} (GeV); Counts", 120, 110, 230);
	TH1F* bDiscriminant = new TH1F("bDiscriminant", ";b-discriminant; Counts", 50, 0, 1);
	TH1F* WbDiscriminant = new TH1F("WbDiscriminant", "; b-discriminant; Counts", 50, 0, 1);
	TH1F* WDeltaEta = new TH1F("WDeltaEta", "; #Delta#eta; Counts", 50, -3, 3);
	TH1F* RandomDeltaEta = new TH1F("RandomDeltaEta", "; #Delta#eta; Counts", 50, -3, 3);
	TH1F* WDeltaPhi = new TH1F("WDeltaPhi", "; #Delta#varphi; Counts", 50, 0, 3.15);
	TH1F* RandomDeltaPhi = new TH1F("RandomDeltaPhi", "; #Delta#varphi; Counts", 50, 0, 3.15);
	TH1F* WDeltaR= new TH1F("WDeltaR", "; #DeltaR; Counts", 50, 0, 5);
	TH1F* RandomDeltaR = new TH1F("RandomDeltaR", "; #DeltaR; Counts", 50, 0, 5);
	TH1F* WPt= new TH1F("WPt", "; p_{T}; Counts", 100, 0, 500);
	TH1F* RandomPt = new TH1F("RandomPt", "; p_{T}; Counts", 100, 0, 500);
	TH2F* WMassDeltaR = new TH2F("WMassDeltaR", "; mass^{W} (GeV); #DeltaR; Counts", 80, 40, 120, 50, 0, 5);
	TH1F* TopDeltaEta = new TH1F("TopDeltaEta", "; #Delta#eta; Counts", 50, -3, 3);
	TH1F* TopDeltaPhi = new TH1F("TopDeltaPhi", "; #Delta#varphi; Counts", 50, 0, 3.15);
//	TH1F* RandomTopDeltaEta = new TH1F("RandomTopDeltaEta", "; #Delta#eta; Counts", 50, -3, 3);
//	TH1F* RandomTopDeltaPhi = new TH1F("RandomTopDeltaPhi", "; #Delta#varphi; Counts", 50, 0, 3.15);
	TH1F* NTopCandidate = new TH1F("NTopCandidate", "; N^{top}; Counts", 3, -0.5, 2.5);
	TH1F* NTopCandidateCorrect = new TH1F("NTopCandidateCorrect", "; N^{top}; Counts", 3, -0.5, 2.5);
	TH1F* NWCandidate = new TH1F("NWCandidate", "; N^{W}; Counts", 11, -0.5, 10.5);
	TH1F* CorrectWPosition = new TH1F("CorrectWPosition", "; N^{b-jets}; Counts", 4, -0.5, 3.5);
	TH1F* CorrectbCandidate = new TH1F("CorrectbCandidate", "; N^{b}; Counts", 3, -0.5, 2.5);
	TH1F* RecoTopMass = new TH1F("RecoTopMass", ";mass^{top} (GeV); Counts", 120, 110, 230);
	TH1F* RecoTopMassRandom = new TH1F("RecoTopMassRandom", ";mass^{top} (GeV); Counts", 120, 110, 230);
	TH1F* RecoTopPt = new TH1F("RecoTopPt", ";p_{T}^{top} (GeV); Counts", 100, 0, 500);
	TH1F* RecoTopPtRandom = new TH1F("RecoTopPtRandom", ";p_{T}^{top} (GeV); Counts", 100, 0, 500);
	TH1F* RecoTopEta = new TH1F("RecoTopEta", ";#eta^{top}; Counts", 50, -3, 3);
	TH1F* RecoTopEtaRandom = new TH1F("RecoTopEtaRandom", ";#eta^{top}; Counts", 50, -3, 3);
	TH1F* RecoTopDeltaEta = new TH1F("RecoTopDeltaEta", ";#Delta#eta; Counts", 50, -3, 3);
	TH1F* RecoTopDeltaEtaRandom = new TH1F("RecoTopDeltaEtaRandom", ";#Delta#eta; Counts", 50, -3, 3);
	TH1F* RecoTopDeltaPhi = new TH1F("RecoTopDeltaPhi", ";#Delta#varphi; Counts", 50, 0, 3.15);
	TH1F* RecoTopDeltaPhiRandom = new TH1F("RecoTopDeltaPhiRandom", ";#Delta#varphi; Counts", 50, 0, 3.15);

	int nvtx = 0;
	float weight = 0;
	
	float dipho_sumpt = 0;
	float dipho_mass = 0;
	float dipho_leadPt = 0;
	float dipho_leadEta = 0;
	float dipho_leadPhi = 0;
	float dipho_leadR9 = 0;
	float dipho_subleadPt = 0;
	float dipho_subleadEta = 0;
	float dipho_subleadPhi = 0;
	float dipho_subleadR9 = 0;
	float dipho_leadIDMVA = 0;
	float dipho_subleadIDMVA = 0;
	float dipho_mva = 0;

	float jet_pt[9] = {0., 0., 0., 0., 0., 0., 0., 0., 0.};
	float jet_eta[9] = {0., 0., 0., 0., 0., 0., 0., 0., 0.};
	float jet_phi[9] = {0., 0., 0., 0., 0., 0., 0., 0., 0.};
	float jet_bdiscriminant[9] = {0., 0., 0., 0., 0., 0., 0., 0., 0.};
	float jet_energy[9] = {0., 0., 0., 0., 0., 0., 0., 0., 0.};

	float MetPt = 0.;
	float MetPhi = 0.;
	float ttHMVA = 0.;

	float GenLeadPhton_pt = 0.;
	float GenLeadPhton_eta = 0.;
	float GenLeadPhton_phi = 0.;
	float GenLeadPhton_energy = 0.;
	float GenSubleadPhton_pt = 0.;
	float GenSubleadPhton_eta = 0.;
	float GenSubleadPhton_phi = 0.;
	float GenSubleadPhton_energy = 0.;

	float top1_pt[3] = {0., 0., 0.};
	float top1_eta[3] = {0., 0., 0.};
	float top1_phi[3] = {0., 0., 0.};
	float top1_energy[3] = {0., 0., 0.};
	float top1_pdgId[3] = {0., 0., 0.};
	float top2_pt[3] = {0., 0., 0.};
	float top2_eta[3] = {0., 0., 0.};
	float top2_phi[3] = {0., 0., 0.};
	float top2_energy[3] = {0., 0., 0.};
	float top2_pdgId[3] = {0., 0., 0.};


	ttH -> SetBranchAddress("nvtx", &nvtx);
	ttH -> SetBranchAddress("weight", &weight);
	ttH -> SetBranchAddress("dipho_sumpt", &dipho_sumpt);
	ttH -> SetBranchAddress("dipho_mass", &dipho_mass);
	ttH -> SetBranchAddress("dipho_leadPt", &dipho_leadPt);
	ttH -> SetBranchAddress("dipho_leadEta", &dipho_leadEta);
	ttH -> SetBranchAddress("dipho_leadPhi", &dipho_leadPhi);
	ttH -> SetBranchAddress("dipho_leadR9", &dipho_leadR9);
	ttH -> SetBranchAddress("dipho_subleadPt", &dipho_subleadPt);
	ttH -> SetBranchAddress("dipho_subleadEta", &dipho_subleadEta);
	ttH -> SetBranchAddress("dipho_subleadPhi", &dipho_subleadPhi);
	ttH -> SetBranchAddress("dipho_subleadR9", &dipho_subleadR9);
	ttH -> SetBranchAddress("dipho_leadIDMVA", &dipho_leadIDMVA);
	ttH -> SetBranchAddress("dipho_subleadIDMVA", &dipho_subleadIDMVA);
	ttH -> SetBranchAddress("dipho_mva", &dipho_mva);


	ttH -> SetBranchAddress("GenLeadPhton_pt", &GenLeadPhton_pt);
	ttH -> SetBranchAddress("GenLeadPhton_eta", &GenLeadPhton_eta);
	ttH -> SetBranchAddress("GenLeadPhton_phi", &GenLeadPhton_phi);
	ttH -> SetBranchAddress("GenLeadPhton_energy", &GenLeadPhton_energy);
	ttH -> SetBranchAddress("GenSubleadPhton_pt", &GenSubleadPhton_pt);
	ttH -> SetBranchAddress("GenSubleadPhton_eta", &GenSubleadPhton_eta);
	ttH -> SetBranchAddress("GenSubleadPhton_phi", &GenSubleadPhton_phi);
	ttH -> SetBranchAddress("GenSubleadPhton_energy", &GenSubleadPhton_energy);

	for(int i=1; i<10; i++)
	{	ttH -> SetBranchAddress(("jet_pt"+ std::to_string(i)).c_str(), &jet_pt[i-1]);
		ttH -> SetBranchAddress(("jet_eta"+ std::to_string(i)).c_str(), &jet_eta[i-1]);
		ttH -> SetBranchAddress(("jet_phi"+ std::to_string(i)).c_str(), &jet_phi[i-1]);
		ttH -> SetBranchAddress(("jet_bdiscriminant"+ std::to_string(i)).c_str(), &jet_bdiscriminant[i-1]);
		ttH -> SetBranchAddress(("jet_energy"+ std::to_string(i)).c_str(), &jet_energy[i-1]);
	}

	ttH -> SetBranchAddress("ttHMVA", &ttHMVA);
	ttH -> SetBranchAddress("MetPt", &MetPt);
	ttH -> SetBranchAddress("MetPhi", &MetPhi);

	for(int i=1; i<4; i++)
	{
		ttH -> SetBranchAddress(("Top1_pt"+ std::to_string(i)).c_str(), &top1_pt[i-1]);
		ttH -> SetBranchAddress(("Top1_eta"+ std::to_string(i)).c_str(), &top1_eta[i-1]);
		ttH -> SetBranchAddress(("Top1_phi"+ std::to_string(i)).c_str(), &top1_phi[i-1]);
		ttH -> SetBranchAddress(("Top1_energy"+ std::to_string(i)).c_str(), &top1_energy[i-1]);
		ttH -> SetBranchAddress(("Top1_pdgId"+ std::to_string(i)).c_str(), &top1_pdgId[i-1]);

		ttH -> SetBranchAddress(("Top2_pt"+ std::to_string(i)).c_str(), &top2_pt[i-1]);
		ttH -> SetBranchAddress(("Top2_eta"+ std::to_string(i)).c_str(), &top2_eta[i-1]);
		ttH -> SetBranchAddress(("Top2_phi"+ std::to_string(i)).c_str(), &top2_phi[i-1]);
		ttH -> SetBranchAddress(("Top2_energy"+ std::to_string(i)).c_str(), &top2_energy[i-1]);
		ttH -> SetBranchAddress(("Top2_pdgId"+ std::to_string(i)).c_str(), &top2_pdgId[i-1]);
	}


	int nentries = ttH -> GetEntries();
	int ProbeW = 0;
	int OneTop = 0;
	int TwoTops = 0;
//	int nentries = 100;

	for(int i=0; i<nentries; i++)
	{	
		ttH -> GetEntry(i);
		if(i%10000==0) cout << "Processing event " << i << " out of " << nentries << "\r" << flush;
//		cout << "Processing event " << i << " out of " << nentries << endl;

		std::vector<int> MatchedIndex = genMatching(top1_eta, top1_phi, top1_pdgId, top2_eta, top2_phi, top2_pdgId, jet_eta, jet_phi, jet_pt, Dr);

		if(MatchedIndex.size()!= 6) continue;

		ProbeW++;

		TLorentzVector* b1 = new TLorentzVector();
		b1 -> SetPtEtaPhiE(jet_pt[MatchedIndex[0]], jet_eta[MatchedIndex[0]], jet_phi[MatchedIndex[0]], jet_energy[MatchedIndex[0]]);
		TLorentzVector* q1 = new TLorentzVector();
		q1 -> SetPtEtaPhiE(jet_pt[MatchedIndex[1]], jet_eta[MatchedIndex[1]], jet_phi[MatchedIndex[1]], jet_energy[MatchedIndex[1]]);
		TLorentzVector* q2 = new TLorentzVector();
		q2 -> SetPtEtaPhiE(jet_pt[MatchedIndex[2]], jet_eta[MatchedIndex[2]], jet_phi[MatchedIndex[2]], jet_energy[MatchedIndex[2]]);
		TLorentzVector* b2 = new TLorentzVector();
		b2 -> SetPtEtaPhiE(jet_pt[MatchedIndex[3]], jet_eta[MatchedIndex[3]], jet_phi[MatchedIndex[3]], jet_energy[MatchedIndex[3]]);
		TLorentzVector* q3 = new TLorentzVector();
		q3 -> SetPtEtaPhiE(jet_pt[MatchedIndex[4]], jet_eta[MatchedIndex[4]], jet_phi[MatchedIndex[4]], jet_energy[MatchedIndex[4]]);
		TLorentzVector* q4 = new TLorentzVector();
		q4 -> SetPtEtaPhiE(jet_pt[MatchedIndex[5]], jet_eta[MatchedIndex[5]], jet_phi[MatchedIndex[5]], jet_energy[MatchedIndex[5]]);


		float DeltaEta1 = jet_eta[MatchedIndex[1]] - jet_eta[MatchedIndex[2]];
		float DeltaEta2 = jet_eta[MatchedIndex[4]] - jet_eta[MatchedIndex[5]];

		for(int i=0; i<6; i++)
		{	for(int j=i+1; j<6; j++)
			{
				if(i==1 && j==2) continue;
				if(i==4 && j==5) continue;
				float DeltaEta = jet_eta[MatchedIndex[i]] - jet_eta[MatchedIndex[j]];
				RandomDeltaEta -> Fill(DeltaEta, weight*lumiFactor);
				float R = compute_R(jet_eta[MatchedIndex[i]], jet_eta[MatchedIndex[j]], jet_phi[MatchedIndex[i]], jet_phi[MatchedIndex[j]]);
				RandomDeltaR -> Fill(R, weight*lumiFactor);
				float delta_phi = fabs(jet_phi[MatchedIndex[i]] - jet_phi[MatchedIndex[j]]);
				if(delta_phi > 3.14159265359)
					delta_phi =2*3.14159265359 - delta_phi;
				RandomDeltaPhi -> Fill(delta_phi, weight*lumiFactor);
	
				TLorentzVector* q1a = new TLorentzVector();
				TLorentzVector* q2a = new TLorentzVector();

				q1a -> SetPtEtaPhiE(jet_pt[MatchedIndex[i]], jet_eta[MatchedIndex[i]], jet_phi[MatchedIndex[i]], jet_energy[MatchedIndex[i]]);
				q2a -> SetPtEtaPhiE(jet_pt[MatchedIndex[j]], jet_eta[MatchedIndex[j]], jet_phi[MatchedIndex[j]], jet_energy[MatchedIndex[j]]);

				RandomPt -> Fill((*q1a + *q2a).Pt(), weight*lumiFactor);
			}
		}


		float DeltaEtaTop1 = (*q1 + *q2).Eta() - b1->Eta(); 
		TopDeltaEta -> Fill(DeltaEtaTop1,  weight*lumiFactor);
		float DeltaEtaTop2 = (*q3 + *q4).Eta() - b2->Eta(); 
		TopDeltaEta -> Fill(DeltaEtaTop2,  weight*lumiFactor);


		float DeltaPhiTop1 = (*q1 + *q2).Phi() - b1->Phi(); 
		if(DeltaPhiTop1 > 3.14159265359)
			DeltaPhiTop1 =2*3.14159265359 - DeltaPhiTop1;
		TopDeltaPhi -> Fill(DeltaPhiTop1,  weight*lumiFactor);

		float DeltaPhiTop2 = (*q3 + *q4).Phi() - b2->Phi(); 
		if(DeltaPhiTop2 > 3.14159265359)
			DeltaPhiTop2 =2*3.14159265359 - DeltaPhiTop2;
		TopDeltaPhi -> Fill(DeltaPhiTop2,  weight*lumiFactor);

		WMass1 -> Fill((*q1 + *q2).M(), weight*lumiFactor);	
		WMass2 -> Fill((*q3 + *q4).M(), weight*lumiFactor);
		TopMass1 -> Fill((*b1 + *q1 + *q2).M(), weight*lumiFactor);
		TopMass2 -> Fill((*b2 + *q3 + *q4).M(), weight*lumiFactor);
		bDiscriminant -> Fill(jet_bdiscriminant[MatchedIndex[0]], weight*lumiFactor);
		bDiscriminant -> Fill(jet_bdiscriminant[MatchedIndex[3]], weight*lumiFactor);
		WbDiscriminant -> Fill(jet_bdiscriminant[MatchedIndex[1]], weight*lumiFactor);
		WbDiscriminant -> Fill(jet_bdiscriminant[MatchedIndex[2]], weight*lumiFactor);
		WbDiscriminant -> Fill(jet_bdiscriminant[MatchedIndex[4]], weight*lumiFactor);
		WbDiscriminant -> Fill(jet_bdiscriminant[MatchedIndex[5]], weight*lumiFactor);
		WDeltaEta -> Fill(DeltaEta1, weight*lumiFactor);
		WDeltaEta -> Fill(DeltaEta2, weight*lumiFactor);
		float R1 = compute_R(jet_eta[MatchedIndex[1]], jet_eta[MatchedIndex[2]], jet_phi[MatchedIndex[1]], jet_phi[MatchedIndex[2]]);
		float R2 = compute_R(jet_eta[MatchedIndex[4]], jet_eta[MatchedIndex[5]], jet_phi[MatchedIndex[4]], jet_phi[MatchedIndex[5]]);
		WDeltaR -> Fill(R1, weight*lumiFactor);
		WDeltaR -> Fill(R2, weight*lumiFactor);
		WMassDeltaR -> Fill((*q1 + *q2).M(), R1, weight*lumiFactor);
		WMassDeltaR -> Fill((*q3 + *q4).M(), R2, weight*lumiFactor);
		float delta_phi1 = fabs(jet_phi[MatchedIndex[1]] - jet_phi[MatchedIndex[2]]);
		if(delta_phi1 > 3.14159265359)
			delta_phi1 =2*3.14159265359 - delta_phi1;
		WDeltaPhi -> Fill(delta_phi1, weight*lumiFactor);
		float delta_phi2 = fabs(jet_phi[MatchedIndex[4]] - jet_phi[MatchedIndex[5]]);
		if(delta_phi2 > 3.14159265359)
			delta_phi2 =2*3.14159265359 - delta_phi2;
		WDeltaPhi -> Fill(delta_phi2, weight*lumiFactor);
		WPt -> Fill((*q1 + *q2).Pt(), weight*lumiFactor);
		WPt -> Fill((*q3 + *q4).Pt(), weight*lumiFactor);

		bool MatchedTop1 = 0;
		bool MatchedTop2 = 0;
		int NtopCandidate = 0;
		int NCorrectTop = 0;

	// Sort jets by bDiscriminant
		std::vector<int> SortedBdiscriminantIndexes;
		std::vector<int> tmp;
		for(int i=0; i<9; i++)
		{
			if(jet_pt[i]>0.)
				tmp.push_back(i);
		}

		while(tmp.size()>0)
		{
			float maxB = -100.;
			int maxIdx = -1;

			for(unsigned int j=0; j<tmp.size(); j++)
			{
				if(jet_bdiscriminant[tmp[j]]>maxB)
				{
					maxB = jet_bdiscriminant[tmp[j]];
					maxIdx = tmp[j];
				}
			
			}

			SortedBdiscriminantIndexes.push_back(maxIdx);
			std::vector<int> tmp2;

			for(unsigned int idx=0; idx<tmp.size(); idx++)
			{
				if(tmp[idx]!=maxIdx)
					tmp2.push_back(tmp[idx]);
			}

			tmp.clear();
			tmp = tmp2;
			tmp2.clear();
		}


		int isB1 = 0;
		int isB2 = 0;

		if(SortedBdiscriminantIndexes.size()>0 && ((SortedBdiscriminantIndexes[0]==MatchedIndex[0]) || (SortedBdiscriminantIndexes[0]==MatchedIndex[3])))
			isB1 = 1;
		if(SortedBdiscriminantIndexes.size()>1 && ((SortedBdiscriminantIndexes[1]==MatchedIndex[0]) || (SortedBdiscriminantIndexes[1]==MatchedIndex[3])))
			isB2 = 1;

		CorrectbCandidate -> Fill(isB1 + isB2);

		std::vector<int> Wcandidates;



		for(int TopCounter=0; TopCounter<2; TopCounter++)
		{
			for(unsigned int i=1; i<SortedBdiscriminantIndexes.size(); i++)
			{
				for(unsigned int j=i+1; j<SortedBdiscriminantIndexes.size(); j++)
				{
					float DeltaEta = jet_eta[SortedBdiscriminantIndexes[i]] - jet_eta[SortedBdiscriminantIndexes[j]];
					if(abs(DeltaEta)>WDeltaEtaThreshold) continue;
					float delta_phi = fabs(jet_phi[SortedBdiscriminantIndexes[i]] - jet_phi[SortedBdiscriminantIndexes[j]]);
					if(delta_phi > 3.14159265359)
						delta_phi =2*3.14159265359 - delta_phi;
					if(delta_phi>WDeltaPhiThreshold) continue;

					TLorentzVector* q1 = new TLorentzVector();
					TLorentzVector* q2 = new TLorentzVector();

					q1 -> SetPtEtaPhiE(jet_pt[SortedBdiscriminantIndexes[i]], jet_eta[SortedBdiscriminantIndexes[i]], jet_phi[SortedBdiscriminantIndexes[i]], jet_energy[SortedBdiscriminantIndexes[i]]);
		                        q2 -> SetPtEtaPhiE(jet_pt[SortedBdiscriminantIndexes[j]], jet_eta[SortedBdiscriminantIndexes[j]], jet_phi[SortedBdiscriminantIndexes[j]], jet_energy[SortedBdiscriminantIndexes[j]]);

					float mW = fabs((*q1 + *q2).M());
					float pTW = (*q1 + *q2).Pt();

					if(mW<WMassLowThreshold || mW>WMassHighThreshold) continue;
					if(pTW<WPtThreshold) continue;

					Wcandidates.push_back(SortedBdiscriminantIndexes[i]);
					Wcandidates.push_back(SortedBdiscriminantIndexes[j]);
				}
			}



			NWCandidate -> Fill(Wcandidates.size()/2);

			std::vector<int> SortedWCandidates;
			tmp.clear();
			tmp = Wcandidates;
			while(tmp.size()>0 && SortedWCandidates.size()<6)
			{
				float maxPt = -100.;
				int maxIdx = -1;

				for(unsigned int j=0; j<tmp.size(); j=j+2)
				{
					TLorentzVector* q1 = new TLorentzVector();
					TLorentzVector* q2 = new TLorentzVector();

					q1 -> SetPtEtaPhiE(jet_pt[tmp[j]], jet_eta[tmp[j]], jet_phi[tmp[j]], jet_energy[tmp[j]]);
					q2 -> SetPtEtaPhiE(jet_pt[tmp[j+1]], jet_eta[tmp[j+1]], jet_phi[tmp[j+1]], jet_energy[tmp[j+1]]);

					if((*q1 + *q2).Pt()>maxPt)
					{
						maxPt = (*q1 + *q2).Pt();
						maxIdx = j;
					}
			
				}

				SortedWCandidates.push_back(tmp[maxIdx]);
				SortedWCandidates.push_back(tmp[maxIdx+1]);

				std::vector<int> tmp2;

				for(unsigned int idx=0; idx<tmp.size(); idx=idx+2)
				{
					if(idx!=maxIdx)
					{	tmp2.push_back(tmp[idx]);
						tmp2.push_back(tmp[idx +1]);
					}
				}


				tmp.clear();
				tmp = tmp2;
				tmp2.clear();
			}

			int pTCandidate = 0;

			if(pTCandidate==0 && SortedWCandidates.size()>0 && ((SortedWCandidates[0]==MatchedIndex[1] && SortedWCandidates[1]==MatchedIndex[2]) || (SortedWCandidates[0]==MatchedIndex[2] && SortedWCandidates[1]==MatchedIndex[1])))
				pTCandidate = 1;
			if(pTCandidate==0 && SortedWCandidates.size()>0 && ((SortedWCandidates[0]==MatchedIndex[4] && SortedWCandidates[1]==MatchedIndex[5]) || (SortedWCandidates[0]==MatchedIndex[5] && SortedWCandidates[1]==MatchedIndex[4])))
				pTCandidate = 1;	
			if(pTCandidate==0 && SortedWCandidates.size()>2 && ((SortedWCandidates[2]==MatchedIndex[1] && SortedWCandidates[3]==MatchedIndex[2]) || (SortedWCandidates[2]==MatchedIndex[2] && SortedWCandidates[3]==MatchedIndex[1])))
				pTCandidate = 2;
			if(pTCandidate==0 && SortedWCandidates.size()>2 && ((SortedWCandidates[2]==MatchedIndex[4] && SortedWCandidates[3]==MatchedIndex[5]) || (SortedWCandidates[2]==MatchedIndex[5] && SortedWCandidates[3]==MatchedIndex[4])))
				pTCandidate = 2;	
			if(pTCandidate==0 && SortedWCandidates.size()>4 && ((SortedWCandidates[4]==MatchedIndex[1] && SortedWCandidates[5]==MatchedIndex[2]) || (SortedWCandidates[4]==MatchedIndex[2] && SortedWCandidates[5]==MatchedIndex[1])))
				pTCandidate = 3;
			if(pTCandidate==0 && SortedWCandidates.size()>4 && ((SortedWCandidates[4]==MatchedIndex[4] && SortedWCandidates[5]==MatchedIndex[5]) || (SortedWCandidates[4]==MatchedIndex[5] && SortedWCandidates[5]==MatchedIndex[4])))
				pTCandidate = 3;	

			CorrectWPosition -> Fill(pTCandidate);


			bool isTop = 0;
			bool isCorrectTop = 0;
			std::vector<int> TopIndex;
			for(unsigned int i=0; i< SortedWCandidates.size(); i=i+2)
			{
				if(jet_bdiscriminant[SortedBdiscriminantIndexes[0]]<bDiscriminantThreshold) continue;

				TLorentzVector* b = new TLorentzVector();
				b -> SetPtEtaPhiE(jet_pt[SortedBdiscriminantIndexes[0]], jet_eta[SortedBdiscriminantIndexes[0]], jet_phi[SortedBdiscriminantIndexes[0]], jet_energy[SortedBdiscriminantIndexes[0]]);

				TLorentzVector* q1 = new TLorentzVector();
				TLorentzVector* q2 = new TLorentzVector();

				q1 -> SetPtEtaPhiE(jet_pt[SortedWCandidates[i]], jet_eta[SortedWCandidates[i]], jet_phi[SortedWCandidates[i]], jet_energy[SortedWCandidates[i]]);
				q2 -> SetPtEtaPhiE(jet_pt[SortedWCandidates[i+1]], jet_eta[SortedWCandidates[i+1]], jet_phi[SortedWCandidates[i+1]], jet_energy[SortedWCandidates[i+1]]);

				TLorentzVector W = (*q1 + *q2);
				float mT = (*b + W).M();
				float pTT = (*b + W).Pt();
				float deltaEtaT = b->Eta() - W.Eta();
				float deltaPhiT = fabs(b->Phi() - W.Phi());
				if(deltaPhiT > 3.14159265359)
					deltaPhiT =2*3.14159265359 - deltaPhiT;
			
				if(mT < topMassLowThreshold) continue;
				if(mT > topMassHighThreshold) continue;
				if(fabs(deltaEtaT)>topDeltaEtaThreshold) continue;
				if(deltaPhiT>topDeltaPhiThreshold) continue;
				if(pTT<topPtThreshold) continue;

				isTop = 1;
				if(SortedBdiscriminantIndexes[0]==MatchedIndex[0] && ((SortedWCandidates[i]==MatchedIndex[1] && SortedWCandidates[i+1]==MatchedIndex[2]) || (SortedWCandidates[i]==MatchedIndex[2] && SortedWCandidates[i+1]==MatchedIndex[1])))
				       isCorrectTop = 1;
			
				if(SortedBdiscriminantIndexes[0]==MatchedIndex[3] && ((SortedWCandidates[i]==MatchedIndex[4] && SortedWCandidates[i+1]==MatchedIndex[5]) || (SortedWCandidates[i]==MatchedIndex[5] && SortedWCandidates[i+1]==MatchedIndex[4])))
				       isCorrectTop = 1;

				TopIndex.push_back(SortedBdiscriminantIndexes[0]);
				TopIndex.push_back(SortedWCandidates[i]);
				TopIndex.push_back(SortedWCandidates[i+1]);

				break;
			}

			if(isTop)
			{
				NtopCandidate++;
				//remove used jets from vector

				tmp = SortedBdiscriminantIndexes;
				SortedBdiscriminantIndexes.clear();

				for(unsigned int i=1; i<tmp.size(); i++)
				{
					bool isUsed = 0;
					for(unsigned int j=1; j<TopIndex.size(); j++)
					{
						if(tmp[i] == TopIndex[j])
							isUsed = 1;
					}
					if(!isUsed)
						SortedBdiscriminantIndexes.push_back(tmp[i]);
				}
			}

			else
			{
				//discrad the first b-jet

				tmp = SortedBdiscriminantIndexes;
				SortedBdiscriminantIndexes.clear();
					
				for(unsigned int i=1; i<tmp.size(); i++)
					SortedBdiscriminantIndexes.push_back(tmp[i]);
			}

			if(isCorrectTop)
				NCorrectTop++;

		}

/*

		for(unsigned int i=2; i<SortedBdiscriminantIndexes.size(); i++)
		{
			for(unsigned int j=i+1; j<SortedBdiscriminantIndexes.size(); j++)
			{
				float DeltaEta = jet_eta[SortedBdiscriminantIndexes[i]] - jet_eta[SortedBdiscriminantIndexes[j]];
				if(abs(DeltaEta)>1.5) continue;
				float delta_phi = fabs(jet_phi[SortedBdiscriminantIndexes[i]] - jet_phi[SortedBdiscriminantIndexes[j]]);
				if(delta_phi > 3.14159265359)
					delta_phi =2*3.14159265359 - delta_phi;
				if(delta_phi>2.) continue;

				//float R = DeltaEta*DeltaEta+delta_phi*delta_phi;
				//if(R>2.) continue;				

				TLorentzVector* q1 = new TLorentzVector();
				TLorentzVector* q2 = new TLorentzVector();

				q1 -> SetPtEtaPhiE(jet_pt[SortedBdiscriminantIndexes[i]], jet_eta[SortedBdiscriminantIndexes[i]], jet_phi[SortedBdiscriminantIndexes[i]], jet_energy[SortedBdiscriminantIndexes[i]]);
                                q2 -> SetPtEtaPhiE(jet_pt[SortedBdiscriminantIndexes[j]], jet_eta[SortedBdiscriminantIndexes[j]], jet_phi[SortedBdiscriminantIndexes[j]], jet_energy[SortedBdiscriminantIndexes[j]]);

				float mW = fabs((*q1 + *q2).M());

				if(mW<60 || mW>110) continue;

				Wcandidates.push_back(SortedBdiscriminantIndexes[i]);
				Wcandidates.push_back(SortedBdiscriminantIndexes[j]);
			}
		}


		NWCandidate -> Fill(Wcandidates.size()/2);

		std::vector<int> SortedWCandidates;
		tmp.clear();
		tmp = Wcandidates;
		while(tmp.size()>0 && SortedWCandidates.size()<6)
		{
			float maxPt = -100.;
			int maxIdx = -1;

			for(unsigned int j=0; j<tmp.size(); j=j+2)
			{
				TLorentzVector* q1 = new TLorentzVector();
				TLorentzVector* q2 = new TLorentzVector();

				q1 -> SetPtEtaPhiE(jet_pt[tmp[j]], jet_eta[tmp[j]], jet_phi[tmp[j]], jet_energy[tmp[j]]);
				q2 -> SetPtEtaPhiE(jet_pt[tmp[j+1]], jet_eta[tmp[j+1]], jet_phi[tmp[j+1]], jet_energy[tmp[j+1]]);

				if((*q1 + *q2).Pt()>maxPt)
				{
					maxPt = (*q1 + *q2).Pt();
					maxIdx = j;
				}
			
			}

			SortedWCandidates.push_back(tmp[maxIdx]);
			SortedWCandidates.push_back(tmp[maxIdx+1]);

			std::vector<int> tmp2;

			for(unsigned int idx=0; idx<tmp.size(); idx=idx+2)
			{
				if(idx!=maxIdx)
				{	tmp2.push_back(tmp[idx]);
					tmp2.push_back(tmp[idx +1]);
				}
			}


			tmp.clear();
			tmp = tmp2;
			tmp2.clear();
		}

		int pTCandidate = 0;

		if(pTCandidate==0 && SortedWCandidates.size()>0 && ((SortedWCandidates[0]==MatchedIndex[1] && SortedWCandidates[1]==MatchedIndex[2]) || (SortedWCandidates[0]==MatchedIndex[2] && SortedWCandidates[1]==MatchedIndex[1])))
			pTCandidate = 1;
		if(pTCandidate==0 && SortedWCandidates.size()>0 && ((SortedWCandidates[0]==MatchedIndex[4] && SortedWCandidates[1]==MatchedIndex[5]) || (SortedWCandidates[0]==MatchedIndex[5] && SortedWCandidates[1]==MatchedIndex[4])))
			pTCandidate = 1;	
		if(pTCandidate==0 && SortedWCandidates.size()>2 && ((SortedWCandidates[2]==MatchedIndex[1] && SortedWCandidates[3]==MatchedIndex[2]) || (SortedWCandidates[2]==MatchedIndex[2] && SortedWCandidates[3]==MatchedIndex[1])))
			pTCandidate = 2;
		if(pTCandidate==0 && SortedWCandidates.size()>2 && ((SortedWCandidates[2]==MatchedIndex[4] && SortedWCandidates[3]==MatchedIndex[5]) || (SortedWCandidates[2]==MatchedIndex[5] && SortedWCandidates[3]==MatchedIndex[4])))
			pTCandidate = 2;	
		if(pTCandidate==0 && SortedWCandidates.size()>4 && ((SortedWCandidates[4]==MatchedIndex[1] && SortedWCandidates[5]==MatchedIndex[2]) || (SortedWCandidates[4]==MatchedIndex[2] && SortedWCandidates[5]==MatchedIndex[1])))
			pTCandidate = 3;
		if(pTCandidate==0 && SortedWCandidates.size()>4 && ((SortedWCandidates[4]==MatchedIndex[4] && SortedWCandidates[5]==MatchedIndex[5]) || (SortedWCandidates[4]==MatchedIndex[5] && SortedWCandidates[5]==MatchedIndex[4])))
			pTCandidate = 3;	

		CorrectWPosition -> Fill(pTCandidate);


		for(unsigned int i=0; i< SortedWCandidates.size(); i=i+2)
		{
			for(unsigned int j=0; j< 2; j++)
			{
				TLorentzVector* b = new TLorentzVector();
				b -> SetPtEtaPhiE(jet_pt[SortedBdiscriminantIndexes[j]], jet_eta[SortedBdiscriminantIndexes[j]], jet_phi[SortedBdiscriminantIndexes[j]], jet_energy[SortedBdiscriminantIndexes[j]]);

				TLorentzVector* q1 = new TLorentzVector();
				TLorentzVector* q2 = new TLorentzVector();

				q1 -> SetPtEtaPhiE(jet_pt[SortedWCandidates[i]], jet_eta[SortedWCandidates[i]], jet_phi[SortedWCandidates[i]], jet_energy[SortedWCandidates[i]]);
				q2 -> SetPtEtaPhiE(jet_pt[SortedWCandidates[i+1]], jet_eta[SortedWCandidates[i+1]], jet_phi[SortedWCandidates[i+1]], jet_energy[SortedWCandidates[i+1]]);

				TLorentzVector W = (*q1 + *q2);
				float mT = (*b + W).M();
				float pTT = (*b + W).Pt();
				float etaT = (*b + W).Eta();
				float delta_phi = fabs(b->Phi() - W.Phi());
				if(delta_phi > 3.14159265359)
					delta_phi =2*3.14159265359 - delta_phi;


				if(SortedBdiscriminantIndexes[j]==MatchedIndex[0] && ((SortedWCandidates[i]==MatchedIndex[1] && SortedWCandidates[i+1]==MatchedIndex[2]) || (SortedWCandidates[i]==MatchedIndex[2] && SortedWCandidates[i+1]==MatchedIndex[1])))
				{	RecoTopMass -> Fill(mT, weight*lumiFactor);
				        RecoTopPt -> Fill(pTT, weight*lumiFactor);
				        RecoTopEta -> Fill(etaT, weight*lumiFactor);
					RecoTopDeltaEta -> Fill(b->Eta() - W.Eta(), weight*lumiFactor);
					RecoTopDeltaPhi -> Fill(delta_phi, weight*lumiFactor);
				        MatchedTop1 = 1;
				}
				else if(SortedBdiscriminantIndexes[j]==MatchedIndex[3] && ((SortedWCandidates[i]==MatchedIndex[4] && SortedWCandidates[i+1]==MatchedIndex[5]) || (SortedWCandidates[i]==MatchedIndex[5] && SortedWCandidates[i+1]==MatchedIndex[4])))
				{	RecoTopMass -> Fill(mT, weight*lumiFactor);
				        RecoTopPt -> Fill(pTT, weight*lumiFactor);
				        RecoTopEta -> Fill(etaT, weight*lumiFactor);
					RecoTopDeltaEta -> Fill(b->Eta() - W.Eta(), weight*lumiFactor);
					RecoTopDeltaPhi -> Fill(delta_phi, weight*lumiFactor);
				        MatchedTop2 = 1;
				}
				else
				{	RecoTopMassRandom -> Fill(mT, weight*lumiFactor);
				        RecoTopPtRandom -> Fill(pTT, weight*lumiFactor);
				        RecoTopEtaRandom -> Fill(etaT, weight*lumiFactor);
					RecoTopDeltaEtaRandom -> Fill(b->Eta() - W.Eta(), weight*lumiFactor);
					RecoTopDeltaPhiRandom -> Fill(delta_phi, weight*lumiFactor);
				}
			}
		}
		if(MatchedTop1 && MatchedTop2)
		  TwoTops++;
		if(!(MatchedTop1 && MatchedTop2) && (MatchedTop1 || MatchedTop2))
		  OneTop++;



		for(unsigned int i=0; i< SortedWCandidates.size(); i=i+2)
		{
			for(unsigned int j=0; j< 2; j++)
			{
			        if(jet_bdiscriminant[SortedBdiscriminantIndexes[j]]<bDiscriminantThreshold) continue;

				TLorentzVector* b = new TLorentzVector();
				b -> SetPtEtaPhiE(jet_pt[SortedBdiscriminantIndexes[j]], jet_eta[SortedBdiscriminantIndexes[j]], jet_phi[SortedBdiscriminantIndexes[j]], jet_energy[SortedBdiscriminantIndexes[j]]);

				TLorentzVector* q1 = new TLorentzVector();
				TLorentzVector* q2 = new TLorentzVector();

				q1 -> SetPtEtaPhiE(jet_pt[SortedWCandidates[i]], jet_eta[SortedWCandidates[i]], jet_phi[SortedWCandidates[i]], jet_energy[SortedWCandidates[i]]);
				q2 -> SetPtEtaPhiE(jet_pt[SortedWCandidates[i+1]], jet_eta[SortedWCandidates[i+1]], jet_phi[SortedWCandidates[i+1]], jet_energy[SortedWCandidates[i+1]]);

				TLorentzVector W = (*q1 + *q2);
				float mT = (*b + W).M();
				float pTT = (*b + W).Pt();
				float deltaEtaT = b->Eta() - W.Eta();
				float delta_phi = fabs(b->Phi() - W.Phi());
				if(delta_phi > 3.14159265359)
					delta_phi =2*3.14159265359 - delta_phi;
				
				if(mT < topMassLowThreshold) continue;
				if(mT > topMassHighThreshold) continue;
				if(fabs(deltaEtaT)>1.5) continue;
				if(delta_phi>2.) continue;
				if(pTT<topPtThreshold) continue;

				NtopCandidate++;
				if(SortedBdiscriminantIndexes[j]==MatchedIndex[0] && ((SortedWCandidates[i]==MatchedIndex[1] && SortedWCandidates[i+1]==MatchedIndex[2]) || (SortedWCandidates[i]==MatchedIndex[2] && SortedWCandidates[i+1]==MatchedIndex[1])))
				       NCorrectTop++;
				
				else if(SortedBdiscriminantIndexes[j]==MatchedIndex[3] && ((SortedWCandidates[i]==MatchedIndex[4] && SortedWCandidates[i+1]==MatchedIndex[5]) || (SortedWCandidates[i]==MatchedIndex[5] && SortedWCandidates[i+1]==MatchedIndex[4])))
				       NCorrectTop++;
			}
		}
*/		
		NTopCandidate -> Fill(NtopCandidate);
		NTopCandidateCorrect -> Fill(NCorrectTop);


	}

	cout << "Processed " << nentries << " events out of " << nentries << endl;

	cout << "Events with both top quarks selected: " << TwoTops << endl;
	cout << "Events with only one top quark selected " << OneTop << endl;

	MakePlot(Dr, "Dr");
	MakePlot(WMass1, "WMass1");
	MakePlot(WMass2, "WMass2");
	MakePlot(TopMass1, "TopMass1");
	MakePlot(TopMass2, "TopMass2");
	MakePlot(TopDeltaEta, "TopDeltaEta");
	MakePlot(TopDeltaPhi, "TopDeltaPhi");
	MakePlot(NWCandidate, "NWCandidate");
	Make2Plot(NTopCandidate, NTopCandidateCorrect, "NTopCandidate", "Top Candidates", "True Top Candidates");
	MakePlot(CorrectbCandidate, "CorrectbCandidate");
	MakePlot(CorrectWPosition, "CorrectWPosition");
	Make2Plot(bDiscriminant, WbDiscriminant, "bDiscriminant", "b-jets", "W-jets");
	Make2Plot(WDeltaEta, RandomDeltaEta, "WDeltaEta", "W-jets", "Random", 1);
	Make2Plot(WDeltaPhi, RandomDeltaPhi, "WDeltaPhi", "W-jets", "Random", 1);
	Make2Plot(WDeltaR, RandomDeltaR, "WDeltaR", "W-jets", "Random", 1);
	Make2Plot(WPt, RandomPt, "WPt", "W-jets", "Random", 1);
	Make2Plot(RecoTopMass, RecoTopMassRandom, "RecoTopMass", "True top", "Random", 1);
//	Make2Plot(RecoTopPt, RecoTopPtRandom, "RecoTopPt", "True top", "Random", 1);
//	Make2Plot(RecoTopEta, RecoTopEtaRandom, "RecoTopEta", "True top", "Random", 1);
//	Make2Plot(RecoTopDeltaEta, RecoTopDeltaEtaRandom, "RecoTopDeltaEta", "True top", "Random", 1);
//	Make2Plot(RecoTopDeltaPhi, RecoTopDeltaPhiRandom, "RecoTopDeltaPhi", "True top", "Random", 1);
	MakePlot2D(WMassDeltaR, "WMassDeltaR");


	;


	system("mv *.png ~/www/ttH/Full2016Dataset/Hadronic/ttH");
	system("mv *.pdf ~/www/ttH/Full2016Dataset/Hadronic/ttH");





	cout << endl << "Lufthansa, partner of Star Alliance, thanks you for choosing our company, we hope to see you again on board of our airctafts" << endl << endl;



























}









     






