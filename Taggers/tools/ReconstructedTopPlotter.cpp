/*
g++ -Wall -o ReconstructedTopPlotter `root-config --cflags --glibs` -L $ROOTSYS/lib -lRooFitCore -lFoam -lMinuit -lMathMore CMS_lumi.C tdrstyle.C ReconstructedTopPlotter.cpp
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

	TH1F* NTopCandidate = new TH1F("NTopCandidate", "; N^{top}; Counts", 3, -0.5, 2.5);

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

	int nentries = ttH -> GetEntries();


	int passCutBasedSelection = 0;
	int passMVASelection = 0;
	int passTopSelection = 0;


	//Cut Based Efficiency:
	for(int ev=0; ev<nentries; ev++)
	{
		ttH -> GetEntry(ev);
		if(ev%10000==0) cout << "Processing event " << ev << " out of " << nentries << "\r" << flush;

		if(dipho_leadPt<dipho_mass/3. || dipho_subleadPt<dipho_mass/4.) continue;

		int nJet = 0;
		int looseBjet = 0;
		int mediumBjet = 0;

		for(int i=0; i<9; i++)
		{
			if(jet_pt[i]<25 || fabs(jet_eta[i])>2.4) continue;
			float deltaR1 = compute_R(jet_eta[i], dipho_leadEta, jet_phi[i], dipho_leadPhi);
			float deltaR2 = compute_R(jet_eta[i], dipho_subleadEta, jet_phi[i], dipho_subleadPhi);
			if(deltaR1<0.4 || deltaR2<0.4) continue;

			nJet++;
			if(jet_bdiscriminant[i]>0.6)
				looseBjet++;
			if(jet_bdiscriminant[i]>0.8)
				mediumBjet++;
		}

		if(mediumBjet>=1 && nJet>=5 && dipho_leadPt>dipho_mass/2.)
			passCutBasedSelection++;

		if(looseBjet>=1 && nJet>=3 && ttHMVA>0.75)
			passMVASelection++;
	}	


	for(int ev=0; ev<nentries; ev++)
	{	
		ttH -> GetEntry(ev);
		if(ev%10000==0) cout << "Processing event " << ev << " out of " << nentries << "\r" << flush;
//		cout << "Processing event " << ev << " out of " << nentries << endl;

		int NtopCandidate = 0;

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

		if(SortedBdiscriminantIndexes.size() < 3) continue;
		if(jet_bdiscriminant[SortedBdiscriminantIndexes[0]]<bDiscriminantThreshold) continue;


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

			std::vector<int> SortedWCandidates;
			tmp.clear();
			tmp = Wcandidates;
			while(tmp.size()>0 && SortedWCandidates.size()<6)
			{
				float maxPt = -100.;
				unsigned int maxIdx = 0;

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


			bool isTop = 0;

			std::vector<int> TopIndex;
			for(unsigned int i=0; i < SortedWCandidates.size(); i=i+2)
			{
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
		}
	
		NTopCandidate -> Fill(NtopCandidate, weight*lumiFactor);
		if(NtopCandidate>0)
			passTopSelection++;
	}

	cout << "Processed " << nentries << " events out of " << nentries << endl;

	MakePlot(NTopCandidate, "NTop_ttH");

	cout << "Cut based efficiency: " <<  (float)passCutBasedSelection/nentries << endl;
	cout << "MVA efficiency: " <<  (float)passMVASelection/nentries << endl;
	cout << "Top selection efficiency: " <<  (float)passTopSelection/nentries << endl;

	system("mv *.png ~/www/ttH/Full2016Dataset/Hadronic/ttH");
	system("mv *.pdf ~/www/ttH/Full2016Dataset/Hadronic/ttH");





	cout << endl << "Lufthansa, partner of Star Alliance, thanks you for choosing our company, we hope to see you again on board of our airctafts" << endl << endl;



























}









     






