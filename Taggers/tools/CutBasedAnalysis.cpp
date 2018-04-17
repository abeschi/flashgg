/*
g++ -Wall -o CutBasedAnalysis `root-config --cflags --glibs` -L $ROOTSYS/lib -lRooFitCore -lFoam -lMinuit -lMathMore CMS_lumi.C tdrstyle.C AnalysisUtilities.cc CutBasedAnalysis.cpp
*/

#ifndef CMS_LUMI_H
#include "CMS_lumi.h"
#endif

#ifndef CMS_STYLE
#include "tdrstyle.h"
#endif

#ifndef ANALYSIS_UTILITIES_H
#include "AnalysisUtilities.h"
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
#include <fstream>

using namespace std;


int main(int argc, char *argv[])
{
	writeExtraText = true;       // if extra text
	extraText  = "Preliminary";  // default extra text is "Preliminary"
	lumi_sqrtS = "35.9 fb^{-1} (13 TeV)";       // used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)

	setTDRStyle();
	gStyle -> SetOptFit(0);
	gStyle -> SetOptStat(0);

        TString SignalFold = "/afs/cern.ch/work/a/abeschi/ttH_MassTop_v5/";
        TString BkgFold = "/afs/cern.ch/work/a/abeschi/ttH_MassTop_Bkg_v3/";
        TString DataFold = "/afs/cern.ch/work/a/abeschi/ttH_MassTop_Data_v2/";
        TString CSFold = "/afs/cern.ch/work/a/abeschi/ttH_CS_v1/";
        TString names[9] = {"ttH", "ggH", "VBF", "VH", "bbH", "tHq", "tHW", "Data", "Control Sample"};


	// Thresholds parsed form Paramenters.txt, check the file for more details
	// nJets, nLoosebjets, nMediumbjets, nTightbjets, Pt/mggLeading, Pt/mggSubleading, ptLeadingJet, ptSubleadingJet, ptSubsubleadingJet, EtaLeadingJet, EtaSubleadingJet, EtaSubsubleadingJet,
	// DeltaEtaLeadingSubleadingJet, DeltaEtaLeadinSubsbuleadingJet, DeltaEtaSubleadingJetSubsubleadingjet

	std::vector<float> thresholds = parse("Parameters.txt");

	//for(unsigned int i=0; i<thresholds.size(); i++)
	//  cout << thresholds.at(i) << endl;
	
        TChain* ttH;
        TChain* ggH;
        TChain* vbf;
        TChain* vh;
        TChain* bbH;
        TChain* tHq;
        TChain* tHW;

        TChain* data;
	TChain* controlSample;

        float lumiFactor = 35.9;
        float bDiscriminantThresholdLoose = 0.6;
        float bDiscriminantThresholdMedium = 0.8;
        float bDiscriminantThresholdTight = 0.97;

	float passCutBasedSelection[9] = {0., 0., 0., 0., 0., 0., 0., 0., 0.};
	float passMVASelection[9] = {0., 0., 0., 0., 0., 0., 0., 0., 0.,};
	float passNewSelection[9] = {0., 0., 0., 0., 0., 0., 0., 0., 0.};

        bool isLept = (argc!=1 ? 0 : 1);

        if(isLept)
	{       cout << "Processing leptonic tag" << endl;

	        ttH = new TChain("TTHLeptonicDumper/trees/tth_13TeV_all");
		ttH -> Add(SignalFold + "ttHJetToGG_M125_13TeV_*.root");
		ggH = new TChain("TTHLeptonicDumper/trees/ggh_13TeV_all");
		ggH -> Add(SignalFold + "GluGluHToGG_M125_13TeV_*.root");
		vbf = new TChain("TTHLeptonicDumper/trees/vbf_13TeV_all");
		vbf -> Add(SignalFold + "VBFHToGG_M125_13TeV*.root");
		vh = new TChain("TTHLeptonicDumper/trees/vh_13TeV_all");
		vh -> Add(SignalFold + "VHToGG_M125_13TeV_*.root");
		bbH = new TChain("TTHLeptonicDumper/trees/bbh_13TeV_all");
		bbH -> Add(SignalFold + "bbHToGG_M125_*.root");
		tHq = new TChain("TTHLeptonicDumper/trees/thq_13TeV_all");
		tHq -> Add(SignalFold + "THQ_HToGG_13TeV*.root");
		tHW = new TChain("TTHLeptonicDumper/trees/thw_13TeV_all");
		tHW -> Add(SignalFold + "THW_HToGG_13TeV*.root");

		data = new TChain("TTHLeptonicDumper/trees/Data_13TeV_all");
		data -> Add(DataFold + "output_DoubleEG_Run2016*.root");
		controlSample = new TChain("TTHLeptonicDumper/trees/Data_13TeV_all");
		controlSample -> Add(CSFold + "output_DoubleEG_Run2016*.root");
	}

        else
	{       cout << "Processing hadronic tag" << endl;

	        ttH = new TChain("TTHHadronicDumper/trees/tth_13TeV_all");
		ttH -> Add(SignalFold + "ttHJetToGG_M125_13TeV_*.root");
		ggH = new TChain("TTHHadronicDumper/trees/ggh_13TeV_all");
		ggH -> Add(SignalFold + "GluGluHToGG_M125_13TeV_*.root");
		vbf = new TChain("TTHHadronicDumper/trees/vbf_13TeV_all");
		vbf -> Add(SignalFold + "VBFHToGG_M125_13TeV*.root");
		vh = new TChain("TTHHadronicDumper/trees/vh_13TeV_all");
		vh -> Add(SignalFold + "VHToGG_M125_13TeV_*.root");
		bbH = new TChain("TTHHadronicDumper/trees/bbh_13TeV_all");
		bbH -> Add(SignalFold + "bbHToGG_M125_*.root");
		tHq = new TChain("TTHHadronicDumper/trees/thq_13TeV_all");
		tHq -> Add(SignalFold + "THQ_HToGG_13TeV*.root");
		tHW = new TChain("TTHHadronicDumper/trees/thw_13TeV_all");
		tHW -> Add(SignalFold + "THW_HToGG_13TeV*.root");

		data = new TChain("TTHHadronicDumper/trees/Data_13TeV_all");
		data -> Add(DataFold + "output_DoubleEG_Run2016*.root");
		controlSample = new TChain("TTHHadronicDumper/trees/Data_13TeV_all");
		controlSample -> Add(CSFold + "output_DoubleEG_Run2016*.root");
	}


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

	float ele_pt[2] = {0., 0.};
        float ele_eta[2] = {0., 0.};
        float ele_phi[2] = {0., 0.};

        float mu_pt[2] = {0., 0.};
        float mu_eta[2] = {0., 0.};
        float mu_phi[2] = {0., 0.};

        float MetPt = 0.;
        float MetPhi = 0.;
        float ttHMVA = 0.;



	TH1F* Signal_histo_CutBased = new TH1F("Signal_histo_CutBased", "; m_{#gamma#gamma} (GeV); Counts/GeV", 640, 100, 180);
	TH1F* Bkg_histo_CutBased = new TH1F("Bkg_histo_CutBased", "; m_{#gamma#gamma} (GeV); Counts/GeV", 80, 100, 180);
	TH1F* CS_histo_CutBased = new TH1F("CS_histo_CutBased", "; m_{#gamma#gamma} (GeV); Counts/GeV", 80, 100, 180);
	TH1F* Signal_histo_MVA = new TH1F("Signal_histo_MVA", "; m_{#gamma#gamma} (GeV); Counts/GeV", 640, 100, 180);
	TH1F* Bkg_histo_MVA = new TH1F("Bkg_histo_MVA", "; m_{#gamma#gamma} (GeV); Counts/GeV", 80, 100, 180);
	TH1F* CS_histo_MVA = new TH1F("CS_histo_MVA", "; m_{#gamma#gamma} (GeV); Counts/GeV", 80, 100, 180);
	TH1F* Signal_histo_New = new TH1F("Signal_histo_New", "; m_{#gamma#gamma} (GeV); Counts/GeV", 640, 100, 180);
	TH1F* Bkg_histo_New = new TH1F("Bkg_histo_New", "; m_{#gamma#gamma} (GeV); Counts/GeV", 80, 100, 180);
	TH1F* CS_histo_New = new TH1F("CS_histo_New", "; m_{#gamma#gamma} (GeV); Counts/GeV", 80, 100, 180);

	for(int n=0; n<9; n++)
	{
	        int nentries;
		TChain* serviceTree;
		switch (n)
		{
	                case(0):
		                  nentries = ttH -> GetEntries();
				  serviceTree = (TChain*)ttH -> Clone();
				  break;

	                case(1):
			          nentries = ggH -> GetEntries();
				  serviceTree = (TChain*)ggH -> Clone();
				  break;

		        case(2):
			          nentries = vbf -> GetEntries();
				  serviceTree = (TChain*)vbf -> Clone();
				  break;

	                case(3):
			          nentries = vh -> GetEntries();
				  serviceTree = (TChain*)vh -> Clone();
				  break;

	                case(4):
			          nentries = bbH -> GetEntries();
				  serviceTree = (TChain*)bbH -> Clone();
				  break;

	                case(5):
			          nentries = tHq -> GetEntries();
				  serviceTree = (TChain*)tHq -> Clone();
				  break;

	                case(6):
			          nentries = tHW -> GetEntries();
				  serviceTree = (TChain*)tHW -> Clone();
				  break;

		        case(7):
			          nentries = data -> GetEntries();
				  serviceTree = (TChain*)data -> Clone();
				  break;

		        case(8):
			          nentries = controlSample -> GetEntries();
				  serviceTree = (TChain*)controlSample -> Clone();
				  break;


		        default:
		                  nentries = 0;
		}

		serviceTree -> SetBranchAddress("nvtx", &nvtx);
		serviceTree -> SetBranchAddress("weight", &weight);
		serviceTree -> SetBranchAddress("dipho_sumpt", &dipho_sumpt);
		serviceTree -> SetBranchAddress("dipho_mass", &dipho_mass);
		serviceTree -> SetBranchAddress("dipho_leadPt", &dipho_leadPt);
		serviceTree -> SetBranchAddress("dipho_leadEta", &dipho_leadEta);
		serviceTree -> SetBranchAddress("dipho_leadPhi", &dipho_leadPhi);
		serviceTree -> SetBranchAddress("dipho_leadR9", &dipho_leadR9);
		serviceTree -> SetBranchAddress("dipho_subleadPt", &dipho_subleadPt);
		serviceTree -> SetBranchAddress("dipho_subleadEta", &dipho_subleadEta);
		serviceTree -> SetBranchAddress("dipho_subleadPhi", &dipho_subleadPhi);
		serviceTree -> SetBranchAddress("dipho_subleadR9", &dipho_subleadR9);
		serviceTree -> SetBranchAddress("dipho_leadIDMVA", &dipho_leadIDMVA);
		serviceTree -> SetBranchAddress("dipho_subleadIDMVA", &dipho_subleadIDMVA);
		serviceTree -> SetBranchAddress("dipho_mva", &dipho_mva);
		
		
		if(isLept)
		{
		          for(int i=1; i<6; i++)
			  {       serviceTree -> SetBranchAddress(("jet_pt"+ std::to_string(i)).c_str(), &jet_pt[i-1]);
			          serviceTree -> SetBranchAddress(("jet_eta"+ std::to_string(i)).c_str(), &jet_eta[i-1]);
				  serviceTree -> SetBranchAddress(("jet_phi"+ std::to_string(i)).c_str(), &jet_phi[i-1]);
				  serviceTree -> SetBranchAddress(("jet_bdiscriminant"+ std::to_string(i)).c_str(), &jet_bdiscriminant[i-1]);
			  }
                                                                                                                                                                                                    
			  for(int i=1; i<3; i++)
			  {       serviceTree -> SetBranchAddress(("ele_pt"+ std::to_string(i)).c_str(), &ele_pt[i-1]);
			          serviceTree -> SetBranchAddress(("ele_eta"+ std::to_string(i)).c_str(), &ele_eta[i-1]);
				  serviceTree -> SetBranchAddress(("ele_phi"+ std::to_string(i)).c_str(), &ele_phi[i-1]);                                                
				  serviceTree -> SetBranchAddress(("mu_pt"+ std::to_string(i)).c_str(), &mu_pt[i-1]);                                    
				  serviceTree -> SetBranchAddress(("mu_eta"+ std::to_string(i)).c_str(), &mu_eta[i-1]);                          
				  serviceTree -> SetBranchAddress(("mu_phi"+ std::to_string(i)).c_str(), &mu_phi[i-1]);
			  } 
			  serviceTree -> SetBranchAddress("MetPt", &MetPt); 
			  serviceTree -> SetBranchAddress("MetPhi", &MetPhi);
		}                                                                                                                     
		else
		{                                                                                                                                                                                   
	                for(int i=1; i<10; i++)
			{       serviceTree -> SetBranchAddress(("jet_pt"+ std::to_string(i)).c_str(), &jet_pt[i-1]);
		                serviceTree -> SetBranchAddress(("jet_eta"+ std::to_string(i)).c_str(), &jet_eta[i-1]);
				serviceTree -> SetBranchAddress(("jet_phi"+ std::to_string(i)).c_str(), &jet_phi[i-1]);
				serviceTree -> SetBranchAddress(("jet_bdiscriminant"+ std::to_string(i)).c_str(), &jet_bdiscriminant[i-1]);
			} 
			serviceTree -> SetBranchAddress("ttHMVA", &ttHMVA);
			serviceTree -> SetBranchAddress("MetPt", &MetPt);
			serviceTree -> SetBranchAddress("MetPhi", &MetPhi);                                                  
		}                                                                       

		for(int ev=0; ev<nentries; ev++)
		{
		        serviceTree -> GetEntry(ev);
			if(ev%10000==0) cout << "Processing tag " << names[n] << ", event " << ev << " out of " << nentries << "\r" << flush;

			if(n==7 && (dipho_mass>115 && dipho_mass<135)) continue;
			if(n==7 || n==8) weight=1./lumiFactor;

			if(dipho_leadPt<dipho_mass/3. || dipho_subleadPt<dipho_mass/4.) continue;

			int nJet = 0;
			int looseBjet = 0;
			int mediumBjet = 0;
			int tightBjet = 0;

			for(int i=0; i<9; i++)
			{
			        if(jet_pt[i]<25 || fabs(jet_eta[i])>2.4) continue;
				float deltaR1 = compute_R(jet_eta[i], dipho_leadEta, jet_phi[i], dipho_leadPhi);
				float deltaR2 = compute_R(jet_eta[i], dipho_subleadEta, jet_phi[i], dipho_subleadPhi);
				if(deltaR1<0.4 || deltaR2<0.4) continue;

				nJet++;
				if(jet_bdiscriminant[i]>bDiscriminantThresholdLoose && jet_bdiscriminant[i]<=bDiscriminantThresholdMedium)
				  looseBjet++;
				if(jet_bdiscriminant[i]>bDiscriminantThresholdMedium && jet_bdiscriminant[i]<=bDiscriminantThresholdTight)
				  mediumBjet++;
				if(jet_bdiscriminant[i]>bDiscriminantThresholdTight)
				  tightBjet++;
			}
			if(nJet<3) continue;

			bool passBjet_CutBased = bjetCut(looseBjet, mediumBjet, tightBjet, 0, 1, 0);
			bool passBjet_MVA = bjetCut(looseBjet, mediumBjet, tightBjet, 1, 0, 0);
			bool passBjet_New = bjetCut(looseBjet, mediumBjet, tightBjet, thresholds[1], thresholds[2], thresholds[3]);

			if(passBjet_CutBased && nJet>=5 && dipho_leadPt>dipho_mass/2.)
			{	passCutBasedSelection[n] += weight*lumiFactor;
				if(n<6)	
					Signal_histo_CutBased -> Fill(dipho_mass, weight*lumiFactor);
				if(n==7)
					Bkg_histo_CutBased -> Fill(dipho_mass, weight*lumiFactor);
				if(n==8)
					CS_histo_CutBased -> Fill(dipho_mass, weight*lumiFactor);
			}

			if(passBjet_MVA && nJet>=3 && ttHMVA>0.75)
			{       passMVASelection[n] += weight*lumiFactor;
				if(n<6)
					Signal_histo_MVA -> Fill(dipho_mass, weight*lumiFactor);
				if(n==7)
					Bkg_histo_MVA -> Fill(dipho_mass, weight*lumiFactor);
				if(n==8)
					CS_histo_MVA -> Fill(dipho_mass, weight*lumiFactor);
			}


			if(nJet>=thresholds[0] && passBjet_New && dipho_leadPt>dipho_mass*thresholds[4] && dipho_subleadPt>dipho_mass*thresholds[5] && jet_pt[0]>thresholds[6] && jet_pt[1]>thresholds[7] && jet_pt[2]>thresholds[8] && fabs(jet_eta[0])<thresholds[9] && fabs(jet_eta[1])<thresholds[10] && fabs(jet_eta[2])<thresholds[11] && fabs(jet_eta[0] - jet_eta[1])<thresholds[12] && fabs(jet_eta[0] - jet_eta[2])<thresholds[13] && fabs(jet_eta[1] - jet_eta[2])<thresholds[14])
			{
			        passNewSelection[n] += weight*lumiFactor;
				if(n<6)
					Signal_histo_New -> Fill(dipho_mass, weight*lumiFactor);
				if(n==7)
					Bkg_histo_New -> Fill(dipho_mass, weight*lumiFactor);
				if(n==8)
					CS_histo_New -> Fill(dipho_mass, weight*lumiFactor);
			}

		}

		cout << "Processing tag " << names[n] << ", processed " << nentries << " events out of " << nentries << endl;
	}
        
	

        std::vector<float> effectiveSigma_CutBased = FindSmallestInterval(Signal_histo_CutBased);
        Signal_histo_CutBased -> Rebin(8);      
        TF1* bkgFunc_CutBased = FitBkg(Bkg_histo_CutBased);

        std::vector<float> effectiveSigma_MVA = FindSmallestInterval(Signal_histo_MVA);
        Signal_histo_MVA -> Rebin(8);      
        TF1* bkgFunc_MVA = FitBkg(Bkg_histo_MVA);

        std::vector<float> effectiveSigma_New = FindSmallestInterval(Signal_histo_New);
        Signal_histo_New -> Rebin(8);      
        TF1* bkgFunc_New = FitBkg(Bkg_histo_New);


	//Read the progressive index for the name
        std::ifstream file("ProgressiveNaming.txt");
	std::string line;
	std::getline(file, line);
	int ProgressiveNumber = stoi(line);
	ProgressiveNumber++;
	file.close();

	system("rm ProgressiveNaming.txt");
	std::ofstream outfile("ProgressiveNaming.txt");
	outfile << ProgressiveNumber;
	outfile.close();

	TString name = ("NewSelections" + to_string(ProgressiveNumber)).c_str();

        packageResults(Signal_histo_CutBased, Bkg_histo_CutBased, Bkg_histo_CutBased, bkgFunc_CutBased, effectiveSigma_CutBased[2], effectiveSigma_CutBased[3], thresholds, passCutBasedSelection, lumiFactor, "2016CutBased");
        packageResults(Signal_histo_MVA, Bkg_histo_MVA, CS_histo_MVA, bkgFunc_MVA, effectiveSigma_MVA[2], effectiveSigma_MVA[3], thresholds, passMVASelection, lumiFactor, "2016MVA");
        packageResults(Signal_histo_New, Bkg_histo_New, CS_histo_New, bkgFunc_New, effectiveSigma_New[2], effectiveSigma_New[3], thresholds, passNewSelection, lumiFactor, name);




	//	MakePlot(NTopCandidate, "NTop_ttH");

        cout << endl << "|***********************************************************************************************************************|" << endl << endl;
	for(int i=0; i<7; i++)
        {        cout << "Cut based analysis on tag " <<  names[i] << ": number of events " << (float)passCutBasedSelection[i] << endl;
	         cout << "MVA analysis on tag " <<  names[i] << ": number of events " << (float)passMVASelection[i] << endl;
		 cout << "New cut based analysis on tag " <<  names[i] << ": number of events " << (float)passNewSelection[i] << endl;
	}
        cout << endl << "|***********************************************************************************************************************|" << endl;



	//	system("mv *.png ~/www/ttH/Full2016Dataset/Hadronic/ttH");
	//	system("mv *.pdf ~/www/ttH/Full2016Dataset/Hadronic/ttH");



	cout << endl << "Lufthansa, partner of Star Alliance, thanks you for choosing our company, we hope to see you again on board of our airctafts" << endl << endl;



























}









     






