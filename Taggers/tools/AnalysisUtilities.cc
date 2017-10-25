#ifndef ANALYSYS_UTILITIES_H
#include "AnalysisUtilities.h"
#endif


float compute_R(float eta_1, float eta_2, float phi_1, float phi_2)
{
	float delta_eta_squared = (eta_1 - eta_2)*(eta_1 - eta_2);
	float delta_phi = fabs(phi_1 - phi_2);
		if(delta_phi > 3.14159265359)
			delta_phi =2*3.14159265359 - delta_phi;

	return sqrt(delta_eta_squared + delta_phi*delta_phi);
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


std::vector<float>  parse (TString fileName)
{
        std::vector<float> values;
        std::ifstream file(fileName);

	std::string line;
	while( std::getline( file, line ) )   
	{
	        std::istringstream iss( line );
		std::istringstream iss2( line );
		std::string value;
		if( iss2.get()=='#' ) continue;

		std::getline(iss, value);
		if(value== "END") break;
		values.push_back(stof(value));
	}
	file.close();


  return values;
}




bool bjetCut(int loose, int medium, int tight, int loose_cut, int medium_cut, int tight_cut) 
{
	if(tight<tight_cut)
		return 0;
	else
	{
		if((medium + tight - tight_cut) < medium_cut)
			return 0;
		else
		{
			if((loose + medium + tight - medium_cut - tight_cut) < loose_cut)
				return 0;
			else
				return 1;
		}
	}
}




std::vector<float> FindSmallestInterval(TH1F* histo, const float& fraction, const bool& verbosity)
{

  std::vector<float> output;
  float integralMax = fraction * histo->Integral();
  
  int N = histo -> GetNbinsX();
  int M1 = 0;
  int M2 = 0;
  for(int bin1 = 0; bin1 < N; ++bin1)
  {
    if( histo->GetBinContent(bin1+1) > 0. && M1 == 0 ) M1 = bin1-1;
    if( histo->GetBinContent(bin1+1) > 0. ) M2 = bin1+2;
  }
  
  std::map<int,float> binCenters;
  std::map<int,float> binContents;
  std::map<int,float> binIntegrals;
  for(int bin1 = M1; bin1 < M2; ++bin1)
  {
    binCenters[bin1] = histo->GetBinCenter(bin1+1);
    binContents[bin1] = histo->GetBinContent(bin1+1);
    
    for(int bin2 = M1; bin2 <= bin1; ++bin2)
      binIntegrals[bin1] += binContents[bin2];
  }
  
  float min = 0.;
  float max = 0.;
  float delta = 999999.;
  for(int bin1 = M1; bin1 < M2; ++bin1)
  {
    for(int bin2 = bin1+1; bin2 < M2; ++bin2)
    {
      if( (binIntegrals[bin2]-binIntegrals[bin1]) < integralMax ) continue;
      
      float tmpMin = histo -> GetBinCenter(bin1+1);
      float tmpMax = histo -> GetBinCenter(bin2+1);
      
      if( (tmpMax-tmpMin) < delta )
      {
        delta = (tmpMax - tmpMin);
        min = tmpMin;
        max = tmpMax;
      }
      
      break;
    }
  }
  
  TH1F* smallHisto = (TH1F*)( histo->Clone("smallHisto") );
  for(int bin = 1; bin <= smallHisto->GetNbinsX(); ++bin)
  {
    if( smallHisto->GetBinCenter(bin) < min )
      smallHisto -> SetBinContent(bin,0);
    
    if( smallHisto->GetBinCenter(bin) > max )
      smallHisto -> SetBinContent(bin,0);
  }
  smallHisto -> SetFillColor(kYellow);
  
  float mean = smallHisto -> GetMean();
  float meanErr = smallHisto -> GetMeanError();  
  
  output.push_back(mean);
  output.push_back(meanErr);
  output.push_back(min);
  output.push_back(max);

  return output;
}





TF1* FitBkg(TH1F* histo)
{

	TF1* fitFunc = new TF1("fitFunc", "expo", 110, 180);
	fitFunc -> SetParameter(0, 5);
	fitFunc -> SetParameter(1, -1);

	TFitResultPtr rp = histo -> Fit("fitFunc", "QENRS+");
	int fStatus = rp;
	int nTrials = 0;
	while( (fStatus != 0) && (nTrials < 10) )
	{
		rp = histo -> Fit("fitFunc", "QENRS+");
		fStatus = rp;
		if( fStatus == 0 ) break;
		++nTrials;
	}

	return fitFunc;
}









void packageResults(TH1F* Signal, TH1F* Bkg, TH1F* Cs, TF1* fitFunc, float min, float max, std::vector<float> thresholds, float* purity, float lumiFactor, TString title)
{
	TCanvas* c = new TCanvas();
	//int wtopx, wtopy;
	//unsigned int ww, wh;
	//c -> GetCanvasPar(wtopx, wtopy, ww, wh);
	unsigned int ww = c -> GetWw();	
	unsigned int wh = c -> GetWh();	
	c -> SetCanvasSize(2*ww, 1.2*wh);

	TPad *pad1 = new TPad("pad1", "pad1", 0., 0.3, 0.5, 0.97);
	pad1->SetBottomMargin(-0.1);
	pad1->Draw();
	pad1->cd();

	Bkg -> SetMarkerStyle(20);
	Bkg -> SetMarkerColor(kBlack);
	Bkg -> SetLineWidth(1);
	Bkg -> GetXaxis() -> SetRangeUser(100, 180);
	Bkg -> GetYaxis() -> SetRangeUser(0, 1.2*Bkg->GetMaximum());

	Signal -> SetLineColor(kAzure);
	Signal -> SetLineWidth(3);
	Signal -> SetFillStyle(0);
	//Signal -> Scale(100);

	Cs -> SetFillColor(kTeal-5);
	Cs -> SetFillStyle(3003);
	Cs -> SetLineColor(kTeal-6);
	Cs -> SetLineWidth(3);

	fitFunc -> SetLineWidth(2);
	fitFunc -> SetLineColor(kRed);

	Bkg -> Draw("EP");
	Signal -> Draw("histo SAME");
	Cs -> Draw("histo SAME");
	fitFunc -> Draw("SAME");	
	
	TLine* l1 = new TLine(min, 0, min, 0.5*Bkg->GetMaximum());
	l1 -> SetLineStyle(9);
	l1 -> SetLineWidth(3);
	l1 -> SetLineColor(kAzure+10);

	TLine* l2 = new TLine(max, 0, max, 0.5*Bkg->GetMaximum());
	l2 -> SetLineStyle(9);
	l2 -> SetLineWidth(3);
	l2 -> SetLineColor(kAzure+10);

	l1 -> Draw("SAME");
	l2 -> Draw("SAME");


	TLegend* leg = new TLegend(0.6, 0.7, 0.9, 0.9);
	leg -> AddEntry(Signal, "ttH", "l");
	leg -> AddEntry(Bkg, "Sidebands", "p");
	leg -> AddEntry(Bkg, "Control sample", "l");
	leg -> AddEntry(l1, "1 #sigma^{eff}", "l");
	leg -> Draw("SAME");

	stringstream ss;
	string texString = "#scale[1.5]{#font[22]{s/#sqrt{B}: ";
	ss << fixed << setprecision(2) << Signal->Integral(Signal->FindBin(min), Signal->FindBin(max))/sqrt(fitFunc ->Integral(min, max));
	texString += ss.str();
	texString += "}}";
	TString texString2 = texString;

	TLatex* tex = new TLatex(0.6, 0.5, texString2);
	tex -> SetNDC();
	tex -> Draw("SAME");


	CMS_lumi(pad1, 0, 0);

	c -> cd();

	TPad *pad2 = new TPad("pad2", "pad2", 0.5, 0.3, 1., 0.97);
	pad2->SetBottomMargin(-0.05);
	pad2->Draw();
	pad2->cd();

	TString infoTex[15] = {"Number of jets #geq ", "Number of loose b-jets #geq ", "Number of medium b-jets #geq ", "Number of tight b-jets #geq ", "Leading #gamma p_{T}/m_{#gamma#gamma} > ", "Subleading #gamma p_{T}/m_{#gamma#gamma} > ", "Leading jet p_{T} > ", "Subleading jet p_{T} > ", "Subsubleading jet p_{T} > ", "Leading jet |#eta| < ", "Subleading jet |#eta| < ", "Subsubleading jet |#eta| < ", "Leading-Subleading jet |#Delta#eta| < ", "Leading-Subsubleading jet |#Delta#eta| < ", "Subleading-Subsubleading jet |#Delta#eta| < "};


	for(int i=0; i<15; i++)
	{	
		TLatex* tex2a = new TLatex(0.1, 0.9-(0.06*i), infoTex[i]);
		tex2a -> SetNDC();
		tex2a -> Draw("SAME");

		stringstream ss;
		TString tmp2;

		if(i==0 || i==1 || i==2 || i==3 )
		{
			ss << fixed << setprecision(0) << thresholds[i];	
			string tmp = ss.str();
			tmp2 = tmp;
		}

		else
		{
			ss << fixed << setprecision(2) << thresholds[i];	
			string tmp = ss.str();
			tmp2 = tmp;
		}

		TLatex* tex2b = new TLatex(0.75, 0.9-(0.06*i), tmp2);
		tex2b -> SetNDC();
		tex2b -> Draw("SAME");
	}



	c -> cd();

	TPad *pad3 = new TPad("pad3", "pad3", 0.0, 0.05, 0.5, 0.25);
	pad3->SetBottomMargin(0.05);
	pad3->Draw();
	pad3->cd();

	float totalEvents = purity[0] + purity[1] + purity[2] + purity[3] + purity[4] + purity[5] + purity[6];
	float ttHpurity = purity[0]/totalEvents;

	TH1F* tth = new TH1F("tth", "", 1, 0, 1);
	TH1F* ggh = new TH1F("ggh", "", 1, 0, 1);
	TH1F* vbf = new TH1F("vbf", "", 1, 0, 1);
	TH1F* vh = new TH1F("vh", "", 1, 0, 1);
	TH1F* bbh = new TH1F("bbh", "", 1, 0, 1);
	TH1F* thq = new TH1F("thq", "", 1, 0, 1);
	TH1F* thw = new TH1F("thw", "", 1, 0, 1);

	tth -> SetBinContent(1, purity[0]);
	ggh -> SetBinContent(1, purity[1]);
	vbf -> SetBinContent(1, purity[2]);
	vh  -> SetBinContent(1, purity[3]);
	bbh -> SetBinContent(1, purity[4]);
	thq -> SetBinContent(1, purity[5]);
	thw -> SetBinContent(1, purity[6]);

	tth -> Scale(1./totalEvents);
	ggh -> Scale(1./totalEvents);
	vbf -> Scale(1./totalEvents);
	vh  -> Scale(1./totalEvents);
	bbh -> Scale(1./totalEvents);
	thq -> Scale(1./totalEvents);
	thw -> Scale(1./totalEvents);

	tth -> SetFillColor(kAzure);
	ggh -> SetFillColor(kGreen+2);
	vbf -> SetFillColor(kRed+1);
	vh  -> SetFillColor(kViolet-2);
	bbh -> SetFillColor(kOrange);
	thq -> SetFillColor(kAzure+8);
	thw -> SetFillColor(kViolet+2);

	ggh -> Add(tth);
	vbf  -> Add(ggh);
	vh  -> Add(vbf);
	bbh -> Add(vh);
	thq -> Add(bbh);
	thw -> Add(thq);

	TH1F* axis = new TH1F("axis", "", 1, 0, 1);

	axis -> GetYaxis() -> SetLabelSize(0);
	axis -> GetYaxis() -> SetLabelColor(kWhite);
	axis -> GetYaxis() -> SetTickLength(0);
	axis -> GetYaxis() -> SetRangeUser(0., 1.);
	axis -> GetYaxis() -> SetTickLength(0.03);
	axis -> GetYaxis() -> SetTickLength(0.03);
	axis -> GetXaxis() -> SetLabelSize(0.1);
	axis -> GetXaxis() -> SetRangeUser(0., 1.);
	axis -> GetXaxis() -> SetNdivisions(-10);

	
	axis -> Draw("histo");
	thw -> Draw("hbar same");
	thq -> Draw("hbar same");
	bbh -> Draw("hbar same");
	vh  -> Draw("hbar same");
	vbf -> Draw("hbar same");
	ggh -> Draw("hbar same");
	tth -> Draw("hbar same");

	c -> cd();

	TPad *pad4 = new TPad("pad4", "pad4", 0.5, 0.1, 0.62, 0.25);
	pad4->SetBottomMargin(0.05);
	pad4->Draw();
	pad4->cd();

	TLegend* leg4_1 = new TLegend(0., 0., 1., 1.);
	leg4_1 -> SetLineColor(kWhite);
	leg4_1 -> SetNColumns(2);
	leg4_1 -> AddEntry(tth, "ttH", "f"); 
	leg4_1 -> AddEntry(ggh, "ggH", "f"); 
	leg4_1 -> AddEntry(vbf, "VBF", "f"); 
	leg4_1 -> AddEntry(vh,  "VH",  "f"); 
	leg4_1 -> AddEntry(bbh, "bbH", "f"); 
	leg4_1 -> AddEntry(thq, "tHq", "f"); 
	leg4_1 -> AddEntry(thw, "tHW", "f"); 

	leg4_1 -> Draw("");


	c -> cd();

	TPad *pad5 = new TPad("pad5", "pad5", 0.62, 0.05, 0.95, 0.25);
	pad5->SetBottomMargin(0.05);
	pad5->Draw();
	pad5->cd();


	TColor *color=gROOT->GetColor(20);
	color->SetRGB(1.-ttHpurity, ttHpurity, 0.);

	TBox* box = new TBox(0., 0., 1., 1.);
	box -> SetLineColor(20);
	box -> SetLineStyle(1);
	box -> SetLineWidth(2);

	box -> Draw("l");

	TLatex* tex5_1 = new TLatex(0.37, 0.85, "#scale[4.]{#font[22]{Summary}}");
	tex5_1 -> SetTextColor(20);
	tex5_1 -> SetNDC();
	tex5_1 -> Draw("SAME");

	TLatex* tex5_2 = new TLatex(0.03, 0.67, "#scale[3.]{#font[22]{Expected signal events: }}");
	tex5_2 -> SetTextColor(20);
	tex5_2 -> SetNDC();
	tex5_2 -> Draw("SAME");

	TLatex* tex5_3 = new TLatex(0.03, 0.47, "#scale[3.]{#font[22]{Expected background events: }}");
	tex5_3 -> SetTextColor(20);
	tex5_3 -> SetNDC();
	tex5_3 -> Draw("SAME");

	TLatex* tex5_4 = new TLatex(0.03, 0.27, "#scale[3.]{#font[22]{Tag purity: }}");
	tex5_4 -> SetTextColor(20);
	tex5_4 -> SetNDC();
	tex5_4 -> Draw("SAME");

	TLatex* tex5_5 = new TLatex(0.03, 0.07, "#scale[3.]{#font[22]{Efficiency on ttH: }}");
	tex5_5-> SetTextColor(20);
	tex5_5 -> SetNDC();
	tex5_5 -> Draw("SAME");



	stringstream ss5_1;
	ss5_1 << "#scale[3.2]{#font[22]{" << fixed << setprecision(2) << totalEvents;
	string texString5_1 = ss5_1.str();
	texString5_1 += "}}";
	TString texString5_11 = texString5_1;

	TLatex* tex5_6 = new TLatex(0.8, 0.67, texString5_11);
	tex5_6-> SetTextColor(20);
	tex5_6 -> SetNDC();
	tex5_6 -> Draw("SAME");


	stringstream ss5_2;
	ss5_2 << "#scale[3.2]{#font[22]{" << fixed << setprecision(2) << fitFunc ->Integral(min, max);
	string texString5_2 = ss5_2.str();
	texString5_2 += "}}";
	TString texString5_22 = texString5_2;

	TLatex* tex5_7 = new TLatex(0.8, 0.47, texString5_22);
	tex5_7 -> SetTextColor(20);
	tex5_7 -> SetNDC();
	tex5_7 -> Draw("SAME");


	stringstream ss5_3;
	ss5_3 << "#scale[3.2]{#font[22]{" << fixed << setprecision(2) << ttHpurity;
	string texString5_3 = ss5_3.str();
	texString5_3 += "}}";
	TString texString5_33 = texString5_3;

	TLatex* tex5_8 = new TLatex(0.8, 0.27, texString5_33);
	tex5_8 -> SetTextColor(20);
	tex5_8 -> SetNDC();
	tex5_8 -> Draw("SAME");


	stringstream ss5_4;
	ss5_4 << "#scale[3.2]{#font[22]{" << fixed << setprecision(2) << purity[0]/(508.5*lumiFactor*0.002)*100.;
	string texString5_4 = ss5_4.str();
	texString5_4 +=  "% }}";
	TString texString5_44 = texString5_4;

	TLatex* tex5_9 = new TLatex(0.8, 0.07, texString5_44);
	tex5_9 -> SetTextColor(20);
	tex5_9 -> SetNDC();
	tex5_9 -> Draw("SAME");


	c -> cd();

	c -> SaveAs(title + ".pdf");
	c -> SaveAs(title + ".png");

	system("mv *.png *.pdf ~/www/ttH/Full2016Dataset/Hadronic/CutBasedOptimization");

	delete tth;
	delete ggh;
	delete vbf;
	delete vh;
	delete bbh;
	delete thq;
	delete thw;
	delete axis;
	delete c;

	return;
}
























