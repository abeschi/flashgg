#ifndef __IsolationCorrection__
#define __IsolationCorrection__

#include "TFile.h"
#include "TVectorD.h"
#include "TH1.h"
#include "TRandom.h"
    
#include <iostream>


class IsolationCorrection
{
	public:
		IsolationCorrection(const char * fname);
		~IsolationCorrection();
	
		float getExtra(float eta, float rho, float extraMult=0);
	
		int index(int ieta, int irho)  { return ieta*n_rho_centers_ + irho;}
		int findIndex(float val, std::vector<float> & vec);
		std::vector<int> convertToStdInt(TVectorD * input);
		std::vector<float> convertToStd(TVectorD * input);

	private:    
		std::vector<TH1 *> histograms_; 
		size_t n_rho_centers_;
		std::vector<float>  eta_centers_, rho_centers_eb_, rho_centers_ee_, extra_multiplicity_,  extra_multiplicity_slope_;
		std::vector<int> multiplicity_offset_;
    
};

#endif // __IsolationCorrection__
