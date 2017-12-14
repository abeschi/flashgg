#include "flashgg/DataFormats/interface/TTHLeptonicEfficiency.h"
#include "flashgg/DataFormats/interface/Jet.h"

using namespace flashgg;

TTHLeptonicEfficiency::TTHLeptonicEfficiency() : DiPhotonTagBase::DiPhotonTagBase() {}

TTHLeptonicEfficiency::~TTHLeptonicEfficiency() {}


TTHLeptonicEfficiency::TTHLeptonicEfficiency(std::vector<reco::GenJet> GenJets, std::vector<reco::GenParticle> GenPhotons, std::vector<reco::GenParticle> GenLeptons, std::vector<flashgg::TTHLeptonicTag> ttH)
{
    ttH_ = ttH;
    Photons_ = GenPhotons;
    GenLeptons_ = GenLeptons;
    GenJets_ = GenJets;
}




const float TTHLeptonicEfficiency::LeadingPhotonPt() const
{
	if(Photons_.size()!=2)
		return -100.;
	else
	{	if(Photons_[0].pt()>Photons_[1].pt())
			return Photons_[0].pt();
		else
			return Photons_[1].pt();
	}
}

const float TTHLeptonicEfficiency::LeadingPhotonEta() const
{
	if(Photons_.size()!=2)
		return -100.;
	else
	{	if(Photons_[0].pt()>Photons_[1].pt())
			return Photons_[0].eta();
		else
			return Photons_[1].eta();
	}
}

const float TTHLeptonicEfficiency::LeadingPhotonPhi() const
{
	if(Photons_.size()!=2)
		return -100.;
	else
	{	if(Photons_[0].pt()>Photons_[1].pt())
			return Photons_[0].phi();
		else
			return Photons_[1].phi();
	}
}

const float TTHLeptonicEfficiency::LeadingPhotonEnergy() const
{
	if(Photons_.size()!=2)
		return -100.;
	else
	{	if(Photons_[0].pt()>Photons_[1].pt())
			return Photons_[0].energy();
		else
			return Photons_[1].energy();
	}
}



const float TTHLeptonicEfficiency::SubleadingPhotonPt() const
{
	if(Photons_.size()!=2)
		return -100.;
	else
	{	if(Photons_[0].pt()>Photons_[1].pt())
			return Photons_[1].pt();
		else
			return Photons_[0].pt();
	}
}


const float TTHLeptonicEfficiency::SubleadingPhotonEta() const
{
	if(Photons_.size()!=2)
		return -100.;
	else
	{	if(Photons_[0].pt()>Photons_[1].pt())
			return Photons_[1].eta();
		else
			return Photons_[0].eta();
	}
}


const float TTHLeptonicEfficiency::SubleadingPhotonPhi() const
{
	if(Photons_.size()!=2)
		return -100.;
	else
	{	if(Photons_[0].pt()>Photons_[1].pt())
			return Photons_[1].phi();
		else
			return Photons_[0].phi();
	}
}



const float TTHLeptonicEfficiency::SubleadingPhotonEnergy() const
{
	if(Photons_.size()!=2)
		return -100.;
	else
	{	if(Photons_[0].pt()>Photons_[1].pt())
			return Photons_[1].energy();
		else
			return Photons_[0].energy();
	}
}



const int TTHLeptonicEfficiency::NGenJets() const
{
	int n = 0;
	for(unsigned int i=0; i<GenJets_.size(); i++)
	{
		if(GenJets_[i].pt()>20)
			n++;
	}
	
	return n;
}


const int TTHLeptonicEfficiency::NGenJetsEta2p4() const
{
	int n = 0;
	for(unsigned int i=0; i<GenJets_.size(); i++)
	{
		if(GenJets_[i].pt()>20 && abs(GenJets_[i].eta())<2.4)
			n++;
	}
	
	return n;
}


// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4

