#include "flashgg/DataFormats/interface/TTHHadronicEfficiency.h"
#include "flashgg/DataFormats/interface/Jet.h"

using namespace flashgg;

TTHHadronicEfficiency::TTHHadronicEfficiency() : DiPhotonTagBase::DiPhotonTagBase() {}

TTHHadronicEfficiency::~TTHHadronicEfficiency() {}


TTHHadronicEfficiency::TTHHadronicEfficiency(std::vector<reco::GenJet> GenJets, std::vector<reco::GenParticle> GenPhotons, std::vector<flashgg::TTHHadronicTag> ttH)
{
    ttH_ = ttH;
    Photons_ = GenPhotons;
    GenJets_ = GenJets;

}




const float TTHHadronicEfficiency::LeadingPhotonPt() const
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

const float TTHHadronicEfficiency::LeadingPhotonEta() const
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

const float TTHHadronicEfficiency::LeadingPhotonPhi() const
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

const float TTHHadronicEfficiency::LeadingPhotonEnergy() const
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



const float TTHHadronicEfficiency::SubleadingPhotonPt() const
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


const float TTHHadronicEfficiency::SubleadingPhotonEta() const
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


const float TTHHadronicEfficiency::SubleadingPhotonPhi() const
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



const float TTHHadronicEfficiency::SubleadingPhotonEnergy() const
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



const int TTHHadronicEfficiency::NGenJets() const
{
	int n = 0;
	for(unsigned int i=0; i<GenJets_.size(); i++)
	{
		if(GenJets_[i].pt()>20)
			n++;
	}
	
	return n;
}


const int TTHHadronicEfficiency::NGenJetsEta2p4() const
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

