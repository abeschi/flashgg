#include "flashgg/DataFormats/interface/TTHHadronicTagTruth.h"
#include <iostream>

using namespace flashgg;

TTHHadronicTagTruth::TTHHadronicTagTruth() {}

TTHHadronicTagTruth::~TTHHadronicTagTruth() {}

// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4


TTHHadronicTagTruth *TTHHadronicTagTruth::clone()
{
    TTHHadronicTagTruth *result = new TTHHadronicTagTruth;
    result->setGenJets( GenJets() );
    result->setTop1( GenTop1() );
    result->setTop2( GenTop2() );
    result->setPhotons( GenPhotons() );
    //result->copyBaseInfo( *this );
    return result;

}


const float TTHHadronicTagTruth::LeadingPhotonPt() const
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

const float TTHHadronicTagTruth::LeadingPhotonEta() const
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

const float TTHHadronicTagTruth::LeadingPhotonPhi() const
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

const float TTHHadronicTagTruth::LeadingPhotonEnergy() const
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



const float TTHHadronicTagTruth::SubleadingPhotonPt() const
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


const float TTHHadronicTagTruth::SubleadingPhotonEta() const
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


const float TTHHadronicTagTruth::SubleadingPhotonPhi() const
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



const float TTHHadronicTagTruth::SubleadingPhotonEnergy() const
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



const int TTHHadronicTagTruth::NGenJets() const
{
	int n = 0;
	for(unsigned int i=0; i<GenJets_.size(); i++)
	{
		if(GenJets_[i].pt()>20)
			n++;
	}
	
	return n;
}


const int TTHHadronicTagTruth::NGenJetsEta2p4() const
{
	int n = 0;
	for(unsigned int i=0; i<GenJets_.size(); i++)
	{
		if(GenJets_[i].pt()>20 && abs(GenJets_[i].eta())<2.4)
			n++;
	}
	
	return n;
}




