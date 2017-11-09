#ifndef FLASHgg_TTHHadronicEfficiency_h
#define FLASHgg_TTHHadronicEfficiency_h

#include "flashgg/DataFormats/interface/DiPhotonTagBase.h"
#include "flashgg/DataFormats/interface/TTHHadronicTag.h"
#include "flashgg/DataFormats/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/VertexReco/interface/Vertex.h"


namespace flashgg {

    class TTHHadronicEfficiency: public DiPhotonTagBase
    {
    public:

       typedef math::XYZPoint Point;

        TTHHadronicEfficiency();
        TTHHadronicEfficiency(std::vector<reco::GenJet>, std::vector<reco::GenParticle>, std::vector<flashgg::TTHHadronicTag> );
        TTHHadronicEfficiency *clone() const override { return ( new TTHHadronicEfficiency( *this ) ); }
        ~TTHHadronicEfficiency();

	const std::vector<flashgg::TTHHadronicTag> GetTTH() const {return ttH_;} 

        const std::vector<reco::GenJet> GenJets() const {return GenJets_;}

        const float LeadingPhotonPt() const;
        const float LeadingPhotonEta() const;
        const float LeadingPhotonPhi() const;
        const float LeadingPhotonEnergy() const;

        const float SubleadingPhotonPt() const;
        const float SubleadingPhotonEta() const;
        const float SubleadingPhotonPhi() const;
        const float SubleadingPhotonEnergy() const;

        const int NGenJets() const;
        const int NGenJetsEta2p4() const;

	void setHiggsVertex(Point higgsVtx) {higgsVtx_ = higgsVtx;}
	void setVertex0(reco::Vertex vertex0) {vertex0_ = vertex0;}

	const Point higgsVertex() const {return higgsVtx_;}
	const reco::Vertex vertex0() const {return vertex0_;}

        DiPhotonTagBase::tag_t tagEnum() const override {return DiPhotonTagBase::kTTHHadronic; }


    private:
        std::vector<flashgg::TTHHadronicTag> ttH_;
        std::vector<reco::GenJet> GenJets_;
        std::vector<reco::GenParticle> Photons_;
        Point higgsVtx_;
        reco::Vertex vertex0_;


    };
}

#endif





