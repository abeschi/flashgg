#ifndef FLASHgg_TTHLeptonicEfficiency_h
#define FLASHgg_TTHLeptonicEfficiency_h

#include "flashgg/DataFormats/interface/DiPhotonTagBase.h"
#include "flashgg/DataFormats/interface/TTHLeptonicTag.h"
#include "flashgg/DataFormats/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/VertexReco/interface/Vertex.h"


namespace flashgg {

    class TTHLeptonicEfficiency: public DiPhotonTagBase
    {
    public:

       typedef math::XYZPoint Point;

        TTHLeptonicEfficiency();
        TTHLeptonicEfficiency(std::vector<reco::GenJet>, std::vector<reco::GenParticle>, std::vector<reco::GenParticle>, std::vector<flashgg::TTHLeptonicTag> );
        TTHLeptonicEfficiency *clone() const override { return ( new TTHLeptonicEfficiency( *this ) ); }
        ~TTHLeptonicEfficiency();

	const std::vector<flashgg::TTHLeptonicTag> GetTTH() const {return ttH_;} 

        const std::vector<reco::GenJet> GenJets() const {return GenJets_;}
        const std::vector<reco::GenParticle> GenLeptons() const {return GenLeptons_;}

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

	void setNDiphotons(int ndiphotons) {nDiphotons_ = ndiphotons;}
	const int NDiphotons() const {return nDiphotons_;}


	void setAllMuons(std::vector<flashgg::Muon> AllMuons) {AllMuons_ = AllMuons;}
	void setAllElectrons(std::vector<flashgg::Electron> AllElectrons) {AllElectrons_ = AllElectrons;}

	const std::vector<flashgg::Muon> Muons() const {return AllMuons_;}
	const std::vector<flashgg::Electron> Electrons() const {return AllElectrons_;}


        DiPhotonTagBase::tag_t tagEnum() const override {return DiPhotonTagBase::kTTHHadronic; }


    private:
        std::vector<flashgg::TTHLeptonicTag> ttH_;
        std::vector<reco::GenJet> GenJets_;
        std::vector<reco::GenParticle> Photons_;
        std::vector<reco::GenParticle> GenLeptons_;
	std::vector<flashgg::Muon> AllMuons_;
	std::vector<flashgg::Electron> AllElectrons_;
        Point higgsVtx_;
        reco::Vertex vertex0_;
	int nDiphotons_;


    };
}

#endif





