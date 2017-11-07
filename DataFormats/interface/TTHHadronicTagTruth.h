#ifndef FLASHgg_TTHHadronicTagTruth_h
#define FLASHgg_TTHHadronicTagTruth_h

#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Common/interface/Ptr.h"
#include "flashgg/DataFormats/interface/TagTruthBase.h"
#include "flashgg/DataFormats/interface/Jet.h"
#include "DataFormats/Math/interface/deltaR.h"

namespace flashgg {

    class TTHHadronicTagTruth : public TagTruthBase
    {

    public:

        TTHHadronicTagTruth();
        ~TTHHadronicTagTruth();
        TTHHadronicTagTruth *clone();


        void setGenJets(std::vector<reco::GenJet> jets) {GenJets_ = jets;}
        void setTop1(std::vector<reco::GenParticle> top) {GenTop1_ = top;}
        void setTop2(std::vector<reco::GenParticle> top) {GenTop2_ = top;}
        void setPhotons(std::vector<reco::GenParticle> ph) {Photons_ = ph;}

        const std::vector<reco::GenJet> GenJets() const {return GenJets_;}
        const std::vector<reco::GenParticle> GenTop1() const {return GenTop1_;}
        const std::vector<reco::GenParticle> GenTop2() const {return GenTop2_;}
        const std::vector<reco::GenParticle> GenPhotons() const {return Photons_;}

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

     private:
        std::vector<reco::GenJet> GenJets_;
        std::vector<reco::GenParticle> GenTop1_;
        std::vector<reco::GenParticle> GenTop2_;
        std::vector<reco::GenParticle> Photons_;
    };
}

#endif
// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
