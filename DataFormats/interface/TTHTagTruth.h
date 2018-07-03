#ifndef FLASHgg_TTHTagTruth_h
#define FLASHgg_TTHTagTruth_h

#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Common/interface/Ptr.h"
#include "flashgg/DataFormats/interface/TagTruthBase.h"
#include "flashgg/DataFormats/interface/Jet.h"
#include "DataFormats/Math/interface/deltaR.h"

namespace flashgg {

    class TTHTagTruth : public TagTruthBase
    {

    public:

        TTHTagTruth();
        ~TTHTagTruth();
        //TTHTagTruth(const TTHTagTruth &b);
        
        
        edm::Ptr<reco::GenParticle> H() const { return H_; }
        edm::Ptr<reco::GenParticle> leadPhoton() const { return leadPhoton_; }
        edm::Ptr<reco::GenParticle> subleadPhoton() const { return subleadPhoton_; }
        edm::Ptr<reco::GenParticle> t() const { return t_; }
        edm::Ptr<reco::GenParticle> b() const { return b_; }
        edm::Ptr<reco::GenParticle> Wplus1() const { return Wplus1_; }
        edm::Ptr<reco::GenParticle> Wplus2() const { return Wplus2_; }
        edm::Ptr<reco::GenParticle> tbar() const { return tbar_; }
        edm::Ptr<reco::GenParticle> bbar() const { return bbar_; }
        edm::Ptr<reco::GenParticle> Wminus1() const { return Wminus1_; }
        edm::Ptr<reco::GenParticle> Wminus2() const { return Wminus2_; }
        //Setter methods
        void setH( edm::Ptr<reco::GenParticle> val ) { H_ = val; }
        void setLeadPhoton( edm::Ptr<reco::GenParticle> val ) { leadPhoton_ = val; }
        void setSubleadPhoton( edm::Ptr<reco::GenParticle> val ) { subleadPhoton_ = val; }
        void setT( edm::Ptr<reco::GenParticle> val ) { t_ = val; }
        void setB( edm::Ptr<reco::GenParticle> val ) { b_ = val; }
        void setWplus1( edm::Ptr<reco::GenParticle> val ) { Wplus1_ = val; }
        void setWplus2( edm::Ptr<reco::GenParticle> val ) { Wplus2_ = val; }
        void setTbar( edm::Ptr<reco::GenParticle> val ) { tbar_ = val; }
        void setBbar( edm::Ptr<reco::GenParticle> val ) { bbar_ = val; }
        void setWminus1( edm::Ptr<reco::GenParticle> val ) { Wminus1_ = val; }
        void setWminus2( edm::Ptr<reco::GenParticle> val ) { Wminus2_ = val; }
        
        //Counts
        //Clone
        TTHTagTruth *clone() const;

    private:
        edm::Ptr<reco::GenParticle> H_;
        edm::Ptr<reco::GenParticle> leadPhoton_;
        edm::Ptr<reco::GenParticle> subleadPhoton_;
        edm::Ptr<reco::GenParticle> t_;
        edm::Ptr<reco::GenParticle> b_;
        edm::Ptr<reco::GenParticle> Wplus1_;
        edm::Ptr<reco::GenParticle> Wplus2_;
        edm::Ptr<reco::GenParticle> tbar_;
        edm::Ptr<reco::GenParticle> bbar_;
        edm::Ptr<reco::GenParticle> Wminus1_;
        edm::Ptr<reco::GenParticle> Wminus2_;
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
