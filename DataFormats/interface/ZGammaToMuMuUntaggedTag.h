#ifndef flashgg_ZGammaToMuMuUntaggedTag
#define flashgg_ZGammaToMuMuUntaggedTag

#include "flashgg/DataFormats/interface/DiPhotonTagBase.h"
#include "flashgg/DataFormats/interface/Muon.h"
#include "flashgg/DataFormats/interface/Electron.h"
#include "flashgg/DataFormats/interface/Jet.h"
#include "flashgg/DataFormats/interface/Met.h"

#include "flashgg/DataFormats/interface/MuMuGammaCandidate.h"
#include "flashgg/DataFormats/interface/DiMuonCandidate.h"

namespace flashgg {

    class ZGammaToMuMuUntaggedTag: public DiPhotonTagBase
    {
    public:
        ZGammaToMuMuUntaggedTag();
        ~ZGammaToMuMuUntaggedTag();

        ZGammaToMuMuUntaggedTag( edm::Ptr<flashgg::MuMuGammaCandidate>);
        ZGammaToMuMuUntaggedTag( flashgg::MuMuGammaCandidate);

        const flashgg::MuMuGammaCandidate getMMG() const {return MMG_;}

        void setMuons( std::vector<edm::Ptr<Muon> > Muons ) {Muons_ = Muons;}
        const std::vector<edm::Ptr<Muon> > getMuons() const { return Muons_;}

        void setElectrons( std::vector<edm::Ptr<Electron> > Electrons ) {Electrons_ = Electrons;}
        const std::vector<edm::Ptr<flashgg::Electron> > getElectrons() const {return Electrons_;}

        void setJets( std::vector<edm::Ptr<Jet>> Jets ) { Jets_ = Jets; }
        const std::vector<edm::Ptr<Jet> > getJets() const { return Jets_;}

        void setMet( std::vector<edm::Ptr<Met>> Met ) { Met_ = Met; }
        const std::vector<edm::Ptr<Met> > getMet() const { return Met_;}

        void  setMva( float mva ) { mvaScore_ = mva; }
        float getMva() const { return mvaScore_;}

    private:

        flashgg::MuMuGammaCandidate               MMG_;
        std::vector<edm::Ptr<flashgg::Muon>>      Muons_;
        std::vector<edm::Ptr<flashgg::Electron>>  Electrons_;
        std::vector<edm::Ptr<flashgg::Jet>>       Jets_;
        std::vector<edm::Ptr<flashgg::Met>>       Met_;
        float mvaScore_;

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

