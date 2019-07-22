#include "flashgg/DataFormats/interface/DiElectronCandidate.h"
#include "CommonTools/CandUtils/interface/AddFourMomenta.h"

using namespace flashgg;

DiElectronCandidate::DiElectronCandidate() {}

DiElectronCandidate::~DiElectronCandidate() {}

DiElectronCandidate::DiElectronCandidate( edm::Ptr<flashgg::Electron> electron1, edm::Ptr<flashgg::Electron> electron2 )
{
    addDaughter( *electron1 );
    addDaughter( *electron2 );

    // Adding momenta
    // Needs its own object - but why?
    // Copied from example
    AddFourMomenta addP4;
    addP4.set( *this );
}

DiElectronCandidate::DiElectronCandidate( const flashgg::Electron &electron1, const flashgg::Electron &electron2 )
{
    addDaughter( electron1 );
    addDaughter( electron2 );

    // Adding momenta
    // Needs its own object - but why?
    // Copied from example
    AddFourMomenta addP4;
    addP4.set( *this );
}

const flashgg::Electron *DiElectronCandidate::leadingElectron() const
{
    if( daughter( 0 )->pt() > daughter( 1 )->pt() ) {
        return dynamic_cast<const flashgg::Electron *>( daughter( 0 ) );
    } else {
        return dynamic_cast<const flashgg::Electron *>( daughter( 1 ) );
    }
}

const flashgg::Electron *DiElectronCandidate::subleadingElectron() const
{
    if( daughter( 0 )->pt() > daughter( 1 )->pt() ) {
        return dynamic_cast<const flashgg::Electron *>( daughter( 1 ) );
    } else {
        return dynamic_cast<const flashgg::Electron *>( daughter( 0 ) );
    }
}

double DiElectronCandidate::getMass() const
{
    TLorentzVector v1, v2;
    v1.SetPtEtaPhiE(daughter(0)->pt(), daughter(0)->eta(), daughter(0)->phi(), daughter(0)->energy());
    v2.SetPtEtaPhiE(daughter(1)->pt(), daughter(1)->eta(), daughter(1)->phi(), daughter(1)->energy());

    return (v1+v2).M();
}

double DiElectronCandidate::getPt() const
{
    TLorentzVector v1, v2;
    v1.SetPtEtaPhiE(daughter(0)->pt(), daughter(0)->eta(), daughter(0)->phi(), daughter(0)->energy());
    v2.SetPtEtaPhiE(daughter(1)->pt(), daughter(1)->eta(), daughter(1)->phi(), daughter(1)->energy());

    return (v1+v2).Pt();
}

// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
