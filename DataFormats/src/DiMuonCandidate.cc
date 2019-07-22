#include "flashgg/DataFormats/interface/DiMuonCandidate.h"
#include "CommonTools/CandUtils/interface/AddFourMomenta.h"

//-----------J. Tao from IHEP-Beijing--------------

using namespace flashgg;

DiMuonCandidate::DiMuonCandidate() {}

DiMuonCandidate::~DiMuonCandidate() {}

DiMuonCandidate::DiMuonCandidate( edm::Ptr<flashgg::Muon> muon1, edm::Ptr<flashgg::Muon> muon2 )
{
    addDaughter( *muon1 );
    addDaughter( *muon2 );

    // Adding momenta
    // Needs its own object - but why?
    // Copied from example
    AddFourMomenta addP4;
    addP4.set( *this );
}

DiMuonCandidate::DiMuonCandidate( const flashgg::Muon &muon1, const flashgg::Muon &muon2 )
{
    addDaughter( muon1 );
    addDaughter( muon2 );

    // Adding momenta
    // Needs its own object - but why?
    // Copied from example
    AddFourMomenta addP4;
    addP4.set( *this );
}

const flashgg::Muon *DiMuonCandidate::leadingMuon() const
{
    if( daughter( 0 )->pt() > daughter( 1 )->pt() ) {
        return dynamic_cast<const flashgg::Muon *>( daughter( 0 ) );
    } else {
        return dynamic_cast<const flashgg::Muon *>( daughter( 1 ) );
    }
}

const flashgg::Muon *DiMuonCandidate::subleadingMuon() const
{
    if( daughter( 0 )->pt() > daughter( 1 )->pt() ) {
        return dynamic_cast<const flashgg::Muon *>( daughter( 1 ) );
    } else {
        return dynamic_cast<const flashgg::Muon *>( daughter( 0 ) );
    }
}

double DiMuonCandidate::getMass() const
{
    TLorentzVector v1, v2;
    v1.SetPtEtaPhiE(daughter(0)->pt(), daughter(0)->eta(), daughter(0)->phi(), daughter(0)->energy());
    v2.SetPtEtaPhiE(daughter(1)->pt(), daughter(1)->eta(), daughter(1)->phi(), daughter(1)->energy());

    return (v1+v2).M();
}

double DiMuonCandidate::getPt() const
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
