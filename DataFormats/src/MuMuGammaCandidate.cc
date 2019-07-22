#include "flashgg/DataFormats/interface/MuMuGammaCandidate.h"
#include "CommonTools/CandUtils/interface/AddFourMomenta.h"

//-----------J. Tao from IHEP-Beijing--------------

using namespace flashgg;

MuMuGammaCandidate::MuMuGammaCandidate() {}

MuMuGammaCandidate::~MuMuGammaCandidate() {}

MuMuGammaCandidate::MuMuGammaCandidate( edm::Ptr<flashgg::DiMuonCandidate> dimuon, edm::Ptr<flashgg::Photon> photon, edm::Ptr<reco::Vertex> vertex )
{
    addDaughter( *dimuon );
    addDaughter( *photon );
    vertex_ = vertex;
    dimuptr_ = dimuon;
    // Adding momenta
    // Needs its own object - but why?
    // Copied from example
    AddFourMomenta addP4;
    addP4.set( *this );
}

MuMuGammaCandidate::MuMuGammaCandidate( edm::Ptr<flashgg::DiMuonCandidate> dimuon,  const flashgg::Photon &photon, edm::Ptr<reco::Vertex> vertex )
{
    addDaughter( *dimuon );
    addDaughter( photon );
    vertex_ = vertex;
    dimuptr_ = dimuon;

    AddFourMomenta addP4;
    addP4.set( *this );
}

MuMuGammaCandidate::MuMuGammaCandidate( flashgg::DiMuonCandidate dimuon,  const flashgg::Photon &photon, edm::Ptr<reco::Vertex> vertex )
{
    addDaughter( dimuon );
    addDaughter( photon );
    vertex_ = vertex;

    AddFourMomenta addP4;
    addP4.set( *this );
}

const flashgg::DiMuonCandidate *MuMuGammaCandidate::MMG_DiMu() const
{
    return dynamic_cast<const flashgg::DiMuonCandidate *>( daughter( 0 ) );
}


const flashgg::Photon *MuMuGammaCandidate::MMG_Photon() const
{
    return dynamic_cast<const flashgg::Photon *>( daughter( 1 ) );
}


const flashgg::Muon* MuMuGammaCandidate::leadingMuon() const
{
    double pt1 = daughter(0)->daughter(0)->pt();
    double pt2 = daughter(0)->daughter(1)->pt();

    if(pt1>pt2)
        return dynamic_cast<const flashgg::Muon *> (daughter(0)->daughter(0));
    else
        return dynamic_cast<const flashgg::Muon *> (daughter(0)->daughter(1));
}

const flashgg::Muon* MuMuGammaCandidate::subleadingMuon() const
{
    double pt1 = daughter(0)->daughter(0)->pt();
    double pt2 = daughter(0)->daughter(1)->pt();

    if(pt1<pt2)
        return dynamic_cast<const flashgg::Muon *> (daughter(0)->daughter(0));
    else
        return dynamic_cast<const flashgg::Muon *> (daughter(0)->daughter(1));
}


double MuMuGammaCandidate::getMass() const
{
    TLorentzVector v1, v2, v3;
    v1.SetPtEtaPhiE(daughter(0)->daughter(0)->pt(), daughter(0)->daughter(0)->eta(), daughter(0)->daughter(0)->phi(), daughter(0)->daughter(0)->energy());
    v2.SetPtEtaPhiE(daughter(0)->daughter(1)->pt(), daughter(0)->daughter(1)->eta(), daughter(0)->daughter(1)->phi(), daughter(0)->daughter(1)->energy());
    v3.SetPtEtaPhiE(daughter(1)->pt(), daughter(1)->eta(), daughter(1)->phi(), daughter(1)->energy());
    
    return (v1+v2+v3).M();
}

double MuMuGammaCandidate::getPt() const
{
    TLorentzVector v1, v2, v3;
    v1.SetPtEtaPhiE(daughter(0)->daughter(0)->pt(), daughter(0)->daughter(0)->eta(), daughter(0)->daughter(0)->phi(), daughter(0)->daughter(0)->energy());
    v2.SetPtEtaPhiE(daughter(0)->daughter(1)->pt(), daughter(0)->daughter(1)->eta(), daughter(0)->daughter(1)->phi(), daughter(0)->daughter(1)->energy());
    v3.SetPtEtaPhiE(daughter(1)->pt(), daughter(1)->eta(), daughter(1)->phi(), daughter(1)->energy());
    
    return (v1+v2+v3).Pt();
}

double MuMuGammaCandidate::getDiMuonMass() const
{
    TLorentzVector v1, v2;
    v1.SetPtEtaPhiE(daughter(0)->daughter(0)->pt(), daughter(0)->daughter(0)->eta(), daughter(0)->daughter(0)->phi(), daughter(0)->daughter(0)->energy());
    v2.SetPtEtaPhiE(daughter(0)->daughter(1)->pt(), daughter(0)->daughter(1)->eta(), daughter(0)->daughter(1)->phi(), daughter(0)->daughter(1)->energy());
    
    return (v1+v2).M();
}

double MuMuGammaCandidate::getDiMuonPt() const
{
    TLorentzVector v1, v2;
    v1.SetPtEtaPhiE(daughter(0)->daughter(0)->pt(), daughter(0)->daughter(0)->eta(), daughter(0)->daughter(0)->phi(), daughter(0)->daughter(0)->energy());
    v2.SetPtEtaPhiE(daughter(0)->daughter(1)->pt(), daughter(0)->daughter(1)->eta(), daughter(0)->daughter(1)->phi(), daughter(0)->daughter(1)->energy());
    
    return (v1+v2).Pt();
}

// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
