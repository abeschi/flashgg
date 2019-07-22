#include "flashgg/DataFormats/interface/EleEleGammaCandidate.h"
#include "CommonTools/CandUtils/interface/AddFourMomenta.h"

using namespace flashgg;

EleEleGammaCandidate::EleEleGammaCandidate() {}

EleEleGammaCandidate::~EleEleGammaCandidate() {}

EleEleGammaCandidate::EleEleGammaCandidate( edm::Ptr<flashgg::DiElectronCandidate> dielectron, edm::Ptr<flashgg::Photon> photon, edm::Ptr<reco::Vertex> vertex )
{
    addDaughter( *dielectron );
    addDaughter( *photon );
    vertex_ = vertex;
    dieleptr_ = dielectron;

    AddFourMomenta addP4;
    addP4.set( *this );
}


EleEleGammaCandidate::EleEleGammaCandidate( edm::Ptr<flashgg::DiElectronCandidate> dielectron,  const flashgg::Photon &photon, edm::Ptr<reco::Vertex> vertex )
{
    addDaughter( *dielectron );
    addDaughter( photon );
    vertex_ = vertex;
    dieleptr_ = dielectron;

    AddFourMomenta addP4;
    addP4.set( *this );
}

EleEleGammaCandidate::EleEleGammaCandidate( flashgg::DiElectronCandidate dielectron, const flashgg::Photon &photon, edm::Ptr<reco::Vertex> vertex )
{
    addDaughter( dielectron );
    addDaughter( photon );
    vertex_ = vertex;

    AddFourMomenta addP4;
    addP4.set( *this );
}

const flashgg::DiElectronCandidate *EleEleGammaCandidate::EEG_DiEle() const
{
    return dynamic_cast<const flashgg::DiElectronCandidate *>( daughter( 0 ) );
}


const flashgg::Photon *EleEleGammaCandidate::EEG_Photon() const
{
    return dynamic_cast<const flashgg::Photon *>( daughter( 1 ) );
}

const flashgg::Electron* EleEleGammaCandidate::leadingElectron() const
{
    double pt1 = daughter(0)->daughter(0)->pt();
    double pt2 = daughter(0)->daughter(1)->pt();

    if(pt1>pt2)
        return dynamic_cast<const flashgg::Electron *> (daughter(0)->daughter(0));
    else
        return dynamic_cast<const flashgg::Electron *> (daughter(0)->daughter(1));
}

const flashgg::Electron* EleEleGammaCandidate::subleadingElectron() const
{
    double pt1 = daughter(0)->daughter(0)->pt();
    double pt2 = daughter(0)->daughter(1)->pt();

    if(pt1<pt2)
        return dynamic_cast<const flashgg::Electron *> (daughter(0)->daughter(0));
    else
        return dynamic_cast<const flashgg::Electron *> (daughter(0)->daughter(1));
}


double EleEleGammaCandidate::getMass() const
{
    TLorentzVector v1, v2, v3;
    v1.SetPtEtaPhiE(daughter(0)->daughter(0)->pt(), daughter(0)->daughter(0)->eta(), daughter(0)->daughter(0)->phi(), daughter(0)->daughter(0)->energy());
    v2.SetPtEtaPhiE(daughter(0)->daughter(1)->pt(), daughter(0)->daughter(1)->eta(), daughter(0)->daughter(1)->phi(), daughter(0)->daughter(1)->energy());
    v3.SetPtEtaPhiE(daughter(1)->pt(), daughter(1)->eta(), daughter(1)->phi(), daughter(1)->energy());
    
    return (v1+v2+v3).M();
}

double EleEleGammaCandidate::getPt() const
{
    TLorentzVector v1, v2, v3;
    v1.SetPtEtaPhiE(daughter(0)->daughter(0)->pt(), daughter(0)->daughter(0)->eta(), daughter(0)->daughter(0)->phi(), daughter(0)->daughter(0)->energy());
    v2.SetPtEtaPhiE(daughter(0)->daughter(1)->pt(), daughter(0)->daughter(1)->eta(), daughter(0)->daughter(1)->phi(), daughter(0)->daughter(1)->energy());
    v3.SetPtEtaPhiE(daughter(1)->pt(), daughter(1)->eta(), daughter(1)->phi(), daughter(1)->energy());
    
    return (v1+v2+v3).Pt();
}

double EleEleGammaCandidate::getDiEleMass() const
{
    TLorentzVector v1, v2;
    v1.SetPtEtaPhiE(daughter(0)->daughter(0)->pt(), daughter(0)->daughter(0)->eta(), daughter(0)->daughter(0)->phi(), daughter(0)->daughter(0)->energy());
    v2.SetPtEtaPhiE(daughter(0)->daughter(1)->pt(), daughter(0)->daughter(1)->eta(), daughter(0)->daughter(1)->phi(), daughter(0)->daughter(1)->energy());
    
    return (v1+v2).M();
}

double EleEleGammaCandidate::getDiElePt() const
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
