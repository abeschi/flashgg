#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "flashgg/DataFormats/interface/MuMuGammaCandidate.h"
#include "flashgg/MicroAOD/interface/PhotonIdUtils.h"

//-----------J. Tao from IHEP-Beijing--------------

using namespace edm;
using namespace std;

namespace flashgg {

    class MuMuGammaRandomizedPhotonProducer : public EDProducer
    {

    public:
        MuMuGammaRandomizedPhotonProducer( const ParameterSet & );
    private:
        void produce( Event &, const EventSetup & ) override;
        EDGetTokenT<View<flashgg::DiMuonCandidate> > dimuToken_;
        EDGetTokenT<View<flashgg::Photon> > photonToken_;
        EDGetTokenT<View<reco::Vertex> > vertexToken_;

        double minPhotonPT_;
        double leadingMuonMinPt_;
        double subleadingMuonMinPt_;
        double MuonMinSumPt_;
        double leadingMuonIsoOverPt_;
        double subleadingMuonIsoOverPt_;
        double MinDimuonMass_;
        double MinMuMuGammaMass_;
        double MaxMuMuGammaMass_;
        double MaxMuMuPlusMuMuGammaMass_;
        double MaxDeltaR_;
        double MinClosestMuonPt_;
        //double maxPhotonEta_;

    };

    MuMuGammaRandomizedPhotonProducer::MuMuGammaRandomizedPhotonProducer( const ParameterSet &iConfig ) :
        dimuToken_( consumes<View<flashgg::DiMuonCandidate> >( iConfig.getParameter<InputTag> ( "DiMuonTag" ) ) ),
        photonToken_( consumes<View<flashgg::Photon> >( iConfig.getParameter<InputTag> ( "PhotonTag" ) ) ),
        vertexToken_( consumes<View<reco::Vertex> >( iConfig.getParameter<InputTag> ( "VertexTag" ) ) )
    {
        minPhotonPT_ = iConfig.getParameter<double>( "minPhotonPT" );
        leadingMuonMinPt_ = iConfig.getParameter<double>( "leadingMuonMinPt" );
        subleadingMuonMinPt_ = iConfig.getParameter<double>( "subleadingMuonMinPt" );
        MuonMinSumPt_ = iConfig.getParameter<double>( "MuonMinSumPt" );
        leadingMuonIsoOverPt_ = iConfig.getParameter<double>( "leadingMuonIsoOverPt" );
        subleadingMuonIsoOverPt_ = iConfig.getParameter<double>( "subleadingMuonIsoOverPt" );
        MinDimuonMass_ = iConfig.getParameter<double>( "MinDimuonMass" );
        MinMuMuGammaMass_ = iConfig.getParameter<double>( "MinMuMuGammaMass" );
        MinDimuonMass_ = iConfig.getParameter<double>( "MinDimuonMass" );
        MaxMuMuGammaMass_ = iConfig.getParameter<double>( "MaxMuMuGammaMass" );
        MaxDeltaR_ = iConfig.getParameter<double>( "MaxDeltaR" );
        MinClosestMuonPt_ = iConfig.getParameter<double>( "MinClosestMuonPt" );
        produces<vector<flashgg::MuMuGammaCandidate> >();
    }

    void MuMuGammaRandomizedPhotonProducer::produce( Event &evt, const EventSetup & )
    {

        Handle<View<reco::Vertex> > primaryVertices;
        evt.getByToken( vertexToken_, primaryVertices );
  
        const std::vector<edm::Ptr<reco::Vertex>> &pvPointers = primaryVertices->ptrs();
        edm::Ptr<reco::Vertex> pvx = pvPointers[0]; //selected vertex 0

        Handle<View<flashgg::DiMuonCandidate> > dimuons;
        evt.getByToken( dimuToken_, dimuons );

        Handle<View<flashgg::Photon> > photons;
        evt.getByToken( photonToken_, photons );

        auto_ptr<vector<flashgg::MuMuGammaCandidate> > MuMuGammaColl( new vector<flashgg::MuMuGammaCandidate> );
 
        for( unsigned int i = 0 ; i < dimuons->size() ; i++ )
        {
            edm::Ptr<flashgg::DiMuonCandidate> dimuon = dimuons->ptrAt(i);
 
            for( unsigned int j = 0; j < photons->size() ; j++ )
            {
                edm::Ptr<flashgg::Photon> photon = photons->ptrAt(j);

                // A number of things need to be done once the vertex is chosen
                // recomputing photon 4-momenta accordingly
                flashgg::Photon photon_corr = PhotonIdUtils::pho4MomCorrection( photon, pvx );
                photon_corr.addUserFloat("rnd_g_E", photon->userFloat("rnd_g_E"));

                float PhotonET =  photon_corr.pt();
                if( PhotonET < minPhotonPT_ )  continue; 

                float PhotonSCEta = photon_corr.superCluster()->position().Eta();
                if( fabs( PhotonSCEta ) > 2.5 || ( fabs( PhotonSCEta ) > 1.4442 && fabs( PhotonSCEta ) < 1.566 ) ) continue;


                const pat::Muon *muon_lead = dimuon->leadingMuon();
                const pat::Muon *muon_sublead = dimuon->subleadingMuon();

                reco::MuonPFIsolation  leadmuIso04 = muon_lead->pfIsolationR04();
                reco::MuonPFIsolation  subleadmuIso04 = muon_sublead->pfIsolationR04();
                double DeltaR1 = reco::deltaR( photon_corr.eta(), photon_corr.phi(), muon_lead->eta(), muon_lead->phi() );
                double DeltaR2 = reco::deltaR( photon_corr.eta(), photon_corr.phi(), muon_sublead->eta(), muon_sublead->phi() );

               //MuMuGammaCandidate mumugamma(dimu, photon_corr);
                flashgg::MuMuGammaCandidate* mumugamma = new flashgg::MuMuGammaCandidate(dimuon, photon_corr, pvx );
                 mumugamma->setVertex( pvx );
                //====================

                if( muon_lead->pt() < leadingMuonMinPt_ || muon_sublead->pt() < subleadingMuonMinPt_) continue;

                if((muon_lead->pt() + muon_sublead->pt()) < MuonMinSumPt_) continue;
                if(!dimuon->IsOSDiMuPair() || !dimuon->IfBothTightMu() || !dimuon->IfBothGlobalAndPF()) continue;
                if(leadmuIso04.sumChargedHadronPt/muon_lead->pt()>leadingMuonIsoOverPt_  || subleadmuIso04.sumChargedHadronPt/muon_sublead->pt()>subleadingMuonIsoOverPt_) continue;
                if(dimuon->mass() < MinDimuonMass_ || mumugamma->mass() < MinMuMuGammaMass_ || mumugamma->mass() > MaxMuMuGammaMass_ || ( dimuon->mass() + mumugamma->mass() ) > MaxMuMuGammaMass_ ) continue;
                if(min( DeltaR1, DeltaR2 ) > 0.8) continue;

                bool closestMuonPt;
                if(DeltaR1 > DeltaR2)
                    closestMuonPt = muon_lead->pt() > MinClosestMuonPt_;
                else
                    closestMuonPt = muon_sublead->pt() > MinClosestMuonPt_;
                
                if(!closestMuonPt) continue;

                //=====================
                double PhotonTrkIsoHollow03 = photon_corr.trkSumPtHollowConeDR03();
                mumugamma->setPhotonTrkIsoHollow03( PhotonTrkIsoHollow03 );
                //===Near mu within ISO cone: first near then far==
                double PhotonTrkIsoHollow03_MuCorr = PhotonTrkIsoHollow03;
                double LeadMuTrackPT = muon_lead->pt(), SubLeadMuTrackPT = muon_sublead->pt();
                if( ! muon_lead->track().isNull() ) { LeadMuTrackPT = muon_lead->track()->pt(); }
                if( ! muon_sublead->track().isNull() ) { SubLeadMuTrackPT = muon_sublead->track()->pt(); }
                double DeltaR_PhotonNearMu = DeltaR1 < DeltaR2 ? DeltaR1 : DeltaR2;
                double DeltaR_PhotonFarMu = DeltaR2 > DeltaR1 ? DeltaR2 : DeltaR1;
                double MuNear_TrkPT =  DeltaR1 < DeltaR2 ? LeadMuTrackPT : SubLeadMuTrackPT;
                double MuFar_TrkPT =  DeltaR1 < DeltaR2 ? SubLeadMuTrackPT : LeadMuTrackPT;
                if( DeltaR_PhotonNearMu < 0.3 && PhotonTrkIsoHollow03 > 0.99 * MuNear_TrkPT ) { PhotonTrkIsoHollow03_MuCorr -= MuNear_TrkPT; }
                if( DeltaR_PhotonFarMu < 0.3 && PhotonTrkIsoHollow03 > 0.99 * MuFar_TrkPT ) { PhotonTrkIsoHollow03_MuCorr -= MuFar_TrkPT; }
                mumugamma->setPhotonTrkIsoHollow03MuCorr( PhotonTrkIsoHollow03_MuCorr );


 
                // store the dimuon into the collection
                MuMuGammaColl->push_back( *mumugamma );
            }
        }

        evt.put( MuMuGammaColl );   
     }
}

typedef flashgg::MuMuGammaRandomizedPhotonProducer FlashggMuMuGammaRandomizedPhotonProducer;
DEFINE_FWK_MODULE( FlashggMuMuGammaRandomizedPhotonProducer );

// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
