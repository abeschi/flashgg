#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDMException.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "flashgg/DataFormats/interface/Jet.h"
#include "flashgg/DataFormats/interface/DiPhotonCandidate.h"
#include "flashgg/DataFormats/interface/TTHLeptonicTag.h"
#include "flashgg/DataFormats/interface/TTHLeptonicEfficiency.h"
#include "flashgg/DataFormats/interface/DiPhotonMVAResult.h"
#include "flashgg/DataFormats/interface/Electron.h"
#include "flashgg/DataFormats/interface/Muon.h"
#include "flashgg/DataFormats/interface/Met.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "flashgg/Taggers/interface/LeptonSelection.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/Common/interface/RefToPtr.h"

#include <vector>
#include <algorithm>
#include <string>
#include <utility>

#include "TMVA/Reader.h"

using namespace std;
using namespace edm;

namespace flashgg {

    class TTHFullyLeptonicEfficiencyProducer : public EDProducer
    {

    public:
        typedef math::XYZPoint Point;
        typedef math::Error<3>::type Error;

        TTHFullyLeptonicEfficiencyProducer( const ParameterSet & );
    private:
        void produce( Event &, const EventSetup & ) override;

        std::vector<edm::EDGetTokenT<View<flashgg::Jet> > > tokenJets_;
        EDGetTokenT<View<DiPhotonCandidate> > diPhotonToken_;
        std::vector<edm::InputTag> inputTagJets_;
        EDGetTokenT<View<Electron> > electronToken_;
        EDGetTokenT<View<flashgg::Muon> > muonToken_;
        EDGetTokenT<View<flashgg::Met> > METToken_;
        EDGetTokenT<edm::TriggerResults> HLTToken_;
        EDGetTokenT<View<reco::Vertex> > vertexToken_;
        EDGetTokenT<View<DiPhotonMVAResult> > mvaResultToken_;
        EDGetTokenT<View<reco::GenParticle> > genParticleToken_;
        EDGetTokenT<View<reco::GenJet> > genJetToken_;
        EDGetTokenT<int> stage0catToken_, stage1catToken_, njetsToken_;
        EDGetTokenT<float> pTHToken_,pTVToken_;
        EDGetTokenT<double> rhoTag_;
        string systLabel_;


        typedef std::vector<edm::Handle<edm::View<flashgg::Jet> > > JetCollectionVector;
        bool useTTHHadronicMVA_;
        double matchingGenPhotons_;
        //---thresholds---
        //---jets
        double jetPtThreshold_;
        double jetEtaThreshold_;
        double dRJetPhoLeadCut_;
        double dRJetPhoSubleadCut_;
        vector<double> bDiscriminator_;
         string bTag_;
        //leptons
        double leptonPtThreshold_;
        double muonEtaThreshold_;
        vector<double>  electronEtaThresholds_;
        double muPFIsoSumRelThreshold_;
        double muMiniIsoSumRelThreshold_;

        bool useStdLeptonID_;
        bool useElectronMVARecipe_;
        bool useElectronLooseID_;

    };

    TTHFullyLeptonicEfficiencyProducer::TTHFullyLeptonicEfficiencyProducer( const ParameterSet &iConfig ) :
        diPhotonToken_( consumes<View<flashgg::DiPhotonCandidate> >( iConfig.getParameter<InputTag> ( "DiPhotonTag" ) ) ),
        inputTagJets_( iConfig.getParameter<std::vector<edm::InputTag> >( "inputTagJets" ) ),
        electronToken_( consumes<View<flashgg::Electron> >( iConfig.getParameter<InputTag>( "ElectronTag" ) ) ),
        muonToken_( consumes<View<flashgg::Muon> >( iConfig.getParameter<InputTag>( "MuonTag" ) ) ),
        METToken_( consumes<View<flashgg::Met> >( iConfig.getParameter<InputTag>( "MetTag" ) ) ),
        HLTToken_( consumes<edm::TriggerResults>( iConfig.getParameter<InputTag>( "HLTTag" ) ) ),
        vertexToken_( consumes<View<reco::Vertex> >( iConfig.getParameter<InputTag> ( "VertexTag" ) ) ),
        mvaResultToken_( consumes<View<flashgg::DiPhotonMVAResult> >( iConfig.getParameter<InputTag>( "MVAResultTag" ) ) ),
        genParticleToken_( consumes<View<reco::GenParticle> >( iConfig.getParameter<InputTag> ( "GenParticleTag" ) ) ),
        genJetToken_( consumes<View<reco::GenJet> >( iConfig.getParameter<InputTag> ( "GenJetTag" ) ) ),
        rhoTag_( consumes<double>( iConfig.getParameter<InputTag>( "rhoTag" ) ) ),
        systLabel_( iConfig.getParameter<string> ( "SystLabel" ) )
    {
        ParameterSet HTXSps = iConfig.getParameterSet( "HTXSTags" );
        stage0catToken_ = consumes<int>( HTXSps.getParameter<InputTag>("stage0cat") );
        stage1catToken_ = consumes<int>( HTXSps.getParameter<InputTag>("stage1cat") );
        njetsToken_ = consumes<int>( HTXSps.getParameter<InputTag>("njets") );
        pTHToken_ = consumes<float>( HTXSps.getParameter<InputTag>("pTH") );
        pTVToken_ = consumes<float>( HTXSps.getParameter<InputTag>("pTV") );

        matchingGenPhotons_ = iConfig.getParameter<double>( "matchingGenPhotons" );

        leptonPtThreshold_ = iConfig.getParameter<double>( "leptonPtThreshold");
        muonEtaThreshold_ = iConfig.getParameter<double>( "muonEtaThreshold");
        jetPtThreshold_ = iConfig.getParameter<double>( "jetPtThreshold");
        jetEtaThreshold_ = iConfig.getParameter<double>( "jetEtaThreshold");
        dRJetPhoLeadCut_ = iConfig.getParameter<double>( "dRJetPhoLeadCut");
        dRJetPhoSubleadCut_ = iConfig.getParameter<double>( "dRJetPhoSubleadCut");
        bDiscriminator_ = iConfig.getParameter<vector<double > >( "bDiscriminator");
        bTag_ = iConfig.getParameter<string> ( "bTag");
 
        muPFIsoSumRelThreshold_ = iConfig.getParameter<double>( "muPFIsoSumRelThreshold");
        muMiniIsoSumRelThreshold_ = iConfig.getParameter<double>( "muMiniIsoSumRelThreshold");
        electronEtaThresholds_ = iConfig.getParameter<vector<double > >( "electronEtaThresholds");
        useStdLeptonID_=iConfig.getParameter<bool>("useStdLeptonID");
        useElectronMVARecipe_=iConfig.getParameter<bool>("useElectronMVARecipe");
        useElectronLooseID_=iConfig.getParameter<bool>("useElectronLooseID");
        

         for (unsigned i = 0 ; i < inputTagJets_.size() ; i++) {
            auto token = consumes<View<flashgg::Jet> >(inputTagJets_[i]);
            tokenJets_.push_back(token);
        }

        produces<vector<TTHLeptonicEfficiency> >();
    }

    void TTHFullyLeptonicEfficiencyProducer::produce( Event &evt, const EventSetup & )
    {

        Handle<int> stage0cat, stage1cat, njets;
        Handle<float> pTH, pTV;
        evt.getByToken(stage0catToken_, stage0cat);
        evt.getByToken(stage1catToken_,stage1cat);
        evt.getByToken(njetsToken_,njets);
        evt.getByToken(pTHToken_,pTH);
        evt.getByToken(pTVToken_,pTV);

         JetCollectionVector Jets( inputTagJets_.size() );
        for( size_t j = 0; j < inputTagJets_.size(); ++j ) {
            evt.getByToken( tokenJets_[j], Jets[j] );
        }

        Handle<View<flashgg::DiPhotonCandidate> > diPhotons;
        evt.getByToken( diPhotonToken_, diPhotons );

        Handle<View<flashgg::Muon> > theMuons;
        evt.getByToken( muonToken_, theMuons );

        Handle<View<flashgg::Electron> > theElectrons;
        evt.getByToken( electronToken_, theElectrons );

        edm::Handle<double>  rho;
        evt.getByToken(rhoTag_,rho);
        double rho_    = *rho;

        Handle<View<reco::Vertex> > vertices;
        evt.getByToken( vertexToken_, vertices );

        Handle<View<flashgg::DiPhotonMVAResult> > mvaResults;
        evt.getByToken( mvaResultToken_, mvaResults );

        Handle<View<flashgg::Met> > theMet_;
        evt.getByToken( METToken_, theMet_ );
 
       Handle<edm::TriggerResults> triggerResults_;
        evt.getByToken( HLTToken_, triggerResults_ );

        Handle<View<reco::GenParticle> > genParticles;
        Handle<View<reco::GenJet> > genJets;

        std::auto_ptr<vector<TTHLeptonicEfficiency> > ttheff  ( new vector<TTHLeptonicEfficiency> );

        Point higgsVtx;

        if( ! evt.isRealData() )
        {

            evt.getByToken( genParticleToken_, genParticles );
            evt.getByToken( genJetToken_, genJets );
 
            bool SavedHiggs = 0;
            bool  SavedTop1 = 0;
            bool  SavedTop2 = 0;

            for( unsigned int genLoop = 0 ; genLoop < genParticles->size(); genLoop++ )
            {
                int pdgid = genParticles->ptrAt( genLoop )->pdgId();
                if( pdgid == 25 || pdgid == 22 )
                {
                    
                    break;
                }
            }
    

            vector<reco::GenJet> genJetsVector;
            vector<reco::GenParticle> genPhotons;
            vector<reco::GenParticle> genLeptons;

            for(unsigned int i=0; i< genJets->size(); i++)
            {
                edm::Ptr<reco::GenJet> j = genJets->ptrAt(i);
                genJetsVector.push_back(*j);
            }

            for(unsigned int i=0; i<genParticles->size(); i++)
            {
                edm::Ptr<reco::GenParticle> p = genParticles->ptrAt(i);
                if(p->pdgId() == 6 && p->status()==62 && !SavedTop1)
                {   if(p->numberOfDaughters()!=2)
                    cout << "Warning: top quark is not decayed as expected" << endl;
                    else
                    {
                        for(unsigned int topDaugthers = 0; topDaugthers<p->numberOfDaughters(); topDaugthers++)
                        {
                            if(abs(p->daughter(topDaugthers)->pdgId())==24)
                            {
                                const reco::Candidate* Wboson = p->daughter(topDaugthers);
                                while(Wboson -> numberOfDaughters()!=2)
                                    Wboson = Wboson->daughter(0);
                                for(unsigned int WDaugthers = 0; WDaugthers<Wboson->numberOfDaughters(); WDaugthers++)
                                {
                                    if( abs((Wboson->daughter(WDaugthers))->pdgId())==11 || abs((Wboson->daughter(WDaugthers))->pdgId())==13 || abs((Wboson->daughter(WDaugthers))->pdgId())==15 )
                                       genLeptons.push_back(*((reco::GenParticle*)(Wboson->daughter(WDaugthers))));
                                }
                            }      
                        }
                    }
                    SavedTop1 = 1;
                }

                if(p->pdgId() == -6 && p->status()==62 && !SavedTop2)
                {   if(p->numberOfDaughters()!=2)
                    cout << "Warning: antitop quark is not decayed as expected" << endl;
                    else
                    {
                        for(unsigned int topDaugthers = 0; topDaugthers<p->numberOfDaughters(); topDaugthers++)
                        {
                            if(abs(p->daughter(topDaugthers)->pdgId())==24)
                            {
                                const reco::Candidate* Wboson = p->daughter(topDaugthers);
                                while(Wboson -> numberOfDaughters()!=2)
                                    Wboson = Wboson->daughter(0);
                                for(unsigned int WDaugthers = 0; WDaugthers<Wboson->numberOfDaughters(); WDaugthers++)
                                {
                                   if( abs((Wboson->daughter(WDaugthers))->pdgId())==11 || abs((Wboson->daughter(WDaugthers))->pdgId())==13 || abs((Wboson->daughter(WDaugthers))->pdgId())==15 )
                                       genLeptons.push_back(*((reco::GenParticle*)(Wboson->daughter(WDaugthers))));
                                }
                            }
                        }
                    }
                    SavedTop2 = 1;
                }

                if(p->pdgId() == 25 && !SavedHiggs)
                {   
                    const reco::Candidate* Hboson = (reco::Candidate*)(&(*p));
                    while(Hboson -> numberOfDaughters()!=2)
                        Hboson = Hboson->daughter(0);

                    higgsVtx = Hboson->vertex();

                     for(unsigned int HDaugthers = 0; HDaugthers<Hboson->numberOfDaughters(); HDaugthers++)
                         genPhotons.push_back(*((reco::GenParticle*)(Hboson->daughter(HDaugthers))));
                     SavedHiggs = 1;
                 }
            }
       

            reco::Vertex vtx0;
            
            
            if(vertices->size()!=0)
                vtx0 = *(vertices->ptrAt(0));
            else
            {   reco::Vertex::Error err;
                Point p(-137., -137., -137.);
                reco::Vertex vtx_tmp(p, err);
                vtx0 = vtx_tmp;
            }
            
            std::vector<edm::Ptr<flashgg::Muon> > goodMuons;
            if( !useStdLeptonID_)
            {
                goodMuons = selectAllMuonsSum16( theMuons->ptrs(), vertices->ptrs(), muonEtaThreshold_ , leptonPtThreshold_, muMiniIsoSumRelThreshold_ );
            } 
            else
            {
                goodMuons = selectAllMuons( theMuons->ptrs(), vertices->ptrs(), muonEtaThreshold_ , leptonPtThreshold_, muPFIsoSumRelThreshold_ );
            }
            
            std::vector<edm::Ptr<Electron> > goodElectrons ;
            goodElectrons = selectStdAllElectrons(theElectrons->ptrs(), vertices->ptrs(), leptonPtThreshold_, electronEtaThresholds_, useElectronMVARecipe_, useElectronLooseID_, rho_, evt.isRealData() );

            int nEle = (int)goodElectrons.size();
            int nMu = (int)goodMuons.size();
            int diphoSize = (int)diPhotons->size();

            std::vector<flashgg::TTHLeptonicTag> ttHtags;

            for( unsigned int diphoIndex = 0; diphoIndex < (double)diPhotons->size(); diphoIndex++ )
            {
                unsigned int jetCollectionIndex = diPhotons->ptrAt( diphoIndex )->jetCollectionIndex();
                edm::Ptr<flashgg::DiPhotonMVAResult> mvares = mvaResults->ptrAt( diphoIndex );
                bool passPreselection = 0;
                bool passHLT = 0;          
                
                int jetcount_ = 0;
                int njets_btagloose_ = 0;
                int njets_btagmedium_ = 0;
                int njets_btagtight_ = 0;

                std::vector<edm::Ptr<Muon> > tagMuons;
                std::vector<edm::Ptr<Electron> > tagElectrons;
                std::vector<edm::Ptr<Jet> > JetVect;
                std::vector<edm::Ptr<Jet> > BJetVect;

                edm::Ptr<flashgg::DiPhotonCandidate> dipho = diPhotons->ptrAt( diphoIndex );

                //require matching with genPhotons
                float dr11 = deltaR( dipho->leadingPhoton()->superCluster()->eta(), dipho->leadingPhoton()->superCluster()->phi(), genPhotons[0].eta(), genPhotons[0].phi() );
                float dr12 = deltaR( dipho->leadingPhoton()->superCluster()->eta(), dipho->leadingPhoton()->superCluster()->phi(), genPhotons[1].eta(), genPhotons[1].phi() );
                float dr21 = deltaR( dipho->subLeadingPhoton()->superCluster()->eta(), dipho->subLeadingPhoton()->superCluster()->phi(), genPhotons[0].eta(), genPhotons[0].phi() );
                float dr22 = deltaR( dipho->subLeadingPhoton()->superCluster()->eta(), dipho->subLeadingPhoton()->superCluster()->phi(), genPhotons[1].eta(), genPhotons[1].phi() );

                if(!((dr11<matchingGenPhotons_ && dr22<matchingGenPhotons_) || (dr12<matchingGenPhotons_ && dr21<matchingGenPhotons_))) continue;

                for( unsigned int jetIndex = 0; jetIndex < Jets[jetCollectionIndex]->size() ; jetIndex++ )
                {
                    edm::Ptr<flashgg::Jet> thejet = Jets[jetCollectionIndex]->ptrAt( jetIndex );
                    if( fabs( thejet->eta() ) > jetEtaThreshold_ ) { continue; }
                    if(!thejet->passesJetID  ( flashgg::Loose ) ) { continue; }
                    float dRPhoLeadJet = deltaR( thejet->eta(), thejet->phi(), dipho->leadingPhoton()->superCluster()->eta(), dipho->leadingPhoton()->superCluster()->phi() ) ;
                    float dRPhoSubLeadJet = deltaR( thejet->eta(), thejet->phi(), dipho->subLeadingPhoton()->superCluster()->eta(),
                    dipho->subLeadingPhoton()->superCluster()->phi() );

                    if( dRPhoLeadJet < dRJetPhoLeadCut_ || dRPhoSubLeadJet < dRJetPhoSubleadCut_ ) { continue; }
                    if( thejet->pt() < jetPtThreshold_ ) { continue; }

                    jetcount_++;
                    JetVect.push_back( thejet );
                                                
                    float bDiscriminatorValue = thejet->bDiscriminator( bTag_ );
                    
                    if( bDiscriminatorValue > bDiscriminator_[0] ) njets_btagloose_++;
                    if( bDiscriminatorValue > bDiscriminator_[1] )
                    {
                        njets_btagmedium_++;
                        BJetVect.push_back( thejet );
                    }
                    if( bDiscriminatorValue > bDiscriminator_[2] ) njets_btagtight_++;
                }
        
                TTHLeptonicTag tthhtags_obj( dipho, mvares);

                tthhtags_obj.setNjet( jetcount_ );
                tthhtags_obj.setNBLoose( njets_btagloose_ );
                tthhtags_obj.setNBMedium( njets_btagmedium_ );
                tthhtags_obj.setNBTight( njets_btagtight_ );
                tthhtags_obj.setSystLabel( systLabel_ );
                tthhtags_obj.setJets( JetVect );
                tthhtags_obj.setBJets( BJetVect );
                tthhtags_obj.setMuons( goodMuons );
                tthhtags_obj.setElectrons( goodElectrons );


                for( unsigned num = 0; num < JetVect.size(); num++ )
                    tthhtags_obj.includeWeightsByLabel( *JetVect[num] , "JetBTagCutWeight");

                if( nMu>0 && nEle>0)
                {
                    if( goodMuons.at(0)->pt() > goodElectrons.at(0)->pt() ) 
                    {
                        tthhtags_obj.includeWeightsByLabel( *goodMuons.at(0), "MuonMiniIsoWeight");
                    } 
                    else
                    {
                        tthhtags_obj.includeWeights( *goodElectrons.at(0) );
                    }
                } 
                else if( nMu>0 && nEle==0)
                {
                    tthhtags_obj.includeWeightsByLabel( *goodMuons.at(0), "MuonMiniIsoWeight" );
                }
                else if( nMu==0 && nEle>0)
                {
                    tthhtags_obj.includeWeights( *goodElectrons.at(0) );
                }            


                if( theMet_ -> size() != 1 )
                std::cout << "WARNING number of MET is not equal to 1" << std::endl;
                Ptr<flashgg::Met> Met = theMet_->ptrAt( 0 );
                tthhtags_obj.setMetPt(Met->pt());
                tthhtags_obj.setMetPhi(Met->phi());

                tthhtags_obj.includeWeights( *dipho );

                const edm::TriggerNames& trigNames = evt.triggerNames(*triggerResults_);   
                std::string pathName="HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90_v";
                for(unsigned int j=0; j<triggerResults_->size(); j++)
                {
                    for(int i=1; i<8; i++)
                    {
                        std::string pathName2 = pathName + to_string(i);
 
                        if(trigNames.triggerName(j)==pathName2)
                        {   passHLT =  triggerResults_->accept(j);
                           break;
                        }
                    }

                }

                if( (dipho->leadingPhoton()->full5x5_r9()>0.8 || dipho->leadingPhoton()->egChargedHadronIso()> 20 || dipho->leadingPhoton()->egChargedHadronIso()/dipho->leadingPhoton()->pt()<0.3) && (dipho->subLeadingPhoton()->full5x5_r9()>0.8 || dipho->subLeadingPhoton()->egChargedHadronIso()> 20 || dipho->subLeadingPhoton()->egChargedHadronIso()/dipho->subLeadingPhoton()->pt()<0.3) && (dipho->leadingPhoton()->pt()>30 && dipho->subLeadingPhoton()->pt()>20 && dipho->leadingPhoton()->passElectronVeto() && dipho->subLeadingPhoton()->passElectronVeto()) && (dipho->leadingPhoton()->phoIdMvaDWrtVtx(dipho->vtx())>-0.9 && dipho->subLeadingPhoton()->phoIdMvaDWrtVtx(dipho->vtx())>-0.9) && (dipho->leadingPhoton()->hadronicOverEm() < 0.08 && dipho->subLeadingPhoton()->hadronicOverEm()<0.08))
                    passPreselection = 1;

                tthhtags_obj.setPassPreselection(passPreselection);
                tthhtags_obj.setPassHLT(passHLT);

                ttHtags.push_back(tthhtags_obj);
            }

            if(genLeptons.size()==2)
            {   
                TTHLeptonicEfficiency tthhe(genJetsVector, genPhotons, genLeptons, ttHtags);
                tthhe.setHiggsVertex(higgsVtx);
                tthhe.setVertex0(vtx0);
                tthhe.setNDiphotons(diphoSize);
                ttheff->push_back( tthhe );
             }
        }

        evt.put( ttheff );


       // cout << "tagged events = " << count << endl;

    }
}
typedef flashgg::TTHFullyLeptonicEfficiencyProducer FlashggTTHFullyLeptonicEfficiencyProducer;
DEFINE_FWK_MODULE( FlashggTTHFullyLeptonicEfficiencyProducer );

// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4

