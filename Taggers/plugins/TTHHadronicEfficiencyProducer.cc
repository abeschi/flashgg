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
#include "flashgg/DataFormats/interface/TTHHadronicTag.h"
#include "flashgg/DataFormats/interface/TTHHadronicEfficiency.h"
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

    class TTHHadronicEfficiencyProducer : public EDProducer
    {

    public:
        typedef math::XYZPoint Point;
        typedef math::Error<3>::type Error;

        TTHHadronicEfficiencyProducer( const ParameterSet & );
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
        //---photons
        double MVAThreshold_;
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

        unique_ptr<TMVA::Reader>TThMva_;
        FileInPath tthMVAweightfile_;
        string _MVAMethod;

        int jetcount_;
        float nJets_;
        int njets_btagloose_;
        int njets_btagmedium_;
        int njets_btagtight_;
        double idmva1_;
        double idmva2_;
        float leadJetPt_;
        float subLeadJetPt_;
        float sumJetPt_;
        float maxBTagVal_;
        float secondMaxBTagVal_;
        float tthMvaVal_;

    };

    TTHHadronicEfficiencyProducer::TTHHadronicEfficiencyProducer( const ParameterSet &iConfig ) :
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
        systLabel_( iConfig.getParameter<string> ( "SystLabel" ) ),
         _MVAMethod( iConfig.getParameter<string> ( "MVAMethod" ) )
    {
        ParameterSet HTXSps = iConfig.getParameterSet( "HTXSTags" );
        stage0catToken_ = consumes<int>( HTXSps.getParameter<InputTag>("stage0cat") );
        stage1catToken_ = consumes<int>( HTXSps.getParameter<InputTag>("stage1cat") );
        njetsToken_ = consumes<int>( HTXSps.getParameter<InputTag>("njets") );
        pTHToken_ = consumes<float>( HTXSps.getParameter<InputTag>("pTH") );
        pTVToken_ = consumes<float>( HTXSps.getParameter<InputTag>("pTV") );

        MVAThreshold_ = iConfig.getParameter<double>( "MVAThreshold");
        matchingGenPhotons_ = iConfig.getParameter<double>( "matchingGenPhotons" );

        leptonPtThreshold_ = iConfig.getParameter<double>( "leptonPtThreshold");
        muonEtaThreshold_ = iConfig.getParameter<double>( "muonEtaThreshold");
        jetPtThreshold_ = iConfig.getParameter<double>( "jetPtThreshold");
        jetEtaThreshold_ = iConfig.getParameter<double>( "jetEtaThreshold");
        dRJetPhoLeadCut_ = iConfig.getParameter<double>( "dRJetPhoLeadCut");
        dRJetPhoSubleadCut_ = iConfig.getParameter<double>( "dRJetPhoSubleadCut");
        bDiscriminator_ = iConfig.getParameter<vector<double > >( "bDiscriminator");
        bTag_ = iConfig.getParameter<string> ( "bTag");

        useTTHHadronicMVA_ = iConfig.getParameter<bool>( "useTTHHadronicMVA");
 
        muPFIsoSumRelThreshold_ = iConfig.getParameter<double>( "muPFIsoSumRelThreshold");
        muMiniIsoSumRelThreshold_ = iConfig.getParameter<double>( "muMiniIsoSumRelThreshold");
        electronEtaThresholds_ = iConfig.getParameter<vector<double > >( "electronEtaThresholds");
        useStdLeptonID_=iConfig.getParameter<bool>("useStdLeptonID");
        useElectronMVARecipe_=iConfig.getParameter<bool>("useElectronMVARecipe");
        useElectronLooseID_=iConfig.getParameter<bool>("useElectronLooseID");
        

        tthMVAweightfile_ = iConfig.getParameter<edm::FileInPath>( "tthMVAweightfile" ); 

        nJets_ = 0;
        leadJetPt_ = 0.;
        subLeadJetPt_ = 0.;
        sumJetPt_ = 0.;
        maxBTagVal_ = -999.;
        secondMaxBTagVal_ = -999.;

        if (_MVAMethod != ""){
            TThMva_.reset( new TMVA::Reader( "!Color:Silent" ) );
            TThMva_->AddVariable( "nJets", &nJets_);
            TThMva_->AddVariable( "maxBTagVal",&maxBTagVal_);
            TThMva_->AddVariable( "secondMaxBTagVal", &secondMaxBTagVal_);
            TThMva_->AddVariable( "leadJetPt", &leadJetPt_);
        
            TThMva_->BookMVA( _MVAMethod.c_str() , tthMVAweightfile_.fullPath() );
        
        }       

        for (unsigned i = 0 ; i < inputTagJets_.size() ; i++) {
            auto token = consumes<View<flashgg::Jet> >(inputTagJets_[i]);
            tokenJets_.push_back(token);
        }

        produces<vector<TTHHadronicEfficiency> >();
    }

    void TTHHadronicEfficiencyProducer::produce( Event &evt, const EventSetup & )
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

        std::auto_ptr<vector<TTHHadronicEfficiency> > ttheff  ( new vector<TTHHadronicEfficiency> );

        Point higgsVtx;
        bool isLeptonic = 0;

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
                                    if( abs((Wboson->daughter(WDaugthers))->pdgId())>6) isLeptonic = 1;
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
                                    if( abs((Wboson->daughter(WDaugthers))->pdgId())>6) isLeptonic = 1;
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
            
            bool hasRecoLeptons = 0;

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

            if( goodElectrons.size() > 0 ||  goodMuons.size() > 0 )  hasRecoLeptons = 1;

            jetcount_ = 0;
            nJets_ = 0;
            njets_btagloose_ = 0;
            njets_btagmedium_ = 0;
            njets_btagtight_ = 0;
            idmva1_ = -999.;
            idmva2_ = -999.;
            leadJetPt_ = 0.;
            subLeadJetPt_ = 0.;
            sumJetPt_ = 0.;
            maxBTagVal_ = -999.;
            secondMaxBTagVal_ = -999.;
            tthMvaVal_ = -999.;        

            std::vector<flashgg::TTHHadronicTag> ttHtags;

            for( unsigned int diphoIndex = 0; diphoIndex < (double)diPhotons->size(); diphoIndex++ )
            {
                unsigned int jetCollectionIndex = diPhotons->ptrAt( diphoIndex )->jetCollectionIndex();
                edm::Ptr<flashgg::DiPhotonMVAResult> mvares = mvaResults->ptrAt( diphoIndex );
                bool passPreselection = 0;
                bool passHLT = 0;          
                

                std::vector<edm::Ptr<flashgg::Jet> > JetVect;
                JetVect.clear();
                std::vector<edm::Ptr<flashgg::Jet> > BJetVect;
                BJetVect.clear();
                std::vector<edm::Ptr<flashgg::Jet> > BJetTTHHMVAVect;
                BJetTTHHMVAVect.clear();
                std::vector<float> JetBTagVal;
                JetBTagVal.clear();

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
                    nJets_ = jetcount_;
                    JetVect.push_back( thejet );
                        
                    float jetPt = thejet->pt();
                    if(jetPt > leadJetPt_)
                    {
                        if(leadJetPt_ > subLeadJetPt_) { subLeadJetPt_ = leadJetPt_; }
                        leadJetPt_ = jetPt;
                    }
                    else if(jetPt > subLeadJetPt_)
                    {
                            subLeadJetPt_ = jetPt;
                    }
                    sumJetPt_ += jetPt;
                        
                    float bDiscriminatorValue = -2.;
                    bDiscriminatorValue = thejet->bDiscriminator( bTag_ );
                    
                    if(bDiscriminatorValue > maxBTagVal_)
                    {
                        BJetTTHHMVAVect.insert( BJetTTHHMVAVect.begin(), thejet );
                        if(BJetTTHHMVAVect.size() >= 3){ BJetTTHHMVAVect.pop_back(); }
                        if(maxBTagVal_ > secondMaxBTagVal_) { secondMaxBTagVal_ = maxBTagVal_; }
                        maxBTagVal_ = bDiscriminatorValue;

                    } 
                    else if(bDiscriminatorValue > secondMaxBTagVal_)
                    {
                        secondMaxBTagVal_ = bDiscriminatorValue;
                        if(BJetTTHHMVAVect.size() >= 2){BJetTTHHMVAVect.pop_back();} 
                        BJetTTHHMVAVect.push_back( thejet );
                    }
                        
                    JetBTagVal.push_back( bDiscriminatorValue );
                    if( bDiscriminatorValue > bDiscriminator_[0] ) njets_btagloose_++;
                    if( bDiscriminatorValue > bDiscriminator_[1] )
                    {
                        njets_btagmedium_++;
                        //JetVect.pop_back();
                        BJetVect.push_back( thejet );
                    }
                    if( bDiscriminatorValue > bDiscriminator_[2] ) njets_btagtight_++;


                }
                tthMvaVal_ = TThMva_->EvaluateMVA( _MVAMethod.c_str() );
         
                if(useTTHHadronicMVA_)
                {
                    BJetVect.clear();
                    BJetVect = BJetTTHHMVAVect;
                }

                TTHHadronicTag tthhtags_obj( dipho, mvares, JetVect, BJetVect );
                tthhtags_obj.setNjet( jetcount_ );
                tthhtags_obj.setNBLoose( njets_btagloose_ );
                tthhtags_obj.setNBMedium( njets_btagmedium_ );
                tthhtags_obj.setNBTight( njets_btagtight_ );
                tthhtags_obj.setDiPhotonIndex( diphoIndex );
                tthhtags_obj.setLeadJetPt( leadJetPt_ );
                tthhtags_obj.setSubLeadJetPt( subLeadJetPt_ );
                tthhtags_obj.setSumJetPt( sumJetPt_ );
                tthhtags_obj.setMaxBTagVal( maxBTagVal_ );
                tthhtags_obj.setSecondMaxBTagVal( secondMaxBTagVal_ );
                tthhtags_obj.setSystLabel( systLabel_ );
                tthhtags_obj.setMVAres(tthMvaVal_);
                tthhtags_obj.setHasRecoLeptons(hasRecoLeptons);
                tthhtags_obj.setDiphoIndex(diphoIndex);

                if(!useTTHHadronicMVA_)
                {
                    for( unsigned num = 0; num < JetVect.size(); num++ )
                    {
                        tthhtags_obj.includeWeightsByLabel( *JetVect[num] , "JetBTagCutWeight");
                    }
                }
                else
                {
                    for( unsigned num = 0; num < JetVect.size(); num++ )
                    {
                        tthhtags_obj.includeWeightsByLabel( *JetVect[num] , "JetBTagReshapeWeight");
                    }                    
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

            if(!isLeptonic)
            {   
                TTHHadronicEfficiency tthhe(genJetsVector, genPhotons, ttHtags);
                tthhe.setHiggsVertex(higgsVtx);
                tthhe.setVertex0(vtx0);
                ttheff->push_back( tthhe );
             }

        }

        evt.put( ttheff );


       // cout << "tagged events = " << count << endl;

    }
}
typedef flashgg::TTHHadronicEfficiencyProducer FlashggTTHHadronicEfficiencyProducer;
DEFINE_FWK_MODULE( FlashggTTHHadronicEfficiencyProducer );

// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4

