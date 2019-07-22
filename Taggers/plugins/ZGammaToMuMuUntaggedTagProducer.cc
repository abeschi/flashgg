#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDMException.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "flashgg/DataFormats/interface/Jet.h"
#include "flashgg/DataFormats/interface/DiMuonCandidate.h"
#include "flashgg/DataFormats/interface/MuMuGammaCandidate.h"
#include "flashgg/DataFormats/interface/ZGammaToMuMuUntaggedTag.h"
#include "flashgg/DataFormats/interface/Electron.h"
#include "flashgg/DataFormats/interface/Muon.h"
#include "flashgg/DataFormats/interface/Met.h"
#include "flashgg/DataFormats/interface/Photon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "flashgg/DataFormats/interface/TagTruthBase.h"
#include "DataFormats/Common/interface/RefToPtr.h"
#include "SimDataFormats/HTXS/interface/HiggsTemplateCrossSections.h"

#include <vector>
#include <algorithm>
#include <string>
#include <utility>
#include "TLorentzVector.h"
#include "TMath.h"
#include "TRandom.h"
#include "TMVA/Reader.h"

#define MZ 91.1876


using namespace std;
using namespace edm;


namespace flashgg {
    class ZGammaToMuMuUntaggedTagProducer : public EDProducer
    {

    public:
        typedef math::XYZPoint Point;

        ZGammaToMuMuUntaggedTagProducer( const ParameterSet & );
    private:
        void produce( Event &, const EventSetup & ) override;
        int  chooseCategory( float MVA , bool debug_ );
        bool checkMuId(flashgg::Muon mu, std::string WP, edm::Ptr<reco::Vertex>);
        void makeMuMuGammaCandidate(std::vector<flashgg::MuMuGammaCandidate>* myMMG, Handle<View<flashgg::Muon>> Muons, Handle<View<flashgg::Photon>> Photons, edm::Ptr<reco::Vertex> pvx, bool debug_);

        std::vector<edm::EDGetTokenT<View<flashgg::Jet> > > tokenJets_;
        std::vector<edm::InputTag> inputTagJets_;
        EDGetTokenT<View<flashgg::MuMuGammaCandidate> > MuMuGammaToken_;
        EDGetTokenT<View<flashgg::Photon> > photonToken_;
        EDGetTokenT<View<Electron> > electronToken_;
        EDGetTokenT<View<flashgg::Muon> > muonToken_;
        EDGetTokenT<View<flashgg::Met> > METToken_;
        EDGetTokenT<View<reco::Vertex> > vertexToken_;
        EDGetTokenT<View<reco::GenParticle> > genParticleToken_;
        string systLabel_;

        typedef std::vector<edm::Handle<edm::View<flashgg::Jet> > > JetCollectionVector;

//        unique_ptr<TMVA::Reader> MVA_;
//        FileInPath MVAweightfile_;

        bool debug_;
        bool useMuMuGammaCandidate_;
        bool chooseByZMass_;

        //Thresholds
        double preselectedMuonPt_;
        double preselectedMuonEta_;
        double preselectedPhoPt_;
        double MinLeadMuonPt_;
        double MinSubleadMuonPt_;
        double MaxMuonEta_;
        double MuonIso_;
        double minDiMuonMass_;
        double minLeptDR_;
        string MuIdWP_;
        double minMMGMass_;
        double maxMMGMass_;
        double minMMGPlusDiMuonMass_;
        double minPhoId_;
        double minPtGammaOverMass_;
        double minLeptPhoDR_;
        std::vector<double> mvaBoundaries_;
        
    };


    ZGammaToMuMuUntaggedTagProducer::ZGammaToMuMuUntaggedTagProducer( const ParameterSet &iConfig ) :
        inputTagJets_( iConfig.getParameter<std::vector<edm::InputTag> >( "inputTagJets" ) ),
        MuMuGammaToken_( consumes<View<flashgg::MuMuGammaCandidate> >( iConfig.getParameter<InputTag>( "MuMuGammaTag" ) ) ),
        photonToken_( consumes<View<flashgg::Photon> >( iConfig.getParameter<InputTag>( "PhotonTag" ) ) ),
        electronToken_( consumes<View<flashgg::Electron> >( iConfig.getParameter<InputTag>( "ElectronTag" ) ) ),
        muonToken_( consumes<View<flashgg::Muon> >( iConfig.getParameter<InputTag>( "MuonTag" ) ) ),
        METToken_( consumes<View<flashgg::Met> >( iConfig.getParameter<InputTag>( "MetTag" ) ) ),
        vertexToken_( consumes<View<reco::Vertex> >( iConfig.getParameter<InputTag> ( "VertexTag" ) ) ),
        genParticleToken_( consumes<View<reco::GenParticle> >( iConfig.getParameter<InputTag> ( "GenParticleTag" ) ) ),
        systLabel_( iConfig.getParameter<string> ( "SystLabel" ) )
    {
        debug_ = iConfig.getParameter<bool>( "debug");
        useMuMuGammaCandidate_ = iConfig.getParameter<bool>( "useMuMuGammaCandidate");
        chooseByZMass_ = iConfig.getParameter<bool>( "chooseByZMass");
        preselectedMuonPt_ = iConfig.getParameter<double>( "preselectedMuonPt" );
        preselectedMuonEta_ = iConfig.getParameter<double>( "preselectedMuonEta" );
        MinLeadMuonPt_ = iConfig.getParameter<double>( "MinLeadMuonPt");
        MinSubleadMuonPt_ = iConfig.getParameter<double>( "MinSubleadMuonPt");
        MaxMuonEta_ = iConfig.getParameter<double>( "MaxMuonEta");
        MuonIso_ = iConfig.getParameter<double>( "MuonIso");
        minDiMuonMass_ = iConfig.getParameter<double>( "minDiMuonMass");
        minLeptDR_ = iConfig.getParameter<double>( "minLeptDR");
        MuIdWP_ = iConfig.getParameter<string>( "MuIdWP");
        minMMGMass_ = iConfig.getParameter<double>( "minMMGMass");
        maxMMGMass_ = iConfig.getParameter<double>( "maxMMGMass");
        minMMGPlusDiMuonMass_ = iConfig.getParameter<double>( "minMMGPlusDiMuonMass");
        minPhoId_ = iConfig.getParameter<double>( "minPhoId");
        minPtGammaOverMass_ = iConfig.getParameter<double>( "minPtGammaOverMass");
        minLeptPhoDR_ = iConfig.getParameter<double>( "minLeptPhoDR");
        mvaBoundaries_ = iConfig.getParameter<std::vector<double>>( "mvaBoundaries" );

/*
        MVA_.reset( new TMVA::Reader( "!Color:Silent" ) );
        MVA_->AddVariable( "dipho_leadEta", &leadeta_ );
        MVA_->AddVariable( "dipho_subleadEta", &subleadeta_ );
        MVA_->AddVariable( "dipho_lead_ptoM", &leadptom_ );
        MVA_->AddVariable( "dipho_sublead_ptoM", &subleadptom_ );
        MVA_->AddVariable( "dipho_leadIDMVA", &leadIDMVA_ );
        MVA_->AddVariable( "dipho_subleadIDMVA", &subleadIDMVA_ );
        MVA_->AddVariable( "dipho_deltaphi", &deltaphi_ );
        MVA_->AddVariable( "dipho_lead_PSV", &leadPSV_ );
        MVA_->AddVariable( "dipho_sublead_PSV", &subleadPSV_ );
        MVA_->AddVariable( "nJets", &nJets_ );
        MVA_->AddVariable( "nJets_bTagMedium", &nJets_bTagMedium_ );
        MVA_->AddVariable( "jet1_pt", &jet_pt1_ );
        MVA_->AddVariable( "jet2_pt", &jet_pt2_ );
        MVA_->AddVariable( "jet3_pt", &jet_pt3_ );
        MVA_->AddVariable( "jet1_eta", &jet_eta1_ );
        MVA_->AddVariable( "jet2_eta", &jet_eta2_ );
        MVA_->AddVariable( "jet3_eta", &jet_eta3_ );
        MVA_->AddVariable( "bTag1", &bTag1_ );
        MVA_->AddVariable( "bTag2", &bTag2_ );
        MVA_->AddVariable( "MetPt", &MetPt_ );
        MVA_->AddVariable( "lepton_leadPt", &lepton_leadPt_ );
        MVA_->AddVariable( "lepton_leadEta", &lepton_leadEta_ );

        MVA_->BookMVA( "BDT", MVAweightfile_.fullPath() );
*/

        for (unsigned i = 0 ; i < inputTagJets_.size() ; i++) {
            auto token = consumes<View<flashgg::Jet> >(inputTagJets_[i]);
            tokenJets_.push_back(token);
        }

        produces<vector<ZGammaToMuMuUntaggedTag> >();
        produces<vector<TagTruthBase> >();
    }


    bool ZGammaToMuMuUntaggedTagProducer::checkMuId(flashgg::Muon mu, string WP, edm::Ptr<reco::Vertex> pvx)
    {
        if(WP=="Loose" || WP=="loose" || WP=="LOOSE")
            return mu.isLooseMuon();
        else if(WP=="Medium" || WP=="medium" || WP=="MEDIUM")
            return mu.isMediumMuon();
        else if(WP=="Tight" || WP=="tight" || WP=="TIGHT")
            return mu.isTightMuon(*pvx);
        else
        {
            cout << "WP provided for Muon Id does not exist. Using Medium as default" << endl;
            return mu.isMediumMuon();
        }
    }

    int ZGammaToMuMuUntaggedTagProducer::chooseCategory( float MVA , bool debug_ )
    {
        // should return 0 if mva above all the numbers, 1 if below the first, ..., boundaries.size()-N if below the Nth, ...

        if(debug_)
            cout << "Selecting the category" << endl;

        std::vector<double> bound =mvaBoundaries_;
        //bound.push_back(1.);
        bound.push_back(MVA);
        std::sort(bound.begin(), bound.end());
        std::vector<double>::iterator it = std::find(bound.begin(), bound.end(), MVA);
        int index = std::distance(bound.begin(), it);

        int cat = ((int)bound.size() -1 - index)== (int)mvaBoundaries_.size() ? -1 : (int)bound.size() -1 - index;

        if(debug_)
        {   cout << "MVA boundaries: ";
            for(unsigned int i=0; i<mvaBoundaries_.size(); ++i)
                cout << mvaBoundaries_[i] << " ";
            cout << endl;
            cout << "MVA value: " << MVA << endl;
            cout << "Category: " << cat << endl;
        }

        return cat;
    }

    void ZGammaToMuMuUntaggedTagProducer::makeMuMuGammaCandidate(std::vector<flashgg::MuMuGammaCandidate>* myMMG, Handle<View<flashgg::Muon>> Muons, Handle<View<flashgg::Photon>> Photons, edm::Ptr<reco::Vertex> pvx, bool debug_)
    {
        std::vector<flashgg::DiMuonCandidate> diMuons;

        for( unsigned int i=0; i < Muons->size(); ++i )
        {
            if(debug_)
            {
                cout << "Checking muons preselections, muon candidate 1 " << i << " out of " << Muons->size() << endl;
                cout << "Muon pT: " << Muons->ptrAt( i )->pt() << ", eta: " << Muons->ptrAt( i )->eta() << endl;
            }

            edm::Ptr<flashgg::Muon> Mu1 = Muons->ptrAt( i );
            if(Mu1->pt()<preselectedMuonPt_ || fabs(Mu1->eta()>preselectedMuonEta_)) continue;
    
            for( unsigned int j = i+1; j < Muons->size(); ++j )
            {
                if(debug_)
                {
                    cout << "Checking muons preselections, muon candidate 2 " << j << " out of " << Muons->size() << endl;
                    cout << "Muon pT: " << Muons->ptrAt( j )->pt() << ", eta: " << Muons->ptrAt( j )->eta() << endl;
                }

                edm::Ptr<flashgg::Muon> Mu2 = Muons->ptrAt( j );

                if(Mu2->pt()<preselectedMuonPt_ || fabs(Mu2->eta()>preselectedMuonEta_)) continue;

                if(Mu1->pt()>Mu2->pt())
                {    DiMuonCandidate diMu( Mu1, Mu2 );
                     diMuons.push_back(diMu);

                }
                else
                {    DiMuonCandidate diMu( Mu2, Mu1 );
                     diMuons.push_back(diMu);
                }
            }
        }


        if(diMuons.size()==0) return;

        std::vector<flashgg::Photon> phos;
        for( unsigned int i=0; i < Photons->size(); ++i )
        {
            if(debug_)
            {
                cout << "Checking photon preselections, photon candidate " << i << " out of " << Photons->size() << endl;
                cout << "Photon pT: " << Photons->ptrAt( i )->pt() << endl;
            }

            edm::Ptr<flashgg::Photon> pho = Photons->ptrAt( i );
            if(pho->pt()<preselectedPhoPt_) continue;
            phos.push_back(*pho);
        }

        if(phos.size()==0) return;

        for( unsigned int i=0; i<diMuons.size(); ++i )
        {
             for( unsigned int j=0; j<phos.size(); ++j )
            {
                flashgg::MuMuGammaCandidate mmg(diMuons.at(i), phos.at(j), pvx);
                myMMG->push_back(mmg);
            }
        }
        
        cout << myMMG->size() << endl;
        return;

    }   

    void ZGammaToMuMuUntaggedTagProducer::produce( Event &evt, const EventSetup & )
    {
        JetCollectionVector Jets( inputTagJets_.size() );
        for( size_t j = 0; j < inputTagJets_.size(); ++j ) {
            evt.getByToken( tokenJets_[j], Jets[j] );
        }

        Handle<View<flashgg::MuMuGammaCandidate> > MuMuGamma;
        if(useMuMuGammaCandidate_)
                    evt.getByToken( MuMuGammaToken_, MuMuGamma );

        Handle<View<flashgg::Photon> > Photons;
        evt.getByToken( photonToken_, Photons );

        Handle<View<flashgg::Electron> > Electrons;
        evt.getByToken( electronToken_, Electrons );

        Handle<View<flashgg::Muon> > Muons;
        evt.getByToken( muonToken_, Muons );

        Handle<View<flashgg::Met> > Met;
        evt.getByToken( METToken_, Met );

        Handle<View<reco::GenParticle> > genParticles;

        Handle<View<reco::Vertex> > vertices;
        evt.getByToken( vertexToken_, vertices );
        const std::vector<edm::Ptr<reco::Vertex>> &pvPointers = vertices->ptrs();
        edm::Ptr<reco::Vertex> pvx = pvPointers[0];

        std::unique_ptr<vector<ZGammaToMuMuUntaggedTag> > ZGToMuMuTag( new vector<ZGammaToMuMuUntaggedTag> );
        std::unique_ptr<vector<TagTruthBase> > truths( new vector<TagTruthBase> );
        edm::RefProd<vector<TagTruthBase> > rTagTruth = evt.getRefBeforePut<vector<TagTruthBase> >();
        unsigned int idx = 0;

        Point higgsVtx;

        if( ! evt.isRealData() )
        {
            evt.getByToken( genParticleToken_, genParticles );
            for( unsigned int genLoop = 0 ; genLoop < genParticles->size(); genLoop++ )
            {
                int pdgid = genParticles->ptrAt( genLoop )->pdgId();
                if( pdgid == 25 || pdgid == 22 )
                {
                    higgsVtx = genParticles -> ptrAt( genLoop )-> vertex();
                    break;
                }
            }
        }

        //If the MuMuGamma does not exist first create the collection
        std::vector<flashgg::MuMuGammaCandidate> myMMG;
        if(!useMuMuGammaCandidate_)
        {    
            if(debug_)
                cout << "Creating the collection." << endl;

            makeMuMuGammaCandidate(&myMMG, Muons, Photons, pvx, debug_);

           if(debug_)
                cout << "Collection created, number of candidates: " << myMMG.size() << endl;

        }
        else
        {
            for( unsigned int index = 0; index < MuMuGamma->size(); ++index )
            {
                edm::Ptr<flashgg::MuMuGammaCandidate> MMG = MuMuGamma->ptrAt( index );
                myMMG.push_back(*MMG);
            }
            if(debug_)
                cout << "Using MuMuGammaCandidate collection. Number of candidates " << myMMG.size() << endl;
        }

        //Loop on MuMuGamma candidates to select the events
        std::vector<flashgg::MuMuGammaCandidate> selectedMMG;

        for(unsigned int MMGIndex=0; MMGIndex<myMMG.size(); ++MMGIndex)
        {
            if(debug_)
                cout << "Analysing candidate " << MMGIndex << " out of " << myMMG.size() << endl;

            const flashgg::Muon* Mu1 = myMMG.at(MMGIndex).leadingMuon();
            const flashgg::Muon* Mu2 = myMMG.at(MMGIndex).subleadingMuon();

            if(debug_)
            {
                cout << "Checking DiMuon selections" << endl;
                cout << "Muon 1: pT=" << Mu1->pt() << ", eta=" << Mu1->eta() << ", phi=" << Mu1->phi() << ", MediumId=" << Mu1->isMediumMuon() << endl;
                cout << "Muon 2: pT=" << Mu2->pt() << ", eta=" << Mu2->eta() << ", phi=" << Mu2->phi() << ", MediumId=" << Mu2->isMediumMuon() << endl;
                cout << "DeltaR between muons: " << deltaR(Mu1->eta(), Mu1->phi(), Mu2->eta(), Mu2->phi()) << endl;
                cout << "DiMuon invariant mass: " << myMMG.at(MMGIndex).getDiMuonMass() << endl;
                cout << "Three body system invariant mass: " << myMMG.at(MMGIndex).getMass() << endl;
            }

            if(Mu1->pt()<MinLeadMuonPt_ || Mu2->pt()<MinSubleadMuonPt_) continue;
            if(fabs(Mu1->eta())>MaxMuonEta_ || fabs(Mu2->eta())>MaxMuonEta_) continue;
            if(myMMG.at(MMGIndex).getDiMuonMass() < minDiMuonMass_) continue;
            if( deltaR(Mu1->eta(), Mu1->phi(), Mu2->eta(), Mu2->phi()) < minLeptDR_) continue;
            if(!checkMuId(*Mu1, MuIdWP_, pvx) || !checkMuId(*Mu2, MuIdWP_, pvx)) continue;
            if(myMMG.at(MMGIndex).getMass()<minMMGMass_ || myMMG.at(MMGIndex).getMass()>maxMMGMass_) continue;
            if(myMMG.at(MMGIndex).getDiMuonMass() + myMMG.at(MMGIndex).getMass()<minMMGPlusDiMuonMass_) continue;

            //Check muon isolation

             double muPFIsoSumRel1 = ( Mu1->pfIsolationR04().sumChargedHadronPt + max( 0., Mu1->pfIsolationR04().sumNeutralHadronEt + Mu1->pfIsolationR04().sumPhotonEt - 0.5 * Mu1->pfIsolationR04().sumPUPt ) ) / ( Mu1->pt() );
             double muPFIsoSumRel2 = ( Mu2->pfIsolationR04().sumChargedHadronPt + max( 0., Mu2->pfIsolationR04().sumNeutralHadronEt + Mu2->pfIsolationR04().sumPhotonEt - 0.5 * Mu2->pfIsolationR04().sumPUPt ) ) / ( Mu2->pt() );

            if( muPFIsoSumRel1 > MuonIso_ || muPFIsoSumRel2 > MuonIso_ ) continue; 


            const flashgg::Photon* pho = myMMG.at(MMGIndex).MMG_Photon();

            if(debug_)
            {
                cout << "DiMuon selection passed, now checking photon ones" << endl;
                cout << "Photon 1: pT=" << pho->pt() << ", eta=" << pho->eta() << ", phi=" <<  pho->phi() << ", PhoId=" << pho->phoIdMvaDWrtVtx( myMMG.at(MMGIndex).Vertex()) << endl;
                cout << "DeltaR between photon and Muon 1: " << deltaR(pho->eta(), pho->phi(), Mu1->eta(), Mu1->phi()) << endl;
                cout << "DeltaR between photon and Muon 2: " << deltaR(pho->eta(), pho->phi(), Mu2->eta(), Mu2->phi()) << endl;
                cout << "Photon pT over mass  " << pho->pt()/myMMG.at(MMGIndex).getMass() << endl;
            }

            if(pho->phoIdMvaDWrtVtx( myMMG.at(MMGIndex).Vertex()) < minPhoId_) continue;
            if(pho->pt()/myMMG.at(MMGIndex).getMass()<minPtGammaOverMass_) continue;
            if( deltaR(pho->eta(), pho->phi(), Mu1->eta(), Mu1->phi()) < minLeptPhoDR_ || deltaR(pho->eta(), pho->phi(), Mu2->eta(), Mu2->phi()) < minLeptPhoDR_ ) continue;

            if(debug_)
                cout << "Candidate " << MMGIndex << " passed the selections" << endl;

            selectedMMG.push_back(myMMG.at(MMGIndex));
        }

        if(debug_)
            cout << "Selected a total of " << selectedMMG.size() << " candidates" << endl;

        //Choose candidates if more than one...
        //This is just an idea, the logic should be improved!
        std::vector<flashgg::MuMuGammaCandidate> finalMMGCandidate;
        if(selectedMMG.size()>1)
        {
            if(chooseByZMass_)
            {
                if(debug_)
                    cout << "Selecting candidate with diMuon mass closest to Z boson one" << endl;

                int id= -1;
                double mass = 1e6;
                for(unsigned int index=0; index<selectedMMG.size(); ++index)
                {   double deltaM = fabs(selectedMMG.at(index).getDiMuonMass() - MZ);
                    if(debug_)
                        cout << "Candidate " << index << " DiElectron mass: " << selectedMMG.at(index).getDiMuonMass() << ", difference from Z boson mass: " << deltaM << endl; 
                    if(deltaM<mass)
                    {
                        mass = deltaM;
                        id = index;
                    }   
                }

                finalMMGCandidate.push_back(selectedMMG.at(id));
                if(debug_)
                    cout << "Selected candidate " << id << endl;
            }

            else
            {
                if(debug_)
                    cout << "Selecting candidate with highest three body pT" << endl;

                int id= -1;
                double maxPt = -1.;
                for(unsigned int index=0; index<selectedMMG.size(); ++index)
                {   double pt = fabs(selectedMMG.at(index).getPt());
                    if(debug_)
                        cout << "Candidate " << index << " pT: " << pt << endl; 
                    if(pt>maxPt)
                    {
                        maxPt = pt;
                        id = index;
                    }   
                }

                finalMMGCandidate.push_back(selectedMMG.at(id));
                if(debug_)
                    cout << "Selected candidate " << id << endl;
            }
        }
        else if(selectedMMG.size()==1)
            finalMMGCandidate.push_back(selectedMMG.at(0));

        //MVA categories
        if(finalMMGCandidate.size()>0)
        {    if(debug_)
                    cout << "Checking category of the candidate" << endl; 

            // At this stage the MVA to categorize the event should be evluated
            // For now the mva is just a random number to test the code

            double mvaValue = gRandom -> Uniform();
            int cat = chooseCategory(mvaValue, debug_);
            if(debug_)
                    cout << "Value of the MVA: " << mvaValue << ", category: " << cat << endl;

            //Create the final object
            if(cat!=-1)
            {
                flashgg::MuMuGammaCandidate mmg = finalMMGCandidate.at(0);

                ZGammaToMuMuUntaggedTag ZGToMuMuTag_obj( mmg );
                ZGToMuMuTag_obj.setCategoryNumber(cat);

                //ZGToEETag_obj.includeWeights( *dipho ); HERE ADD WEIGHTS FOR SYST

                ZGToMuMuTag_obj.setSystLabel( systLabel_ );
                ZGToMuMuTag_obj.setMva(mvaValue);
                ZGToMuMuTag -> push_back( ZGToMuMuTag_obj );

                if( ! evt.isRealData() )
                {
                    TagTruthBase truth_obj;
                    truth_obj.setGenPV( higgsVtx );

                    truths -> push_back( truth_obj );
                    ZGToMuMuTag -> back().setTagTruth( edm::refToPtr( edm::Ref<vector<TagTruthBase> >( rTagTruth, idx++ ) ) );
                }
            }
        }
        assert(ZGToMuMuTag->size()==1 || ZGToMuMuTag->size()==0); // Produce at most one tag per event

        evt.put( std::move( ZGToMuMuTag ) );
        evt.put( std::move( truths ) );
    }

}
typedef flashgg::ZGammaToMuMuUntaggedTagProducer FlashggZGammaToMuMuUntaggedTagProducer;
DEFINE_FWK_MODULE( FlashggZGammaToMuMuUntaggedTagProducer );
// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4

