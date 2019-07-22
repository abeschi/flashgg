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
#include "flashgg/DataFormats/interface/DiElectronCandidate.h"
#include "flashgg/DataFormats/interface/EleEleGammaCandidate.h"
#include "flashgg/DataFormats/interface/ZGammaToEEUntaggedTag.h"
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
    class ZGammaToEEUntaggedTagProducer : public EDProducer
    {

    public:
        typedef math::XYZPoint Point;

        ZGammaToEEUntaggedTagProducer( const ParameterSet & );
    private:
        void produce( Event &, const EventSetup & ) override;
        int  chooseCategory( float MVA , bool debug_ );
        bool checkEleId(flashgg::Electron ele, std::string WP);
        void makeEleEleGammaCandidate(std::vector<flashgg::EleEleGammaCandidate>* myEEG, Handle<View<flashgg::Electron>> Electrons, Handle<View<flashgg::Photon>> Photons, edm::Ptr<reco::Vertex> pvx, bool debug_);

        std::vector<edm::EDGetTokenT<View<flashgg::Jet> > > tokenJets_;
        std::vector<edm::InputTag> inputTagJets_;
        EDGetTokenT<View<flashgg::EleEleGammaCandidate> > EleEleGammaToken_;
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
        bool useEleEleGammaCandidate_;
        bool chooseByZMass_;

        //Thresholds
        double preselectedElePt_;
        double preselectedEleEta_;
        double preselectedPhoPt_;
        double MinLeadElePt_;
        double MinSubleadElePt_;
        double MaxEleEta_;
        double minDiEleMass_;
        double minLeptDR_;
        string EleIdWP_;
        double minEEGMass_;
        double maxEEGMass_;
        double minEEGPlusDiEleMass_;
        double minPhoId_;
        double minPtGammaOverMass_;
        double minLeptPhoDR_;
        std::vector<double> mvaBoundaries_;
        
    };


    ZGammaToEEUntaggedTagProducer::ZGammaToEEUntaggedTagProducer( const ParameterSet &iConfig ) :
        inputTagJets_( iConfig.getParameter<std::vector<edm::InputTag> >( "inputTagJets" ) ),
        EleEleGammaToken_( consumes<View<flashgg::EleEleGammaCandidate> >( iConfig.getParameter<InputTag>( "EleEleGammaTag" ) ) ),
        photonToken_( consumes<View<flashgg::Photon> >( iConfig.getParameter<InputTag>( "PhotonTag" ) ) ),
        electronToken_( consumes<View<flashgg::Electron> >( iConfig.getParameter<InputTag>( "ElectronTag" ) ) ),
        muonToken_( consumes<View<flashgg::Muon> >( iConfig.getParameter<InputTag>( "MuonTag" ) ) ),
        METToken_( consumes<View<flashgg::Met> >( iConfig.getParameter<InputTag>( "MetTag" ) ) ),
        vertexToken_( consumes<View<reco::Vertex> >( iConfig.getParameter<InputTag> ( "VertexTag" ) ) ),
        genParticleToken_( consumes<View<reco::GenParticle> >( iConfig.getParameter<InputTag> ( "GenParticleTag" ) ) ),
        systLabel_( iConfig.getParameter<string> ( "SystLabel" ) )
    {
        debug_ = iConfig.getParameter<bool>( "debug");
        useEleEleGammaCandidate_ = iConfig.getParameter<bool>( "useEleEleGammaCandidate");
        chooseByZMass_ = iConfig.getParameter<bool>( "chooseByZMass");
        preselectedElePt_ = iConfig.getParameter<double>( "preselectedElePt" );
        preselectedEleEta_ = iConfig.getParameter<double>( "preselectedEleEta" );
        MinLeadElePt_ = iConfig.getParameter<double>( "MinLeadElePt");
        MinSubleadElePt_ = iConfig.getParameter<double>( "MinSubleadElePt");
        MaxEleEta_ = iConfig.getParameter<double>( "MaxEleEta");
        minDiEleMass_ = iConfig.getParameter<double>( "minDiEleMass");
        minLeptDR_ = iConfig.getParameter<double>( "minLeptDR");
        EleIdWP_ = iConfig.getParameter<string>( "EleIdWP");
        minEEGMass_ = iConfig.getParameter<double>( "minEEGMass");
        maxEEGMass_ = iConfig.getParameter<double>( "maxEEGMass");
        minEEGPlusDiEleMass_ = iConfig.getParameter<double>( "minEEGPlusDiEleMass");
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

        produces<vector<ZGammaToEEUntaggedTag> >();
        produces<vector<TagTruthBase> >();
    }


    bool ZGammaToEEUntaggedTagProducer::checkEleId(flashgg::Electron ele, string WP)
    {
        if(WP=="Veto" || WP=="veto" || WP=="VETO")
            return ele.passVetoId();
        else if(WP=="Loose" || WP=="loose" || WP=="LOOSE")
            return ele.passLooseId();
        else if(WP=="Medium" || WP=="medium" || WP=="MEDIUM")
            return ele.passMediumId();
        else if(WP=="Tight" || WP=="tight" || WP=="TIGHT")
            return ele.passTightId();
        else if(WP=="MVALoose" || WP=="mvaloose" || WP=="mvaLoose" || WP=="MVAloose" || WP=="MVALOOSE")
            return ele.passMVALooseId();
        else if(WP=="MVAMedium" || WP=="mvamedium" || WP=="mvaMedium" || WP=="MVAmedium" || WP=="MVAMedium")
            return ele.passMVAMediumId();
        else if(WP=="MVATight" || WP=="mvatight" || WP=="mvaTight" || WP=="MVAtight" || WP=="MVATIGHT")
            return ele.passMVATightId();
        else if(WP=="MVALooseNoIso" || WP=="mvaloosenoiso" || WP=="mvaLooseNoIso" || WP=="MVAlooseNoIso" || WP=="MVALOOSENOISO")
            return ele.passMVALooseNoIsoId();
        else if(WP=="MVAMediumNoIso" || WP=="mvamediumnoiso" || WP=="mvaMediumNoIso" || WP=="MVAmediumNoIso" || WP=="MVAMEDIUMNOISO")
            return ele.passMVAMediumNoIsoId();
        else if(WP=="MVATightNoIso" || WP=="mvatightnoiso" || WP=="mvaTightNoIso" || WP=="MVAtightNoIso" || WP=="MVATIGHTNOISO")
            return ele.passMVATightNoIsoId();
        else
        {
            cout << "WP provided for Ele Id does not exist. Using MediumMVAId as default" << endl;
            return ele.passMVAMediumId();
        }


    }

    int ZGammaToEEUntaggedTagProducer::chooseCategory( float MVA , bool debug_ )
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

    void ZGammaToEEUntaggedTagProducer::makeEleEleGammaCandidate(std::vector<flashgg::EleEleGammaCandidate>* myEEG, Handle<View<flashgg::Electron>> Electrons, Handle<View<flashgg::Photon>> Photons, edm::Ptr<reco::Vertex> pvx, bool debug_)
    {
        std::vector<flashgg::DiElectronCandidate> diEles;

        for( unsigned int i=0; i < Electrons->size(); ++i )
        {
            if(debug_)
            {
                cout << "Checking electrons preselections, electron candidate 1 " << i << " out of " << Electrons->size() << endl;
                cout << "Electron pT: " << Electrons->ptrAt( i )->pt() << ", eta: " << Electrons->ptrAt( i )->eta() << endl;
            }

            edm::Ptr<flashgg::Electron> Ele1 = Electrons->ptrAt( i );
            if(Ele1->pt()<preselectedElePt_ || fabs(Ele1->eta()>preselectedEleEta_)) continue;
    
            for( unsigned int j = i+1; j < Electrons->size(); ++j )
            {
                if(debug_)
                {
                    cout << "Checking electrons preselections, electron candidate 2 " << j << " out of " << Electrons->size() << endl;
                    cout << "Electron pT: " << Electrons->ptrAt( j )->pt() << ", eta: " << Electrons->ptrAt( j )->eta() << endl;
                }

                edm::Ptr<flashgg::Electron> Ele2 = Electrons->ptrAt( j );

                if(Ele2->pt()<preselectedElePt_ || fabs(Ele2->eta()>preselectedEleEta_)) continue;

                if(Ele1->pt()>Ele2->pt())
                {    DiElectronCandidate diEle( Ele1, Ele2 );
                     diEles.push_back(diEle);

                }
                else
                {    DiElectronCandidate diEle( Ele2, Ele1 );
                     diEles.push_back(diEle);
                }
            }
        }


        if(diEles.size()==0) return;

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

        for( unsigned int i=0; i<diEles.size(); ++i )
        {
             for( unsigned int j=0; j<phos.size(); ++j )
            {
                flashgg::EleEleGammaCandidate eeg(diEles.at(i), phos.at(j), pvx);
                myEEG->push_back(eeg);
            }
        }
        
cout << myEEG->size() << endl;
        return;

    }   

    void ZGammaToEEUntaggedTagProducer::produce( Event &evt, const EventSetup & )
    {
        JetCollectionVector Jets( inputTagJets_.size() );
        for( size_t j = 0; j < inputTagJets_.size(); ++j ) {
            evt.getByToken( tokenJets_[j], Jets[j] );
        }

        Handle<View<flashgg::EleEleGammaCandidate> > EleEleGamma;
        if(useEleEleGammaCandidate_)
                    evt.getByToken( EleEleGammaToken_, EleEleGamma );

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

        std::unique_ptr<vector<ZGammaToEEUntaggedTag> > ZGToEETag( new vector<ZGammaToEEUntaggedTag> );
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

        //If the EleEleGamma does not exist first create the collection
        std::vector<flashgg::EleEleGammaCandidate> myEEG;
        if(!useEleEleGammaCandidate_)
        {    
            if(debug_)
                cout << "Creating the collection." << endl;

            makeEleEleGammaCandidate(&myEEG, Electrons, Photons, pvx, debug_);

           if(debug_)
                cout << "Collection created, number of candidates: " << myEEG.size() << endl;

        }
        else
        {
            for( unsigned int index = 0; index < EleEleGamma->size(); ++index )
            {
                edm::Ptr<flashgg::EleEleGammaCandidate> EEG = EleEleGamma->ptrAt( index );
                myEEG.push_back(*EEG);
            }
            if(debug_)
                cout << "Using EleEleGammaCandidate collection. Number of candidates " << myEEG.size() << endl;
        }

        //Loop on EleEleGamma candidates to select the events
        std::vector<flashgg::EleEleGammaCandidate> selectedEEG;

        for(unsigned int EEGIndex=0; EEGIndex<myEEG.size(); ++EEGIndex)
        {
            if(debug_)
                cout << "Analysing candidate " << EEGIndex << " out of " << myEEG.size() << endl;

            const flashgg::Electron* Ele1 = myEEG.at(EEGIndex).leadingElectron();
            const flashgg::Electron* Ele2 = myEEG.at(EEGIndex).subleadingElectron();

            if(debug_)
            {
                cout << "Checking DiElectron selections" << endl;
                cout << "Electron 1: pT=" << Ele1->pt() << ", eta=" << Ele1->eta() << ", phi=" << Ele1->phi() << ", MVAMediumId=" << Ele1->passMVAMediumId() << endl;
                cout << "Electron 2: pT=" << Ele2->pt() << ", eta=" << Ele2->eta() << ", phi=" << Ele2->phi() << ", MVAMediumId=" << Ele2->passMVAMediumId() << endl;
                cout << "DeltaR between electrons: " << deltaR(Ele1->eta(), Ele1->phi(), Ele2->eta(), Ele2->phi()) << endl;
                cout << "DiElectron invariant mass: " << myEEG.at(EEGIndex).getDiEleMass() << endl;
                cout << "Three body system invariant mass: " << myEEG.at(EEGIndex).getMass() << endl;
            }

            if(Ele1->pt()<MinLeadElePt_ || Ele2->pt()<MinSubleadElePt_) continue;
            if(fabs(Ele1->eta())>MaxEleEta_ || fabs(Ele2->eta())>MaxEleEta_) continue;
            if(myEEG.at(EEGIndex).getDiEleMass() < minDiEleMass_) continue;
            if( deltaR(Ele1->eta(), Ele1->phi(), Ele2->eta(), Ele2->phi()) < minLeptDR_) continue;
            if(!checkEleId(*Ele1, EleIdWP_) || !checkEleId(*Ele2, EleIdWP_)) continue;
            if(myEEG.at(EEGIndex).getMass()<minEEGMass_ || myEEG.at(EEGIndex).getMass()>maxEEGMass_) continue;
            if(myEEG.at(EEGIndex).getDiEleMass() + myEEG.at(EEGIndex).getMass()<minEEGPlusDiEleMass_) continue;

            const flashgg::Photon* pho = myEEG.at(EEGIndex).EEG_Photon();


            if(debug_)
            {
                cout << "DiElectron selection passed, now checking photon ones" << endl;
                cout << "Photon 1: pT=" << pho->pt() << ", eta=" << pho->eta() << ", phi=" <<  pho->phi() << ", PhoId=" << pho->phoIdMvaDWrtVtx( myEEG.at(EEGIndex).Vertex()) << endl;
                cout << "DeltaR between photon and Electron 1: " << deltaR(pho->eta(), pho->phi(), Ele1->eta(), Ele1->phi()) << endl;
                cout << "DeltaR between photon and Electron 2: " << deltaR(pho->eta(), pho->phi(), Ele2->eta(), Ele2->phi()) << endl;
                cout << "Photon pT over mass  " << pho->pt()/myEEG.at(EEGIndex).getMass() << endl;
            }

            if(pho->phoIdMvaDWrtVtx( myEEG.at(EEGIndex).Vertex()) < minPhoId_) continue;
            if(pho->pt()/myEEG.at(EEGIndex).getMass()<minPtGammaOverMass_) continue;
            if( deltaR(pho->eta(), pho->phi(), Ele1->eta(), Ele1->phi()) < minLeptPhoDR_ || deltaR(pho->eta(), pho->phi(), Ele2->eta(), Ele2->phi()) < minLeptPhoDR_ ) continue;

            if(debug_)
                cout << "Candidate " << EEGIndex << " passed the selections" << endl;

            selectedEEG.push_back(myEEG.at(EEGIndex));
        }

        if(debug_)
            cout << "Selected a total of " << selectedEEG.size() << " candidates" << endl;

        //Choose candidates if more than one...
        //This is just an idea, the logic should be improved!
        std::vector<flashgg::EleEleGammaCandidate> finalEEGCandidate;
        if(selectedEEG.size()>1)
        {
            if(chooseByZMass_)
            {
                if(debug_)
                    cout << "Selecting candidate with diElectron mass closest to Z boson one" << endl;

                int id= -1;
                double mass = 1e6;
                for(unsigned int index=0; index<selectedEEG.size(); ++index)
                {   double deltaM = fabs(selectedEEG.at(index).getDiEleMass() - MZ);
                    if(debug_)
                        cout << "Candidate " << index << " DiElectron mass: " << selectedEEG.at(index).getDiEleMass() << ", difference from Z boson mass: " << deltaM << endl; 
                    if(deltaM<mass)
                    {
                        mass = deltaM;
                        id = index;
                    }   
                }

                finalEEGCandidate.push_back(selectedEEG.at(id));
                if(debug_)
                    cout << "Selected candidate " << id << endl;
            }

            else
            {
                if(debug_)
                    cout << "Selecting candidate with highest three body pT" << endl;

                int id= -1;
                double maxPt = -1.;
                for(unsigned int index=0; index<selectedEEG.size(); ++index)
                {   double pt = fabs(selectedEEG.at(index).getPt());
                    if(debug_)
                        cout << "Candidate " << index << " pT: " << pt << endl; 
                    if(pt>maxPt)
                    {
                        maxPt = pt;
                        id = index;
                    }   
                }

                finalEEGCandidate.push_back(selectedEEG.at(id));
                if(debug_)
                    cout << "Selected candidate " << id << endl;
            }
        }
        else if(selectedEEG.size()==1)
            finalEEGCandidate.push_back(selectedEEG.at(0));

        //MVA categories
        if(finalEEGCandidate.size()>0)
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
                flashgg::EleEleGammaCandidate eeg = finalEEGCandidate.at(0);

                ZGammaToEEUntaggedTag ZGToEETag_obj( eeg );
                ZGToEETag_obj.setCategoryNumber(cat);

                //ZGToEETag_obj.includeWeights( *dipho ); HERE ADD WEIGHTS FOR SYST

                ZGToEETag_obj.setSystLabel( systLabel_ );
                ZGToEETag_obj.setMva(mvaValue);
                ZGToEETag -> push_back( ZGToEETag_obj );

                if( ! evt.isRealData() )
                {
                    TagTruthBase truth_obj;
                    truth_obj.setGenPV( higgsVtx );

                    truths -> push_back( truth_obj );
                    ZGToEETag -> back().setTagTruth( edm::refToPtr( edm::Ref<vector<TagTruthBase> >( rTagTruth, idx++ ) ) );
                }
            }
        }
        assert(ZGToEETag->size()==1 || ZGToEETag->size()==0); // Produce at most one tag per event

        evt.put( std::move( ZGToEETag ) );
        evt.put( std::move( truths ) );
    }

}
typedef flashgg::ZGammaToEEUntaggedTagProducer FlashggZGammaToEEUntaggedTagProducer;
DEFINE_FWK_MODULE( FlashggZGammaToEEUntaggedTagProducer );
// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4

