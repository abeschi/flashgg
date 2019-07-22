#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDMException.h"

#include "flashgg/DataFormats/interface/DiPhotonCandidate.h"
#include "flashgg/DataFormats/interface/Jet.h"
#include "flashgg/DataFormats/interface/VBFMVAResult.h"
#include "flashgg/DataFormats/interface/DiPhotonTagBase.h"
#include "flashgg/DataFormats/interface/UntaggedTag.h"
#include "flashgg/DataFormats/interface/TagTruthBase.h"
#include "DataFormats/Common/interface/RefToPtr.h"
#include "flashgg/DataFormats/interface/VBFTag.h"
#include "flashgg/DataFormats/interface/NoTag.h"

#include "SimDataFormats/HTXS/interface/HiggsTemplateCrossSections.h"

#include "TMVA/Reader.h"
#include "TMath.h"
//#include <typeinfo>

#include <algorithm>


using namespace std;
using namespace edm;

namespace flashgg {

    struct TagPriorityRange
    {
        string name;
        int minCat;
        int maxCat;
        unsigned int collIndex;

        TagPriorityRange( string s, int c1, int c2, unsigned int i )
        {
            name = s;
            minCat = c1;
            maxCat = c2;
            collIndex = i;
        }
    };
//Assumes that the producers create a single candidate per event, thus the sorter only assign the priority between the different tags.

    class ZGammaTagSorter : public EDProducer
    {

    public:
        ZGammaTagSorter( const ParameterSet & );
    private:
        void produce( Event &, const EventSetup & ) override;

        std::vector<edm::EDGetTokenT<View<flashgg::DiPhotonTagBase> > > TagList_;
        std::vector<TagPriorityRange> TagPriorityRanges;

        double minObjectWeightException;
        double maxObjectWeightException;
        double minObjectWeightWarning;
        double maxObjectWeightWarning;

        bool debug_;

        bool createNoTag_;

//        string tagName(DiPhotonTagBase::tag_t) const;
    };

    ZGammaTagSorter::ZGammaTagSorter( const ParameterSet &iConfig )
    {
        minObjectWeightException = iConfig.getParameter<double>( "MinObjectWeightException" );
        maxObjectWeightException = iConfig.getParameter<double>( "MaxObjectWeightException" );
        minObjectWeightWarning = iConfig.getParameter<double>( "MinObjectWeightWarning" );
        maxObjectWeightWarning = iConfig.getParameter<double>( "MaxObjectWeightWarning" );


        debug_ = iConfig.getUntrackedParameter<bool>( "Debug", false );

        createNoTag_ = iConfig.getParameter<bool>("CreateNoTag");

        const auto &vpset = iConfig.getParameterSetVector( "TagPriorityRanges" );

        vector<string> labels;

        for( const auto &pset : vpset )
        {
            InputTag tag = pset.getParameter<InputTag>( "TagName" );
            int c1 = pset.getUntrackedParameter<int>( "MinCategory", -999 );
            int c2 = pset.getUntrackedParameter<int>( "MaxCategory", 999 );
            unsigned int i = 0;
            for( ; i < labels.size() ; i++ )
            {
                if( labels[i] == tag.label() ) { break; }
            }
            if( i == TagList_.size() )
            {
                labels.push_back( tag.label() );
                TagList_.push_back( consumes<View<flashgg::DiPhotonTagBase> >( tag ) );
            }
            TagPriorityRanges.emplace_back( tag.label(), c1, c2, i );
        }

        produces<edm::OwnVector<flashgg::DiPhotonTagBase> >();
        produces<edm::OwnVector<flashgg::TagTruthBase> >();
    }

    void ZGammaTagSorter::produce( Event &evt, const EventSetup & )
    {
        unique_ptr<edm::OwnVector<flashgg::DiPhotonTagBase> > SelectedTag( new edm::OwnVector<flashgg::DiPhotonTagBase> );
        unique_ptr<edm::OwnVector<flashgg::TagTruthBase> > SelectedTagTruth( new edm::OwnVector<flashgg::TagTruthBase> );

        for( auto tpr = TagPriorityRanges.begin() ; tpr != TagPriorityRanges.end() ; tpr++ )
        {
            Handle<View<flashgg::DiPhotonTagBase> > TagVectorEntry;
            evt.getByToken( TagList_[tpr->collIndex], TagVectorEntry );

            edm::RefProd<edm::OwnVector<TagTruthBase> > rTagTruth = evt.getRefBeforePut<edm::OwnVector<TagTruthBase> >();

            assert(TagVectorEntry->size()<=1); //It should be one or zero, i.e. max one event per tag

            if(TagVectorEntry->size()==0) continue;

            float centralObjectWeight = TagVectorEntry->ptrAt( 0 )->centralWeight();
            if (centralObjectWeight < minObjectWeightException || centralObjectWeight > maxObjectWeightException)
            {
                throw cms::Exception( "TagObjectWeight" ) << " Tag centralWeight=" << centralObjectWeight << " outside of bound ["
                                                          << minObjectWeightException << "," << maxObjectWeightException
                                                          << "] - " << tpr->name << " - change bounds or debug tag";
            }
            if (centralObjectWeight < minObjectWeightWarning || centralObjectWeight > maxObjectWeightWarning)
            {
                std::cout << "WARNING Tag centralWeight=" << centralObjectWeight << " outside of bound ["
                          << minObjectWeightWarning << "," << maxObjectWeightWarning
                          << "] - " << tpr->name << " - consider investigating!" << std::endl;
            }

            SelectedTag->push_back( *TagVectorEntry->ptrAt( 0 ) );
            edm::Ptr<TagTruthBase> truth = TagVectorEntry->ptrAt( 0 )->tagTruth();
            if( truth.isNonnull() )
            {
                SelectedTagTruth->push_back( *truth );
                SelectedTag->back().setTagTruth( edm::refToPtr( edm::Ref<edm::OwnVector<TagTruthBase> >( rTagTruth, 0 ) ) ); // Normally this 0 would be the index number
            }
        } 


        assert( SelectedTag->size() == 1 || SelectedTag->size() == 0 );

        if (createNoTag_ && SelectedTag->size() == 0)
        {
            SelectedTag->push_back(NoTag());
            edm::RefProd<edm::OwnVector<TagTruthBase> > rTagTruth = evt.getRefBeforePut<edm::OwnVector<TagTruthBase> >();
            TagTruthBase truth_obj;

            SelectedTagTruth->push_back(truth_obj);
            SelectedTag->back().setTagTruth( edm::refToPtr( edm::Ref<edm::OwnVector<TagTruthBase> >( rTagTruth, 0 ) ) );

        }

        evt.put( std::move( SelectedTag ) );
        evt.put( std::move( SelectedTagTruth ) );
    }

/*    string TagSorter::tagName(DiPhotonTagBase::tag_t tagEnumVal) const
    {
        switch(tagEnumVal)
        {
        case DiPhotonTagBase::tag_t::kUndefined:
            return string("UNDEFINED");
        case DiPhotonTagBase::tag_t::kUntagged: 
            return string("Untagged");
        case DiPhotonTagBase::tag_t::kVBF:
            return string("VBF");
        case DiPhotonTagBase::tag_t::kTTHHadronic:
            return string("TTHHadronic");
        case DiPhotonTagBase::tag_t::kTTHLeptonic:
            return string("TTHLeptonic");
        case DiPhotonTagBase::tag_t::kTTHDiLepton:
            return string("TTHDiLepton");
        case DiPhotonTagBase::tag_t::kVHTight:
            return string("VHTight");
        case DiPhotonTagBase::tag_t::kVHLoose:
            return string("VHLoose");
        case DiPhotonTagBase::tag_t::kVHHadronic:
            return string("VHHadronic");
        case DiPhotonTagBase::tag_t::kVHEt:
            return string("VHEt");
        case DiPhotonTagBase::tag_t::kZHLeptonic:
            return string("ZHLeptonic");
        case DiPhotonTagBase::tag_t::kWHLeptonic:
            return string("WHLeptonic");
        case DiPhotonTagBase::tag_t::kVHLeptonicLoose:
            return string("VHLeptonicLoose");
        case DiPhotonTagBase::tag_t::kVHMet:
            return string("VHMet");
        }
        return string("TAG NOT ON LIST");
    }
*/

}

typedef flashgg::ZGammaTagSorter FlashggZGammaTagSorter;
DEFINE_FWK_MODULE( FlashggZGammaTagSorter );
// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4

