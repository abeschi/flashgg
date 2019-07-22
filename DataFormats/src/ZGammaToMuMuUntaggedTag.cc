#include "flashgg/DataFormats/interface/ZGammaToMuMuUntaggedTag.h"

using namespace flashgg;

ZGammaToMuMuUntaggedTag::ZGammaToMuMuUntaggedTag() {}

ZGammaToMuMuUntaggedTag::~ZGammaToMuMuUntaggedTag() {}

ZGammaToMuMuUntaggedTag::ZGammaToMuMuUntaggedTag( edm::Ptr<flashgg::MuMuGammaCandidate> MuMuG)
{
	MMG_ = *MuMuG;
}

ZGammaToMuMuUntaggedTag::ZGammaToMuMuUntaggedTag( flashgg::MuMuGammaCandidate MuMuG)
{
	MMG_ = MuMuG;
}

