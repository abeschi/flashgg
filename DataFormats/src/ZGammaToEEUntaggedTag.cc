#include "flashgg/DataFormats/interface/ZGammaToEEUntaggedTag.h"

using namespace flashgg;

ZGammaToEEUntaggedTag::ZGammaToEEUntaggedTag() {}

ZGammaToEEUntaggedTag::~ZGammaToEEUntaggedTag() {}

ZGammaToEEUntaggedTag::ZGammaToEEUntaggedTag( edm::Ptr<flashgg::EleEleGammaCandidate> EleEleG)
{
	EEG_=   *EleEleG;
}

ZGammaToEEUntaggedTag::ZGammaToEEUntaggedTag( flashgg::EleEleGammaCandidate EleEleG)
{
	EEG_=   EleEleG;
}

