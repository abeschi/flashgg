#include "flashgg/DataFormats/interface/TTHTagTruth.h"
#include <iostream>

using namespace flashgg;

TTHTagTruth::TTHTagTruth() {}

TTHTagTruth::~TTHTagTruth() {}

TTHTagTruth *TTHTagTruth::clone() const
{
    TTHTagTruth *result = new TTHTagTruth;
    result->setH( H() );
    result->setLeadPhoton( leadPhoton() );
    result->setSubleadPhoton( subleadPhoton() );
    result->setT( t() );
    result->setB( b() );
    result->setWplus1( Wplus1() );
    result->setWplus2( Wplus2() );
    result->setTbar( tbar() );
    result->setBbar( bbar() );
    result->setWminus1( Wminus1() );
    result->setWminus2( Wminus2() );
    
    result->copyBaseInfo( *this );
    
    return result;
}

// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
