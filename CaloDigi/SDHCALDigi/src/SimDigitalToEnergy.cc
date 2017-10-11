#include "SimDigitalToEnergy.h"

#include <limits>
#include <marlin/Exceptions.h>
#include <EVENT/SimCalorimeterHit.h>
#include <EVENT/CalorimeterHit.h>
#include <IMPL/CalorimeterHitImpl.h>
#include <IMPL/LCRelationImpl.h>
#include <IMPL/LCFlagImpl.h>
#include <UTIL/LCRelationNavigator.h>

#include <UTIL/CellIDEncoder.h>
#include <UTIL/CellIDDecoder.h>

#include <DD4hep/Factories.h>

using namespace lcio ;
using namespace marlin ;

SimDigitalToEnergy SimDigitalToEnergy ;

SimDigitalToEnergy::SimDigitalToEnergy()
	: RealisticCaloReco::Processor("SimDigitalToEnergy")
{
	_description = "This processor transforms the threshold value of digitized SDHCAL CalorimeterHits into energy" ;

	std::vector<float> energyCoefficients = {0.4f} ;
	registerProcessorParameter("EnergyCalibration" ,
							   "Threshold to Energy correspondace" ,
							   _energyCoefficients,
							   energyCoefficients) ;
}

void SimDigitalToEnergy::init()
{
	//to avoid crash
	_calibrCoeff.push_back(0) ;
	_calLayers.push_back(0) ;

	RealisticCaloReco::init() ;
	assert( _energyCoefficients.size() > 0 ) ;

	_flag.setBit(LCIO::RCHBIT_ID1) ;
}

float SimDigitalToEnergy::reconstructEnergy(const CalorimeterHit* hit)
{
	float threshold = hit->getEnergy() - 1.0f ;

	if ( threshold < 0 )
		return 0.f ;

	unsigned int iCoeff = static_cast<unsigned int>( threshold ) ;
	if ( iCoeff >= _energyCoefficients.size() )
		return _energyCoefficients.at( _energyCoefficients.size()-1 ) ;
	else
		return _energyCoefficients.at( iCoeff ) ;
}

