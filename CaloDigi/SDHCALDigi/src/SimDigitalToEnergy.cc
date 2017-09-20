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
	: Processor("SimDigitalToEnergy")
{
	_description = "This processor transforms the threshold value of digitized SDHCAL CalorimeterHits into energy" ;


	std::vector<std::string> inputHcalCollections = { "HCALBarrelDigi" , "HCALEndcapDigi" , "HCALOtherDigi" } ;
	registerInputCollections( LCIO::CALORIMETERHIT,
							  "inputHitCollections" ,
							  "Digi CalorimeterHit Collections" ,
							  _inputHcalCollections ,
							  inputHcalCollections) ;


	std::vector <std::string> outputHcalCollections = { "HCALBarrelRec" , "HCALEndcapRec" , "HCALOtherRec" } ;
	registerProcessorParameter( "outputHitCollections",
								"output hit collection names",
								_outputHcalCollections,
								outputHcalCollections) ;


	registerInputCollection( LCIO::LCRELATION,
							 "inputRelationCollection",
							 "input relation collection name",
							 _inputRelCollection,
							 std::string("HCALDigiRelationSim") ) ;


	registerOutputCollection( LCIO::LCRELATION,
							  "outputRelationCollection" ,
							  "CaloHit Relation Collection" ,
							  _outputRelCollection ,
							  std::string("HCALRecRelationSim") ) ;


	std::vector<float> energyCoefficients = {0.4f} ;
	registerProcessorParameter("EnergyCalibration" ,
							   "Threshold to Energy correspondace" ,
							   _energyCoefficients,
							   energyCoefficients) ;

	if ( _energyCoefficients.size() == 0 )
		throw ParseException( "ERROR : No energy Coefficient provided" ) ;
}

void SimDigitalToEnergy::init()
{

}

void SimDigitalToEnergy::processEvent( LCEvent* evt )
{
	// create the output collections
	_relcol = new LCCollectionVec(LCIO::LCRELATION) ;

	for ( unsigned int i = 0 ; i < _inputHcalCollections.size() ; ++i )
	{
		try
		{
			const std::string& colName = _inputHcalCollections.at(i) ;
			LCCollection* digiCol = evt->getCollection( colName.c_str() ) ;
			std::string initString = digiCol->getParameters().getStringVal(LCIO::CellIDEncoding) ;

			inputRelCol = evt->getCollection( _inputRelCollection.c_str() ) ;

			LCCollectionVec* recCol = processCollection(digiCol) ;

			recCol->parameters().setValue(LCIO::CellIDEncoding , initString) ;
			evt->addCollection( recCol , _outputHcalCollections.at(i).c_str() ) ;
		}
		catch ( DataNotAvailableException& )
		{
		}
	}

	evt->addCollection( _relcol , _outputRelCollection.c_str() ) ;
}

LCCollectionVec* SimDigitalToEnergy::processCollection(LCCollection* col)
{
	LCCollectionVec* recCol = new LCCollectionVec(LCIO::CALORIMETERHIT) ;

	CellIDDecoder<CalorimeterHit> decoder(col) ;

	LCFlagImpl flag = col->getFlag() ;
	recCol->setFlag( flag.getFlag() ) ;

	LCRelationNavigator navi(inputRelCol) ;


	unsigned long nThresholds = _energyCoefficients.size() ;

	int nHit = col->getNumberOfElements() ;
	for ( int i = 0 ; i < nHit ; ++i )
	{
		CalorimeterHit* digiHit = dynamic_cast<CalorimeterHit*>( col->getElementAt( i ) ) ;

		//		dd4hep::long64 cellIDvalue = decoder( digiHit ).getValue() ;
		//		encoder.setValue(cellIDvalue) ;

		float threshold = digiHit->getEnergy() - 1.0f ;

		if ( threshold < 0 )
			continue ;

		CalorimeterHitImpl* recHit = new CalorimeterHitImpl ;

		unsigned int iCoeff = static_cast<unsigned int>( threshold ) ;
		if ( iCoeff >= nThresholds )
			recHit->setEnergy( _energyCoefficients.at( nThresholds-1 ) ) ;
		else
			recHit->setEnergy( _energyCoefficients.at( iCoeff ) ) ;

		recHit->setCellID0( digiHit->getCellID0() ) ;
		recHit->setCellID1( digiHit->getCellID1() ) ;
		//		encoder.setCellID( recHit ) ;
		recHit->setPosition( digiHit->getPosition() ) ;
		recHit->setTime( digiHit->getTime() ) ;
		recHit->setRawHit( digiHit->getRawHit() ) ;

		recHit->setType( digiHit->getType() ) ;

		recCol->addElement(recHit) ;



		for ( unsigned int iRel = 0 ; iRel < navi.getRelatedFromObjects( digiHit ).size() ; ++iRel )
		{
			SimCalorimeterHit* simHit = dynamic_cast<SimCalorimeterHit*>( navi.getRelatedFromObjects(digiHit).at(iRel) ) ;
			float weight = navi.getRelatedFromWeights(digiHit).at(iRel) ;
			_relcol->addElement( new LCRelationImpl(simHit , recHit , weight) ) ;
		}
	}

	return recCol ;
}

void SimDigitalToEnergy::end()
{

}

