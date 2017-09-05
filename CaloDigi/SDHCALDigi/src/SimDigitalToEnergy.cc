#include "SimDigitalToEnergy.h"

#include <limits>
#include <marlin/Exceptions.h>
#include <EVENT/CalorimeterHit.h>
#include <IMPL/CalorimeterHitImpl.h>
#include <IMPL/LCRelationImpl.h>
#include <IMPL/LCFlagImpl.h>

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


	std::vector<std::string> hcalCollections = { "HCALBarrelDigi" , "HCALEndcapDigi" , "HCALOtherDigi" } ;
	registerInputCollections( LCIO::SIMCALORIMETERHIT,
							  "HCALInputCollections" ,
							  "Digi CalorimeterHit Collections" ,
							  _hcalCollections ,
							  hcalCollections) ;

	_outputHcalCollections = { "HCALBarrelRec" , "HCALEndcapRec" , "HCALOtherRec" } ;

	registerOutputCollection( LCIO::CALORIMETERHIT,
							  "HCALOutputCollections" ,
							  "HCAL Collection of real Hits" ,
							  _outputHcalCollections[0] ,
			std::string("HCALBarrelRec")  ) ;

	registerOutputCollection( LCIO::CALORIMETERHIT,
							  "HCALOutputCollection1" ,
							  "HCAL Collection of real Hits" ,
							  _outputHcalCollections[1] ,
			std::string("HCALEndcapRec") );

	registerOutputCollection( LCIO::CALORIMETERHIT,
							  "HCALOutputCollection2" ,
							  "HCAL Collection of real Hits" ,
							  _outputHcalCollections[2] ,
			std::string("HCALOtherRec") ) ;

	registerOutputCollection( LCIO::LCRELATION,
							  "RelationOutputCollection" ,
							  "CaloHit Relation Collection" ,
							  _outputRelCollection ,
							  std::string("RelationCaloHit")) ;


	std::vector<float> energyCoefficients = {0.4f} ;
	registerProcessorParameter("HCALThreshold" ,
							   "Threshold for HCAL Hits in pC" ,
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

	for ( unsigned int i = 0 ; i < _hcalCollections.size() ; ++i )
	{
		try
		{
			const std::string& colName = _hcalCollections.at(i) ;
			LCCollection* digiCol = evt->getCollection( colName.c_str() ) ;

			LCCollectionVec* recCol = processCollection(digiCol) ;

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
	CellIDEncoder<CalorimeterHitImpl> encoder( col->getParameters().getStringVal(LCIO::CellIDEncoding) , recCol ) ;

	LCFlagImpl flag = col->getFlag() ;
	recCol->setFlag( flag.getFlag() ) ;


	unsigned long nThresholds = _energyCoefficients.size() ;

	int nHit = col->getNumberOfElements() ;
	for ( int i = 0 ; i < nHit ; ++i )
	{
		CalorimeterHit* digiHit = dynamic_cast<CalorimeterHit*>( col->getElementAt( i ) ) ;

		dd4hep::long64 cellIDvalue = decoder( digiHit ).getValue() ;
		encoder.setValue(cellIDvalue) ;

		CalorimeterHitImpl* recHit = new CalorimeterHitImpl ;

		float threshold = digiHit->getEnergy() - 1.0f ;

		if ( threshold < 0 )
		{
			delete recHit ;
			continue ;
		}

		unsigned int iCoeff = static_cast<unsigned int>( threshold ) ;
		if ( iCoeff > nThresholds )
			recHit->setEnergy( _energyCoefficients.at( nThresholds-1 ) ) ;
		else
			recHit->setEnergy( _energyCoefficients.at( iCoeff ) ) ;

//		recHit->setCellID0( digiHit->getCellID0() ) ;
//		recHit->setCellID1( digiHit->getCellID1() ) ;
		encoder.setCellID( recHit ) ;
		recHit->setPosition( digiHit->getPosition() ) ;
		recHit->setTime( digiHit->getTime() ) ;
		recHit->setRawHit( digiHit ) ;

		recCol->addElement(recHit) ;

		LCRelationImpl* rel = new LCRelationImpl(recHit , digiHit , 1.0) ;
		_relcol->addElement( rel ) ;
	}

	return recCol ;
}

void SimDigitalToEnergy::end()
{

}

