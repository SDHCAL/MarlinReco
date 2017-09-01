#ifndef SimDigital_HHH
#define SimDigital_HHH

#include "marlin/Processor.h"
#include "lcio.h"
#include <vector>
#include <string>
#include <map>
#include <set>
#include <utility>
#include <stdlib.h>
#include <EVENT/LCCollection.h>
#include <EVENT/SimCalorimeterHit.h>
#include <IMPL/CalorimeterHitImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCFlagImpl.h>
#include <UTIL/CellIDEncoder.h>
#include <UTIL/CellIDDecoder.h>
#include <marlin/Global.h>
#include <gear/GearParameters.h>
#include <gear/CalorimeterParameters.h>
#include <gear/LayerLayout.h>
#include <TROOT.h>
#include <TFile.h>
#include <TF2.h>
#include <TH1.h>
#include "TH1F.h"

#include "CalorimeterHitType.h" //in MarlinUtil
#include "marlinutil/LCGeometryTypes.h"

#include "DD4hep/Factories.h"
#include "DD4hep/DD4hepUnits.h"
#include "DDRec/DetectorData.h"
#include "DD4hep/DetType.h"
#include "DDRec/DetectorData.h"
#include "DDRec/DDGear.h"
#include "DDRec/MaterialManager.h"
#include "DDRec/API/Calorimeter.h"
#include "DDRec/DetectorSurfaces.h"

#include "SimDigitalGeom.h"
#include "ChargeSpreader.h"
#include "ChargeInducer.h"
#include "EfficiencyManager.h"

class TH1F;
class TF1;
class TTree;

namespace AIDA
{
class ITuple ;
}

using namespace lcio ;
using namespace marlin ;


//TK   added TH_TOWER to root tuples  Nov 2016
//
// Gerald Grenier :
// start removing the ECAL part whih has no use for SDHCAL
// this can be fully activated when there is an ECAL digitizer 
// that don't do HCAL digitization

/** Digitization for the SDHcal - based on NewLDCCaloDigi. 
 *
 *  @author  G.Grenier, INPL
 *  @author  R.Han, INPL
 *  @version $Id$
 */



class SimDigitalGeomRPCFrame ;

struct StepAndCharge
{
		StepAndCharge()
			: step()
		{}
		StepAndCharge(LCVector3D vec)
			: step(vec)
		{}
		LCVector3D step ;
		double charge = 0 ;
} ;

struct AsicKey
{
		AsicKey(int l , int aI = -1 , int aJ = -1) : layerID(l) , asicI(aI) , asicJ(aJ) {}
		int layerID ;
		int asicI ;
		int asicJ ;

		bool operator<(const AsicKey& b) const
		{
			if ( this->layerID != b.layerID )
				return this->layerID < b.layerID ;
			else if ( this->asicI != b.asicI )
				return this->asicI < b.asicI ;
			else
				return this->asicJ < b.asicJ ;
		}
		bool operator==(const AsicKey& b) const
		{
			return ( this->layerID == b.layerID ) && ( this->asicI == b.asicI ) && ( this->asicJ == b.asicJ ) ;
		}

} ;




class SimDigital : public Processor
{
	public:
		virtual Processor*  newProcessor() { return new SimDigital;}
		SimDigital() ;

		/** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
		virtual void init() ;

		/** Called for every run.
   */
		virtual void processRunHeader( LCRunHeader* run ) ;

		/** Called for every event - the working horse.
   */
		virtual void processEvent( LCEvent * evt ) ;


		/** Called after data processing for clean up.
   */
		virtual void end() ;


		SimDigital(const SimDigital &toCopy) = delete ;
		void operator=(const SimDigital &toCopy) = delete ;


	private :

		std::vector<std::string> _hcalCollections ;
		std::vector<std::string> _outputHcalCollections ;
		std::map<std::string, int> _counters ;
		std::vector<float> _thresholdHcal ;
		std::vector<float> _calibrCoeffHcal ;
		std::string _outputRelCollection = "" ;


		LCCollectionVec* _relcol = nullptr ;

		void processHCAL(LCEvent* evt, LCFlagImpl& flag) ;

		static bool sortStepWithCharge(StepAndCharge s1, StepAndCharge s2) {return s1.charge>s2.charge;}

		//intermediate storage class
		struct hitMemory
		{
				hitMemory()
					: relatedHits()
				{}
				CalorimeterHitImpl* ahit = nullptr ;
				std::set<int> relatedHits ;
				float maxEnergydueToHit = -1 ;
				int rawHit = -1 ;

				hitMemory(const hitMemory& other) = delete ;
				hitMemory& operator=(const hitMemory& other)
				{
					if (this != &other)
						*this = other ;

					return *this ;
				}
		} ;

		typedef std::map<dd4hep::long64, hitMemory> cellIDHitMap ;

		float depositedEnergyInRPC = 0.0f ;

		//charge spreader
		std::string chargeSpreaderOption = "Uniform ";
		ChargeSpreaderParameters chargeSpreaderParameters ;
		ChargeSpreader* chargeSpreader = nullptr ;


		std::string polyaOption = "Uniform" ;
		double polyaQbar = 0 ;
		double polyaTheta = 0 ;
		ChargeInducer* chargeInducer = nullptr ;
		int _polyaRandomSeed = 1 ;


		bool _doThresholds = true ;

		std::string _effMapFileName = "" ;
		float _constEffMapValue = 0.97f ;
		std::string efficiencyOption  = "Uniform" ;
		EfficiencyManager* efficiency = nullptr ;


		float _absZstepFilter  = 0.0005f ;
		bool _keepAtLeastOneStep = true ;
		float _minXYdistanceBetweenStep  = 0.5f ;
		AIDA::ITuple* _debugTupleStepFilter = nullptr ;
		AIDA::ITuple* _tupleStepFilter = nullptr ;
		AIDA::ITuple* _tupleCollection = nullptr ;


		//helper function to remove steps too close in I,J
		void remove_adjacent_step(std::vector<StepAndCharge>& vec);
		void fillTupleStep(std::vector<StepAndCharge>& vec,int level);

		LCCollectionVec* processHCALCollection(LCCollection * col ,CHT::Layout layout, LCFlagImpl& flag) ;
		void createPotentialOutputHits(cellIDHitMap& myHitMap, LCCollection* col, SimDigitalGeomCellId& aGeomCellId) ;
		void removeHitsBelowThreshold(cellIDHitMap& myHitMap, float threshold) ;
		void applyThresholds(cellIDHitMap& myHitMap) ;

		std::string _encodingType  = "LCGEO" ;
		std::string _hcalOption = "VIDEAU" ;

} ;

#endif
