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

/** Digitization for the SDHcal - based on NewLDCCaloDigi. 
 *
 *  @author  G.Grenier, IPNL
 *  @author  R.Han, IPNL
 *	@author  A.Steen, IPNL
 *	@author  G.Garillot, IPNL
 *	@author  B.Li, IPNL
 *  @version $Id$
 */


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
		virtual Processor* newProcessor() { return new SimDigital ; }
		SimDigital() ;

		/** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
		virtual void init() ;

		/** Called for every event - the working horse.
   */
		virtual void processEvent( LCEvent * evt ) ;

		SimDigital(const SimDigital &toCopy) = delete ;
		void operator=(const SimDigital &toCopy) = delete ;


	private :
		//intermediate storage class
		struct hitMemory
		{
				hitMemory() {;}

				CalorimeterHitImpl* ahit = nullptr ;
				std::set<int> relatedHits {} ;
				float maxEnergydueToHit = -1 ;
				int rawHit = -1 ;

				hitMemory(const hitMemory& other) = delete ;
				hitMemory& operator=(const hitMemory& other) = delete ;
		} ;

		typedef std::map<dd4hep::long64, hitMemory> cellIDHitMap ;


		void processCollection(LCCollection* inputCol , LCCollectionVec*& outputCol , LCCollectionVec*& outputRelCol , CHT::Layout layout, LCFlagImpl& flag) ;
		void createPotentialOutputHits(cellIDHitMap& myHitMap , LCCollection* col , SimDigitalGeomCellId* aGeomCellId) ;

		void removeAdjacentStep(std::vector<StepAndCharge>& vec) ;
		void fillTupleStep(std::vector<StepAndCharge>& vec , int level) ;
		void removeHitsBelowThreshold(cellIDHitMap& myHitMap , float threshold) ;
		void applyThresholds(cellIDHitMap& myHitMap) ;

		std::vector<std::string> _inputCollections {} ;
		std::vector<std::string> _outputCollections {} ;
		std::vector<std::string> _outputRelCollections {} ;

		std::map<std::string, int> _counters {} ;
		std::vector<float> _thresholdHcal {} ;

		//charge spreader
		std::string chargeSpreaderOption = "Uniform" ;
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
		float _minXYdistanceBetweenStep  = 0.5f ;
		bool _keepAtLeastOneStep = true ;
		AIDA::ITuple* _debugTupleStepFilter = nullptr ;
		AIDA::ITuple* _tupleStepFilter = nullptr ;
		AIDA::ITuple* _tupleCollection = nullptr ;

		std::string _encodingType  = "LCGEO" ;
		std::string _hcalOption = "VIDEAU" ;
} ;

#endif
