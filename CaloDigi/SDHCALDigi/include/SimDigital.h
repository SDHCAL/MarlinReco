#ifndef SimDigital_HHH
#define SimDigital_HHH

#include "marlin/Processor.h"
#include "lcio.h"
#include <vector>
#include <string>
#include <map>
#include <set>
#include <utility>
#include <limits>
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

#include <AIDA/IHistogramFactory.h>
#include <AIDA/IHistogram1D.h>

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
		StepAndCharge(LCVector3D vec , float _time)
			: step(vec) , time(_time)
		{}
		LCVector3D step ;
		float charge = 0 ;
		float stepLength = 0 ;
		float time = 0 ;
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
				hitMemory()
					: ahit(nullptr) , relatedHits() , maxEnergydueToHit(-1) , rawHit(-1)
				{}

				CalorimeterHitImpl* ahit = nullptr ;
				std::set<int> relatedHits {} ;
				float maxEnergydueToHit = -1 ;
				int rawHit = -1 ;

				hitMemory(const hitMemory& other) = delete ;
				hitMemory& operator=(const hitMemory& other) = delete ;
		} ;

		typedef std::map<dd4hep::long64, hitMemory> cellIDHitMap ;

		void createGenericObjects(LCCollection* col) ;

		void processCollection(LCCollection* inputCol , LCCollectionVec*& outputCol , LCCollectionVec*& outputRelCol , CHT::Layout layout, LCFlagImpl& flag) ;
		cellIDHitMap createPotentialOutputHits(LCCollection* col , SimDigitalGeomCellId* aGeomCellId) ;

		void removeAdjacentStep(std::vector<StepAndCharge>& vec) ;
		void fillTupleStep(std::vector<StepAndCharge>& vec , int level) ;
		void removeHitsBelowThreshold(cellIDHitMap& myHitMap , float threshold) ;
		void applyThresholds(cellIDHitMap& myHitMap) ;

		std::vector<std::string> _inputCollections {} ;
		std::vector<std::string> _inputGenericCollections {} ;

		std::vector<std::string> _outputCollections {} ;
		std::vector<std::string> _outputRelCollections {} ;

		std::map<std::string, int> _counters {} ;
		std::vector<float> _thresholdHcal {} ;

		std::vector<double> _hitCharge = {} ;

		std::map<dd4hep::long64 , std::vector<LCGenericObject*>> geneMap = {} ;

		float _cellSize = 0 ;
		float _gasGapWidth = 1.2f ;

		//charge spreader
		std::string chargeSpreaderOption = "Uniform" ;
		std::string spreaderMapFile = "" ;
		ChargeSpreaderParameters chargeSpreaderParameters ;
		ChargeSpreader* chargeSpreader = nullptr ;

		std::string polyaOption = "Uniform" ;
		std::string polyaMapFile = "" ;
		float polyaQbar = 0.0f ;
		float polyaTheta = 0.0f ;
		ChargeInducer* chargeInducer = nullptr ;
		int _polyaRandomSeed = 1 ;

		double timeCut = std::numeric_limits<double>::max() ;
		double stepLengthCut = -1.0 ;
		float _angleCorrPow = 0.4f ;

		bool _doThresholds = true ;

		std::string efficiencyOption  = "Uniform" ;
		std::string effMapFile = "" ;
		float _constEffMapValue = 0.97f ;
		EfficiencyManager* efficiency = nullptr ;

		float _absZstepFilter  = 0.0005f ;
		float _minXYdistanceBetweenStep  = 0.5f ;
		bool _keepAtLeastOneStep = true ;
		AIDA::ITuple* _debugTupleStepFilter = nullptr ;
		AIDA::ITuple* _tupleStepFilter = nullptr ;
		AIDA::ITuple* _tupleCollection = nullptr ;

		AIDA::IHistogram1D* _histoCellCharge = nullptr ;

		std::string _encodingType  = "LCGEO" ;
		std::string _hcalOption = "VIDEAU" ;
} ;

#endif
