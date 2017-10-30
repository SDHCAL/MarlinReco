#ifndef SimDigitalGeom_h
#define SimDigitalGeom_h

#include <marlin/Processor.h>
#include <IMPL/LCCollectionVec.h>

#include <UTIL/CellIDEncoder.h>
#include <UTIL/CellIDDecoder.h>

#include <EVENT/SimCalorimeterHit.h>
#include <EVENT/LCGenericObject.h>
#include <IMPL/CalorimeterHitImpl.h>

#include "CalorimeterHitType.h" //in MarlinUtil
#include "marlinutil/LCGeometryTypes.h"

#include <gear/GearParameters.h>
#include <gear/CalorimeterParameters.h>
#include <gear/LayerLayout.h>

#include "DD4hep/Factories.h"
#include "DD4hep/DD4hepUnits.h"
#include "DDRec/DetectorData.h"
#include "DD4hep/DetType.h"
#include "DDRec/DetectorData.h"
#include "DDRec/DDGear.h"
#include "DDRec/MaterialManager.h"
#include "DDRec/API/Calorimeter.h"
#include "DDRec/DetectorSurfaces.h"

namespace AIDA
{
class ITuple ;
}

struct StepAndCharge ;

class SimDigitalGeomCellId
{
	public :
		SimDigitalGeomCellId(LCCollection* inputCol, LCCollectionVec* outputCol) ;
		virtual ~SimDigitalGeomCellId() ;

		virtual float getCellSize() = 0 ;
		virtual void setLayerLayout(CHT::Layout layout) = 0 ;

		std::vector<StepAndCharge> decode(SimCalorimeterHit* hit) ;
		std::vector<StepAndCharge> decode(SimCalorimeterHit* hit , const std::map<dd4hep::long64, std::vector<LCGenericObject*> >& map) ;

	protected :
		virtual void processGeometry(SimCalorimeterHit* hit) = 0 ;
		void createStepAndChargeVec(SimCalorimeterHit* hit , std::vector<StepAndCharge>& vec) ;
		void createStepAndChargeVec(SimCalorimeterHit* hit , const std::vector<LCGenericObject*>& genericVec , std::vector<StepAndCharge>& vec) ;

	public :
		virtual void encode(CalorimeterHitImpl *hit , int delta_I , int delta_J) = 0 ;


		int I() const { return _Iy ; }
		int J() const { return _Jz ; }
		int K() const { return _trueLayer ; }
		int stave() const { return _stave ; }
		int module() const { return _module ; }
		int tower() const { return _tower ; }

		const LCVector3D& normalToRPCPlane() const { return _normal ; }
		const LCVector3D& Iaxis() const { return _Iaxis ; }
		const LCVector3D& Jaxis() const { return _Jaxis ; }



		SimDigitalGeomCellId(const SimDigitalGeomCellId &toCopy) = delete ;
		void operator=(const SimDigitalGeomCellId &toCopy) = delete ;

	protected :

		CHT::Layout _currentHCALCollectionCaloLayout = CHT::any ;

		dd4hep::long64 _cellIDvalue = 0 ;
		CellIDDecoder<SimCalorimeterHit> _decoder ;
		CellIDEncoder<CalorimeterHitImpl> _encoder ;

		int _trueLayer = -999 ;
		int _stave = -999 ;
		int _module = -999 ;
		int _tower = -999 ;
		int _Iy = -999 ;
		int _Jz = -999 ;

		LCVector3D _normal ;
		LCVector3D _Iaxis ;
		LCVector3D _Jaxis ;

		const float* _hitPosition = nullptr ;

		std::string _cellIDEncodingString = "" ;



		//geometry debug tuples
	public :
		static void bookTuples(const marlin::Processor* proc) ;
	protected :
		void fillDebugTupleGeometryHit() ;
		void fillDebugTupleGeometryStep(SimCalorimeterHit* hit , const std::vector<StepAndCharge>& stepsInIJZcoord) ;

		static AIDA::ITuple* _tupleHit ;
		enum {TH_CHTLAYOUT,TH_MODULE,TH_TOWER,TH_STAVE,TH_LAYER,TH_I,TH_J,
			  TH_X,TH_Y,TH_Z,
			  TH_NORMALX,TH_NORMALY,TH_NORMALZ,
			  TH_IX,TH_IY,TH_IZ,
			  TH_JX,TH_JY,TH_JZ} ;
		static AIDA::ITuple* _tupleStep ;
		enum {TS_CHTLAYOUT,TS_HITCELLID,TS_NSTEP,
			  TS_HITX,TS_HITY,TS_HITZ,
			  TS_STEPX,TS_STEPY,TS_STEPZ,
			  TS_DELTAI,TS_DELTAJ,TS_DELTALAYER,TS_TIME} ;
} ;

class SimDigitalGeomCellIdLCGEO : public SimDigitalGeomCellId
{
	public :
		SimDigitalGeomCellIdLCGEO(LCCollection* inputCol, LCCollectionVec* outputCol) ;
		virtual ~SimDigitalGeomCellIdLCGEO() ;

		void setCellSize(float size) { _cellSize = size ; }
		virtual float getCellSize() ;
		virtual void setLayerLayout(CHT::Layout layout) ;

		virtual void encode(CalorimeterHitImpl *hit , int delta_I , int delta_J) ;

		SimDigitalGeomCellIdLCGEO(const SimDigitalGeomCellIdLCGEO &toCopy) = delete ;
		void operator=(const SimDigitalGeomCellIdLCGEO &toCopy) = delete ;

	protected :
		virtual void processGeometry(SimCalorimeterHit* hit) ;

		std::vector<std::string> _encodingString = { "layer", "stave", "module", "tower", "x", "y" } ;

		float _cellSize = 0.0f ;

		dd4hep::rec::LayeredCalorimeterData* _caloData = nullptr ;

		//				dd4hep::DetElement theDetector;

} ;

class SimDigitalGeomRPCFrame ;
class SimDigitalGeomCellIdMOKKA : public SimDigitalGeomCellId
{
	public :
		SimDigitalGeomCellIdMOKKA(LCCollection* inputCol, LCCollectionVec* outputCol) ;
		virtual ~SimDigitalGeomCellIdMOKKA() ;


		virtual float getCellSize() ;
		virtual void setLayerLayout(CHT::Layout layout) ;

		enum HCAL_GEOM {VIDEAU,TESLA} ;
		void setGeom(HCAL_GEOM geom) { _geom = geom ; }

		HCAL_GEOM getGeom() const { return _geom ; }


		virtual void encode(CalorimeterHitImpl *hit , int delta_I , int delta_J) ;

		SimDigitalGeomCellIdMOKKA(const SimDigitalGeomCellIdMOKKA &toCopy) = delete ;
		void operator=(const SimDigitalGeomCellIdMOKKA &toCopy) = delete ;

	protected :
		virtual void processGeometry(SimCalorimeterHit* hit) ;

		std::vector<std::string> _encodingString = { "K-1", "S-1", "M", "", "I", "J" } ;


		HCAL_GEOM _geom = TESLA ;

		SimDigitalGeomRPCFrame* _normal_I_J_setter = nullptr ;
		const gear::LayerLayout* _layerLayout = nullptr ;


		friend class SimDigitalGeomRPCFrame ;
} ;


//hierarchy of classes to determine the RPC reference frame
class SimDigitalGeomRPCFrame
{
	public:
		SimDigitalGeomRPCFrame(SimDigitalGeomCellIdMOKKA& h) : _layerInfo(h) {}
		virtual ~SimDigitalGeomRPCFrame() ;
		virtual void setRPCFrame() = 0 ;
	private :
		SimDigitalGeomCellIdMOKKA& _layerInfo ;
	protected :
		int stave() const { return _layerInfo._stave ; }
		int module() const { return _layerInfo._module ; }
		LCVector3D& normal() const { return _layerInfo._normal ; }
		LCVector3D& Iaxis() const { return _layerInfo._Iaxis ; }
		LCVector3D& Jaxis() const { return _layerInfo._Jaxis ; }
};

class SimDigitalGeomRPCFrame_TESLA_BARREL : public SimDigitalGeomRPCFrame
{
	public:
		SimDigitalGeomRPCFrame_TESLA_BARREL(SimDigitalGeomCellIdMOKKA& h) : SimDigitalGeomRPCFrame(h) {}
		void setRPCFrame();
};
class SimDigitalGeomRPCFrame_VIDEAU_BARREL : public SimDigitalGeomRPCFrame
{
	public:
		SimDigitalGeomRPCFrame_VIDEAU_BARREL(SimDigitalGeomCellIdMOKKA& h) : SimDigitalGeomRPCFrame(h) {}
		void setRPCFrame();
};
class SimDigitalGeomRPCFrame_TESLA_ENDCAP : public SimDigitalGeomRPCFrame
{
	public:
		SimDigitalGeomRPCFrame_TESLA_ENDCAP(SimDigitalGeomCellIdMOKKA& h) : SimDigitalGeomRPCFrame(h) {}
		void setRPCFrame();
};
class SimDigitalGeomRPCFrame_VIDEAU_ENDCAP : public SimDigitalGeomRPCFrame
{
	public:
		SimDigitalGeomRPCFrame_VIDEAU_ENDCAP(SimDigitalGeomCellIdMOKKA& h) : SimDigitalGeomRPCFrame(h) {}
		void setRPCFrame();
};

#endif //SimDigitalGeom_h
