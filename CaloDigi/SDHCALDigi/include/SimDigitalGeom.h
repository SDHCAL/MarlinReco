#ifndef SimDigitalGeom_h
#define SimDigitalGeom_h

#include <marlin/Processor.h>
#include <IMPL/LCCollectionVec.h>

#include <UTIL/CellIDEncoder.h>
#include <UTIL/CellIDDecoder.h>

#include <EVENT/SimCalorimeterHit.h>
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

//helper class to manage cellId and local geometry
const int ENCODINGTYPES        = 2 ;
const int ENCODINGSTRINGLENGTH = 6 ;

struct StepAndCharge ;
class SimDigitalGeomRPCFrame ;

class SimDigitalGeomCellId
{

	public:
		static void bookTuples(const marlin::Processor* proc);
		SimDigitalGeomCellId(LCCollection *inputCol, LCCollectionVec *outputCol);
		~SimDigitalGeomCellId();
		//return the list of step positions in coordinates corresponding to 'I' ,'J' and 'layer'
		std::vector<StepAndCharge> decode(SimCalorimeterHit *hit) ;
		void encode(CalorimeterHitImpl *hit,int delta_I, int delta_J);
		void setLayerLayout( CHT::Layout layout);
		static void setEncodingType(std::string type);
		static void setHcalOption(std::string hcalOption);
		float getCellSize();
		const LCVector3D& normalToRPCPlane() {return _normal;}
		const LCVector3D& Iaxis() {return _Iaxis;}
		const LCVector3D& Jaxis() {return _Jaxis;}

		inline int I() const {return _Iy;}
		inline int J() const {return _Jz;}
		inline int K() const {return _trueLayer;}
		inline int stave() const {return _stave;}
		inline int module() const {return _module;}
		inline int tower() const {return _tower;}

		SimDigitalGeomCellId(const SimDigitalGeomCellId &toCopy) = delete ;
		void operator=(const SimDigitalGeomCellId &toCopy) = delete ;

	private :

		enum HCAL_GEOM {VIDEAU,TESLA};
		HCAL_GEOM _geom = TESLA ;
		int _trueLayer = -999 ;
		int _stave = -999 ;
		int _module = -999 ;
		int _tower = -999 ;
		int _Iy = -999 ;
		int _Jz = -999 ;
		dd4hep::long64 _cellIDvalue = 0 ;
		static int _encodingType;
		static std::string _hcalOption;
		const float* _hitPosition = nullptr ;
		CellIDDecoder<SimCalorimeterHit> _decoder;
		CellIDEncoder<CalorimeterHitImpl> _encoder;
		const gear::LayerLayout* _layerLayout = nullptr ;
		dd4hep::rec::LayeredCalorimeterData* _caloData = nullptr ;
		dd4hep::DetElement theDetector;

		SimDigitalGeomRPCFrame* _normal_I_J_setter = nullptr ;
		CHT::Layout _currentHCALCollectionCaloLayout = CHT::any ;
		LCVector3D _normal;
		LCVector3D _Iaxis;
		LCVector3D _Jaxis;
		static AIDA::ITuple* _tupleHit ;
		enum {TH_DETECTOR,TH_CHTLAYOUT,TH_MODULE,TH_TOWER,TH_STAVE,TH_LAYER,TH_I,TH_J,
			  TH_X,TH_Y,TH_Z,
			  TH_NORMALX,TH_NORMALY,TH_NORMALZ,
			  TH_IX,TH_IY,TH_IZ,
			  TH_JX,TH_JY,TH_JZ};
		static AIDA::ITuple* _tupleStep ;
		enum {TS_DETECTOR,TS_CHTLAYOUT,TS_HITCELLID,TS_NSTEP,
			  TS_HITX,TS_HITY,TS_HITZ,
			  TS_STEPX,TS_STEPY,TS_STEPZ,
			  TS_DELTAI,TS_DELTAJ,TS_DELTALAYER};

		static std::string _encodingStrings[ENCODINGTYPES][ENCODINGSTRINGLENGTH];

		std::string _cellIDEncodingString = "" ;

		bool _useGear = false ;

		friend class SimDigitalGeomRPCFrame;
};



//hierarchy of classes to determine the RPC reference frame
class SimDigitalGeomRPCFrame
{
	public:
		SimDigitalGeomRPCFrame(SimDigitalGeomCellId& h) : _layerInfo(h) {}
		virtual ~SimDigitalGeomRPCFrame() ;
		virtual void setRPCFrame() = 0 ;
	private :
		SimDigitalGeomCellId& _layerInfo ;
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
		SimDigitalGeomRPCFrame_TESLA_BARREL(SimDigitalGeomCellId& h) : SimDigitalGeomRPCFrame(h) {}
		void setRPCFrame();
};
class SimDigitalGeomRPCFrame_VIDEAU_BARREL : public SimDigitalGeomRPCFrame
{
	public:
		SimDigitalGeomRPCFrame_VIDEAU_BARREL(SimDigitalGeomCellId& h) : SimDigitalGeomRPCFrame(h) {}
		void setRPCFrame();
};
class SimDigitalGeomRPCFrame_TESLA_ENDCAP : public SimDigitalGeomRPCFrame
{
	public:
		SimDigitalGeomRPCFrame_TESLA_ENDCAP(SimDigitalGeomCellId& h) : SimDigitalGeomRPCFrame(h) {}
		void setRPCFrame();
};
class SimDigitalGeomRPCFrame_VIDEAU_ENDCAP : public SimDigitalGeomRPCFrame
{
	public:
		SimDigitalGeomRPCFrame_VIDEAU_ENDCAP(SimDigitalGeomCellId& h) : SimDigitalGeomRPCFrame(h) {}
		void setRPCFrame();
};

#endif //SimDigitalGeom_h
