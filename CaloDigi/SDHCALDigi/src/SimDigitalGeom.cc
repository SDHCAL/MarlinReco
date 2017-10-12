#include "SimDigitalGeom.h"

#include "SimDigital.h"

#include <marlin/Global.h>
#include <marlin/AIDAProcessor.h>
#include <AIDA/ITupleFactory.h>
#include <AIDA/ITuple.h>

#include <DDSegmentation/BitField64.h>
#include <DDRec/CellIDPositionConverter.h>
#include <DD4hep/DetectorSelector.h>

using namespace lcio ;
using namespace marlin ;


AIDA::ITuple* SimDigitalGeomCellId::_tupleHit = NULL;
AIDA::ITuple* SimDigitalGeomCellId::_tupleStep = NULL;

void SimDigitalGeomCellId::bookTuples(const marlin::Processor* proc)
{
	_tupleHit  = AIDAProcessor::tupleFactory( proc )->create("SimDigitalGeom",
															 "SimDigital_Debug",
															 "int chtlayout,module,tower,stave,layer,I,J, float x,y,z, normalx,normaly,normalz, Ix,Iy,Iz,Jx,Jy,Jz");
	streamlog_out(DEBUG) << "Tuple for Hit has been initialized to " << _tupleHit << std::endl;
	streamlog_out(DEBUG)<< "it has " << _tupleHit->columns() << " columns" <<std::endl;

	_tupleStep = AIDAProcessor::tupleFactory( proc )->create("SimDigitalStep",
															 "SimDigital_DebugStep",
															 "int chtlayout,hitcellid,nstep, float hitx,hity,hitz,stepx,stepy,stepz,deltaI,deltaJ,deltaLayer,time");
	streamlog_out(DEBUG) << "Tuple for Step has been initialized to " << _tupleStep << std::endl;
	streamlog_out(DEBUG) << "it has " << _tupleStep->columns() << " columns" <<std::endl;
}


SimDigitalGeomCellId::SimDigitalGeomCellId(LCCollection* inputCol, LCCollectionVec* outputCol)
	: _decoder(inputCol) , _encoder(inputCol->getParameters().getStringVal(LCIO::CellIDEncoding),outputCol) ,
	  _normal() , _Iaxis() , _Jaxis()
{
	outputCol->parameters().setValue(LCIO::CellIDEncoding,inputCol->getParameters().getStringVal(LCIO::CellIDEncoding)) ;
	_cellIDEncodingString = inputCol->getParameters().getStringVal(LCIO::CellIDEncoding) ;
}

SimDigitalGeomCellIdLCGEO::SimDigitalGeomCellIdLCGEO(LCCollection* inputCol, LCCollectionVec* outputCol)
	: SimDigitalGeomCellId(inputCol , outputCol)
	//	  theDetector()
{
	streamlog_out( DEBUG ) << "we will use lcgeo!" << std::endl ;
}

SimDigitalGeomCellIdMOKKA::SimDigitalGeomCellIdMOKKA(LCCollection* inputCol, LCCollectionVec* outputCol)
	: SimDigitalGeomCellId(inputCol , outputCol)
{
	std::string gearFile = Global::parameters->getStringVal("GearXMLFile") ;

	_geom = VIDEAU ;

	// maybe we can also set the geometry by the hcal option parameter
	try
	{
		Global::GEAR->getHcalBarrelParameters().getIntVal("Hcal_outer_polygon_order") ; //it is VIDEAU geometry if it is OK
	}
	catch (gear::Exception &)
	{
		_geom = TESLA ;
	}

	streamlog_out( DEBUG ) << "we will use gear!" << std::endl ;
	streamlog_out( DEBUG ) << "gear: " << Global::GEAR << std::endl ;
	streamlog_out( DEBUG )<< "!!!!!!(Videau=0, TESLA=1) Geometry is _geom= : "<<_geom << std::endl;
}

SimDigitalGeomCellId::~SimDigitalGeomCellId()
{
}

SimDigitalGeomCellIdLCGEO::~SimDigitalGeomCellIdLCGEO()
{
}

SimDigitalGeomCellIdMOKKA::~SimDigitalGeomCellIdMOKKA()
{
	if (_normal_I_J_setter != nullptr)
		delete _normal_I_J_setter;
}


SimDigitalGeomRPCFrame::~SimDigitalGeomRPCFrame()
{
}

void SimDigitalGeomRPCFrame_TESLA_BARREL::setRPCFrame()
{
	normal().set(1,0,0);
	Iaxis().set(0,1,0);
	Jaxis().set(0,0,1);
	double angle=(stave()/2)*(-45)*CLHEP::degree;
	normal().rotateZ(angle);
	Iaxis().rotateZ(angle);
}

void SimDigitalGeomRPCFrame_VIDEAU_BARREL::setRPCFrame()
{
	normal().set(-1,0,0);
	Iaxis().set(0,1,0);
	Jaxis().set(0,0,1);
	double angle=(stave()-1)*(45)*CLHEP::degree;
	normal().rotateZ(angle);
	Iaxis().rotateZ(angle);
}

void SimDigitalGeomRPCFrame_TESLA_ENDCAP::setRPCFrame()
{
	if (module()==6)
	{
		normal().set(0,0,1);
		Iaxis().set(1,0,0);
		Jaxis().set(0,1,0);
		Iaxis().rotateZ(stave()*(-90)*CLHEP::degree);
		Jaxis().rotateZ(stave()*(-90)*CLHEP::degree);
	}
	else if (module()==4)
	{
		normal().set(0,0,-1);
		Iaxis().set(-1,0,0);
		Jaxis().set(0,1,0);
		Iaxis().rotateZ(stave()*90*CLHEP::degree);
		Jaxis().rotateZ(stave()*90*CLHEP::degree);
	}
	else
	{
		streamlog_out(ERROR) << "ERROR ; TESLA detector : unknown module for endcap " << module() << std::endl;
	}
}

//valid also for all endcap rings (VIDEAU and TESLA)
void SimDigitalGeomRPCFrame_VIDEAU_ENDCAP::setRPCFrame()
{
	if (module()==6)
	{
		//streamlog_out( DEBUG )<< "!!!!!!!!!! VIDEAU_ENDCAP : module=6 "<< std::endl;
		normal().set(0,0,1);
		Iaxis().set(1,0,0);
		Jaxis().set(0,1,0);
		Iaxis().rotateZ(stave()*(90)*CLHEP::degree);
		Jaxis().rotateZ(stave()*(90)*CLHEP::degree);
	}
	else if (module()==0)
	{
		//streamlog_out( DEBUG )<< "!!!!!!!!!! VIDEAU_ENDCAP : module 0"<< std::endl;
		normal().set(0,0,-1);
		Iaxis().set(-1,0,0);
		Jaxis().set(0,1,0);
		Iaxis().rotateZ(stave()*(-90)*CLHEP::degree);
		Jaxis().rotateZ(stave()*(-90)*CLHEP::degree);
	}
	else
	{
		streamlog_out(ERROR) << "ERROR : unknown module for endcap or endcapring " << module() << std::endl;
	}
}

void SimDigitalGeomCellId::createStepAndChargeVec(SimCalorimeterHit* hit , std::vector<StepAndCharge>& vec)
{
	LCVector3D hitpos ;
	if (NULL != _hitPosition)
		hitpos.set(_hitPosition[0],_hitPosition[1],_hitPosition[2]) ;
	for (int imcp = 0 ; imcp < hit->getNMCContributions() ; imcp++)
	{
		LCVector3D stepposvec;
		const float* steppos = hit->getStepPosition(imcp) ;
		if (NULL != steppos)
			stepposvec.set(steppos[0],steppos[1],steppos[2]) ;
		if (stepposvec.mag2() != 0)
		{
			stepposvec -= hitpos ;
			vec.push_back( StepAndCharge(LCVector3D(stepposvec*_Iaxis,stepposvec*_Jaxis,stepposvec*_normal) , hit->getTimeCont(imcp)) );
		}
		else
			streamlog_out(WARNING) << "DIGITISATION : STEP POSITION IS (0,0,0)" << std::endl;
	}

	//if no steps have been found, then put one step at the center of the cell :
	if (vec.size() == 0)
	{
		streamlog_out(MESSAGE) << "no Steps in hit" << std::endl ;
		//		vec.push_back(StepAndCharge(LCVector3D(0,0,0))) ;
	}
}


std::vector<StepAndCharge> SimDigitalGeomCellIdLCGEO::decode(SimCalorimeterHit* hit)
{
	_cellIDvalue = _decoder( hit ).getValue() ;

	_trueLayer = _decoder( hit )[_encodingString.at(0)] - 1 ; // -1;
	_stave     = _decoder( hit )[_encodingString.at(1)] ;     // +1
	_module    = _decoder( hit )[_encodingString.at(2)] ;

	if( _encodingString.at(3).size() != 0)
		_tower     = _decoder( hit )[_encodingString.at(3)] ;

	_Iy        = _decoder( hit )[_encodingString.at(4)] ;
	try
	{
		_Jz = _decoder( hit )[_encodingString.at(5)] ;
	}
	catch (lcio::Exception &)
	{
		_encodingString.at(5) = "z" ;

		try
		{
			_Jz = _decoder( hit )[_encodingString.at(5)] ;
		}
		catch(lcio::Exception &)
		{
			_encodingString.at(5) = "y" ;
			_Jz = _decoder( hit )[_encodingString.at(5)] ;
		}
	}

	// _slice     = _decoder( hit )["slice"];
	_hitPosition = hit->getPosition() ;
	if(abs(_Iy)<1 && abs(_Iy)!=0.0)
		streamlog_out(DEBUG) << "_Iy, _Jz:"<<_Iy <<" "<<_Jz<< std::endl;
	//if(_module==0||_module==6) streamlog_out( DEBUG )<<"tower "<<_tower<<" layer "<<_trueLayer<<" stave "<<_stave<<" module "<<_module<<std::endl;
	//<<" Iy " << _Iy <<"  Jz "<<_Jz<<" hitPosition "<<_hitPosition<<std::endl
	//<<" _hitPosition[0] "<<_hitPosition[0]<<" _hitPosition[1] "<<_hitPosition[1]<<" _hitPosition[2] "<<_hitPosition[2]<<std::endl;


	dd4hep::Detector& ild = dd4hep::Detector::getInstance() ;
	dd4hep::rec::CellIDPositionConverter idposConv( ild )  ;

	dd4hep::BitField64 idDecoder( _cellIDEncodingString ) ;

	dd4hep::long64 id0 = hit->getCellID0() ;
	dd4hep::long64 id1 = hit->getCellID1() ;

	idDecoder.setValue( id0 , id1 ) ;

	dd4hep::long64 id = idDecoder.getValue() ;
	dd4hep::Position pos_0 = idposConv.position( id ) ;

#if 0
	const float* hitPos = hit->getPosition();

	streamlog_out( DEBUG ) << "hit pos: " << hitPos[0] << " " << hitPos[1] << " " << hitPos[2] << std::endl;
	streamlog_out( DEBUG ) << "cell pos: " << pos_0.X() << " " << pos_0.Y() << " " << pos_0.Z() << std::endl;

	streamlog_out( DEBUG ) << "layer: "    << idDecoder[_encodingStrings[_encodingType][0]]
			<< ", stave: "  << idDecoder[_encodingStrings[_encodingType][1]]
			<< ", module: " << idDecoder[_encodingStrings[_encodingType][2]]
			<< ", tower: "  << idDecoder[_encodingStrings[_encodingType][3]]
			<< ", x: "      << idDecoder[_encodingStrings[_encodingType][4]]
			<< ", y: "      << idDecoder[_encodingStrings[_encodingType][5]] << std::endl;
#endif

	double const epsilon = 1.e-3;

	///// for direction x
	dd4hep::Position pos_i_plus_1;
	dd4hep::Position dir_i;

	std::string xEncoding = _encodingString.at(4) ;
	std::string yEncoding = _encodingString.at(5) ;

	idDecoder[ xEncoding ]   =  idDecoder[ xEncoding ]  + 1 ;
	pos_i_plus_1 = idposConv.position(   idDecoder.getValue()   ) ;

	dir_i =  pos_i_plus_1  - pos_0  ;

	if(dir_i.R() < epsilon)
	{
		idDecoder[ xEncoding ]   =  idDecoder[ xEncoding ] - 2 ;
		pos_i_plus_1 = idposConv.position(   idDecoder.getValue()   ) ;
		dir_i =  - ( pos_i_plus_1  - pos_0 )  ;
	}

	// reset
	idDecoder.setValue( id0 , id1 ) ;

	////// for direction y
	dd4hep::Position pos_j_plus_1;
	dd4hep::Position dir_j;

	idDecoder[ yEncoding ]   =  idDecoder[ yEncoding ]  + 1 ;
	pos_j_plus_1 = idposConv.position(   idDecoder.getValue()   ) ;

	dir_j =  pos_j_plus_1  - pos_0  ;

	if(dir_j.R() < epsilon)
	{
		idDecoder[ yEncoding ]   =  idDecoder[ yEncoding ] - 2 ;
		pos_j_plus_1 = idposConv.position(   idDecoder.getValue()   ) ;
		dir_j =  - ( pos_j_plus_1  - pos_0 )  ;
	}

	dd4hep::Position dir_layer = dir_i.Cross( dir_j ) ;

	dir_layer = - dir_layer.Unit();

	//streamlog_out( DEBUG ) << "layer dir: " << dir_layer.X() << " " << dir_layer.Y() << " " << dir_layer.Z() << std::endl;

	dir_i = dir_i.Unit();
	dir_j = dir_j.Unit();

	_normal.set(dir_layer.X(), dir_layer.Y(), dir_layer.Z()) ;
	_Iaxis.set(dir_i.X(), dir_i.Y(), dir_i.Z()) ;
	_Jaxis.set(dir_j.X(), dir_j.Y(), dir_j.Z()) ;


	std::vector<StepAndCharge> stepsInIJZcoord ;
	createStepAndChargeVec(hit , stepsInIJZcoord) ;

	fillDebugTupleGeometryHit() ;
	fillDebugTupleGeometryStep(hit , stepsInIJZcoord) ;

	return stepsInIJZcoord ;
}

std::vector<StepAndCharge> SimDigitalGeomCellIdMOKKA::decode(SimCalorimeterHit* hit)
{
	_cellIDvalue = _decoder( hit ).getValue() ;

	_trueLayer = _decoder( hit )[_encodingString.at(0)] - 1; // -1;
	_stave     = _decoder( hit )[_encodingString.at(1)] ;     // +1
	_module    = _decoder( hit )[_encodingString.at(2)] ;

	if( _encodingString.at(3).size() != 0)
		_tower     = _decoder( hit )[_encodingString.at(3)] ;

	_Iy        = _decoder( hit )[_encodingString.at(4)] ;
	try
	{
		_Jz = _decoder( hit )[_encodingString.at(5)] ;
	}
	catch (lcio::Exception &)
	{
		_encodingString.at(5) = "z" ;

		try
		{
			_Jz = _decoder( hit )[_encodingString.at(5)] ;
		}
		catch(lcio::Exception &)
		{
			_encodingString.at(5) = "y" ;
			_Jz = _decoder( hit )[_encodingString.at(5)] ;
		}
	}

	// _slice     = _decoder( hit )["slice"];
	_hitPosition = hit->getPosition();
	if(abs(_Iy)<1 && abs(_Iy)!=0.0)
		streamlog_out(DEBUG) << "_Iy, _Jz:"<<_Iy <<" "<<_Jz<< std::endl;
	//if(_module==0||_module==6) streamlog_out( DEBUG )<<"tower "<<_tower<<" layer "<<_trueLayer<<" stave "<<_stave<<" module "<<_module<<std::endl;
	//<<" Iy " << _Iy <<"  Jz "<<_Jz<<" hitPosition "<<_hitPosition<<std::endl
	//<<" _hitPosition[0] "<<_hitPosition[0]<<" _hitPosition[1] "<<_hitPosition[1]<<" _hitPosition[2] "<<_hitPosition[2]<<std::endl;

	_normal_I_J_setter->setRPCFrame() ;


	std::vector<StepAndCharge> stepsInIJZcoord ;
	createStepAndChargeVec(hit , stepsInIJZcoord) ;

	fillDebugTupleGeometryHit() ;
	fillDebugTupleGeometryStep(hit , stepsInIJZcoord) ;

	return stepsInIJZcoord ;
}

void SimDigitalGeomCellId::fillDebugTupleGeometryHit()
{
	//these tuples are for debugging geometry aspects
	if (_tupleHit != nullptr)
	{
		_tupleHit->fill(TH_CHTLAYOUT,int(_currentHCALCollectionCaloLayout));
		_tupleHit->fill(TH_MODULE,_module);
		_tupleHit->fill(TH_TOWER,_tower);
		_tupleHit->fill(TH_STAVE,_stave);
		_tupleHit->fill(TH_LAYER,_trueLayer);
		_tupleHit->fill(TH_I,_Iy);
		_tupleHit->fill(TH_J,_Jz);
		if (_hitPosition != NULL)
		{
			_tupleHit->fill(TH_X,_hitPosition[0]); //x
			_tupleHit->fill(TH_Y,_hitPosition[1]); //y
			_tupleHit->fill(TH_Z,_hitPosition[2]); //z
		}
		else
		{
			float notset=-88888;
			_tupleHit->fill(TH_X,notset);
			_tupleHit->fill(TH_Y,notset);
			_tupleHit->fill(TH_Z,notset);
		}
		for (int i=0; i<3; i++)
		{
			_tupleHit->fill(TH_NORMALX+i,_normal[i]);
			_tupleHit->fill(TH_IX+i,_Iaxis[i]);
			_tupleHit->fill(TH_JX+i,_Jaxis[i]);
		}
		_tupleHit->addRow() ;
	}
}

void SimDigitalGeomCellId::fillDebugTupleGeometryStep(SimCalorimeterHit* hit , const std::vector<StepAndCharge>& stepsInIJZcoord)
{
	if (_tupleStep != nullptr )
	{
		int nsteps = hit->getNMCContributions() ;
		float notset=-88888;
		for (int imcp = 0 ; imcp<nsteps ; imcp++)
		{
			_tupleStep->fill(TS_CHTLAYOUT,int(_currentHCALCollectionCaloLayout));
			_tupleStep->fill(TS_HITCELLID,hit->getCellID0());
			_tupleStep->fill(TS_NSTEP,nsteps);
			const float* steppos = hit->getStepPosition(imcp) ;
			for (int i=0; i<3; i++)
			{
				if (_hitPosition != NULL )
					_tupleStep->fill(TS_HITX+i,_hitPosition[i]);
				else
					_tupleStep->fill(TS_HITX+i,notset);


				if (steppos != NULL)
					_tupleStep->fill(TS_STEPX+i,steppos[i]);
				else
					_tupleStep->fill(TS_STEPX+i,notset);


				if (imcp < (int)stepsInIJZcoord.size() )
					_tupleStep->fill(TS_DELTAI+i,stepsInIJZcoord[imcp].step[i]);
				else
					_tupleStep->fill(TS_DELTAI+i,notset);
			}
			if (imcp < (int)stepsInIJZcoord.size() )
				_tupleStep->fill(TS_TIME,hit->getTimeCont(imcp)) ;
			else
				_tupleStep->fill(TS_TIME,notset) ;

			_tupleStep->addRow() ;
		}
	}
}


void SimDigitalGeomCellIdLCGEO::encode(CalorimeterHitImpl* hit, int delta_I, int delta_J)
{
	_encoder.setValue(_cellIDvalue) ;

	int RealIy = _Iy+delta_I ;
	//  streamlog_out( DEBUG ) << "RealIy, _Iy, delta_I" << std::endl <<RealIy <<" "<<_Iy <<" " <<delta_I<<std::endl;

	_encoder[_encodingString.at(4)] = RealIy ;
	int RealJz = _Jz+delta_J ;
	//  streamlog_out( DEBUG ) << "RealIy, _Iy, delta_I" <<RealIy<<" "<<_Iy<<" "<<delta_I<< std::endl;
	//  streamlog_out( DEBUG ) << "RealJz, _Jz, delta_J" <<RealJz<<" "<<_Jz<<" "<<delta_J<< std::endl;
	//       <<RealJz <<" "<<_Jz <<" " <<delta_J<<std::endl;

	// Depending on the segmentation type:   Barrel,EndcapRing - CartesianGridXY;  EndCaps- CartesianGridXZ !!!

	_encoder[_encodingString.at(5)] = RealJz ;

	_encoder.setCellID( hit );
	//streamlog_out( DEBUG ) << "CellID0: " << hit->getCellID0() << ", CellID1: " << hit->getCellID1() << " --> " << _encoder.valueString() << std::endl;
	hit->setType( CHT( CHT::had, CHT::hcal , _currentHCALCollectionCaloLayout,  _trueLayer ) );
	//  streamlog_out( DEBUG )   <<"getCellSize() "<<getCellSize()<<" layer " <<_trueLayer<< std::endl;

	float posB[3] ;
	posB[0]=_hitPosition[0]+getCellSize()*(delta_I*_Iaxis.x()+delta_J*_Jaxis.x()) ;
	posB[1]=_hitPosition[1]+getCellSize()*(delta_I*_Iaxis.y()+delta_J*_Jaxis.y()) ;
	posB[2]=_hitPosition[2]+getCellSize()*(delta_I*_Iaxis.z()+delta_J*_Jaxis.z()) ;

	hit->setPosition(posB) ;
}

void SimDigitalGeomCellIdMOKKA::encode(CalorimeterHitImpl* hit, int delta_I, int delta_J)
{
	_encoder.setValue(_cellIDvalue);

	int RealIy=_Iy+delta_I;
	//  streamlog_out( DEBUG ) << "RealIy, _Iy, delta_I" << std::endl <<RealIy <<" "<<_Iy <<" " <<delta_I<<std::endl;

	if (abs(RealIy)>330)
		RealIy=0; //FIXME the 330 value should depend on the cellSize and on the Layer
	_encoder[_encodingString.at(4)] = RealIy ;
	int RealJz=_Jz+delta_J;
	//  streamlog_out( DEBUG ) << "RealIy, _Iy, delta_I" <<RealIy<<" "<<_Iy<<" "<<delta_I<< std::endl;
	//  streamlog_out( DEBUG ) << "RealJz, _Jz, delta_J" <<RealJz<<" "<<_Jz<<" "<<delta_J<< std::endl;
	//       <<RealJz <<" "<<_Jz <<" " <<delta_J<<std::endl;

	if (abs(RealJz)>330)
		RealJz=0; //FIXME the 330 value should depend on the cellSize and on the Layer

	if (abs(RealIy)>330||abs(RealJz)>330)
		streamlog_out( DEBUG ) << "RealIy, RealJz" << std::endl <<RealIy <<"   "<<RealJz <<std::endl;

	// Depending on the segmentation type:   Barrel,EndcapRing - CartesianGridXY;  EndCaps- CartesianGridXZ !!!

	_encoder[_encodingString.at(5)] = RealJz ;

	_encoder.setCellID( hit );
	//streamlog_out( DEBUG ) << "CellID0: " << hit->getCellID0() << ", CellID1: " << hit->getCellID1() << " --> " << _encoder.valueString() << std::endl;
	hit->setType( CHT( CHT::had, CHT::hcal , _currentHCALCollectionCaloLayout,  _trueLayer ) );
	//  streamlog_out( DEBUG )   <<"getCellSize() "<<getCellSize()<<" layer " <<_trueLayer<< std::endl;

	float posB[3] ;

	posB[0]=_hitPosition[0]+getCellSize()*(delta_I*_Iaxis.x()+delta_J*_Jaxis.x()) ;
	posB[1]=_hitPosition[1]+getCellSize()*(delta_I*_Iaxis.y()+delta_J*_Jaxis.y()) ;
	posB[2]=_hitPosition[2]+getCellSize()*(delta_I*_Iaxis.z()+delta_J*_Jaxis.z()) ;
	hit->setPosition(posB) ;
}



void SimDigitalGeomCellIdLCGEO::setLayerLayout(CHT::Layout layout)
{
	_currentHCALCollectionCaloLayout = layout ;

	dd4hep::Detector & ild = dd4hep::Detector::getInstance() ;

	if(_currentHCALCollectionCaloLayout == CHT::barrel)
	{
		const std::vector< dd4hep::DetElement>& det = dd4hep::DetectorSelector(ild).detectors(
														  (dd4hep::DetType::CALORIMETER | dd4hep::DetType::HADRONIC | dd4hep::DetType::BARREL),
														  (dd4hep::DetType::AUXILIARY  |  dd4hep::DetType::FORWARD) ) ;

		_caloData = det.at(0).extension<dd4hep::rec::LayeredCalorimeterData>();
		//streamlog_out( DEBUG ) << "det size: " << det.size() << ", type: " << det.at(0).type() << endl;
	}

	if(_currentHCALCollectionCaloLayout == CHT::ring)
	{
		const std::vector< dd4hep::DetElement>& det = dd4hep::DetectorSelector(ild).detectors(
														  (dd4hep::DetType::CALORIMETER | dd4hep::DetType::HADRONIC | dd4hep::DetType::AUXILIARY),
														  dd4hep::DetType::FORWARD) ;

		_caloData = det.at(0).extension<dd4hep::rec::LayeredCalorimeterData>();

		//streamlog_out( DEBUG ) << "det size: " << det.size() << ", type: " << det.at(0).type() << endl;
	}

	if(_currentHCALCollectionCaloLayout == CHT::endcap)
	{
		const std::vector< dd4hep::DetElement>& det = dd4hep::DetectorSelector(ild).detectors(
														  (dd4hep::DetType::CALORIMETER | dd4hep::DetType::HADRONIC | dd4hep::DetType::ENDCAP),
														  (dd4hep::DetType::AUXILIARY  |  dd4hep::DetType::FORWARD) ) ;

		_caloData = det.at(0).extension<dd4hep::rec::LayeredCalorimeterData>();
		//streamlog_out( DEBUG ) << "det size: " << det.size() << ", type: " << det.at(0).type() << endl;
	}

}

void SimDigitalGeomCellIdMOKKA::setLayerLayout(CHT::Layout layout)
{
	if(_normal_I_J_setter != nullptr)
		delete _normal_I_J_setter ;

	_currentHCALCollectionCaloLayout=layout;


	if (_currentHCALCollectionCaloLayout == CHT::endcap)
	{
		_layerLayout = & Global::GEAR->getHcalEndcapParameters().getLayerLayout();

		if (_geom == TESLA)
			_normal_I_J_setter= new SimDigitalGeomRPCFrame_TESLA_ENDCAP(*this);
		else
			_normal_I_J_setter= new SimDigitalGeomRPCFrame_VIDEAU_ENDCAP(*this);
	}
	else if (_currentHCALCollectionCaloLayout == CHT::ring)
	{
		_layerLayout = & Global::GEAR->getHcalRingParameters().getLayerLayout();

		// no TESLA ring ?
		_normal_I_J_setter= new SimDigitalGeomRPCFrame_VIDEAU_ENDCAP(*this);
	}
	else
	{
		_layerLayout = & Global::GEAR->getHcalBarrelParameters().getLayerLayout();

		if (_geom == TESLA)
			_normal_I_J_setter= new SimDigitalGeomRPCFrame_TESLA_BARREL(*this);
		else
			_normal_I_J_setter= new SimDigitalGeomRPCFrame_VIDEAU_BARREL(*this);
	}
}


float SimDigitalGeomCellIdLCGEO::getCellSize()
{
	float cellSize = 0.f ;
	const double CM2MM = 10.0 ;

	if ( _caloData != nullptr )
	{
		if ( _cellSize > 0.0f )
		{
			cellSize = _cellSize ;
		}
		else
		{
			const std::vector<dd4hep::rec::LayeredCalorimeterStruct::Layer>& hcalBarrelLayers = _caloData->layers ;
			cellSize = hcalBarrelLayers[_trueLayer].cellSize0 * CM2MM ;
		}
	}

	//streamlog_out( MESSAGE ) << "cellSize: " << cellSize << endl;

	return cellSize ;
}

float SimDigitalGeomCellIdMOKKA::getCellSize()
{
	float cellSize = 0.f ;

	if ( _layerLayout != nullptr )
		cellSize = _layerLayout->getCellSize0(_trueLayer) ;

	//streamlog_out( MESSAGE ) << "cellSize: " << cellSize << endl;

	return cellSize ;
}
