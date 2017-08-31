#ifndef ChargeSpreader_h
#define ChargeSpreader_h

#include <marlin/Global.h>

#include <vector>
#include <map>
#include <utility>
#include <set>

struct AsicKey ;
class SimDigitalGeomCellId ;

struct ChargeSpreaderParameters
{
		double cellSize = 10 ;
		double range = 30 ;
		double padSeparation = 0 ;

		//erf
		std::vector<float> erfWidth = {2} ;
		std::vector<float> erfWeigth = {1} ;

		//exact
		double d = 1 ;
} ;


class ChargeSpreader
{
	public :
		ChargeSpreader() ;
		virtual ~ChargeSpreader() ;

		virtual void setParameters(ChargeSpreaderParameters param) { parameters = param ; }

		virtual void init() = 0 ;

		typedef std::pair<int,int> I_J_Coordinates ;
		virtual void addCharge( double charge , double posI , double posJ , const SimDigitalGeomCellId& ) ;
		void newHit(float cellSize_) { chargeMap.clear() ; parameters.cellSize = cellSize_ ; }

		const std::map<I_J_Coordinates,double>& getChargeMap() const { return chargeMap ; }

	protected :
		virtual double computeIntegral(double x1 , double x2 , double y1 , double y2) const = 0 ;

		std::map<I_J_Coordinates,double> chargeMap ;
		ChargeSpreaderParameters parameters ;

		float normalisation = 0 ;
} ;


class GaussianSpreader : public ChargeSpreader
{
	public :
		GaussianSpreader() ;
		virtual ~GaussianSpreader() ;
		virtual void init() ;

	protected :
		virtual double computeIntegral(double x1 , double x2 , double y1 , double y2) const ;

	private :
		friend class SimDigital ;
} ;

class ExactSpreader : public ChargeSpreader
{
	public :
		ExactSpreader() ;
		virtual ~ExactSpreader() ;
		virtual void init() ;

	protected :
		double computeIntegral(double x1 , double x2 , double y1 , double y2) const ;

	private :
		friend class SimDigital ;
} ;

class ExactSpreaderPerAsic : public ExactSpreader
{
	public :
		ExactSpreaderPerAsic(std::string fileName) ;
		virtual ~ExactSpreaderPerAsic() ;

		virtual void setParameters(ChargeSpreaderParameters param) { parameters = param ; dGlobal = parameters.d ; }

		virtual void addCharge(double charge, double posI, double posJ , const SimDigitalGeomCellId& cellID) ;

	protected :

		double dGlobal = 1 ;
		void readFile(std::string fileName) ;


		std::map<AsicKey,double> dMap ;

} ;


#endif //ChargeSpreader_h
