#ifndef ChargeInducer_h
#define ChargeInducer_h

#include <string>
#include <map>


#include <boost/random/mersenne_twister.hpp>
#include <boost/random/gamma_distribution.hpp>


class SimDigitalGeomCellId ;
struct AsicKey ;

class ChargeInducer
{
	public :
		ChargeInducer() ;
		virtual ~ChargeInducer() ;
		virtual double getCharge(const SimDigitalGeomCellId& cellID) = 0 ;

		void setSeed(unsigned int value) ;

	protected :
		boost::mt19937 generator ;
} ;

class UniformPolya : public ChargeInducer
{
	public :
		UniformPolya(double _qbar , double _theta) ;
		~UniformPolya() ;

		virtual double getCharge(const SimDigitalGeomCellId& cellID) ;


	protected :
		boost::gamma_distribution<double> gammadist ;

} ;

class AsicPolya : public UniformPolya
{
	public :
		AsicPolya(double _qbar , double _theta , std::string fileName) ;
		~AsicPolya() ;

		virtual double getCharge(const SimDigitalGeomCellId& cellID) ;


	protected :
		void readFile(std::string fileName) ;

		std::map<AsicKey , boost::gamma_distribution<double> > polyaMap ;
} ;

#endif //ChargeInducer_h
