#ifndef SimDigitalToEnergy_h
#define SimDigitalToEnergy_h

#include <marlin/Processor.h>
#include <IMPL/LCCollectionVec.h>

class SimDigitalToEnergy : public marlin::Processor
{
	public :
		virtual Processor* newProcessor() { return new SimDigitalToEnergy ; }
		SimDigitalToEnergy() ;

		/** Called at the begin of the job before anything is read.
		* Use to initialize the processor, e.g. book histograms.*/
		virtual void init() ;

		/** Called for every event - the working horse. */
		virtual void processEvent(LCEvent* evt) ;

		/** Called after data processing for clean up. */
		virtual void end() ;

		SimDigitalToEnergy(const SimDigitalToEnergy &toCopy) = delete ;
		void operator=(const SimDigitalToEnergy &toCopy) = delete ;

	protected :
		LCCollectionVec* processCollection(LCCollection* col) ;

	private :

		std::vector<std::string> _inputHcalCollections {} ;
		std::vector<std::string> _outputHcalCollections {} ;
		std::string _inputRelCollection = "" ;
		std::string _outputRelCollection = "" ;

		LCCollection* inputRelCol = nullptr ;
		LCCollectionVec* _relcol = nullptr ;

		std::vector<float> _energyCoefficients {} ;


} ;

#endif //SimDigitalToEnergy_h
