#ifndef HdecayMode_h
#define HdecayMode_h 1

#include "marlin/Processor.h"
#include "IMPL/LCCollectionVec.h"
#include "lcio.h"
#include <string>
#include <TFile.h>
#include <TTree.h>
#include <vector>
class TFile;
class TH1F;
class TH1I;
class TH2I;
class TTree;

using namespace lcio ;
using namespace marlin ;
class HdecayMode : public Processor
{
	public:

		virtual Processor*  newProcessor()
		{
			return new HdecayMode;
		}
		HdecayMode();
		virtual ~HdecayMode() = default;
		HdecayMode(const HdecayMode&) = delete;
		HdecayMode& operator=(const HdecayMode&) = delete;
		virtual void init();
		virtual void processRunHeader();
		virtual void processEvent( EVENT::LCEvent *pLCEvent );
		virtual void check();
		virtual void end();
		void Clear();

	private:

		typedef std::vector<int>		IntVector;
		typedef std::vector<double>		DoubleVector;
		typedef std::vector<float>		FloatVector;

		std::string				m_mcParticleCollection{};
		std::string				m_HdecayModeCollection{};
		LCCollectionVec				*m_col_HDecayMode{};
		int                                     _nhiggs{};
		int					m_isDecayedTob;
		int					m_isDecayedToc;
		int					m_isDecayedToother;
		int                                     m_isBothDecayedToHeavy;
		float					m_ISR1Energy;
		float					m_ISR1Px;
		float					m_ISR1Py;
		float					m_ISR1Pz;
		float					m_ISR2Energy;
		float					m_ISR2Px;
		float					m_ISR2Py;
		float					m_ISR2Pz;
};

#endif
