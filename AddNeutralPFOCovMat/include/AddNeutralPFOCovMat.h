#ifndef AddNeutralPFOCovMat_h
#define AddNeutralPFOCovMat_h 1
#include <marlin/Processor.h>
#include <marlin/Global.h>
#include "UTIL/LCRelationNavigator.h"
#include <EVENT/MCParticle.h>
#include "IMPL/LCCollectionVec.h"
#include <IMPL/ReconstructedParticleImpl.h>
#include "TLorentzVector.h"
#include "TVector3.h"
#include "lcio.h"
#include <string>
#include <vector>
#include <math.h>
#include <set>
#include <vector>
#include "TMatrixD.h"
#include <TFile.h>
#include <TTree.h>

class TFile;
class TDirectory;
class TH1F;
class TH1I;
class TH2I;
class TH2F;
class TTree;

using namespace lcio ;
using namespace marlin ;

class AddNeutralPFOCovMat : public Processor
{
	public:

		virtual Processor*  newProcessor()
		{
			return new AddNeutralPFOCovMat;
		}
		AddNeutralPFOCovMat();
		virtual ~AddNeutralPFOCovMat() = default;
		AddNeutralPFOCovMat(const AddNeutralPFOCovMat&) = delete;
		AddNeutralPFOCovMat& operator=(const AddNeutralPFOCovMat&) = delete;
		virtual void init();
		virtual void Clear();
		virtual void processRunHeader();
		virtual void processEvent( EVENT::LCEvent *pLCEvent );
		std::vector<float> getNeutralCovMat( TVector3 clusterPosition , float pfoEc , float pfoMass , std::vector<float> clusterPositionError , float clusterEnergyError );
		virtual void check( EVENT::LCEvent *pLCEvent );
		virtual void end();

	private:
		std::string				m_inputPfoCollection{};
		std::string				m_outputPfoCollection{};

		bool					m_AssumeNeutralPFOMassive = true;
		bool					m_isClusterEnergyKinEnergy = false;
		bool					m_updatePFO4Momentum = false;
		bool					m_useTrueJacobian = false;

};
#endif
