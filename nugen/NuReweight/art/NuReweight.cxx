////////////////////////////////////////////////////////////////////////
/// \file  NuReweight.cxx
/// \brief Wrapper for reweightings neutrino interactions within the ART framework
///
/// \author  nathan.mayer@tufts.edu
////////////////////////////////////////////////////////////////////////

// C/C++ includes
#include <math.h>
#include <map>
#include <fstream>
#include <memory>

//ROOT includes
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TSystem.h"

//GENIE includes
#ifdef GENIE_PRE_R3
  #include "Conventions/Units.h"
  #include "EVGCore/EventRecord.h"
  #include "GHEP/GHepUtils.h"

#include "ReWeight/GReWeightI.h"
#include "ReWeight/GSystSet.h"
#include "ReWeight/GSyst.h"
#include "ReWeight/GReWeight.h"
#include "ReWeight/GReWeightNuXSecNCEL.h"
#include "ReWeight/GReWeightNuXSecCCQE.h"
#include "ReWeight/GReWeightNuXSecCCRES.h"
#include "ReWeight/GReWeightNuXSecCOH.h"
#include "ReWeight/GReWeightNonResonanceBkg.h"
#include "ReWeight/GReWeightFGM.h"
#include "ReWeight/GReWeightDISNuclMod.h"
#include "ReWeight/GReWeightResonanceDecay.h"
#include "ReWeight/GReWeightFZone.h"
#include "ReWeight/GReWeightINuke.h"
#include "ReWeight/GReWeightAGKY.h"
#include "ReWeight/GReWeightNuXSecCCQEvec.h"
#include "ReWeight/GReWeightNuXSecNCRES.h"
#include "ReWeight/GReWeightNuXSecDIS.h"
#include "ReWeight/GReWeightNuXSecNC.h"
#include "ReWeight/GSystUncertainty.h"
#include "ReWeight/GReWeightUtils.h"

  #include "Interaction/InitialState.h"
  #include "Interaction/Interaction.h"
  #include "Interaction/Kinematics.h"
  #include "Interaction/KPhaseSpace.h"
  #include "Interaction/ProcessInfo.h"
  #include "Interaction/XclsTag.h"
  #include "GHEP/GHepParticle.h"
  #include "PDG/PDGCodeList.h"
  #include "Conventions/Constants.h" //for calculating event kinematics

#else
  // GENIE includes R-3 and beyond
  #include "GENIE/Framework/Messenger/Messenger.h"
  #include "GENIE/Framework/Conventions/Units.h"
  #include "GENIE/Framework/Conventions/Constants.h"
  #include "GENIE/Framework/GHEP/GHepUtils.h"
  #include "GENIE/Framework/EventGen/EventRecord.h"

  #include "GENIE/Framework/GHEP/GHepParticle.h"

  #include "RwFramework/GReWeightI.h"
  #include "RwFramework/GSystSet.h"
  #include "RwFramework/GSyst.h"
  #include "RwFramework/GReWeight.h"
  #include "RwFramework/GSystUncertainty.h"
  #include "RwCalculators/GReWeightNuXSecNCEL.h"
  #include "RwCalculators/GReWeightNuXSecCCQE.h"
  #include "RwCalculators/GReWeightNuXSecCCRES.h"
  #include "RwCalculators/GReWeightNuXSecCOH.h"
  #include "RwCalculators/GReWeightNonResonanceBkg.h"
  #include "RwCalculators/GReWeightFGM.h"
  #include "RwCalculators/GReWeightDISNuclMod.h"
  #include "RwCalculators/GReWeightResonanceDecay.h"
  #include "RwCalculators/GReWeightFZone.h"
  #include "RwCalculators/GReWeightINuke.h"
  #include "RwCalculators/GReWeightAGKY.h"
  #include "RwCalculators/GReWeightNuXSecCCQEvec.h"
  #include "RwCalculators/GReWeightNuXSecNCRES.h"
  #include "RwCalculators/GReWeightNuXSecDIS.h"
  #include "RwCalculators/GReWeightNuXSecNC.h"
  #include "RwCalculators/GReWeightUtils.h"

#endif

#include "nugen/EventGeneratorBase/GENIE/GENIE2ART.h"

// NuGen includes
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "nusimdata/SimulationBase/GTruth.h"
#include "nugen/NuReweight/art/NuReweight.h"

namespace rwgt
{
  bool operator!=(const simb::GTruth & lhs, const simb::GTruth & rhs)
  {
    bool eq = lhs.fVertex == rhs.fVertex
           && lhs.fweight == rhs.fweight
           && lhs.fprobability == rhs.fprobability
           && lhs.fXsec == rhs.fXsec
           && lhs.fDiffXsec == rhs.fDiffXsec
           && lhs.fGPhaseSpace == rhs.fGPhaseSpace
           && lhs.fProbePDG == rhs.fProbePDG
           && lhs.fProbeP4 == rhs.fProbeP4
           && lhs.fTgtP4 == rhs.fTgtP4
           && lhs.ftgtZ == rhs.ftgtZ
           && lhs.ftgtA == rhs.ftgtA
           && lhs.ftgtPDG == rhs.ftgtPDG
           && lhs.fHitNucPDG == rhs.fHitNucPDG
           && lhs.fHitQrkPDG == rhs.fHitQrkPDG
           && lhs.fIsSeaQuark == rhs.fIsSeaQuark
           && lhs.fHitNucP4 == rhs.fHitNucP4
           && lhs.fHitNucPos == rhs.fHitNucPos
           && lhs.fGscatter == rhs.fGscatter
           && lhs.fGint == rhs.fGint
           && lhs.fgQ2 == rhs.fgQ2
           && lhs.fgq2 == rhs.fgq2
           && lhs.fgW == rhs.fgW
           && lhs.fgT == rhs.fgT
           && lhs.fgX == rhs.fgX
           && lhs.fgY == rhs.fgY
           && lhs.fgWrun == rhs.fgWrun
           && lhs.fFSleptonP4 == rhs.fFSleptonP4
           && lhs.fFShadSystP4 == rhs.fFShadSystP4
           && lhs.fIsCharm == rhs.fIsCharm
           && lhs.fCharmHadronPdg == rhs.fCharmHadronPdg
           && lhs.fStrangeHadronPdg == rhs.fStrangeHadronPdg
           && lhs.fNumProton == rhs.fNumProton
           && lhs.fNumNeutron == rhs.fNumNeutron
           && lhs.fNumPi0 == rhs.fNumPi0
           && lhs.fNumPiPlus == rhs.fNumPiPlus
           && lhs.fNumPiMinus == rhs.fNumPiMinus
           && lhs.fNumSingleGammas == rhs.fNumSingleGammas
           && lhs.fNumRho0 == rhs.fNumRho0
           && lhs.fNumRhoPlus == rhs.fNumRhoPlus
           && lhs.fNumRhoMinus == rhs.fNumRhoMinus
           && lhs.fResNum == rhs.fResNum
           && lhs.fDecayMode == rhs.fDecayMode
           && lhs.fFinalQuarkPdg == rhs.fFinalQuarkPdg
           && lhs.fFinalLeptonPdg == rhs.fFinalLeptonPdg;
   return !eq;
  }


  ///<constructor
  NuReweight::NuReweight() {
    //mf::LogVerbatim("GENIEReweight") << "Create GENIEReweight object";
  }

  ///<destructor
  NuReweight::~NuReweight() {
    // Don't delete fWcalc here. The GENIEReweight parent class handles it.
  }

  double NuReweight::CalcWeight(const simb::MCTruth & truth, const simb::GTruth & gtruth) const {

  // do some rudimentary caching to prevent repeated requests for the same event
  // from creating new genie::EventRecords over and over
  static simb::GTruth cached_gtruth = simb::GTruth();
  static std::unique_ptr<genie::EventRecord> evr = nullptr;
  if (cached_gtruth != gtruth)
  {
    cached_gtruth = gtruth;
    evr.reset(evgb::RetrieveGHEP(truth, gtruth));
  }

    double wgt = this->CalculateWeight(*evr);

    //mf::LogVerbatim("GENIEReweight") << "New Event Weight is: " << wgt;
    return wgt;
  }


}
