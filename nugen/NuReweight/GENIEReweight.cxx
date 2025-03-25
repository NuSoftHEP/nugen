////////////////////////////////////////////////////////////////////////
/// \file  GENIEReweight.cxx
/// \brief Wrapper for reweight neutrino interactions with GENIE base class
///
/// \author  nathan.mayer@tufts.edu
////////////////////////////////////////////////////////////////////////

// C/C++ includes
#include <math.h>
#include <map>
#include <fstream>

//ROOT includes
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TSystem.h"

//GENIE includes
#ifdef GENIE_PRE_R3
  #include "Conventions/Units.h"
  #include "EVGCore/EventRecord.h"
  #include "GHEP/GHepUtils.h"
  #include "Messenger/Messenger.h"

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

//#include "Geo/ROOTGeomAnalyzer.h"
//#include "Geo/GeomVolSelectorFiducial.h"
//#include "Geo/GeomVolSelectorRockBox.h"
  #include "Utils/StringUtils.h"
  #include "Utils/XmlParserUtils.h"
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

//  #include "GENIE/Framework/Interaction/InitialState.h"
//  #include "GENIE/Framework/Interaction/Interaction.h"
//  #include "GENIE/Framework/Interaction/Kinematics.h"
//  #include "GENIE/Framework/Interaction/KPhaseSpace.h"
//  #include "GENIE/Framework/Interaction/ProcessInfo.h"
//  #include "GENIE/Framework/Interaction/XclsTag.h"

//  #include "GENIE/Framework/ParticleData/PDGCodes.h"
//  #include "GENIE/Framework/ParticleData/PDGCodeList.h"
//  #include "GENIE/Framework/ParticleData/PDGLibrary.h"
//  #include "GENIE/Framework/GHEP/GHepUtils.h"
//  #include "GENIE/Framework/GHEP/GHepParticle.h"

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

//#include "Geo/ROOTGeomAnalyzer.h"
//#include "Geo/GeomVolSelectorFiducial.h"
//#include "Geo/GeomVolSelectorRockBox.h"
//#include "Utils/StringUtils.h"
//#include "Utils/XmlParserUtils.h"

#endif
// Necessary because the GENIE LOG_* macros don't fully qualify Messenger
using genie::Messenger;


// NuGen includes
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "nusimdata/SimulationBase/GTruth.h"
#include "nugen/NuReweight/GENIEReweight.h"

// Framework includes
//#include "messagefacility/MessageLogger/MessageLogger.h"



namespace rwgt {
  ///<constructor
  GENIEReweight::GENIEReweight() :
        fUseSigmaDef(true) {

    LOG_INFO("GENIEReweight") << "Create GENIEReweight object";

    fWcalc = new genie::rew::GReWeight();
    this->SetNominalValues();
  }

  ///<destructor
  GENIEReweight::~GENIEReweight() {
    delete fWcalc;
  }

  ///<Set the nominal values for the reweight parameters.
  void GENIEReweight::SetNominalValues() {
    //CCQE Nominal Values
    fNominalParameters[(int)rwgt::fReweightMaNCEL] = 0.99;

    fNominalParameters[(int)rwgt::fReweightEtaNCEL] = 0.12;

    //CCQE Nominal Values
    fNominalParameters[(int)rwgt::fReweightNormCCQE] = 1.0;

    fNominalParameters[(int)rwgt::fReweightNormCCQEenu] = 1.0;

    fNominalParameters[(int)rwgt::fReweightMaCCQEshape] = 0.99;

    fNominalParameters[(int)rwgt::fReweightMaCCQE] = 0.99;
    fNominalParameters[(int)rwgt::fReweightVecCCQEshape] = 0.84;

    //Resonance Nominal Values
    fNominalParameters[(int)rwgt::fReweightNormCCRES] = 1.0;
    fNominalParameters[(int)rwgt::fReweightMaCCRESshape] = 1.12;
    fNominalParameters[(int)rwgt::fReweightMvCCRESshape]  = 0.84;
    fNominalParameters[(int)rwgt::fReweightMaCCRES] = 1.12;
    fNominalParameters[(int)rwgt::fReweightMvCCRES] = 0.84;

    fNominalParameters[(int)rwgt::fReweightNormNCRES] = 1.0;
    fNominalParameters[(int)rwgt::fReweightMaNCRESshape] = 1.12;
    fNominalParameters[(int)rwgt::fReweightMvNCRESshape] = 0.84;
    fNominalParameters[(int)rwgt::fReweightMaNCRES] = 1.12;
    fNominalParameters[(int)rwgt::fReweightMvNCRES] = 0.84;

    //Coherent pion nominal values
    fNominalParameters[(int)rwgt::fReweightMaCOHpi] = 1.0;
    fNominalParameters[(int)rwgt::fReweightR0COHpi] = 1.0;


    // Non-resonance background tweaking parameters:
    fNominalParameters[(int)rwgt::fReweightRvpCC1pi] = 1.0;
    fNominalParameters[(int)rwgt::fReweightRvpCC2pi] = 1.0;
    fNominalParameters[(int)rwgt::fReweightRvpNC1pi] = 1.0;
    fNominalParameters[(int)rwgt::fReweightRvpNC2pi] = 1.0;
    fNominalParameters[(int)rwgt::fReweightRvnCC1pi] = 1.0;
    fNominalParameters[(int)rwgt::fReweightRvnCC2pi] = 1.0;
    fNominalParameters[(int)rwgt::fReweightRvnNC1pi] = 1.0;
    fNominalParameters[(int)rwgt::fReweightRvnNC2pi] = 1.0;
    fNominalParameters[(int)rwgt::fReweightRvbarpCC1pi] = 1.0;
    fNominalParameters[(int)rwgt::fReweightRvbarpCC2pi] = 1.0;
    fNominalParameters[(int)rwgt::fReweightRvbarpNC1pi] = 1.0;
    fNominalParameters[(int)rwgt::fReweightRvbarpNC2pi] = 1.0;
    fNominalParameters[(int)rwgt::fReweightRvbarnCC1pi] = 1.0;
    fNominalParameters[(int)rwgt::fReweightRvbarnCC2pi] = 1.0;
    fNominalParameters[(int)rwgt::fReweightRvbarnNC1pi] = 1.0;
    fNominalParameters[(int)rwgt::fReweightRvbarnNC2pi] = 1.0;


    // DIS tweaking parameters - applied for DIS events with [Q2>Q2o, W>Wo], typically Q2okReweight =1GeV^2, WokReweight =1.7-2.0GeV
    fNominalParameters[(int)rwgt::fReweightAhtBY] = 0.538;
    fNominalParameters[(int)rwgt::fReweightBhtBY] = 0.305;
    fNominalParameters[(int)rwgt::fReweightCV1uBY] = 0.291;
    fNominalParameters[(int)rwgt::fReweightCV2uBY] = 0.189;

    fNominalParameters[(int)rwgt::fReweightAhtBYshape] = 0.538;
    fNominalParameters[(int)rwgt::fReweightBhtBYshape] = 0.305;
    fNominalParameters[(int)rwgt::fReweightCV1uBYshape] = 0.291;
    fNominalParameters[(int)rwgt::fReweightCV2uBYshape] = 0.189;

    fNominalParameters[(int)rwgt::fReweightNormDISCC] = 1.0;


    fNominalParameters[(int)rwgt::fReweightRnubarnuCC] = 0.0;//v to vbar ratio reweighting is not working in GENIE at the moment
    fNominalParameters[(int)rwgt::fReweightDISNuclMod] = 0.0;//The DIS nuclear model switch is not working in GENIE at the moment
    //

    fNominalParameters[(int)rwgt::fReweightNC] = 1.0;

    //
    // Hadronization [free nucleon target]
    //
    fNominalParameters[(int)rwgt::fReweightAGKY_xF1pi] = 0.385;
    fNominalParameters[(int)rwgt::fReweightAGKY_pT1pi] = 1./6.625;

    //
    // Medium-effects to hadronization
    //
    fNominalParameters[(int)rwgt::fReweightFormZone] = 1.0;

    //
    // Intranuclear rescattering systematics.
    fNominalParameters[(int)rwgt::fReweightMFP_pi] = 1.0;
    fNominalParameters[(int)rwgt::fReweightMFP_N] = 1.0;
    fNominalParameters[(int)rwgt::fReweightFrCEx_pi] = 1.0;
#ifdef GENIE_PRE_R3
    fNominalParameters[(int)rwgt::fReweightFrElas_pi] = 1.0;
#endif
    fNominalParameters[(int)rwgt::fReweightFrInel_pi] = 1.0;
    fNominalParameters[(int)rwgt::fReweightFrAbs_pi] = 1.0;
    fNominalParameters[(int)rwgt::fReweightFrPiProd_pi] = 1.0;
    fNominalParameters[(int)rwgt::fReweightFrCEx_N] = 1.0;
#ifdef GENIE_PRE_R3
    fNominalParameters[(int)rwgt::fReweightFrElas_N] = 1.0;
#endif
    fNominalParameters[(int)rwgt::fReweightFrInel_N] = 1.0;
    fNominalParameters[(int)rwgt::fReweightFrAbs_N] = 1.0;
    fNominalParameters[(int)rwgt::fReweightFrPiProd_N] = 1.0;


    //
    //RFG Nuclear model
    //
    fNominalParameters[(int)rwgt::fReweightCCQEPauliSupViaKF] = 1.0;
    fNominalParameters[(int)rwgt::fReweightCCQEMomDistroFGtoSF] = 0.0;
    //Continous "switch" at 0.0 full FG model.  At 1.0 full spectral function model.
    //From genie code it looks like weird stuff may happen for <0 and >1.
    //This parameter does not have an "uncertainty" value associated with it.
    //The tweaked dial value gets passed all the way through unchanged to the weight calculator

    //
    // Resonance decays
    //
    fNominalParameters[(int)rwgt::fReweightBR1gamma] = 1.0;
    fNominalParameters[(int)rwgt::fReweightBR1eta] = 1.0;
    fNominalParameters[(int)rwgt::fReweightTheta_Delta2Npi] = 0.0;
    //Continous "switch" at 0.0 full isotropic pion angular distribution.  At 1.0 full R/S pion angular distribtuion.
    //This parameter does not have an "uncertainty" value associated with it.
    //The tweaked dial value gets passed all the way through unchanged to the weight calculator
  }

  ///<Return the nominal value for the given parameter.
  double GENIEReweight::NominalParameterValue(ReweightLabel_t rLabel) {
    double nominal_value;
    nominal_value = fNominalParameters[rLabel];
    return nominal_value;
  }

  ///<Return the configured value of the given parameter.
  double GENIEReweight::ReweightParameterValue(ReweightLabel_t rLabel) {
    int label = (int)rLabel;
    bool in_loop = false;
    bool found_par = false;
    double parameter = -10000;
    for(unsigned int i = 0; i < fReWgtParameterName.size(); i++) {
      in_loop = true;
      if(label == fReWgtParameterName[i]) {
        parameter = fReWgtParameterValue[i];
        found_par = true;
      }
    }
    if(in_loop) {
      if(found_par) return parameter;
      else {
        LOG_WARN("GENIEReweight") << "Parameter " << label << " not set yet or does not exist";
        return parameter;
      }
    }
    else {
      LOG_WARN("GENIEReweight") << "Vector of parameters has not been initialized yet (Size is 0).";
      return parameter;
    }
  }

  ///<Add reweight parameters to the list
  void GENIEReweight::AddReweightValue(ReweightLabel_t rLabel, double value) {
    int label = (int)rLabel;
    LOG_INFO("GENIEReweight") << "Adding parameter: " <<  genie::rew::GSyst::AsString(genie::rew::EGSyst(label)) << ".  With value: " << value;
    fReWgtParameterName.push_back(label);
    fReWgtParameterValue.push_back(value);

  }

  ///<Change a reweight parameter.  If it hasn't been added yet add it.
  void GENIEReweight::ChangeParameterValue(ReweightLabel_t rLabel, double value) {
    int label = (int)rLabel;
    bool found = false;
    for(unsigned int i = 0; i < fReWgtParameterName.size(); i++) {
      if(fReWgtParameterName[i] == label) {
        fReWgtParameterValue[i] = value;
        found = true;
      }
    }
    if(!found) {
      this->AddReweightValue(rLabel, value);
    }
  }

  ///<Configure the weight calculators.
  void GENIEReweight::Configure() {
    LOG_INFO("GENIEReweight") << "Configure weight calculator";

    // these used to decide which calculators need to be configured
    bool configure_NCEL_calc = false,
         configure_CCQE_calc = false,
         configure_QEVec_calc = false,
         configure_ZexpVector_calc = false,
         configure_CCRes_calc = false,
         configure_NCRes_calc = false,
         configure_ResBkg_calc = false,
         configure_ResDecay_calcs = false,
         configure_NC_calc = false,
         configure_DIS_calc = false,
         configure_Coh_calc = false,
         configure_AGKY_calc = false,
         configure_DISNucMod_calc = false,
         configure_FGM_calc = false,
         configure_FZone_calc = false,
         configure_INuke_calc = false,
         configure_MEC_calc = false;

    bool have_MEC_Empirical = false,
         have_MEC_Theory = false;

    for(unsigned int i = 0; i < fReWgtParameterName.size(); i++) {

      switch (fReWgtParameterName[i])
      {
        //NC Elastic parameters
        case rwgt::fReweightMaNCEL:
        case rwgt::fReweightEtaNCEL:
          fReweightNCEL = true;
          break;

        //CC QE Axial parameters
        case rwgt::fReweightNormCCQE:
        case rwgt::fReweightNormCCQEenu:
        case rwgt::fReweightMaCCQEshape:
          configure_CCQE_calc = true;
        case rwgt::fReweightMaCCQE:
          fReweightQEMA = true;
          break;

        //CC QE Vector parameters
          configure_CCQE_calc = true;
        case rwgt::fReweightVecCCQEshape:
          fReweightQEVec = true;
          break;

        //CC Resonance parameters
          configure_ZexpVector_calc = true;
          configure_MEC_calc = true;
          configure_MEC_calc = true;
        case rwgt::fReweightNormCCRES:
        case rwgt::fReweightMaCCRESshape:
        case rwgt::fReweightMvCCRESshape:
          configure_CCRes_calc = true;
        case rwgt::fReweightMaCCRES:
        case rwgt::fReweightMvCCRES:
          fReweightCCRes = true;
          break;

        //NC Resonance parameters
        case rwgt::fReweightNormNCRES:
        case rwgt::fReweightMaNCRESshape:
        case rwgt::fReweightMvNCRESshape:
          configure_NCRes_calc = true;
        case rwgt::fReweightMaNCRES:
        case rwgt::fReweightMvNCRES:
          fReweightNCRes = true;
          configure_NCRes_calc = true;
          break;

        //Coherent parameters
        case rwgt::fReweightMaCOHpi:
        case rwgt::fReweightR0COHpi:
          fReweightCoh = true;
          configure_Coh_calc = true;
          break;

        //Non-resonance background KNO parameters
        case rwgt::fReweightRvpCC1pi:
        case rwgt::fReweightRvpCC2pi:
        case rwgt::fReweightRvpNC1pi:
        case rwgt::fReweightRvpNC2pi:
        case rwgt::fReweightRvnCC1pi:
        case rwgt::fReweightRvnCC2pi:
        case rwgt::fReweightRvnNC1pi:
        case rwgt::fReweightRvnNC2pi:
        case rwgt::fReweightRvbarpCC1pi:
        case rwgt::fReweightRvbarpCC2pi:
        case rwgt::fReweightRvbarpNC1pi:
        case rwgt::fReweightRvbarpNC2pi:
        case rwgt::fReweightRvbarnCC1pi:
        case rwgt::fReweightRvbarnCC2pi:
        case rwgt::fReweightRvbarnNC1pi:
        case rwgt::fReweightRvbarnNC2pi:
          fReweightResBkg = true;
          break;

        //DIS parameters
        case rwgt::fReweightAhtBY:
        case rwgt::fReweightBhtBY:
        case rwgt::fReweightCV1uBY:
        case rwgt::fReweightCV2uBY:
        case rwgt::fReweightAhtBYshape:
        case rwgt::fReweightBhtBYshape:
        case rwgt::fReweightCV1uBYshape:
        case rwgt::fReweightCV2uBYshape:
        case rwgt::fReweightNormDISCC:
          configure_DIS_calc = true;
        case rwgt::fReweightRnubarnuCC:
          configure_DIS_calc = true;
          break;

        //DIS nuclear model parameters
        case rwgt::fReweightDISNuclMod:
          configure_DISNucMod_calc = true;
          break;

        //NC cross section
        case rwgt::fReweightNC:
          configure_NC_calc = true;
          break;

        //Hadronization parameters
        case rwgt::fReweightAGKY_xF1pi:
        case rwgt::fReweightAGKY_pT1pi:
          configure_AGKY_calc = true;
          break;

        //Elastic (and QE) nuclear model parameters
        case rwgt::fReweightCCQEPauliSupViaKF:
        case rwgt::fReweightCCQEMomDistroFGtoSF:
          configure_FGM_calc = true;
          break;

        //Formation Zone
        case rwgt::fReweightFormZone:
          configure_FZone_calc = true;
          break;

        //Intranuke Parameters
        case rwgt::fReweightMFP_pi:
        case rwgt::fReweightMFP_N:
        case rwgt::fReweightFrCEx_pi:
#ifdef GENIE_PRE_R3
        case rwgt::fReweightFrElas_pi:
#endif
        case rwgt::fReweightFrInel_pi:
        case rwgt::fReweightFrAbs_pi:
        case rwgt::fReweightFrPiProd_pi:
        case rwgt::fReweightFrCEx_N:
#ifdef GENIE_PRE_R3
        case rwgt::fReweightFrElas_N:
#endif
        case rwgt::fReweightFrInel_N :
        case rwgt::fReweightFrAbs_N :
        case rwgt::fReweightFrPiProd_N:
          configure_INuke_calc = true;
          break;

        //Resonance Decay parameters
        case rwgt::fReweightBR1gamma:
        case rwgt::fReweightBR1eta:
        case rwgt::fReweightTheta_Delta2Npi:
          fReweightResDecay = true;
          break;

        //Z-expansion parameters
        case rwgt::fReweightZNormCCQE:
        case rwgt::fReweightZExpA1CCQE:
        case rwgt::fReweightZExpA2CCQE:
        case rwgt::fReweightZExpA3CCQE:
        case rwgt::fReweightZExpA4CCQE:
        case rwgt::fReweightAxFFCCQEshape:
          fReweightZexp = true;
          configure_ResDecay_calcs = true;
          break;

      }  // switch(fReWgtParameterName[i])

    } //end for loop

    //configure the individual weight calculators
    if(configure_NCEL_calc) this->ConfigureNCEL();
    if(configure_CCQE_calc) this->ConfigureCCQE();
    if(configure_QEVec_calc) this->ConfigureQEVec();
    if(configure_ZexpVector_calc) this->ConfigureZexpVector();
    if(configure_MEC_calc) this->ConfigureMEC();
    if(configure_CCRes_calc) this->ConfigureCCRes();
    if(configure_NCRes_calc) this->ConfigureNCRes();
    if(configure_ResBkg_calc) this->ConfigureResBkg();
    if(configure_ResDecay_calcs) this->ConfgureResDecay();
    if(configure_NC_calc) this->ConfigureNC();
    if(configure_DIS_calc) this->ConfigureDIS();
    if(configure_Coh_calc) this->ConfigureCoh();
    if(configure_AGKY_calc) this->ConfigureAGKY();
    if(configure_DISNucMod_calc) this->ConfigureDISNucMod();
    if(configure_FGM_calc) this->ConfigureFGM();
    if(configure_FZone_calc) this->ConfigureFZone();
    if(configure_INuke_calc) this->ConfigureINuke();
    this->ConfigureParameters();

  }

  ///<Reconfigure the weight calculators
  void GENIEReweight::Reconfigure() {
    delete fWcalc;
    fWcalc = new genie::rew::GReWeight();
    this->Configure();
  }

  ///<Simple Configuration functions for configuring a single weight calculator

  ///<Simple Configuraiton of the NC elastic weight calculator
  void GENIEReweight::ReweightNCEL(double ma, double eta) {
    LOG_INFO("GENIEReweight") << "Configuring GENIEReweight for NC Elastic Reweighting";
    if(ma!=0.0) {
      this->AddReweightValue(rwgt::fReweightMaNCEL, ma);
    }
    if(eta!=0.0) {
      this->AddReweightValue(rwgt::fReweightEtaNCEL, eta);
    }
    this->Configure();
  }

  ///<Simple Configurtion of the CCQE axial weight calculator
  void GENIEReweight::ReweightQEMA(double ma) {
    LOG_INFO("GENIEReweight") << "Configuring GENIEReweight for QE Axial Mass Reweighting";
    fMaQEshape=false;
    this->AddReweightValue(rwgt::fReweightMaCCQE, ma);
    this->Configure();
  }

  ///<Simple Configuration of the CCQE vector weight calculator
  void GENIEReweight::ReweightQEVec(double mv) {
    LOG_INFO("GENIEReweight") << "Configuring GENIEReweight for QE Vector Mass Reweighting";
    this->AddReweightValue(rwgt::fReweightVecCCQEshape, mv);
    this->Configure();
  }

  void GENIEReweight::ReweightQEZExp(double norm, double a1, double a2, double a3, double a4)
  {
    LOG_INFO("GENIEReweight") << "Configuring GENIEReweight for Z-expansion QE Reweighting";
    this->AddReweightValue(rwgt::fReweightZNormCCQE, norm);
    this->AddReweightValue(rwgt::fReweightZExpA1CCQE, a1);
    this->AddReweightValue(rwgt::fReweightZExpA2CCQE, a2);
    this->AddReweightValue(rwgt::fReweightZExpA3CCQE, a3);
    this->AddReweightValue(rwgt::fReweightZExpA4CCQE, a4);
    this->Configure();
  }

  ///<Simple Configuration of the CC Resonance weight calculator
  void GENIEReweight::ReweightCCRes(double ma, double mv) {
    LOG_INFO("GENIEReweight") << "Configuring GENIEReweight for CC Resonance Reweighting";
    fMaCCResShape=false;
    this->AddReweightValue(rwgt::fReweightMaCCRES, ma);
    if(mv!=0.0) {
      this->AddReweightValue(rwgt::fReweightMvCCRES, mv);
    }
    this->Configure();
  }

  ///<Simple Configurtion of the NC Resonance weight calculator
  void GENIEReweight::ReweightNCRes(double ma, double mv) {
    LOG_INFO("GENIEReweight") << "Configuring GENIEReweight for NC Resonance Reweighting";
    fMaNCResShape=false;
    this->AddReweightValue(rwgt::fReweightMaNCRES, ma);
    if(mv!=0.0) {
      this->AddReweightValue(rwgt::fReweightMvNCRES, mv);
    }
    this->Configure();
  }

  ///<Simple Configuration of the NC and CC Resonance weight calculator with the axial mass parameter for NC/CC ganged together
  void GENIEReweight::ReweightResGanged(double ma, double mv) {
    LOG_INFO("GENIEReweight") << "Configuring GENIEReweight for CC and NC Resonance Reweighting";
    fMaCCResShape=false;
    fMaNCResShape=false;
    this->AddReweightValue(rwgt::fReweightMaCCRES, ma);
    this->AddReweightValue(rwgt::fReweightMaNCRES, ma);
    if(mv!=0.0) {
      this->AddReweightValue(rwgt::fReweightMvCCRES, mv);
      this->AddReweightValue(rwgt::fReweightMvNCRES, mv);
    }
    this->Configure();
  }

  ///<Simple Configuration of the Coherant weight calculator
  void GENIEReweight::ReweightCoh(double ma, double r0) {
    LOG_INFO("GENIEReweight") << "Configuring GENIEReweight for Coherant Reweighting";
    this->AddReweightValue(rwgt::fReweightMaCOHpi, ma);
    this->AddReweightValue(rwgt::fReweightR0COHpi, r0);
    this->Configure();
  }

  ///<Simple Configuration of the Non-Resonance Background weight calculator.
  //Here it is being configured for v+p and vbar + n (1 pi) type interactions
  void GENIEReweight::ReweightNonResRvp1pi(double sigma) {
    LOG_INFO("GENIEReweight") << "Configuring GENIEReweight for  Non-Resonance Background Reweighting (Neutrino Single Pion)";
    this->AddReweightValue(rwgt::fReweightRvpCC1pi, sigma);
    this->AddReweightValue(rwgt::fReweightRvbarnCC1pi, sigma);
    this->AddReweightValue(rwgt::fReweightRvpNC1pi, sigma);
    this->AddReweightValue(rwgt::fReweightRvbarnNC1pi, sigma);
    this->Configure();
  }

  ///<Simple Configuration of the Non-Resonance Background weight calculator.
  //Here it is being configured for v+n and vbar + p (1 pi) type interactions
  void GENIEReweight::ReweightNonResRvbarp1pi(double sigma) {
    LOG_INFO("GENIEReweight") << "Configuring GENIEReweight for  Non-Resonance Background Reweighting (Anti-Neutrino Single Pion)";
    this->AddReweightValue(rwgt::fReweightRvnCC1pi, sigma);
    this->AddReweightValue(rwgt::fReweightRvbarpCC1pi, sigma);
    this->AddReweightValue(rwgt::fReweightRvnNC1pi, sigma);
    this->AddReweightValue(rwgt::fReweightRvbarpNC1pi, sigma);
    this->Configure();
  }

  ///<Simple Configuration of the Non-Resonance Background weight calculator.  Here it is being configured for v+p and vbar + n (2 pi) type interactions
  void GENIEReweight::ReweightNonResRvp2pi(double sigma) {
    LOG_INFO("GENIEReweight") << "Configuring GENIEReweight for  Non-Resonance Background Reweighting (Neutrino Two Pion)";
    this->AddReweightValue(rwgt::fReweightRvpCC2pi, sigma);
    this->AddReweightValue(rwgt::fReweightRvbarnCC2pi, sigma);
    this->AddReweightValue(rwgt::fReweightRvpNC2pi, sigma);
    this->AddReweightValue(rwgt::fReweightRvbarnNC2pi, sigma);
    this->Configure();
  }

  ///<Simple Configuration of the Non-Resonance Background weight calculator.
  // Here it is being configured for v+n and vbar + p (2 pi) type interactions
  void GENIEReweight::ReweightNonResRvbarp2pi(double sigma) {
    LOG_INFO("GENIEReweight") << "Configuring GENIEReweight for  Non-Resonance Background Reweighting (Anti-Neutrino Two Pion)";
    this->AddReweightValue(rwgt::fReweightRvnCC2pi, sigma);
    this->AddReweightValue(rwgt::fReweightRvbarpCC2pi, sigma);
    this->AddReweightValue(rwgt::fReweightRvnNC2pi, sigma);
    this->AddReweightValue(rwgt::fReweightRvbarpNC2pi, sigma);
    this->Configure();
  }

  ///<Simple Configuration of the Resonance decay model weight calculator
  void GENIEReweight::ReweightResDecay(double gamma, double eta, double theta) {
    LOG_INFO("GENIEReweight") << "Configuring GENIEReweight for Resoncance Decay Parameters";
    if(gamma!=0.0) {
      this->AddReweightValue(rwgt::fReweightBR1gamma, gamma);
    }
    if(eta!=0.0) {
      this->AddReweightValue(rwgt::fReweightBR1eta, eta);
    }
    if(theta!=0.0) {
      this->AddReweightValue(rwgt::fReweightTheta_Delta2Npi, theta);
    }
    this->Configure();
  }

  ///<Simple Configuration of the Total NC cross section
  void GENIEReweight::ReweightNC(double norm) {
    LOG_INFO("GENIEReweight") << "Configuring GENIEReweight for NC Cross Section Scale";
    this->AddReweightValue(rwgt::fReweightNC, norm);
    this->Configure();
  }

  ///<Simple Configuration of the DIS FF model weight calculator
  void GENIEReweight::ReweightDIS(double aht, double bht, double cv1u, double cv2u) {
    LOG_INFO("GENIEReweight") << "Configuring GENIEReweight for DIS Form Factor Model Reweighting";
    fDISshape = false;
    if(aht != 0.0) {
      this->AddReweightValue(rwgt::fReweightAhtBY, aht);
    }
    if(bht != 0.0) {
      this->AddReweightValue(rwgt::fReweightBhtBY, bht);
    }
    if(cv1u != 0.0) {
      this->AddReweightValue(rwgt:: fReweightCV1uBY, cv1u);
    }
    if(cv2u != 0.0) {
      this->AddReweightValue(rwgt::fReweightCV2uBY, cv2u);
    }
    this->Configure();
  }

  ///<Simple Configuration of the DIS nuclear model
  void GENIEReweight::ReweightDISnucl(bool mode) {
    LOG_INFO("GENIEReweight") << "Configuring GENIEReweight for DIS Nuclear Model";
    this->AddReweightValue(rwgt::fReweightDISNuclMod, mode);
    this->Configure();
  }

  ///<Simple Configuration of the DIS AGKY hadronization model
  void GENIEReweight::ReweightAGKY(double xF, double pT) {
    LOG_INFO("GENIEReweight") << "Configuring GENIEReweight for DIS AGKY Hadronization Model Reweighting";
    if(xF==0.0) {
      this->AddReweightValue(rwgt::fReweightAGKY_xF1pi, xF);
    }
    if(pT==0.0) {
      this->AddReweightValue(rwgt::fReweightAGKY_pT1pi, pT);
    }
    this->Configure();
  }

  ///<Simple Configuration of the Intranuke Nuclear model
  void GENIEReweight::ReweightIntraNuke(ReweightLabel_t name, double sigma) {
    LOG_INFO("GENIEReweight") << "Configuring GENIEReweight for Intranuke Model Reweighting";
    if ( (name==rwgt::fReweightMFP_pi) ||
        (name==rwgt::fReweightMFP_N) ||
        (name==rwgt::fReweightFrCEx_pi) ||
#ifdef GENIE_PRE_R3
        (name==rwgt::fReweightFrElas_pi) ||
#endif
        (name==rwgt::fReweightFrInel_pi) ||
        (name==rwgt::fReweightFrAbs_pi) ||
        (name==rwgt::fReweightFrPiProd_pi) ||
        (name==rwgt::fReweightFrCEx_N) ||
#ifdef GENIE_PRE_R3
        (name==rwgt::fReweightFrElas_N) ||
#endif
        (name==rwgt::fReweightFrInel_N ) ||
        (name==rwgt::fReweightFrAbs_N ) ||
        (name==rwgt::fReweightFrPiProd_N) ) {
      this->AddReweightValue(name, sigma);
    }
    else {
      LOG_WARN("GENIEReweight") << "That is not a valid Intranuke parameter Intranuke not configured";
    }
    this->Configure();
  }

  ///<Simple Configuration of the Formation Zone reweight calculator
  void GENIEReweight::ReweightFormZone(double sigma) {
    LOG_INFO("GENIEReweight") << "Configuring GENIEReweight for Formation Zone Reweighting";
    this->AddReweightValue(rwgt::fReweightFormZone, sigma);
    this->Configure();
  }

  ///<Simple Configuration of the Fermigas model reweight calculator
  void GENIEReweight::ReweightFGM(double kF, double sf) {
    LOG_INFO("GENIEReweight") << "Configuring GENIEReweight for Fermi Gas Model Reweighting";
    this->AddReweightValue(rwgt::fReweightCCQEPauliSupViaKF, kF);
    this->AddReweightValue(rwgt::fReweightCCQEMomDistroFGtoSF, sf);
    this->Configure();
  }

  ///<End of Simple Reweight Configurations.

  ///<Private Member functions to configure individual weight calculators.
  ///<Configure the NCEL weight calculator
  void GENIEReweight::ConfigureNCEL() {
    LOG_INFO("GENIEReweight") << "Adding NC elastic weight calculator";
    fWcalc->AdoptWghtCalc( "xsec_ncel",       new GReWeightNuXSecNCEL      );
  }

  ///<Configure the MaQE weight calculator
  void GENIEReweight::ConfigureQEMA() {
    LOG_INFO("GENIEReweight") << "Adding CCQE axial FF weight calculator ";
    fWcalc->AdoptWghtCalc( "xsec_ccqe",       new GReWeightNuXSecCCQE      );
  }

  ///<Configure the QE vector FF weight calculator
  void GENIEReweight::ConfigureQEVec() {
    LOG_INFO("GENIEReweight") << "Adding CCQE vector FF weight calculator";
    fWcalc->AdoptWghtCalc( "xsec_ccqe_vec",   new GReWeightNuXSecCCQEvec   );
  }

  ///<Configure the CCRES calculator
  void GENIEReweight::ConfigureCCRes() {
    LOG_INFO("GENIEReweight") << "Adding CC resonance weight calculator";
    fWcalc->AdoptWghtCalc( "xsec_ccres",      new GReWeightNuXSecCCRES     );
    if(!fMaCCResShape) {
      LOG_INFO("GENIEReweight") << "in axial mass (Res) rate+shape mode";
      GReWeightNuXSecCCRES * rwccres = dynamic_cast<GReWeightNuXSecCCRES *> (fWcalc->WghtCalc("xsec_ccres"));
      rwccres->SetMode(GReWeightNuXSecCCRES::kModeMaMv);
    }
    else {
      LOG_INFO("GENIEReweight") << "in axial mass (Res) shape only mode";
    }
  }

  ///<Configure the NCRES calculator
  void GENIEReweight::ConfigureNCRes() {
    LOG_INFO("GENIEReweight") << "Adding NC resonance weight calculator";
    fWcalc->AdoptWghtCalc( "xsec_ncres",      new GReWeightNuXSecNCRES     );
    if(!fMaNCResShape) {
      LOG_INFO("GENIEReweight") << "in axial mass (Res) rate+shape mode";
      GReWeightNuXSecNCRES * rwncres = dynamic_cast<GReWeightNuXSecNCRES *> (fWcalc->WghtCalc("xsec_ncres"));
      rwncres->SetMode(GReWeightNuXSecNCRES::kModeMaMv);
    }
    else {
      LOG_INFO("GENIEReweight") << "in axial mass (Res) shape only mode";
    }
  }

  ///<Configure the ResBkg (kno) weight calculator
  void GENIEReweight::ConfigureResBkg() {
    LOG_INFO("GENIEReweight") << "Adding low Q^2 DIS (KNO) weight calculator";
    fWcalc->AdoptWghtCalc( "xsec_nonresbkg",  new GReWeightNonResonanceBkg );
  }

  ///<Configure the ResDecay weight calculator
  void GENIEReweight::ConfgureResDecay() {
    LOG_INFO("GENIEReweight") << "Adding resonance decay weight calculator";
    fWcalc->AdoptWghtCalc( "hadro_res_decay", new GReWeightResonanceDecay  );
  }

  ///<Configure the NC weight calculator
  void GENIEReweight::ConfigureNC() {
    LOG_INFO("GENIEReweight") << "Adding NC total cross section weight calculator";
    fWcalc->AdoptWghtCalc( "xsec_nc", new GReWeightNuXSecNC );
  }

  ///<Configure the DIS (Bodek-Yang) weight calculator
  void GENIEReweight::ConfigureDIS() {
    LOG_INFO("GENIEReweight") << "Adding DIS (Bodek-Yang) weight calculator";
    fWcalc->AdoptWghtCalc( "xsec_dis",        new GReWeightNuXSecDIS       );
    if(!fDISshape) {
      LOG_INFO("GENIEReweight") << "in shape+rate mode";
      GReWeightNuXSecDIS * rwdis = dynamic_cast<GReWeightNuXSecDIS *> (fWcalc->WghtCalc("xsec_dis"));
      rwdis->SetMode(GReWeightNuXSecDIS::kModeABCV12u);
    }
    else {
      LOG_INFO("GENIEReweight") << "in shape only mode";
    }
  }

  ///<Configure the Coherant model weight calculator
  void GENIEReweight::ConfigureCoh() {
    LOG_INFO("GENIEReweight") << "Adding coherant interaction model weight calculator";
    fWcalc->AdoptWghtCalc( "xsec_coh",        new GReWeightNuXSecCOH       );
  }

  ///<Configure the hadronization (AGKY) weight calculator
  void GENIEReweight::ConfigureAGKY() {
    LOG_INFO("GENIEReweight") << "Adding hadronization (AGKY) model weight calculator";
    fWcalc->AdoptWghtCalc( "hadro_agky",      new GReWeightAGKY            );
  }

  ///<Configure the DIS nuclear model weight calculator
  void GENIEReweight::ConfigureDISNucMod() {
    LOG_INFO("GENIEReweight") << "Adding DIS nuclear model weight calculator";
    fWcalc->AdoptWghtCalc( "nuclear_dis",     new GReWeightDISNuclMod      );
  }

  ///<Configure the FG model weight calculator
  void GENIEReweight::ConfigureFGM() {
    LOG_INFO("GENIEReweight") << "Adding Fermi Gas Model (FGM) weight calculator";
    fWcalc->AdoptWghtCalc( "nuclear_qe",      new GReWeightFGM             );
  }

  ///<Configure the Formation Zone weight calculator
  void GENIEReweight::ConfigureFZone() {
    LOG_INFO("GENIEReweight") << "Adding Formation Zone weight calculator";
    fWcalc->AdoptWghtCalc( "hadro_fzone",     new GReWeightFZone           );
  }

  ///<Configure the intranuke weight calculator
  void GENIEReweight::ConfigureINuke() {
    LOG_INFO("GENIEReweight") << "Adding the Intra-Nuke weight calculator";
    fWcalc->AdoptWghtCalc( "hadro_intranuke", new GReWeightINuke           );
  }

  ///<configure the weight parameters being used
  void GENIEReweight::ConfigureParameters() {
    GSystSet & syst = fWcalc->Systematics();
    for(unsigned int i = 0; i < fReWgtParameterName.size(); i++) {
      LOG_INFO("GENIEReweight") << "Configuring GENIEReweight parameter: " << genie::rew::GSyst::AsString(genie::rew::EGSyst(fReWgtParameterName[i])) << " with value: " << fReWgtParameterValue[i];
      if(fUseSigmaDef) {
        syst.Set( (GSyst_t)fReWgtParameterName[i], fReWgtParameterValue[i]);
      }
      else {
        double parameter = this->CalculateSigma((ReweightLabel_t)fReWgtParameterName[i], fReWgtParameterValue[i]);
        syst.Set( (GSyst_t)fReWgtParameterName[i], parameter);
      }
    }
    fWcalc->Reconfigure();
  }

  ///Used in parameter value mode (instead of parameter sigma mode) Given a user passed parameter value calculate the corresponding sigma value
  ///that needs to be passed to genie to give the same weight.
  double GENIEReweight::CalculateSigma(ReweightLabel_t label, double value) {
    //double GENIEReweight::CalculateSigma(int label, double value) {
    int iLabel = (int) label;
    double nominal_def;
    double frac_err;
    double sigma;
    int sign;
    GSystUncertainty * gsysterr = GSystUncertainty::Instance();
    if(label==rwgt::fReweightCCQEMomDistroFGtoSF ||
        label==rwgt::fReweightTheta_Delta2Npi ||
        label==rwgt::fReweightDISNuclMod) {
      //These parameters don't use the sigma definition just pass them through the function unchanged
      sigma = value;
    }
    else {
      nominal_def = this->NominalParameterValue(label);
      sign = genie::utils::rew::Sign(value-nominal_def);
      frac_err = gsysterr->OneSigmaErr( (GSyst_t)iLabel, sign);
      sigma = (value - nominal_def)/(frac_err*nominal_def);
    }
    return sigma;
  }

  ///<Calculate the weights
  //double GENIEReweight::CalcWeight(simb::MCTruth truth, simb::GTruth gtruth) {
  double GENIEReweight::CalculateWeight(const genie::EventRecord& evr) const {
    //genie::EventRecord evr = this->RetrieveGHEP(truth, gtruth);
    double wgt = fWcalc->CalcWeight(evr);
    //mf::LogVerbatim("GENIEReweight") << "New Event Weight is: " << wgt;
    return wgt;
  }

}
