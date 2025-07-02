////////////////////////////////////////////////////////////////////////
///\file ReweightLabels.h
///\brief typedef defining all of the available GENIE reweighting parameters
///
///\author  nathan.mayer@tufts.edu
////////////////////////////////////////////////////////////////////////
#ifndef RWGT_REWEIGHTLABEL_H
#define RWGT_REWEIGHTLABEL_H

#ifdef GENIE_PRE_R3
  #include "ReWeight/GSyst.h"
#else
  #include "RwFramework/GSyst.h"
#endif

namespace rwgt {

  typedef enum EReweightLabel {

    kReweightNull = genie::rew::kNullSystematic,

    // NCEL tweaking parameters:
    fReweightMaNCEL = genie::rew::kXSecTwkDial_MaNCEL,              ///< tweak Ma NCEL, affects dsigma(NCEL)/dQ2 both in shape and normalization
    fReweightEtaNCEL = genie::rew::kXSecTwkDial_EtaNCEL,            ///< tweak NCEL strange axial form factor eta, affects dsigma(NCEL)/dQ2 both in shape and normalization
    // CCQE tweaking parameters:
    fReweightNormCCQE = genie::rew::kXSecTwkDial_NormCCQE,          ///< tweak CCQE normalization (energy independent)
    fReweightNormCCQEenu = genie::rew::kXSecTwkDial_NormCCQEenu,    ///< tweak CCQE normalization (maintains dependence on neutrino energy)
    fReweightMaCCQEshape = genie::rew::kXSecTwkDial_MaCCQEshape,    ///< tweak Ma CCQE, affects dsigma(CCQE)/dQ2 in shape only (normalized to constant integral)
    fReweightMaCCQE = genie::rew::kXSecTwkDial_MaCCQE,              ///< tweak Ma CCQE, affects dsigma(CCQE)/dQ2 both in shape and normalization
    fReweightVecCCQEshape = genie::rew::kXSecTwkDial_VecFFCCQEshape,///< tweak elastic nucleon form factors (BBA/default -> dipole) - shape only effect of dsigma(CCQE)/dQ2
    // Resonance neutrino-production tweaking parameters:
    fReweightNormCCRES = genie::rew::kXSecTwkDial_NormCCRES,         ///< tweak CCRES normalization
    fReweightMaCCRESshape = genie::rew::kXSecTwkDial_MaCCRESshape,   ///< tweak Ma CCRES, affects d2sigma(CCRES)/dWdQ2 in shape only (normalized to constant integral)
    fReweightMvCCRESshape = genie::rew::kXSecTwkDial_MvCCRESshape,   ///< tweak Mv CCRES, affects d2sigma(CCRES)/dWdQ2 in shape only (normalized to constant integral)
    fReweightMaCCRES = genie::rew::kXSecTwkDial_MaCCRES,             ///< tweak Ma CCRES, affects d2sigma(CCRES)/dWdQ2 both in shape and normalization
    fReweightMvCCRES = genie::rew::kXSecTwkDial_MvCCRES,             ///< tweak Mv CCRES, affects d2sigma(CCRES)/dWdQ2 both in shape and normalization

    fReweightNormNCRES = genie::rew::kXSecTwkDial_NormNCRES,         ///< tweak NCRES normalization
    fReweightMaNCRESshape = genie::rew::kXSecTwkDial_MaNCRESshape,   ///< tweak Ma NCRES, affects d2sigma(NCRES)/dWdQ2 in shape only (normalized to constant integral)
    fReweightMvNCRESshape = genie::rew::kXSecTwkDial_MvNCRESshape,   ///< tweak Mv NCRES, affects d2sigma(NCRES)/dWdQ2 in shape only (normalized to constant integral)
    fReweightMaNCRES = genie::rew::kXSecTwkDial_MaNCRES,             ///< tweak Ma NCRES, affects d2sigma(NCRES)/dWdQ2 both in shape and normalization
    fReweightMvNCRES = genie::rew::kXSecTwkDial_MvNCRES,             ///< tweak Mv NCRES, affects d2sigma(NCRES)/dWdQ2 both in shape and normalization

    // Coherent pion production tweaking parameters:
    fReweightMaCOHpi = genie::rew::kXSecTwkDial_MaCOHpi,             ///< tweak Ma for COH pion production
    fReweightR0COHpi = genie::rew::kXSecTwkDial_R0COHpi,             ///< tweak R0 for COH pion production
    // Non-resonance background tweaking parameters:
    fReweightRvpCC1pi = genie::rew::kXSecTwkDial_RvpCC1pi,           ///< tweak the 1pi non-RES bkg in the RES region, for v+p CC
    fReweightRvpCC2pi = genie::rew::kXSecTwkDial_RvpCC2pi,           ///< tweak the 2pi non-RES bkg in the RES region, for v+p CC
    fReweightRvpNC1pi = genie::rew::kXSecTwkDial_RvpNC1pi,           ///< tweak the 1pi non-RES bkg in the RES region, for v+p NC
    fReweightRvpNC2pi = genie::rew::kXSecTwkDial_RvpNC2pi,           ///< tweak the 2pi non-RES bkg in the RES region, for v+p NC
    fReweightRvnCC1pi = genie::rew::kXSecTwkDial_RvnCC1pi,           ///< tweak the 1pi non-RES bkg in the RES region, for v+n CC
    fReweightRvnCC2pi = genie::rew::kXSecTwkDial_RvnCC2pi,           ///< tweak the 2pi non-RES bkg in the RES region, for v+n CC
    fReweightRvnNC1pi = genie::rew::kXSecTwkDial_RvnNC1pi,           ///< tweak the 1pi non-RES bkg in the RES region, for v+n NC
    fReweightRvnNC2pi = genie::rew::kXSecTwkDial_RvnNC2pi,           ///< tweak the 2pi non-RES bkg in the RES region, for v+n NC
    fReweightRvbarpCC1pi = genie::rew::kXSecTwkDial_RvbarpCC1pi,     ///< tweak the 1pi non-RES bkg in the RES region, for vbar+p CC
    fReweightRvbarpCC2pi = genie::rew::kXSecTwkDial_RvbarpCC2pi,     ///< tweak the 2pi non-RES bkg in the RES region, for vbar+p CC
    fReweightRvbarpNC1pi = genie::rew::kXSecTwkDial_RvbarpNC1pi,     ///< tweak the 1pi non-RES bkg in the RES region, for vbar+p NC
    fReweightRvbarpNC2pi = genie::rew::kXSecTwkDial_RvbarpNC2pi,     ///< tweak the 2pi non-RES bkg in the RES region, for vbar+p NC
    fReweightRvbarnCC1pi = genie::rew::kXSecTwkDial_RvbarnCC1pi,     ///< tweak the 1pi non-RES bkg in the RES region, for vbar+n CC
    fReweightRvbarnCC2pi = genie::rew::kXSecTwkDial_RvbarnCC2pi,     ///< tweak the 2pi non-RES bkg in the RES region, for vbar+n CC
    fReweightRvbarnNC1pi = genie::rew::kXSecTwkDial_RvbarnNC1pi,     ///< tweak the 1pi non-RES bkg in the RES region, for vbar+n NC
    fReweightRvbarnNC2pi = genie::rew::kXSecTwkDial_RvbarnNC2pi,     ///< tweak the 2pi non-RES bkg in the RES region, for vbar+n NC
    // DIS tweaking parameters - applied for DIS events with (Q2>Q2o, W>Wo),
    // typically Q2okReweight =1GeV^2, WokReweight =1.7-2.0GeV
    fReweightAhtBY = genie::rew::kXSecTwkDial_AhtBY,                 ///< tweak the Bodek-Yang model parameter A_{ht} - incl. both shape and normalization effect
    fReweightBhtBY = genie::rew::kXSecTwkDial_BhtBY,                 ///< tweak the Bodek-Yang model parameter B_{ht} - incl. both shape and normalization effect
    fReweightCV1uBY = genie::rew::kXSecTwkDial_CV1uBY,               ///< tweak the Bodek-Yang model parameter CV1u - incl. both shape and normalization effect
    fReweightCV2uBY = genie::rew::kXSecTwkDial_CV2uBY,               ///< tweak the Bodek-Yang model parameter CV2u - incl. both shape and normalization effect
    fReweightAhtBYshape = genie::rew::kXSecTwkDial_AhtBYshape,       ///< tweak the Bodek-Yang model parameter A_{ht} - shape only effect to d2sigma(DIS)/dxdy
    fReweightBhtBYshape = genie::rew::kXSecTwkDial_BhtBYshape,       ///< tweak the Bodek-Yang model parameter B_{ht} - shape only effect to d2sigma(DIS)/dxdy
    fReweightCV1uBYshape = genie::rew::kXSecTwkDial_CV1uBYshape,     ///< tweak the Bodek-Yang model parameter CV1u - shape only effect to d2sigma(DIS)/dxdy
    fReweightCV2uBYshape = genie::rew::kXSecTwkDial_CV2uBYshape,     ///< tweak the Bodek-Yang model parameter CV2u - shape only effect to d2sigma(DIS)/dxdy
    fReweightNormDISCC = genie::rew::kXSecTwkDial_NormDISCC,         ///< tweak the inclusive DIS CC normalization (not currently working in genie)
    fReweightRnubarnuCC = genie::rew::kXSecTwkDial_RnubarnuCC,       ///< tweak the ratio of \sigma(\bar\nu CC) / \sigma(\nu CC) (not currently working in genie)
    fReweightDISNuclMod = genie::rew::kXSecTwkDial_DISNuclMod,       ///< tweak DIS nuclear modification (shadowing, anti-shadowing, EMC).  Does not appear to be working in GENIE at the moment
    //
    fReweightNC = genie::rew::kXSecTwkDial_NC,                ///<


    //
    // Hadronization (free nucleon target)
    //

    fReweightAGKY_xF1pi = genie::rew::kHadrAGKYTwkDial_xF1pi,         ///< tweak xF distribution for low multiplicity (N + pi) DIS f/s produced by AGKY
    fReweightAGKY_pT1pi = genie::rew::kHadrAGKYTwkDial_pT1pi,         ///< tweak pT distribution for low multiplicity (N + pi) DIS f/s produced by AGKY


    //
    // Medium-effects to hadronization
    //

    fReweightFormZone = genie::rew::kHadrNuclTwkDial_FormZone,         ///< tweak formation zone


    //
    // Intranuclear rescattering systematics.
    // There are 2 sets of parameters:
    // - parameters that control the total rescattering probability, P(total)
    // - parameters that control the fraction of each process (`fate'), given a total rescat. prob., P(fate|total)
    // These parameters are considered separately for pions and nucleons.
    //

    fReweightMFP_pi = genie::rew::kINukeTwkDial_MFP_pi,           ///< tweak mean free path for pions
    fReweightMFP_N = genie::rew::kINukeTwkDial_MFP_N,             ///< tweak mean free path for nucleons
    fReweightFrCEx_pi = genie::rew::kINukeTwkDial_FrCEx_pi,       ///< tweak charge exchange probability for pions, for given total rescattering probability
#ifdef GENIE_PRE_R3
    fReweightFrElas_pi = genie::rew::kINukeTwkDial_FrElas_pi,     ///< tweak elastic   probability for pions, for given total rescattering probability
#endif
    fReweightFrInel_pi = genie::rew::kINukeTwkDial_FrInel_pi,     ///< tweak inelastic probability for pions, for given total rescattering probability
    fReweightFrAbs_pi = genie::rew::kINukeTwkDial_FrAbs_pi,       ///< tweak absorption probability for pions, for given total rescattering probability
    fReweightFrPiProd_pi = genie::rew::kINukeTwkDial_FrPiProd_pi, ///< tweak pion production probability for pions, for given total rescattering probability
    fReweightFrCEx_N = genie::rew::kINukeTwkDial_FrCEx_N,         ///< tweak charge exchange probability for nucleons, for given total rescattering probability
#ifdef GENIE_PRE_R3
    fReweightFrElas_N = genie::rew::kINukeTwkDial_FrElas_N,       ///< tweak elastic    probability for nucleons, for given total rescattering probability
#endif
    fReweightFrInel_N = genie::rew::kINukeTwkDial_FrInel_N,       ///< tweak inelastic  probability for nucleons, for given total rescattering probability
    fReweightFrAbs_N = genie::rew::kINukeTwkDial_FrAbs_N,         ///< tweak absorption probability for nucleons, for given total rescattering probability
    fReweightFrPiProd_N = genie::rew::kINukeTwkDial_FrPiProd_N,   ///< tweak pion production probability for nucleons, for given total rescattering probability

    //
    // Nuclear model
    //

    fReweightCCQEPauliSupViaKF = genie::rew::kSystNucl_CCQEPauliSupViaKF,   ///<
    fReweightCCQEMomDistroFGtoSF = genie::rew::kSystNucl_CCQEMomDistroFGtoSF, ///<

    //
    // Resonance decays
    //

    fReweightBR1gamma = genie::rew::kRDcyTwkDial_BR1gamma,               ///< tweak Resonance -> X + gamma branching ratio, eg Delta+(1232) -> p gamma
    fReweightBR1eta = genie::rew::kRDcyTwkDial_BR1eta,                   ///< tweak Resonance -> X + eta   branching ratio, eg N+(1440) -> p eta
    fReweightTheta_Delta2Npi = genie::rew::kRDcyTwkDial_Theta_Delta2Npi,  ///< distort pi angular distribution in Delta -> N + pi

    //
    // Alternative approach to CCQE form factors (z-expansion)
    //

    fReweightZNormCCQE = genie::rew::kXSecTwkDial_ZNormCCQE,         ///< tweak Z-expansion CCQE normalization (energy independent)
    fReweightZExpA1CCQE = genie::rew::kXSecTwkDial_ZExpA1CCQE,       ///< tweak Z-expansion coefficient 1, affects dsigma(CCQE)/dQ2 both in shape and normalization
    fReweightZExpA2CCQE = genie::rew::kXSecTwkDial_ZExpA2CCQE,       ///< tweak Z-expansion coefficient 2, affects dsigma(CCQE)/dQ2 both in shape and normalization
    fReweightZExpA3CCQE = genie::rew::kXSecTwkDial_ZExpA3CCQE,       ///< tweak Z-expansion coefficient 3, affects dsigma(CCQE)/dQ2 both in shape and normalization
    fReweightZExpA4CCQE = genie::rew::kXSecTwkDial_ZExpA4CCQE,       ///< tweak Z-expansion coefficient 4, affects dsigma(CCQE)/dQ2 both in shape and normalization
    fReweightAxFFCCQEshape = genie::rew::kXSecTwkDial_AxFFCCQEshape, ///< tweak axial nucleon form factors (dipole -> z-expansion) - shape only effect of dsigma(CCQE)/dQ2

    //
    //   Alternative approach to CCQE form factors (RunningMA)
    //
    fReweightE0CCQEshape = genie::rew::kXSecTwkDial_E0CCQEshape,  ///< tweak E0 CCQE RunningMA, affects dsigma(CCQE)/dQ2 in shape only (normalized to constant integral)
    fReweightE0CCQE = genie::rew::kXSecTwkDial_E0CCQE,            ///< tweak E0 CCQE RunningMA, affects dsigma(CCQE)/dQ2 both in shape and normalization

    //
    // Empirical MEC dials
    //
    fReweightEmpMEC_Mq2d = genie::rew::kXSecTwkDial_EmpMEC_Mq2d,           ///< tweak the "mass^2" parameter in the Q^2 component of Empirical MEC's form factor
    fReweightEmpMEC_Mass = genie::rew::kXSecTwkDial_EmpMEC_Mass,           ///< tweak the "mass" parameter in the W component of Empirical MEC's form factor
    fReweightEmpMEC_Width = genie::rew::kXSecTwkDial_EmpMEC_Width,         ///< tweak the "width" parameter in the W component of Empirical MEC's form factor
    fReweightEmpMEC_FracPN_NC = genie::rew::kXSecTwkDial_EmpMEC_FracPN_NC, ///< tweak fraction of initial-state NP (vs. NN or PP) pairs in NC reactions
    fReweightEmpMEC_FracPN_CC = genie::rew::kXSecTwkDial_EmpMEC_FracPN_CC, ///< tweak fraction of initial-state NP (vs. NN or PP) pairs in CC reactions
    fReweightEmpMEC_FracCCQE = genie::rew::kXSecTwkDial_EmpMEC_FracCCQE,   ///< tweak the overall rate as a function of the QE rate for CC reactions
    fReweightEmpMEC_FracNCQE = genie::rew::kXSecTwkDial_EmpMEC_FracNCQE,   ///< tweak the overall rate as a function of the QE rate for CC reactions
    fReweightEmpMEC_FracPN_EM = genie::rew::kXSecTwkDial_EmpMEC_FracPN_EM, ///< tweak fraction of initial-state NP (vs. NN or PP) pairs in EM (electron-probe) reactions
    fReweightEmpMEC_FracEMQE = genie::rew::kXSecTwkDial_EmpMEC_FracEMQE,   ///< tweak the overall rate as a function of the QE rate for EM (electron-probe) reactions

    //
	// General MEC dials
	//
    fReweightNormCCMEC = genie::rew::kXSecTwkDial_NormCCMEC,               ///< tweak the overall normalization of MEC events in CC reactions
    fReweightNormNCMEC = genie::rew::kXSecTwkDial_NormNCMEC,               ///< tweak the overall normalization of MEC events in NC reactions
    fReweightNormEMMEC = genie::rew::kXSecTwkDial_NormEMMEC,               ///< tweak the overall normalization of MEC events in EM (electron-probe) reactions
    fReweightDecayAngMEC = genie::rew::kXSecTwkDial_DecayAngMEC,           ///< tweak the (reaction frame) angular distribution between outgoing nucleons in MEC events
    fReweightFracPN_CCMEC = genie::rew::kXSecTwkDial_FracPN_CCMEC,         ///< tweak fraction of initial-state NP (vs. NN or PP) pairs in CC reactions (see fXSecTwkDial_EmpMEC_FracPN_CC for Empirical MEC)
    fReweightFracDelta_CCMEC = genie::rew::kXSecTwkDial_FracDelta_CCMEC,   ///< tweak the fraction of CCMEC events involving an internal delta line (currently only used by Valencia MEC)
  	fReweightXSecShape_CCMEC = genie::rew::kXSecTwkDial_XSecShape_CCMEC,   ///< tweak the shape of CCMEC differential cross section (interpolates between models)

  	fReweightRPA_CCQE = genie::rew::kXSecTwkDial_RPA_CCQE,                 ///< interpolate between the default CCQE model and the same one with RPA off (only gives non-unit weights for Nieves CCQE)

    fReweightTheta_Delta2NRad = genie::rew::kRDcyTwkDial_Theta_Delta2NRad, ///< Distort photon angular distribution in Delta -> N + photon

  	fReweightCoulombCCQE = genie::rew::kXSecTwkDial_CoulombCCQE,           ///< Tweak the value of the EM potential used when computing the Coulomb correction factor in the Nieves CCQE model

    fReweightNormCCCOHpi = genie::rew::kXSecTwkDial_NormCCCOHpi,           ///< Scale the normalization of CC coherent pion production
    fReweightNormNCCOHpi = genie::rew::kXSecTwkDial_NormNCCOHpi,           ///< Scale the normalization of NC coherent pion production

    //
    // Alternative approach to CCQE form factors (z-expansion) vector form factor
    //
    fReweightZExpELFF = genie::rew::kXSecTwkDial_ZExpELFF,
    fReweightZExpELFF_AP1 = genie::rew::kXSecTwkDial_ZExpELFF_AP1,
    fReweightZExpELFF_AP2 = genie::rew::kXSecTwkDial_ZExpELFF_AP2,
    fReweightZExpELFF_AP3 = genie::rew::kXSecTwkDial_ZExpELFF_AP3,
    fReweightZExpELFF_AP4 = genie::rew::kXSecTwkDial_ZExpELFF_AP4,
    fReweightZExpELFF_AN1 = genie::rew::kXSecTwkDial_ZExpELFF_AN1,
    fReweightZExpELFF_AN2 = genie::rew::kXSecTwkDial_ZExpELFF_AN2,
    fReweightZExpELFF_AN3 = genie::rew::kXSecTwkDial_ZExpELFF_AN3,
    fReweightZExpELFF_AN4 = genie::rew::kXSecTwkDial_ZExpELFF_AN4,
    fReweightZExpELFF_BP1 = genie::rew::kXSecTwkDial_ZExpELFF_BP1,
    fReweightZExpELFF_BP2 = genie::rew::kXSecTwkDial_ZExpELFF_BP2,
    fReweightZExpELFF_BP3 = genie::rew::kXSecTwkDial_ZExpELFF_BP3,
    fReweightZExpELFF_BP4 = genie::rew::kXSecTwkDial_ZExpELFF_BP4,
    fReweightZExpELFF_BN1 = genie::rew::kXSecTwkDial_ZExpELFF_BN1,
    fReweightZExpELFF_BN2 = genie::rew::kXSecTwkDial_ZExpELFF_BN2,
    fReweightZExpELFF_BN3 = genie::rew::kXSecTwkDial_ZExpELFF_BN3,
    fReweightZExpELFF_BN4 = genie::rew::kXSecTwkDial_ZExpELFF_BN4,


    kNTwkDials, /// < Not a real dial, just keep as last entry for looping purposes

  } ReweightLabel_t;

} // end rwgt namespace
#endif //RWGT_REWEIGHTLABEL_H
