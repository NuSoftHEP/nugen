//____________________________________________________________________________
/*!

\class   genie::flux::GPowerSpectrumAtmoFlux

\brief   A driver for a power spectrum atmospheric neutrino flux.
		 Elements extensively reused from GAtmoFlux.

\author  Pierre Granger <granger@apc.in2p3.fr>
         APC (CNRS)

\created December 16, 2022

*/
//____________________________________________________________________________

#pragma once

#include <TLorentzVector.h>

#include "GENIE/Tools/Flux/GAtmoFlux.h"
#include "GENIE/Framework/ParticleData/PDGCodeList.h"
#include "TH2D.h"

namespace genie {
namespace flux {

class GPowerSpectrumAtmoFlux: public GFluxI {
public:
	GPowerSpectrumAtmoFlux();
	~GPowerSpectrumAtmoFlux();

	const PDGCodeList &FluxParticles(void) override;  ///< declare list of flux neutrinos that can be generated (for init. purposes)
	double MaxEnergy(void) override; ///< declare the max flux neutrino energy that can be generated (for init. purposes)
	bool GenerateNext(void) override; ///< generate the next flux neutrino (return false in err)
	int PdgCode(void) override; ///< returns the flux neutrino pdg code
	double Weight(void) override; ///< returns the flux neutrino weight (if any)
	const TLorentzVector& Momentum(void) override; ///< returns the flux neutrino 4-momentum
	const TLorentzVector& Position(void) override; ///< returns the flux neutrino 4-position (note: expect SI rather than physical units)
	bool End(void) override; ///< true if no more flux nu's can be thrown (eg reaching end of beam sim ntuples)
	long int Index(void) override; ///< returns corresponding index for current flux neutrino (e.g. for a flux ntuple returns the current entry number)
	void Clear(Option_t *opt) override; ///< reset state variables based on opt
	void GenerateWeighted(bool gen_weighted) override; ///< set whether to generate weighted or unweighted neutrinos

	double MinEnergy(void);
	void SetSpectralIndex(double index);
	void SetUserCoordSystem(TRotation &rotation);
	void SetRadii(double Rlongitudinal, double Rtransverse);
	void SetFlavors(std::vector<int> flavors);
	void SetMinEnergy(double Emin);
	void SetMaxEnergy(double Emax);
	long int NFluxNeutrinos(void) const;
	void ResetNFluxNeutrinos(void);
	void AddFluxFile(int neutrino_pdg, string filename);
	bool LoadFluxData(void);

	double GetFlux(int flavour, double energy, double costh, double phi);
	double ComputeWeight(int flavour, double energy, double costh, double phi);
	void InitializeWeight();


private:
	PDGCodeList *fPdgCList; ///< input list of neutrino pdg-codes
	double fMaxEvCut; ///< user-defined cut: maximum energy
	double fMinEvCut; ///< user-defined cut: minimum energy
	int fgPdgC; ///< current generated nu pdg-code
	TLorentzVector fgP4; ///< current generated nu 4-momentum
	TLorentzVector fgX4; ///< current generated nu 4-position
	double fWeight; ///< current generated nu weight
	double fSpectralIndex; ///< power law function
	double fRl; ///< defining flux neutrino generation surface: longitudinal radius
	double fRt; ///< defining flux neutrino generation surface: transverse radius
	double fGlobalGenWeight; ///< global generation weight to apply to all events
	double fAgen; ///< current generation area
	TRotation fRotTHz2User; ///< coord. system rotation: THZ -> Topocentric user-defined
	long int fNNeutrinos; ///< number of flux neutrinos thrown so far
	vector<int> fFluxFlavour; ///< input flux file for each neutrino species
    vector<string> fFluxFile; ///< input flux file for each neutrino species
  	map<int, TH3D*> fRawFluxHistoMap; ///< flux = f(Ev,cos8,phi) for each neutrino species
  	map<int, TH2D*> fRawFluxHistoMap2D; ///< flux = f(Ev,cos8) for each neutrino species

	bool FillFluxHisto(int nu_pdg, string filename);
	void AddAllFluxes(void);
	TH3D* CreateNormalisedFluxHisto( TH3D* hist);  // normalise flux files
	void ResetSelection(void);
	void Initialize(void);
	void CleanUp(void);
};


} // flux namespace
} // genie namespace