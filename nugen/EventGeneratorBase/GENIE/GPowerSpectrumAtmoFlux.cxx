#include <cmath>
#include <fstream>

#include <TMath.h>
#include <TLorentzVector.h>

#include "nugen/EventGeneratorBase/GENIE/GPowerSpectrumAtmoFlux.h"
#include "GENIE/Framework/Numerical/RandomGen.h"
#include "GENIE/Framework/Conventions/Constants.h"
#include "GENIE/Framework/Utils/PrintUtils.h"
#include "GENIE/Tools/Flux/GFluxDriverFactory.h"
#include "GENIE/Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGCodeList.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "TFile.h"
#include "TH3D.h"
#include "TH2.h"

FLUXDRIVERREG4(genie,flux,GPowerSpectrumAtmoFlux,genie::flux::GPowerSpectrumAtmoFlux)

using namespace genie;
using namespace genie::flux;
using namespace genie::constants;

GPowerSpectrumAtmoFlux::GPowerSpectrumAtmoFlux()
{
	LOG("Flux", pNOTICE)
	<< "Instantiating the GENIE Power Spectrum atmospheric neutrino flux driver";
	this->Initialize();
}

//________________________________________________________________________________________

GPowerSpectrumAtmoFlux::~GPowerSpectrumAtmoFlux()
{

}


//__________________________________________________________________________________________________

const PDGCodeList &GPowerSpectrumAtmoFlux::FluxParticles (void)
{
	return *fPdgCList;
}

//_________________________________________________________________

double GPowerSpectrumAtmoFlux::MaxEnergy(void)
{
  return fMaxEvCut;
}

//__________________________________________________________________

double GPowerSpectrumAtmoFlux::MinEnergy(void)
{
  return fMinEvCut;
}

//_________________________________________________________________________

bool GPowerSpectrumAtmoFlux::GenerateNext(void)
{
	// Reset previously generated neutrino code / 4-p / 4-x
	this->ResetSelection();

	// Get a RandomGen instance
	RandomGen * rnd = RandomGen::Instance();

	// Generate (Ev, costheta, phi)
	double Ev       = 0.;
	double costheta = 0.;
	double phi      = 0;
	double weight   = 0;
	int    nu_pdg   = 0;

	// generate events according to a power law spectrum,
	// then weight events by inverse power law
	// (note: cannot use index alpha=1)
	double alpha = fSpectralIndex;

	double emin = TMath::Power(this->MinEnergy(),1.0-alpha);
	double emax = TMath::Power(this->MaxEnergy(),1.0-alpha);
	Ev          = TMath::Power(emin+(emax-emin)*rnd->RndFlux().Rndm(),1.0/(1.0-alpha));
	costheta    = -1+2*rnd->RndFlux().Rndm();
	phi         = 2.*kPi* rnd->RndFlux().Rndm();

	unsigned int nnu = fPdgCList->size();
	unsigned int inu = rnd->RndFlux().Integer(nnu);
	nu_pdg   = (*fPdgCList)[inu];
	weight = this->ComputeWeight(nu_pdg, Ev, costheta, phi);

	// Compute etc trigonometric numbers
	double sintheta  = TMath::Sqrt(1-costheta*costheta);
	double cosphi    = TMath::Cos(phi);
	double sinphi    = TMath::Sin(phi);

	// Set the neutrino pdg code
	fgPdgC = nu_pdg;

	// Set the neutrino weight
	fWeight = weight;

	// Compute the neutrino momentum
	// The `-1' means it is directed towards the detector.
	double pz = -1.* Ev * costheta;
	double py = -1.* Ev * sintheta * sinphi;
	double px = -1.* Ev * sintheta * cosphi;

	// Default vertex is at the origin
	double z = 0.0;
	double y = 0.0;
	double x = 0.0;

	// Shift the neutrino position onto the flux generation surface.
	// The position is computed at the surface of a sphere with R=fRl
	// at the topocentric horizontal (THZ) coordinate system.
	if( fRl>0.0 ){
	z += fRl * costheta;
	y += fRl * sintheta * sinphi;
	x += fRl * sintheta * cosphi;
	}

	// Apply user-defined rotation from THZ -> user-defined topocentric
	// coordinate system.
	if( !fRotTHz2User.IsIdentity() )
	{
	TVector3 tx3(x, y, z );
	TVector3 tp3(px,py,pz);

	tx3 = fRotTHz2User * tx3;
	tp3 = fRotTHz2User * tp3;

	x  = tx3.X();
	y  = tx3.Y();
	z  = tx3.Z();
	px = tp3.X();
	py = tp3.Y();
	pz = tp3.Z();
	}

	// If the position is left as is, then all generated neutrinos
	// would point towards the origin.
	// Displace the position randomly on the surface that is
	// perpendicular to the selected point P(xo,yo,zo) on the sphere
	if( fRt>0.0 ){
	TVector3 vec(x,y,z);               // vector towards selected point
	TVector3 dvec1 = vec.Orthogonal(); // orthogonal vector
	TVector3 dvec2 = dvec1;            // second orthogonal vector
	dvec2.Rotate(-kPi/2.0,vec);        // rotate second vector by 90deg,
	                                   // now forming a new orthogonal cartesian coordinate system
	double psi = 2.*kPi* rnd->RndFlux().Rndm(); // rndm angle [0,2pi]
	double random = rnd->RndFlux().Rndm();      // rndm number  [0,1]
	dvec1.SetMag(TMath::Sqrt(random)*fRt*TMath::Cos(psi));
	dvec2.SetMag(TMath::Sqrt(random)*fRt*TMath::Sin(psi));
	x += dvec1.X() + dvec2.X();
	y += dvec1.Y() + dvec2.Y();
	z += dvec1.Z() + dvec2.Z();
	}

	// Set the neutrino momentum and position 4-vectors with values
	// calculated at previous steps.
	fgP4.SetPxPyPzE(px, py, pz, Ev);
	fgX4.SetXYZT   (x,  y,  z,  0.);

	// Increment flux neutrino counter used for sample normalization purposes.
	fNNeutrinos++;

	// Report and exit
	LOG("Flux", pINFO)
	   << "Generated neutrino: "
	   << "\n pdg-code: " << fgPdgC
	   << "\n p4: " << utils::print::P4AsShortString(&fgP4)
	   << "\n x4: " << utils::print::X4AsString(&fgX4);

	return true;
}

//________________________________________________________________

void GPowerSpectrumAtmoFlux::ResetSelection(void)
{
// initializing running neutrino pdg-code, 4-position, 4-momentum

  fgPdgC = 0;
  fgP4.SetPxPyPzE (0.,0.,0.,0.);
  fgX4.SetXYZT    (0.,0.,0.,0.);
  fWeight = 0;
}

//________________________________________________________________

int GPowerSpectrumAtmoFlux::PdgCode(void)
{
	return fgPdgC;
}

//________________________________________________________________

double GPowerSpectrumAtmoFlux::Weight(void)
{
	return fWeight;
}

//________________________________________________________________

const TLorentzVector& GPowerSpectrumAtmoFlux::Momentum(void)
{
	return fgP4;
}

//________________________________________________________________

const TLorentzVector& GPowerSpectrumAtmoFlux::Position(void)
{
	return fgX4;
}

//________________________________________________________________

bool GPowerSpectrumAtmoFlux::End(void)
{
	return false;
}

//________________________________________________________________

long int GPowerSpectrumAtmoFlux::Index(void)
{
	return -1;
}

//________________________________________________________________

void GPowerSpectrumAtmoFlux::ResetNFluxNeutrinos(void)
{
	fNNeutrinos = 0;
}

//________________________________________________________________

void GPowerSpectrumAtmoFlux::Clear(Option_t * opt)
{
	LOG("Flux", pWARN) << "GPowerSpectrumAtmoFlux::Clear(" << opt << ") called";
	this->ResetNFluxNeutrinos();
}

//________________________________________________________________________

void GPowerSpectrumAtmoFlux::GenerateWeighted(bool gen_weighted)
{
// Dummy clear method needed to conform to GFluxI interface
//
  LOG("Flux", pERROR) <<
  "The neutrinos are always generated weighted with this implementation!! GenerateWeighted method not implemented.";
}

//______________________________________________________________________________________________________________________

void GPowerSpectrumAtmoFlux::SetSpectralIndex(double index)
{
	if( index != 1.0 ){
		fSpectralIndex = index;
	}
	else {
		LOG("Flux", pWARN) << "Warning: cannot use a spectral index of unity";
	}
	LOG("Flux", pNOTICE) << "Using Spectral Index = " << index;
}

//____________________________________________________________________________________

void GPowerSpectrumAtmoFlux::SetUserCoordSystem(TRotation &rotation)
{
  fRotTHz2User = rotation;
}

//__________________________________________________________________________________

void GPowerSpectrumAtmoFlux::Initialize(void)
{
	LOG("Flux", pNOTICE) << "Initializing power spectrum atmospheric flux driver";

	bool allow_dup = false;
	fPdgCList = new PDGCodeList(allow_dup);

	fSpectralIndex = 2.0;

	fMinEvCut = 0.01;
	fMaxEvCut = 9999999999;

	// Default radii
	fRl = 0.0;
	fRt = 0.0;
	fAgen = 0.0;

	// Default detector coord system: Topocentric Horizontal Coordinate system
	fRotTHz2User.SetToIdentity();

	// Reset `current' selected flux neutrino
	this->ResetSelection();

	// Reset number of neutrinos thrown so far
	fNNeutrinos = 0;

	// Init the global gen weight
	this->InitializeWeight();
}

//_________________________________________________________________________________

void GPowerSpectrumAtmoFlux::SetRadii(double Rlongitudinal, double Rtransverse)
{
  LOG ("Flux", pNOTICE) << "Setting R[longitudinal] = " << Rlongitudinal;
  LOG ("Flux", pNOTICE) << "Setting R[transverse]   = " << Rtransverse;

  fRl = Rlongitudinal;
  fRt = Rtransverse;

  fAgen = kPi*fRt*fRt;
}

//________________________________________________________________________________

void GPowerSpectrumAtmoFlux::SetFlavors(std::vector<int> flavors)
{
	fPdgCList->clear();
	for(int flavor : flavors){
		fPdgCList->push_back(flavor);
	}
}

//________________________________________________________________________________

void GPowerSpectrumAtmoFlux::CleanUp(void)
{
	LOG("Flux", pNOTICE) << "Cleaning up...";

	if (fPdgCList) delete fPdgCList;
}

//________________________________________________________________________________

void GPowerSpectrumAtmoFlux::SetMinEnergy(double Emin)
{
	if(Emin > 0){
		fMinEvCut = Emin;
	}
}

//_________________________________________________________________________________

void GPowerSpectrumAtmoFlux::SetMaxEnergy(double Emax){
	fMaxEvCut = Emax;
}

//_________________________________________________________________________________

long int GPowerSpectrumAtmoFlux::NFluxNeutrinos(void) const
{
  return fNNeutrinos;
}

//_________________________________________________________________________________

bool GPowerSpectrumAtmoFlux::FillFluxHisto(int nu_pdg, string filename){
	LOG("Flux", pNOTICE)
    << "Loading fine grained flux for neutrino: " << nu_pdg
    << " from file: " << filename;

	TH3D* histo = nullptr;

	TFile *f = new TFile(filename.c_str(), "READ");
	f->GetObject("flux", histo);

	if(!histo) {
		LOG("Flux", pERROR) << "Null flux histogram!";
		f->Close();
		return false;
	}

	TH3D* h = static_cast<TH3D*>(histo->Clone());
	f->Close();

	fRawFluxHistoMap.insert(std::make_pair(nu_pdg, h));

	TH2D* h2 = static_cast<TH2D*>(h->Project3D("yx"));

	fRawFluxHistoMap2D.insert(std::make_pair(nu_pdg, h2));

	return true;
}

//_________________________________________________________________________________

void GPowerSpectrumAtmoFlux::AddFluxFile(int nu_pdg, string filename)
{
  // Check file exists
  std::ifstream f(filename.c_str());
  if (!f.good()) {
    LOG("Flux", pFATAL) << "Flux file does not exist "<<filename;
    exit(-1);
  }
  if ( pdg::IsNeutrino(nu_pdg) || pdg::IsAntiNeutrino(nu_pdg) ) {
    fFluxFlavour.push_back(nu_pdg);
    fFluxFile.push_back(filename);
  } else {
    LOG ("Flux", pWARN)
      << "Input particle code: " << nu_pdg << " not a neutrino!";
  }
}

//_________________________________________________________________________________

bool GPowerSpectrumAtmoFlux::LoadFluxData(void)
{
  LOG("Flux", pNOTICE)
        << "Loading atmospheric neutrino flux simulation data";

  fPdgCList->clear();

  bool loading_status = true;

  for( unsigned int n=0; n<fFluxFlavour.size(); n++ ){
    int nu_pdg      = fFluxFlavour.at(n);
    string filename = fFluxFile.at(n);
    string pname = PDGLibrary::Instance()->Find(nu_pdg)->GetName();

    LOG("Flux", pNOTICE) << "Loading data for: " << pname;

    bool loaded = this->FillFluxHisto(nu_pdg, filename);

    loading_status = loading_status && loaded;

    if (!loaded) {
        LOG("Flux", pERROR)
          << "Error loading atmospheric neutrino flux simulation data from " << filename;
        break;
    }
  }

  if(loading_status) {
    map<int,TH3D*>::iterator hist_iter = fRawFluxHistoMap.begin();
    for ( ; hist_iter != fRawFluxHistoMap.end(); ++hist_iter) {
      int   nu_pdg = hist_iter->first;
      fPdgCList->push_back(nu_pdg);
    }

    LOG("Flux", pNOTICE)
          << "Atmospheric neutrino flux simulation data loaded!";
    return true;
  }

  LOG("Flux", pERROR)
    << "Error loading atmospheric neutrino flux simulation data";
  return false;
}


//_________________________________________________________________________

double GPowerSpectrumAtmoFlux::GetFlux(int flavour, double energy, double costh, double phi)
{
  TH3D* flux_hist = nullptr;
  std::map<int,TH3D*>::iterator it = fRawFluxHistoMap.find(flavour);
  if(it != fRawFluxHistoMap.end())
  {
    flux_hist = it->second;
  }

  LOG("Flux", pERROR) << "flux_hist: " << flux_hist;

  if(!flux_hist) return 0.0;

  if(flux_hist->GetZaxis()->GetNbins() == 1){ //no binning in phi, bilinear interpolation only so using the 2D hist
	TH2D* h2 = fRawFluxHistoMap2D[it->first];
	return h2->Interpolate(energy, costh);
  }
  else {
	return flux_hist->Interpolate(energy, costh, phi);
  }
}

//_________________________________________________________________________

double GPowerSpectrumAtmoFlux::ComputeWeight(int flavour, double energy, double costh, double phi)
{
  double flux = this->GetFlux(flavour, energy, costh, phi);

  return flux*fAgen*fGlobalGenWeight*pow(energy, fSpectralIndex);
}

//_________________________________________________________________________

void GPowerSpectrumAtmoFlux::InitializeWeight()
{
	LOG("Flux", pNOTICE) << "Initializing generation weight";

	double IE = 1.;

	if(fSpectralIndex!=1){
		IE = (pow(fMaxEvCut,(1.-fSpectralIndex))-pow(fMinEvCut,(1.-fSpectralIndex)))/(1.-fSpectralIndex);
	}

	double ITheta = 4*kPi;

	fGlobalGenWeight=IE*ITheta;
}