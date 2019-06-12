use strict;

use vars qw(%subdir_list);
use vars qw(%header_list);

# explicit headers to avoid conflicts with experiment code
BEGIN { %header_list = (
"Messenger/Messenger.h" => "Framework/Messenger/Messenger.h",
"ReWeight/GReWeight.h" => "RwFramework/GReWeight.h",
"ReWeight/GReWeightI.h" => "RwFramework/GReWeightI.h",
"ReWeight/GSyst.h" => "RwFramework/GSyst.h",
"ReWeight/GSystSet.h" => "RwFramework/GSystSet.h",
"ReWeight/GSystUncertainty.h" => "RwFramework/GSystUncertainty.h",
"ReWeight/GReWeightAGKY.h" => "RwCalculators/GReWeightAGKY.h",
"ReWeight/GReWeightDISNuclMod.h" => "RwCalculators/GReWeightDISNuclMod.h",
"ReWeight/GReWeightFGM.h" => "RwCalculators/GReWeightFGM.h",
"ReWeight/GReWeightFZone.h" => "RwCalculators/GReWeightFZone.h",
"ReWeight/GReWeightINuke.h" => "RwCalculators/GReWeightINuke.h",
"ReWeight/GReWeightINukeParams.h" => "RwCalculators/GReWeightINukeParams.h",
"ReWeight/GReWeightNonResonanceBkg.h" => "RwCalculators/GReWeightNonResonanceBkg.h",
"ReWeight/GReWeightNuXSecCCQEaxial.h" => "RwCalculators/GReWeightNuXSecCCQEaxial.h",
"ReWeight/GReWeightNuXSecCCQE.h" => "RwCalculators/GReWeightNuXSecCCQE.h",
"ReWeight/GReWeightNuXSecCCQEvec.h" => "RwCalculators/GReWeightNuXSecCCQEvec.h",
"ReWeight/GReWeightNuXSecCCRES.h" => "RwCalculators/GReWeightNuXSecCCRES.h",
"ReWeight/GReWeightNuXSecCOH.h" => "RwCalculators/GReWeightNuXSecCOH.h",
"ReWeight/GReWeightNuXSecDIS.h" => "RwCalculators/GReWeightNuXSecDIS.h",
"ReWeight/GReWeightNuXSecHelper.h" => "RwCalculators/GReWeightNuXSecHelper.h",
"ReWeight/GReWeightNuXSecNCEL.h" => "RwCalculators/GReWeightNuXSecNCEL.h",
"ReWeight/GReWeightNuXSecNC.h" => "RwCalculators/GReWeightNuXSecNC.h",
"ReWeight/GReWeightNuXSecNCRES.h" => "RwCalculators/GReWeightNuXSecNCRES.h",
"ReWeight/GReWeightResonanceDecay.h" => "RwCalculators/GReWeightResonanceDecay.h",
"ReWeight/GReWeightUtils.h" => "RwCalculators/GReWeightUtils.h",
"ReWeight/GReWeightIOBranchDesc.h" => "RwIO/GReWeightIOBranchDesc.h",
"ReWeight/GReWeightIORecord.h" => "RwIO/GReWeightIORecord.h",
"GHEP/GHepFlags.h" => "Framework/GHEP/GHepFlags.h",
"GHEP/GHepParticle.h" => "Framework/GHEP/GHepParticle.h",
"GHEP/GHepRecord.h" => "Framework/GHEP/GHepRecord.h",
"GHEP/GHepRecordHistory.h" => "Framework/GHEP/GHepRecordHistory.h",
"GHEP/GHepStatus.h" => "Framework/GHEP/GHepStatus.h",
"GHEP/GHepUtils.h" => "Framework/GHEP/GHepUtils.h",
"GHEP/GHepVirtualListFolder.h" => "Framework/GHEP/GHepVirtualListFolder.h",
"GHEP/GHepVirtualList.h" => "Framework/GHEP/GHepVirtualList.h",
"Interaction/InitialState.h" => "Framework/Interaction/InitialState.h",
"Interaction/InteractionException.h" => "Framework/Interaction/InteractionException.h",
"Interaction/Interaction.h" => "Framework/Interaction/Interaction.h",
"Interaction/InteractionType.h" => "Framework/Interaction/InteractionType.h",
"Interaction/Kinematics.h" => "Framework/Interaction/Kinematics.h",
"Interaction/KPhaseSpace.h" => "Framework/Interaction/KPhaseSpace.h",
"Interaction/ProcessInfo.h" => "Framework/Interaction/ProcessInfo.h",
"Interaction/ScatteringType.h" => "Framework/Interaction/ScatteringType.h",
"Interaction/SppChannel.h" => "Framework/Interaction/SppChannel.h",
"Interaction/Target.h" => "Framework/Interaction/Target.h",
"Interaction/XclsTag.h" => "Framework/Interaction/XclsTag.h",
"Algorithm/AlgCmp.h" => "Framework/Algorithm/AlgCmp.h",
"Algorithm/AlgConfigPool.h" => "Framework/Algorithm/AlgConfigPool.h",
"Algorithm/AlgFactory.h" => "Framework/Algorithm/AlgFactory.h",
"Algorithm/AlgId.h" => "Framework/Algorithm/AlgId.h",
"Algorithm/Algorithm.h" => "Framework/Algorithm/Algorithm.h",
"Algorithm/AlgStatus.h" => "Framework/Algorithm/AlgStatus.h",
"EVGCore/EventGenerator.h" => "Framework/EventGen/EventGenerator.h",
"EVGCore/EventGeneratorI.h" => "Framework/EventGen/EventGeneratorI.h",
"EVGCore/EventGeneratorListAssembler.h" => "Framework/EventGen/EventGeneratorListAssembler.h",
"EVGCore/EventGeneratorList.h" => "Framework/EventGen/EventGeneratorList.h",
"EVGCore/EventRecord.h" => "Framework/EventGen/EventRecord.h",
"EVGCore/EventRecordVisitorI.h" => "Framework/EventGen/EventRecordVisitorI.h",
"EVGCore/EVGThreadException.h" => "Framework/EventGen/EVGThreadException.h",
"EVGCore/GVldContext.h" => "Framework/EventGen/GVldContext.h",
"EVGCore/InteractionGeneratorMap.h" => "Framework/EventGen/InteractionGeneratorMap.h",
"EVGCore/InteractionListAssembler.h" => "Framework/EventGen/InteractionListAssembler.h",
"EVGCore/InteractionListGeneratorI.h" => "Framework/EventGen/InteractionListGeneratorI.h",
"EVGCore/InteractionList.h" => "Framework/EventGen/InteractionList.h",
"EVGCore/InteractionSelectorI.h" => "Framework/EventGen/InteractionSelectorI.h",
"EVGCore/PhysInteractionSelector.h" => "Framework/EventGen/PhysInteractionSelector.h",
"EVGCore/RunningThreadInfo.h" => "Framework/EventGen/RunningThreadInfo.h",
"EVGCore/ToyInteractionSelector.h" => "Framework/EventGen/ToyInteractionSelector.h",
"EVGCore/XSecAlgorithmMap.h" => "Framework/EventGen/XSecAlgorithmMap.h",
"NucleonDecay/DummyInteractionListGenerator.h" => "Physics/NucleonDecay/DummyInteractionListGenerator.h",
"NucleonDecay/DummyPXSec.h" => "Physics/NucleonDecay/DummyPXSec.h",
"NucleonDecay/NucleonDecayMode.h" => "Physics/NucleonDecay/NucleonDecayMode.h",
"NucleonDecay/NucleonDecayPrimaryVtxGenerator.h" => "Physics/NucleonDecay/NucleonDecayPrimaryVtxGenerator.h",
"NucleonDecay/NucleonDecayUtils.h" => "Physics/NucleonDecay/NucleonDecayUtils.h",
"PDG/PDGCodeList.h" => "Framework/ParticleData/PDGCodeList.h",
"PDG/PDGCodes.h" => "Framework/ParticleData/PDGCodes.h",
"PDG/PDGLibrary.h" => "Framework/ParticleData/PDGLibrary.h",
"PDG/PDGUtils.h" => "Framework/ParticleData/PDGUtils.h",
"Utils/AppInit.h" => "Framework/Utils/AppInit.h",
"Utils/RunOpt.h" => "Framework/Utils/RunOpt.h",
"Utils/BWFunc.h" => "Framework/Utils/BWFunc.h",
"Utils/CacheBranchFx.h" => "Framework/Utils/CacheBranchFx.h",
"Utils/CacheBranchI.h" => "Framework/Utils/CacheBranchI.h",
"Utils/CacheBranchNtp.h" => "Framework/Utils/CacheBranchNtp.h",
"Utils/Cache.h" => "Framework/Utils/Cache.h",
"Utils/CmdLnArgParser.h" => "Framework/Utils/CmdLnArgParser.h",
"Utils/ConfigIsotopeMapUtils.h" => "Framework/Utils/ConfigIsotopeMapUtils.h",
"Utils/GSimFiles.h" => "Framework/Utils/GSimFiles.h",
"Utils/GUIUtils.h" => "Framework/Utils/GUIUtils.h",
"Utils/HadXSUtils.h" => "Framework/Utils/HadXSUtils.h",
"Utils/KineUtils.h" => "Framework/Utils/KineUtils.h",
"Utils/PhysUtils.h" => "Framework/Utils/PhysUtils.h",
"Utils/PREM.h" => "Framework/Utils/PREM.h",
"Utils/PrintUtils.h" => "Framework/Utils/PrintUtils.h",
"Utils/Range1.h" => "Framework/Utils/Range1.h",
"Utils/StringUtils.h" => "Framework/Utils/StringUtils.h",
"Utils/Style.h" => "Framework/Utils/Style.h",
"Utils/SystemUtils.h" => "Framework/Utils/SystemUtils.h",
"Utils/T2KEvGenMetaData.h" => "Framework/Utils/T2KEvGenMetaData.h",
"Utils/UnitUtils.h" => "Framework/Utils/UnitUtils.h",
"Utils/XSecSplineList.h" => "Framework/Utils/XSecSplineList.h",
"Utils/XmlParserUtils.h" => "Framework/Utils/XmlParserUtils.h",
"Geo/FidShape.h" => "Tools/Geometry/FidShape.h",
"Geo/GeomVolSelectorBasic.h" => "Tools/Geometry/GeomVolSelectorBasic.h",
"Geo/GeomVolSelectorFiducial.h" => "Tools/Geometry/GeomVolSelectorFiducial.h",
"Geo/GeomVolSelectorI.h" => "Tools/Geometry/GeomVolSelectorI.h",
"Geo/GeomVolSelectorRockBox.h" => "Tools/Geometry/GeomVolSelectorRockBox.h",
"Geo/GeoUtils.h" => "Tools/Geometry/GeoUtils.h",
"Geo/PathSegmentList.h" => "Tools/Geometry/PathSegmentList.h",
"Geo/PointGeomAnalyzer.h" => "Tools/Geometry/PointGeomAnalyzer.h",
"Geo/ROOTGeomAnalyzer.h" => "Tools/Geometry/ROOTGeomAnalyzer.h",
"Conventions/Constants.h" => "Framework/Conventions/Constants.h",
"Conventions/Controls.h" => "Framework/Conventions/Controls.h",
"Conventions/EnvSnapshot.h" => "Framework/Conventions/EnvSnapshot.h",
"Conventions/GBuild.h" => "Framework/Conventions/GBuild.h",
"Conventions/GMode.h" => "Framework/Conventions/GMode.h",
"Conventions/GVersion.h" => "Framework/Conventions/GVersion.h",
"Conventions/KinePhaseSpace.h" => "Framework/Conventions/KinePhaseSpace.h",
"Conventions/KineVar.h" => "Framework/Conventions/KineVar.h",
"Conventions/RefFrame.h" => "Framework/Conventions/RefFrame.h",
"Conventions/Units.h" => "Framework/Conventions/Units.h",
"FluxDrivers/GAstroFlux.h" => "Tools/Flux/GAstroFlux.h",
"FluxDrivers/GAtmoFlux.h" => "Tools/Flux/GAtmoFlux.h",
"FluxDrivers/GBGLRSAtmoFlux.h" => "Tools/Flux/GBGLRSAtmoFlux.h",
"FluxDrivers/GCylindTH1Flux.h" => "Tools/Flux/GCylindTH1Flux.h",
"FluxDrivers/GFlavorMap.h" => "Tools/Flux/GFlavorMap.h",
"FluxDrivers/GFlavorMixerFactory.h" => "Tools/Flux/GFlavorMixerFactory.h",
"FluxDrivers/GFlavorMixerI.h" => "Tools/Flux/GFlavorMixerI.h",
"FluxDrivers/GFLUKAAtmoFlux.h" => "Tools/Flux/GFLUKAAtmoFlux.h",
"FluxDrivers/GFluxBlender.h" => "Tools/Flux/GFluxBlender.h",
"FluxDrivers/GFluxDriverFactory.h" => "Tools/Flux/GFluxDriverFactory.h",
"FluxDrivers/GFluxExposureI.h" => "Tools/Flux/GFluxExposureI.h",
"FluxDrivers/GFluxFileConfigI.h" => "Tools/Flux/GFluxFileConfigI.h",
"FluxDrivers/GHAKKMAtmoFlux.h" => "Tools/Flux/GHAKKMAtmoFlux.h",
"FluxDrivers/GJPARCNuFlux.h" => "Tools/Flux/GJPARCNuFlux.h",
"FluxDrivers/GMonoEnergeticFlux.h" => "Tools/Flux/GMonoEnergeticFlux.h",
"FluxDrivers/GNuMIFlux.h" => "Tools/Flux/GNuMIFlux.h",
"FluxDrivers/GSimpleNtpFlux.h" => "Tools/Flux/GSimpleNtpFlux.h",
"Conventions/XmlParserStatus.h" => "Framework/Conventions/.h",
"NeutronOsc/NeutronOscMode.h" => "Physics/NNBarOscillation/NNBarOscMode.h"
		       ); }

foreach my $inc (sort keys %header_list) {
  s&^(\s*#include\s+["<])\Q$inc\E(.*)&${1}$header_list{$inc}${2}& and last;
  s&^(\s*#include\s+["<]GENIE/)\Q$inc\E(.*)&${1}$header_list{$inc}${2}& and last;
}
