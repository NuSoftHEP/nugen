
art_make( EXCLUDE EventGeneratorTest_module.cc ### MIGRATE-ACTION-RECOMMENDED: use art_make_library(), art_dictionary(), simple_plugin() with explicit source lists
          MODULE_LIBRARIES nusimdata::SimulationBase
                           nugen::EventGeneratorBase_GENIE
                        art::Framework_Services_Optional_RandomNumberGenerator_service
                        art::Framework_Services_Registry
# for MixGenieEvents_module
#                        art::Framework_Principal
#                        art_root_io_TFileService_service
#                        art::Framework_IO_ProductMix
#                        art_root_io
#                        art_root_io_tfile_support ROOT::Core
# end-of-list for MixGenieEvents_module
                        messagefacility::MF_MessageLogger
                        fhiclcpp::fhiclcpp
                        cetlib::cetlib
                        log4cpp::log4cpp
                        LibXml2::LibXml2
                        Pythia6::Pythia6
                        LHAPDF::LHAPDF
                        ${GENIE_LIB_LIST}
                        GSL::gsl  # FWNUM relies on GSL
                        BLAS::BLAS
                        ROOT::EGPythia6
                        ${ROOT_EVE_LIB_LIST}
                        CLHEP::CLHEP
                        dk2nu::Tree
                        ROOT::Core
                        ROOT::EG
                        ROOT::TreePlayer
                        ROOT::FFTW
                        ROOT::EGPythia6
                        ROOT::Gui
                        )

install_headers()
install_fhicl()
install_source( EXTRAS README )
