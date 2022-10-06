
cet_build_plugin(AddGenieEventsToArt art::module
                 LIBRARIES REG nusimdata::SimulationBase
                        nugen::EventGeneratorBase_GENIE
                        art::Framework_Services_Optional_RandomNumberGenerator_service
                        art::Framework_Services_Registry
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
                        ROOT::MathMore
                        ROOT::EG
                        ROOT::TreePlayer
                        ROOT::FFTW
                        ROOT::EGPythia6
                        ROOT::Gui
                        )
cet_build_plugin(GenieOutput art::module
                 LIBRARIES REG nusimdata::SimulationBase
                        nugen::EventGeneratorBase_GENIE
                        art::Framework_Services_Optional_RandomNumberGenerator_service
                        art::Framework_Services_Registry
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
                        ROOT::MathMore
                        ROOT::EG
                        ROOT::TreePlayer
                        ROOT::FFTW
                        ROOT::EGPythia6
                        ROOT::Gui
                        )
cet_build_plugin(TestGENIEHelper art::module
                 LIBRARIES REG nusimdata::SimulationBase
                        nugen::EventGeneratorBase_GENIE
                        art::Framework_Services_Optional_RandomNumberGenerator_service
                        art::Framework_Services_Registry
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
                        dk2nu::Genie
                        ROOT::Core
                        ROOT::Geom
                        ROOT::GeomPainter
                        ROOT::MathMore
                        ROOT::EG
                        ROOT::TreePlayer
                        ROOT::FFTW
                        ROOT::EGPythia6
                        ROOT::Gui
                        )

install_headers()
install_fhicl()
install_source()
