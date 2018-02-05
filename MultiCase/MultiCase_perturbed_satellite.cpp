/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *	  Written by Joep for the Space Flight Assignment calculating errors in the Tudat software.
 */

#include <Tudat/SimulationSetup/tudatSimulationHeader.h>

#include "applicationOutput.h"

#include <boost/lexical_cast.hpp>
#include <boost/filesystem.hpp>
#include <cmath>
#include <fstream>

#include <json/src/json.hpp>
using json = nlohmann::json;

using namespace tudat;
using namespace tudat::interpolators;
using namespace tudat::propagators;
using namespace tudat::numerical_integrators;
using namespace tudat::orbital_element_conversions;
using namespace tudat::basic_mathematics;
using namespace tudat::gravitation;
using namespace tudat::simulation_setup;

static std::map< std::string, TranslationalPropagatorType > translationalPropagatorTypes =
{
    { "cowell", cowell },
    { "encke", encke },
    { "gaussKeplerian", gauss_keplerian },
    { "gaussModifiedEquinoctial", gauss_modified_equinoctial}
};

void getAndSaveSolution(boost::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings,
	boost::shared_ptr< IntegratorSettings< > > integratorSettings, NamedBodyMap bodyMap, std::vector< std::string > bodiesToCreate,
    int i, std::string propagator, std::string integrator, int step_size_exp, std::string outputDirectory) {
		// Create simulation object and propagate dynamics.
	SingleArcDynamicsSimulator< > dynamicsSimulator(
				bodyMap, integratorSettings, propagatorSettings );
	std::map< double, Eigen::VectorXd > integrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );

	// Write satellite propagation history to file.
	input_output::writeDataMapToTextFile( integrationResult,
        "Case_" + boost::lexical_cast< std::string >(i) + "_" + propagator + "_" + integrator + "_E" + boost::lexical_cast< std::string >(step_size_exp) + ".dat",
		outputDirectory, "", std::numeric_limits< double >::digits10, std::numeric_limits< double >::digits10, "," );
}

int main( )
{
	std::cout << "Setting up.." << std::endl;

    // Load Spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    const boost::filesystem::path cppFilePath( __FILE__ );
	const std::string baseName = cppFilePath.stem().string();
  	const std::string cppFolder = cppFilePath.parent_path().string();
    const std::string outputDirectory = cppFolder + "/" + baseName + "_Out/";

	// Read JSON file
    std::ifstream SettingsFile( cppFolder + "\\" + baseName + ".json");
	json Settings;
	SettingsFile >> Settings;

	std::cout << "Found JSon file.." << std::endl;

	// Read constants from JSON file
	// Epoch settings
	const double simulationStartEpoch = Settings[ "initialEpoch" ];
	const double simulationEndEpoch = Settings[ "finalEpoch" ];
	// An amount to extend the data of the bodies for, to make Lagrange Interpolation happy.
	const double simulationEndEpochExtend = 1e6; // 1e6 should be enough (kind of dependend on the max_step_size)

	// Enviroment and body settings
	const std::string frameOriginBody = Settings[ "globalFrameOrigin" ]; // Probably just Earth
	const std::string orientation = Settings[ "globalFrameOrientation" ]; // Probably J2000

	// Read array with different states/cases to be calculated
	json states = Settings[ "satellite" ][ "states" ];
	const std::string satelliteName = Settings[ "satellite" ][ "name" ];
	const double satelliteMass = Settings["satellite"]["mass"];

	// Read settings for fixed integrators
	const signed int step_size_exp_begin = Settings[ "integratorsettings" ][ "fixed" ][ "step_size_exp_begin" ];
	const signed int step_size_exp_end = Settings[ "integratorsettings" ][ "fixed" ][ "step_size_exp_end" ];
	const signed int step_size_exp_step = Settings[ "integratorsettings" ][ "fixed" ][ "step_size_exp_step" ];

	// Read settings for variable integrators
	const signed int stepSizeInit =  Settings[ "integratorsettings" ][ "variable" ][ "stepSizeInit" ];
	const signed int max_StepN = Settings[ "integratorsettings" ][ "variable" ][ "max_StepN" ];
	const double min_step_size = Settings[ "integratorsettings" ][ "variable" ][ "min_step_size" ];
	const double max_step_size = Settings[ "integratorsettings" ][ "variable" ][ "max_step_size" ];
	const signed int rel_error_tol_exp_begin = Settings[ "integratorsettings" ][ "variable" ][ "rel_error_tol_exp_begin" ];
	const signed int rel_error_tol_exp_end = Settings[ "integratorsettings" ][ "variable" ][ "rel_error_tol_exp_end" ];
	const signed int rel_error_tol_exp_step = Settings[ "integratorsettings" ][ "variable" ][ "rel_error_exp_step" ];

	// Integrators and Propagators to be used (Strings are mostly to give the output data-file the correct name)
	std::vector<std::string> propagators = Settings["propagators"];
	
	std::vector<std::string> integrators = Settings["integrators"]["variable"];
	std::vector<std::string> integrators_fix = Settings["integrators"]["fixed"];
	integrators.insert( integrators.end(), integrators_fix.begin(), integrators_fix.end() );

	std::cout << "Settings read.." << std::endl;

	// Create body objects.
	std::vector< std::string > bodiesToCreate;
	for ( unsigned int b = 0; b < bodies.size(); b++ ) {
		// if( bodies.at( b ) != "Sun")
		bodiesToCreate.push_back( bodies[b] );
	}
	
	std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings =
			getDefaultBodySettings( bodiesToCreate, simulationStartEpoch, simulationEndEpoch + simulationEndEpochExtend);
	
	for( unsigned int i = 0; i < bodiesToCreate.size( ); i++ ) {
	    bodySettings[ bodiesToCreate.at( i ) ]->ephemerisSettings->resetFrameOrientation( "J2000" );
	    bodySettings[ bodiesToCreate.at( i ) ]->rotationModelSettings->resetOriginalFrame( "J2000" );
	}

	NamedBodyMap bodyMap = createBodies( bodySettings );

	bodyMap[ satelliteName ] = boost::make_shared< simulation_setup::Body >( );
	bodyMap[ satelliteName ]->setConstantBodyMass( satelliteMass );

	setGlobalFrameBodyEphemerides( bodyMap, "SSB", orientation );

	std::cout << "Bodies are setup.." << std::endl;

	// Define propagator settings variables.
	SelectedAccelerationMap accelerationMap;
	std::vector< std::string > bodiesToPropagate;
	std::vector< std::string > centralBodies;

	bodiesToPropagate.push_back( satelliteName );
	centralBodies.push_back( frameOriginBody );

	// Define propagation settings.
	std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerations;
	accelerations[ frameOriginBody ].push_back( boost::make_shared< AccelerationSettings >(
													 basic_astrodynamics::central_gravity ) );
	accelerationMap[ satelliteName ] = accelerations;

	// Create acceleration models and propagation settings.
	basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
				bodyMap, accelerationMap, bodiesToPropagate, centralBodies );

	// Convert Body state from Keplerian elements to Cartesian elements.
	double earthGravitationalParameter = bodyMap.at( frameOriginBody )->getGravityFieldModel( )->getGravitationalParameter( );

	// Vector with the Kepler elements to be used per case.
	Eigen::Vector6d KeplerianElements;
	// Integer containing the exponent to be used per Integrator fixed step size or relative error
	int signed step_size_exp;

	std::cout << "Setup completed.." << std::endl;

	// Begin the for loops for cases -> propagators -> integrators -> step size/rel error size
	for(unsigned int i = 0; i <states.size(); i++) { // cases

		std::cout << "Calculating state: " << i + 1 << " out of " << states.size() <<std::endl;

		// Set Keplerian elements for satellite.
		KeplerianElements( semiMajorAxisIndex ) = states[i][ "semiMajorAxis" ];
		KeplerianElements( eccentricityIndex ) = states[i][ "eccentricity" ];
		KeplerianElements( inclinationIndex ) = states[i][ "inclination" ];
		KeplerianElements( argumentOfPeriapsisIndex ) = states[i][ "argumentOfPeriapsis" ];
		KeplerianElements( longitudeOfAscendingNodeIndex ) = states[i][ "longitudeOfAscendingNode" ];
		KeplerianElements( trueAnomalyIndex ) = states[i][ "trueAnomaly" ];

		// convert to Cartesian elements
		Eigen::VectorXd systemInitialState = convertKeplerianToCartesianElements(
					KeplerianElements, earthGravitationalParameter );

		// Loop over the different propagators
		// The propagators have to be of type TranslationalPropagatorType, but become integers when evaluated. 
		// So they have to be cast (which is done when making the propagators)
		for(unsigned int j = 0; j < propagators.size(); j++){ // nBodyStateDerivative -> Enum for propagators with intvalue
			
			std::cout << "\tUsing propagator: " << propagators[j] << std::endl;

			// Set-up propagator
			boost::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
					boost::make_shared< TranslationalStatePropagatorSettings< double > >
					( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState, simulationEndEpoch, translationalPropagatorTypes[propagators[j]]);

			// Loop over integrators
			// Since the integrators have different constructors, a switch-case is used to determine how to handle the integrator
			// Prefereably an overloadable function is used, but this will suffice for now.
			// -> After setting up an integrator with a certain stepsize the propagation is calculated and the results are exported
			for(unsigned int k = 0; k < integrators.size(); k++ ){ // integrators

				std::cout << "\t\tUsing Integrator: " << integrators[k] << std::endl;

				if ( integrators[k].compare( "Euler" ) == 0 ){ // Euler -> fixed
					for(step_size_exp = step_size_exp_begin;
						step_size_exp <= step_size_exp_end;
						step_size_exp += step_size_exp_step)
					{
						boost::shared_ptr< IntegratorSettings< > > integratorSettings
							= boost::make_shared< IntegratorSettings< > >(
								euler,                          // integrator
								simulationStartEpoch,           // simulation start
								std::pow(10.0, step_size_exp)   // step size
						);
						getAndSaveSolution(propagatorSettings, integratorSettings, bodyMap, bodiesToCreate,
					 	 i, propagators[j], integrators[k], step_size_exp, outputDirectory);
					}
				}
				else if ( integrators[k].compare( "Runge-Kutta4" ) == 0 ){ // Runge-Kutta4 -> fixed
					for(step_size_exp = step_size_exp_begin;
						step_size_exp <= step_size_exp_end;
						step_size_exp += step_size_exp_step)
					{
						boost::shared_ptr< IntegratorSettings< > > integratorSettings
							= boost::make_shared< IntegratorSettings< > >(
								rungeKutta4,                    // integrator
								simulationStartEpoch,           // simulation start
								std::pow(10.0, step_size_exp)   // step size
						);
						getAndSaveSolution(propagatorSettings, integratorSettings, bodyMap, bodiesToCreate,
					 	 i, propagators[j], integrators[k], step_size_exp, outputDirectory);
					}
				}
				else if ( integrators[k].compare( "Runge-Kutta78" ) == 0 ){ // Runge-Kutta78 - variable
					for(step_size_exp = rel_error_tol_exp_begin;
						step_size_exp >= rel_error_tol_exp_end;
						step_size_exp += rel_error_tol_exp_step ) 
					{
						boost::shared_ptr< IntegratorSettings< > > integratorSettings =
							boost::make_shared< RungeKuttaVariableStepSizeSettings< > >( 
								rungeKuttaVariableStepSize,         // integrator
								simulationStartEpoch,               // init time
								stepSizeInit,                       // init stepsize
								RungeKuttaCoefficients::CoefficientSets::rungeKuttaFehlberg78,
								min_step_size,                      // min step size
								max_step_size,                      // max stepsize
								std::pow(10.0,step_size_exp),       // rel err tol
								std::pow(10.0,step_size_exp)        // abs err tol
						);
						getAndSaveSolution(propagatorSettings, integratorSettings, bodyMap, bodiesToCreate,
						 i, propagators[j], integrators[k], step_size_exp, outputDirectory);
					}		
				}
				else if ( integrators[k].compare( "Bulirsch-Stoer" ) == 0 ){ // Bulirsch-Stoer - variable
					for(step_size_exp = rel_error_tol_exp_begin;
						step_size_exp >= rel_error_tol_exp_end;
						step_size_exp += rel_error_tol_exp_step ) 
					{
						boost::shared_ptr< IntegratorSettings< > > integratorSettings =
							boost::make_shared< BulirschStoerIntegratorSettings< > >( 
								simulationStartEpoch,           // init time
								stepSizeInit,                   // init stepsize
								bulirsch_stoer_sequence,        // sequence
								max_StepN,                      // max stepN
								min_step_size,                  // min step size
								max_step_size,                  // max stepsize
								std::pow(10.0,step_size_exp),   // rel err tol
								std::pow(10.0,step_size_exp)    // abs err tol
						);
						getAndSaveSolution(propagatorSettings, integratorSettings, bodyMap, bodiesToCreate,
					 	 i, propagators[j], integrators[k], step_size_exp, outputDirectory);
					}
				}
				else {
					std::cout << "integrator defined in JSON file is not implemented." << std::endl;
				}
			}
		}
	}
	
	std::cout << "DONE" << std::endl;

	// Final statement.
	// The exit code EXIT_SUCCESS indicates that the program was successfully executed.
	return EXIT_SUCCESS;
}