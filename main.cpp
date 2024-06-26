#include <iostream>
#include <string>
#include <vector>

#include "src/fem.h"
#include "src/mesh.h"
#include "src/solver.h"
#include "src/tests.h"
#include "src/simu.h"

/* Global variables */
std::vector< std::string > arguments;


/* To parse command line arguments */
bool flag_is_used(
    const std::string& flag,
    const std::vector< std::string >& arguments )
{
    for( int i = 0; i < arguments.size(); ++i ) {
        if( flag == arguments[i] ) {
            return true;
        }
    }
    return false;
}

using namespace FEM2A;

void run_tests()
{
    const bool t_opennl = false;
    const bool t_lmesh = false;
    const bool t_io = false;
    const bool t_quad = false;
    const bool t_map = true;
    const bool t_dir = false;
    const bool t_Fe = false;

    if( t_opennl ) test_opennl();
    if( t_lmesh ) Tests::test_load_mesh();
    if( t_io ) Tests::test_load_save_mesh();
    if(t_quad) Tests::test_quadrature();
    if(t_map) Tests::test_mapping();
    if(t_dir) Tests::test_Ke();
    if(t_Fe) Tests::test_Fe();
}

void run_simu()
{

    const bool simu_pure_dirichlet = false;
    const bool simu_dirichlet_et_terme_source = false;
    const bool simu_dirichlet_sinus_bump = false;
    const bool simu_dirichlet_neumann = true;

    const bool verbose = flag_is_used( "-v", arguments )
        || flag_is_used( "--verbose", arguments );

    if( simu_pure_dirichlet ) {
        Simu::pure_dirichlet_pb("data/square_fine.mesh", verbose);
    }
    if (simu_dirichlet_et_terme_source) {
    	Simu::dirichlet_terme_source_pb("data/square_fine.mesh", verbose);
    }
    if (simu_dirichlet_sinus_bump) {
    	Simu::dirichlet_terme_source_pb_sinus_bump("data/square.mesh", verbose);
    }
    if (simu_dirichlet_neumann) {
    	Simu::dirichlet_neumann_pb("data/square_fine.mesh", verbose);
    }
    
}

int main( int argc, const char * argv[] )
/*argc : argument count, le nombre d'arguments, argv[] : pointeur vers un tableau des arguments (chaines de caractères) */
{
    /* Command line parsing */
    for( int i = 1; i < argc; ++i ) {
        arguments.push_back( std::string(argv[i]) );
    }

    /* Show usage if asked or no arguments */
    if( arguments.size() == 0 || flag_is_used("-h", arguments)
        || flag_is_used("--help", arguments) ) {
        std::cout << "Usage: ./fem2a [options]" << std::endl
            << "Options: " << std::endl;
        std::cout << " -h, --help:        show usage" << std::endl;
        std::cout << " -t, --run-tests:   run the tests" << std::endl;
        std::cout << " -s, --run-simu:    run the simulations" << std::endl;
        std::cout << " -v, --verbose:     print lots of details" << std::endl;
        return 0;
    }

    /* Run the tests if asked */
    if( flag_is_used("-t", arguments)
        || flag_is_used("--run-tests", arguments) ) {
        run_tests();
    }

    /* Run the simulation if asked */
    if( flag_is_used("-s", arguments)
        || flag_is_used("--run-simu", arguments) ) {
        run_simu();
    }

    return 0;
}
