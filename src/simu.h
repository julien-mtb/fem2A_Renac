#pragma once

#include "mesh.h"
#include "fem.h"
#include <math.h>
#include <cmath>
#include <iostream>

namespace FEM2A {
    namespace Simu {

        //#################################
        //  Useful functions
        //#################################

        double unit_fct( vertex v )
        {
            return 1.;
        }

        double zero_fct( vertex v )
        {
            return 0.;
        }

        double xy_fct( vertex v )
        {
            return v.x + v.y;
        }

        //#################################
        //  Simulations
        //#################################

        void pure_dirichlet_pb( const std::string& mesh_filename, bool verbose )
        {
            std::cout << "Solving a pure Dirichlet problem" << std::endl;
            if ( verbose ) {
                std::cout << " with lots of printed details..." << std::endl;
            }
            std::cout << "TO BE IMPLEMENTED !!!" << std::endl;
            Mesh mesh;
            mesh.load(mesh_filename);
            std::vector <bool> attribute_is_dirichlet;
            std::cout << " Attribut max : " << mesh.get_attr_max() << std::endl;
            
            for (int num_edge = 0; num_edge < mesh.nb_edges; num_edge++){
            	std::cout << 
            }
            
            /*
            for (int num_edge = 0; num_edge < mesh.nb_edges; num_edge++){
            	
            	if (
            	attribute_is_dirichlet.push_back(true);
            	int edge_attribute = mesh.get_edge_attribute(num_edge);
            
            }
            */
        }

    }

}
