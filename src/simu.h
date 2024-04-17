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

        void pure_dirichlet_pb_v2( const std::string& mesh_filename, bool verbose )
        {
            std::cout << "Solving a pure Dirichlet problem" << std::endl;
            if ( verbose ) {
                std::cout << " with lots of printed details..." << std::endl;
            }
            std::cout << "TO BE IMPLEMENTED !!!" << std::endl;
            Mesh mesh;
            // On charge le maillage
            mesh.load(mesh_filename);
            
            ShapeFunctions ref_functions(2, 1);
            Quadrature quadrature = quadrature.get_quadrature(2);
            
            std::vector <double> F(mesh.nb_vertices());
            std::vector <double> U(mesh.nb_vertices());
            SparseMatrix K(mesh.nb_vertices());
            std::vector <double> values(mesh.nb_vertices());
            std::vector <bool> attribute_is_dirichlet{true, true, true, true, false, false};
          
            for (int triangle = 0; triangle < mesh.nb_triangles(); triangle++){
            	ElementMapping mapping(mesh, false, triangle);
            	DenseMatrix Ke;
            	assemble_elementary_matrix(mapping, ref_functions, quadrature, unit_fct, Ke);
            	local_to_global_matrix(mesh, triangle, Ke, K);
            }
            // Préparation du vecteur U et du vecteur F
            ///*
            for (int vertice = 0; vertice < mesh.nb_vertices(); vertice++){
            	F[vertice] = 0;
            	// get_bdr_attr_max renvoie la valeur max de l'attribut d'un élément border
            	if(mesh.get_vertex_attribute(vertice) <= mesh.get_bdr_attr_max()){
            		vertex v = mesh.get_vertex(vertice);
            		values[vertice] = v.x+v.y;
            		std::cout << "vertex " << vertice << " , valeur : " << v.x + v.y << std::endl;
            	}
            	else{
            		values[vertice] = 0;
            	}
            }
            std::cout << "Taille de values : " << values.size() << std::endl;
            // application des conditions de dirichlet
            bool solved = false;
            /*
            for (int i = 0; i < values.size(); i++){
            	std::cout << "values[ " << i << " ] : " << values[i] << std::endl;
            }
            */
            apply_dirichlet_boundary_conditions(mesh, attribute_is_dirichlet, values, K, F);
            /*
            std::cout << "Affichage de K : " << std::endl;
            K.print();
            */
            solved = solve(K, F, U);
            /*
            std::cout << "Affichage de U : " << std::endl;
            */
            /*
            for (int i = 0; i < U.size(); i++){
            	std::cout << "U[ " << i << " ] : " << U[i] << std::endl;
            }
            */
            std::string export_name = "pure_dirichlet";
            mesh.save(export_name+".mesh");
            save_solution(U, export_name+".bb");
        }
        
        void pure_dirichlet_pb( const std::string& mesh_filename, bool verbose)
        {
            std::cout << "Solving a pure Dirichlet problem" << std::endl;
            if ( verbose ) {
                std::cout << " with lots of printed details..." << std::endl;
            }
            std::cout << "TO BE IMPLEMENTED !!!" << std::endl;
            Mesh mesh;
            // On charge le maillage
            mesh.load(mesh_filename);
            SparseMatrix K_glob(mesh.nb_vertices());
            std::vector<double> F_glob(mesh.nb_vertices());
            for (int triangle = 0; triangle < mesh.nb_triangles(); triangle++){
            	ElementMapping mapping(mesh, false, triangle);
            	ShapeFunctions shpfct(2,1);
            	Quadrature quad = Quadrature::get_quadrature(2);
            	DenseMatrix Ke;
            	assemble_elementary_matrix(mapping, shpfct, quad, unit_fct, Ke);
            	local_to_global_matrix(mesh, triangle, Ke, K_glob);
            }
            std::vector<bool> att_is_dirichlet(2, false);
            att_is_dirichlet[1] = true;
            mesh.set_attribute(unit_fct, 1, true);
            std::vector<double> imposed_values(mesh.nb_vertices());
            for( int i = 0; i < mesh.nb_vertices(); ++i){
            	imposed_values[i] = xy_fct(mesh.get_vertex(i));
            }
            apply_dirichlet_boundary_conditions(mesh, att_is_dirichlet, imposed_values, K_glob, F_glob);
            std::vector<double> u(mesh.nb_vertices());
            solve(K_glob, F_glob, u);
            std::string export_name = "pure_dirichlet";
            mesh.save(export_name+".mesh");
            save_solution(u, export_name+".bb");
        }

    }

}
