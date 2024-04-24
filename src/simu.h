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
        
        double sinus_bump_fct( vertex v)
        {
        	return 2*M_PI*M_PI*sin(M_PI*v.x)*sin(M_PI*v.y);
        }
        
        double solution_fct( vertex v)
        {
        	return sin(M_PI*v.x)*sin(M_PI*v.y);
        }
        
        double neumann_fct( vertex v)
        {
        	double epsilon = 0.001;
        	if (abs(v.x) < epsilon) {
        		return 1;
        	}
        	else {
        		return -1;
        	}
        }
        
        double neumann_cdt_fct( vertex v)
        {
        	return sin(M_PI*v.y);
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
            std::string export_name = "pure_dirichlet_fine";
            mesh.save(export_name+".mesh");
            save_solution(u, export_name+".bb");
        }
        
        void dirichlet_terme_source_pb( const std::string& mesh_filename, bool verbose)
        {
            std::cout << "Solving Dirichlet problem with source term" << std::endl;
            if ( verbose ) {
                std::cout << " with lots of printed details..." << std::endl;
            }
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
            	std::vector <double> Fe;
            	assemble_elementary_vector(mapping, shpfct, quad, unit_fct, Fe);
            	for (int nb = 0; nb < Fe.size(); ++nb) {
            		std::cout << "Fe[ " << nb << "] : " << Fe[nb] << std::endl;
            	}
            	local_to_global_vector(mesh, false, triangle, Fe, F_glob);
            }
            
            for (int nb = 0; nb < F_glob.size(); nb++) {
            	std::cout << "F_glob[ " << nb << " ] : " << F_glob[nb] << std::endl;
            }
            std::vector<bool> att_is_dirichlet(2, false);
            att_is_dirichlet[1] = true;
            mesh.set_attribute(unit_fct, 1, true);
            std::vector<double> imposed_values(mesh.nb_vertices());
            for( int i = 0; i < mesh.nb_vertices(); ++i){
            	imposed_values[i] = 0;
            }
            apply_dirichlet_boundary_conditions(mesh, att_is_dirichlet, imposed_values, K_glob, F_glob);
            std::vector<double> u(mesh.nb_vertices());
            solve(K_glob, F_glob, u);
            std::string export_name = "dirichlet_et_terme_source_fine";
            mesh.save(export_name+".mesh");
            save_solution(u, export_name+".bb");

    }
    
        void dirichlet_terme_source_pb_sinus_bump( const std::string& mesh_filename, bool verbose)
        {
            std::cout << "Solving Dirichlet problem with source term" << std::endl;
            if ( verbose ) {
                std::cout << " with lots of printed details..." << std::endl;
            }
            vertex v;
            v.x = 0.5;
            v.y = 0.5;
            std::cout << "sinus bump " << sinus_bump_fct(v) << std::endl;
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
            	std::vector <double> Fe;
            	assemble_elementary_vector(mapping, shpfct, quad, sinus_bump_fct, Fe);
            	/*
            	for (int nb = 0; nb < Fe.size(); ++nb) {
            		std::cout << "Fe[ " << nb << "] : " << Fe[nb] << std::endl;
            	}
            	*/
            	local_to_global_vector(mesh, false, triangle, Fe, F_glob);
            }
            /*
            for (int nb = 0; nb < F_glob.size(); nb++) {
            	std::cout << "F_glob[ " << nb << " ] : " << F_glob[nb] << std::endl;
            }
            */
            std::vector<bool> att_is_dirichlet(2, false);
            att_is_dirichlet[1] = true;
            mesh.set_attribute(unit_fct, 1, true);
            std::vector<double> imposed_values(mesh.nb_vertices());
            for( int i = 0; i < mesh.nb_vertices(); ++i){
            	imposed_values[i] = 0;
            }
            apply_dirichlet_boundary_conditions(mesh, att_is_dirichlet, imposed_values, K_glob, F_glob);
            std::vector<double> u(mesh.nb_vertices());
            solve(K_glob, F_glob, u);
            double square_error = 0;
            std::vector <double> ecart;
            for (int vertice = 0; vertice < mesh.nb_vertices(); vertice++) {
            	square_error += pow(u[vertice] - solution_fct(mesh.get_vertex(vertice)), 2);
            	ecart.push_back(u[vertice] - solution_fct(mesh.get_vertex(vertice)));
            }
            
            
            std::cout << "square error : " << square_error << std::endl;
            std::string export_name = "dirichlet_sinus_bump_fine";
            mesh.save(export_name+".mesh");
            save_solution(u, export_name+".bb");
            std::string export_name2 = "dirichlet_sinus_bump_ecart_fine";
            mesh.save(export_name2+".mesh");
            save_solution(ecart, export_name2+".bb");

    }
    
    void dirichlet_neumann_pb( const std::string& mesh_filename, bool verbose)
        {
            std::cout << "Solving Dirichlet problem with neumann conditions" << std::endl;
            if ( verbose ) {
                std::cout << " with lots of printed details..." << std::endl;
            }
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
            	std::vector <double> Fe;
            	assemble_elementary_vector(mapping, shpfct, quad, unit_fct, Fe);
            	/*
            	for (int nb = 0; nb < Fe.size(); ++nb) {
            		std::cout << "Fe[ " << nb << "] : " << Fe[nb] << std::endl;
            	}
            	*/
            	local_to_global_vector(mesh, false, triangle, Fe, F_glob);
           
            }
            for (int edge_vertex = 0; edge_vertex < mesh.nb_edges(); edge_vertex++) {
            	ElementMapping mapping_1D(mesh, true, edge_vertex);
            	ShapeFunctions shpfct_1D(1,1);
            	Quadrature quad_1D = Quadrature::get_quadrature(2);
            	std::cout << "nb : " << edge_vertex << std::endl;
            	for (int local_index = 0; local_index < 2; local_index++) {
            		if (neumann_fct(mesh.get_edge_vertex(edge_vertex, local_index)) > 0) {
            			std::cout << "On a un neumann" << std::endl;
            			std::vector <double> Fe;
            			assemble_elementary_neumann_vector(mapping_1D, shpfct_1D, quad_1D, neumann_cdt_fct, Fe);
            			for (int nb = 0; nb < Fe.size(); nb++) {
            				std::cout << "Fe[ " << nb << " ] : " << Fe[nb] << std::endl;
            			}
            			local_to_global_vector(mesh, true, edge_vertex, Fe, F_glob);
            		}
            	}
            }
            
            /*
            for (int nb = 0; nb < F_glob.size(); nb++) {
            	std::cout << "F_glob[ " << nb << " ] : " << F_glob[nb] << std::endl;
            }
            */
            std::vector<bool> att_is_dirichlet(2, false);
            att_is_dirichlet[1] = true;
            mesh.set_attribute(unit_fct, 1, true);
            std::vector<double> imposed_values(mesh.nb_vertices());
            for( int i = 0; i < mesh.nb_vertices(); ++i){
            	imposed_values[i] = 0;
            }
            apply_dirichlet_boundary_conditions(mesh, att_is_dirichlet, imposed_values, K_glob, F_glob);
            std::vector<double> u(mesh.nb_vertices());
            solve(K_glob, F_glob, u);

            std::string export_name = "dirichlet_neumann_square";
            mesh.save(export_name+".mesh");
            save_solution(u, export_name+".bb");


    }

}
}
