#pragma once

#include "mesh.h"
#include "fem.h"
#include "solver.h"
// #include "simu.h"

#include <assert.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <stdlib.h>


namespace FEM2A {
    namespace Tests {

        bool test_load_mesh()
        {
            Mesh mesh;
            mesh.load("data/square.mesh");
            /*Initialisation du maillage*/

		/*Affichage des données du maillage, vertices*/
            std::cout << "Vertices <x> <y> <att>" << std::endl;
            for( int v = 0; v < mesh.nb_vertices(); v++ ) {
                std::cout << mesh.get_vertex(v).x << " " << mesh.get_vertex(v).y
                    << " " << mesh.get_vertex_attribute(v) << std::endl;
            }

            std::cout << "Edges <v0> <v1> <att>" << std::endl ;
            for( int ed = 0; ed < mesh.nb_edges(); ed++ ) {
                std::cout << mesh.get_edge_vertex_index(ed, 0) << " "
                    << mesh.get_edge_vertex_index(ed, 1) << " "
                    << mesh.get_edge_attribute(ed) << std::endl;
            }

            std::cout << "Triangles <v0> <v1> <v2> <att>" << std::endl ;
            for( int tr = 0; tr < mesh.nb_triangles(); tr++ ) {
                std::cout << mesh.get_triangle_vertex_index(tr, 0) << " "
                    << mesh.get_triangle_vertex_index(tr, 1) << " "
                    << mesh.get_triangle_vertex_index(tr, 2) << " "
                    << mesh.get_triangle_attribute(tr) << std::endl;
            }

            return true;
        }

        bool test_load_save_mesh()
        {
            Mesh mesh;
            mesh.load("data/geothermie_4.mesh");
            mesh.save("data/geothermie_4.mesh");
            return true;
        }
        
        // Test de la quadrature
        bool test_quadrature(){
        	std::cout << "coucou"<< std::endl;
        	
        	Quadrature quad = quad.get_quadrature(0);
        	int nb_pts = quad.nb_points();
        	std::cout << "nombre de points dans la quadrature : " << nb_pts << std::endl;
        	
        	double poids = 0;
        	for (int i=0; i < nb_pts; ++i){
        		poids += quad.weight(i);
        	}
        	std::cout << "poids quadrature : " << poids << std::endl;
        	
        	Quadrature quad2 = quad2.get_quadrature(2);
        	int nb_pts2 = quad2.nb_points();
        	std::cout << "nombre de points dans la quadrature : " << nb_pts2 << std::endl;
        	double poids2 = 0;
        	for (int i=0; i < nb_pts2; ++i){
        		poids2 += quad2.weight(i);
        	}
        	std::cout << "poids quadrature : " << poids2 << std::endl;
        	return true;
        	
        
        }
        
        double unit_fct( vertex v )
        {
            return 1.;
        }
        
	bool test_mapping(){
		Mesh mesh;
            	mesh.load("data/square.mesh");
            	mesh.save("data/square.mesh");
            	
        	std::cout << "Données des vertices du triangle numéro 4 du maillage square.mesh."<< std::endl;
        	
        	
        	vertex v0 = mesh.get_triangle_vertex(4, 0);
        	vertex v1 = mesh.get_triangle_vertex(4, 1);
        	vertex v2 = mesh.get_triangle_vertex(4, 2);
        	
        	std::vector< vertex> vecteur1{v0, v1, v2};
        	for (int i =0; i < vecteur1.size(); ++i){
        		std::cout << "vertice nb " << i << " : x = " << vecteur1[i].x << " y = "<< vecteur1[i].y << std::endl;
		}
        	
        	std::cout << "v0 : x = " << v0.x << " y = " << v0.y << std::endl;
        	std::cout << "v1 : x = " << v1.x << " y = " << v1.y << std::endl;
        	std::cout << "v2 : x = " << v2.x << " y = " << v2.y << std::endl;
        	
        	ElementMapping mapping = ElementMapping(mesh, false, 4);
        	
        	vertex point1{0.2, 0.4};
        	
        	vertex image_point1;
        	
        	image_point1 = mapping.transform(point1);
        	
        	std::cout << "Position x : " << image_point1.x << " ; position y : " << image_point1.y << std::endl;
        	
        	std::cout << "" << std::endl;
        	std::cout << "Tests sur la matrice Jacobienne pour le triangle;" << std::endl;
        	std::cout << "" << std::endl;
        	
        	DenseMatrix J = mapping.jacobian_matrix(point1);
        	std::cout << "Matrice Jacobienne : " << std::endl;
        	std::cout << J.get(0, 0) << " " << J.get(0, 1) << std::endl;
        	std::cout << J.get(1, 0) << " " << J.get(1, 1) << std::endl;
        	
        	std::cout << "" << std::endl;
        	
        	std::cout << "Determinant : " << std::endl;
        	std::cout << mapping.jacobian(point1) << std::endl;
        	
        	std::cout << "" << std::endl;
        	std::cout << "Tests sur la matrice Jacobienne pour un segment;" << std::endl;
        	std::cout << "" << std::endl;
        	///*
        	
        	ElementMapping mapping2 = ElementMapping(mesh, true, 4);
        	vertex v0_ = mesh.get_edge_vertex(4, 0);
        	vertex v1_ = mesh.get_edge_vertex(4, 1);
        	std::cout << "v0 : x = " << v0_.x << " y = " << v0_.y << std::endl;
        	std::cout << "v1 : x = " << v1_.x << " y = " << v1_.y << std::endl;
        	
        	//*/
        	/*
        	int dim = 2;
        	int order = 1;
        	ShapeFunctions functions = ShapeFunctions(dim, order);
        	*/
        	
        	
        	DenseMatrix Ke;
        	Ke.set_size(3, 3);
        	ShapeFunctions ref_functions(2, 1);
        	Quadrature quadrature = quadrature.get_quadrature(6);
        	assemble_elementary_matrix(mapping, ref_functions, quadrature, unit_fct, Ke);
        	
        	std::cout << "calcul de Ke : " << std::endl;
        	

        	
        	for (int i = 0; i<3; i++){
        		for (int j =0; j<3; j++) {
        			std::cout << Ke.get(i, j) << " ";
        		}
        		std::cout << "" << std::endl;
        	}
        	
        	SparseMatrix K(mesh.nb_vertices());
        	int t = 4;
        	local_to_global_matrix(mesh, t, Ke, K);
        	
        	std::vector <double> Fe;
        	assemble_elementary_vector(mapping, ref_functions, quadrature, unit_fct, Fe);
        	for (int i = 0; i < Fe.size() ; i++) {
        		std::cout << "Fe " << i << " : " << Fe[i] << std::endl;
        	}
        	
        	std::vector <double> Fe2;
        	ShapeFunctions ref_functions_1D(1, 1);
        	Quadrature quadrature_1D = quadrature_1D.get_quadrature(2, true);
        	assemble_elementary_neumann_vector(mapping2, ref_functions_1D, quadrature_1D, unit_fct, Fe2);
        	std::cout << "Vecteur conditions de Neumann : " << std::endl;
        	for (int i = 0; i < Fe2.size() ; i++) {
        		std::cout << "Fe " << i << " : " << Fe2[i] << std::endl;
        	}
        	
        	return true;
        }
        
        bool test_Ke(){
        
        	Mesh mesh;
            	mesh.load("data/square.mesh");
            	ElementMapping my_map(mesh, false, 4);
            	ShapeFunctions my_shpfct(2, 1);
            	Quadrature my_quad = Quadrature::get_quadrature(2);
            	DenseMatrix Ke;
            	assemble_elementary_matrix(my_map, my_shpfct, my_quad, unit_fct, Ke);
            	Ke.print();
            	SparseMatrix K_glob(mesh.nb_vertices());
            	local_to_global_matrix(mesh, 4, Ke, K_glob);
            	K_glob.print();
            	
            	/*
            	std::vector <bool> attribute_is_dirichlet;
            	std::vector <double> values;
            	SparseMatrix K;
            	std::vector <double> F;
            	
            	apply_dirichlet_conditions(attribute_is_dirichlet, values, K, F);
        	*/
        	return true;
        }
    }

}

        
