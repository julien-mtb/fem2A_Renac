#pragma once

#include "mesh.h"
#include "fem.h"
#include "solver.h"

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
        	
        	return true;
        }
    }
}
