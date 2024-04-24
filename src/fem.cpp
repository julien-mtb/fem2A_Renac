#include "fem.h"
#include "mesh.h"
// #include "simu.h"

#include <iomanip>
#include <iostream>
#include <cmath>
#include <limits>
#include <stdlib.h>
#include <assert.h>

namespace FEM2A {

    void print( const std::vector<double>& x )
    {
        for ( int i = 0; i < x.size(); ++i ) {
            std::cout << x[i] << " ";
        }
        std::cout << std::endl;
    }

    /****************************************************************/
    /* Implementation of Quadrature */
    /****************************************************************/
    int Quadrature::nb_points() const
    {
    // renvoie le nb de points de la quadrature
        return wxy_.size() / 3 ;
    }

    vertex Quadrature::point( int i ) const
    {
        assert( i < nb_points() ) ;
        vertex v ;
        v.x = wxy_[3 * i + 1] ;
        v.y = wxy_[3 * i + 2] ;
        return v ;
    }

    double Quadrature::weight( int i ) const
    {
    // weight renvoie le poids numéro i de la quadrature
        assert( i < nb_points() ) ;
        return wxy_[3 * i + 0] ;
    }

    const double triangle_P0[3] = {
        0.5, 0.333333333333333, 0.333333333333333
    };

    const double triangle_P2[9] = {
        0.166666666666667, 0.166666666666667, 0.166666666666667,
        0.166666666666667, 0.166666666666667, 0.666666666666667,
        0.166666666666667, 0.666666666666667, 0.166666666666667
    };

    const double triangle_P4[18] = {
        0.0549758718276609, 0.0915762135097707, 0.0915762135097707,
        0.0549758718276609, 0.0915762135097707, 0.816847572980459,
        0.0549758718276609, 0.816847572980459, 0.0915762135097707,
        0.111690794839006, 0.445948490915965, 0.445948490915965,
        0.111690794839006, 0.445948490915965, 0.10810301816807,
        0.111690794839006, 0.10810301816807, 0.445948490915965
    };

    const double triangle_P6[36] = {
        0.0254224531851034, 0.0630890144915022, 0.0630890144915022,
        0.0254224531851034, 0.0630890144915022, 0.873821971016996,
        0.0254224531851034, 0.873821971016996, 0.0630890144915022,
        0.0583931378631897, 0.24928674517091, 0.24928674517091,
        0.0583931378631897, 0.24928674517091, 0.501426509658179,
        0.0583931378631897, 0.501426509658179, 0.24928674517091,
        0.0414255378091868, 0.0531450498448169, 0.310352451033784,
        0.0414255378091868, 0.310352451033784, 0.0531450498448169,
        0.0414255378091868, 0.0531450498448169, 0.636502499121399,
        0.0414255378091868, 0.636502499121399, 0.0531450498448169,
        0.0414255378091868, 0.310352451033784, 0.636502499121399,
        0.0414255378091868, 0.636502499121399, 0.310352451033784
    };

    const double segment_P0[2] = {
        1., 0.5
    };

    const double segment_P2[4] = {
        0.5, 0.21132486540518708,
        0.5, 0.7886751345948129
    };

    Quadrature Quadrature::get_quadrature( int order, bool border )
    {
        double* pts = NULL;
        int nb_pts = 0;
        Quadrature Q;
        if ( order == 0 && !border ) {
            pts = const_cast<double*>(triangle_P0);
            nb_pts = 1;
        } else if ( order == 2 && !border ) {
            pts = const_cast<double*>(triangle_P2);
            nb_pts = 3;
        } else if ( order == 4 && !border ) {
            pts = const_cast<double*>(triangle_P4);
            nb_pts = 6;
        } else if ( order == 6 && !border ) {
            pts = const_cast<double*>(triangle_P6);
            nb_pts = 12;
        } else if ( order == 0 && border ) {
            pts = const_cast<double*>(segment_P0);
            nb_pts = 1;
        } else if ( order == 2 && border ) {
            pts = const_cast<double*>(segment_P2);
            nb_pts = 2;
        } else {
            std::cout << "Quadrature not implemented for order " << order << std::endl;
            assert( false );
        }
        Q.wxy_.resize(nb_pts * 3);
        for ( int i = 0; i < nb_pts; ++i ) {
            if ( !border ) {
                Q.wxy_[3*i+0] = pts[3*i+0];
                Q.wxy_[3*i+1] = pts[3*i+1];
                Q.wxy_[3*i+2] = pts[3*i+2];
            } else {
                Q.wxy_[3*i+0] = pts[2*i+0];
                Q.wxy_[3*i+1] = pts[2*i+1];
                Q.wxy_[3*i+2] = 0.;
            }
        }
        return Q;
    }

    /****************************************************************/
    /* Implementation of ElementMapping */
    /****************************************************************/
    ElementMapping::ElementMapping( const Mesh& M, bool border, int i )
        : border_( border ){
    // si border est true, on est sur un segment
    // sinon, on est dans un triangle
	/*
        std::cout << "[ElementMapping] constructor for element " << i << " ";
        */
        if ( border_ ){
        	// std::cout << "(border)" << std::endl;
                vertex vertice0 = M.get_edge_vertex(i, 0);
        	vertex vertice1 = M.get_edge_vertex(i, 1);
        
        	std::vector< vertex > vertices_all{vertice0, vertice1};
        	vertices_ = vertices_all;
        	// for (int i =0; i < vertices_all.size(); ++i){
        	//	std::cout << "vertice nb " << i << " :  x = " << vertices_all[i].x << " y = "<< vertices_all[i].y << std::endl;
        	//}
        }

        else{
        	// std::cout << "(triangle)" << std::endl;
        	vertex vertice0 = M.get_triangle_vertex(i, 0);
        	vertex vertice1 = M.get_triangle_vertex(i, 1);
        	vertex vertice2 = M.get_triangle_vertex(i, 2);
 
        	std::vector< vertex > vertices_all{vertice0, vertice1, vertice2};
        	vertices_ = vertices_all;
        	// for (int i =0; i < vertices_all.size(); ++i){
        	//	std::cout << "vertice nb " << i << " :  x = " << vertices_all[i].x << " y = "<< vertices_all[i].y << std::endl;
        	//}
        
        // TODO
    }
    }

    vertex ElementMapping::transform( vertex x_r ) const
    {
    	/*
        std::cout << "[ElementMapping] transform reference to world space" << '\n';
        */

        vertex r ;
        double eta = x_r.y;
        double ksi = x_r.x;
        double r_x, r_y;
        
        if (border_){
        	r.x = (1-ksi - eta)*vertices_[0].x + ksi*vertices_[1].x;
        	r.y = (1-ksi - eta)*vertices_[0].y + ksi*vertices_[1].y;
        }
        else{
        	r.x = (1-ksi - eta)*vertices_[0].x + ksi*vertices_[1].x + eta*vertices_[2].x;
        	r.y = (1-ksi - eta)*vertices_[0].y + ksi*vertices_[1].y + eta*vertices_[2].y;
        }

        return r ;
    }

    DenseMatrix ElementMapping::jacobian_matrix( vertex x_r ) const
    {
    	/*
        std::cout << "[ElementMapping] compute jacobian matrix" << '\n';
        */
        // TODO
        DenseMatrix J;
        if (border_){
        	J.set_size(2, 1);
        	J.set(0, 0, vertices_[1].x - vertices_[0].x);
        	J.set(1, 0, vertices_[1].y - vertices_[0].y);
        }
        else{
        	J.set_size(2, 2);
        	J.set(0, 0, vertices_[1].x - vertices_[0].x);
        	J.set(0, 1, vertices_[2].x - vertices_[0].x);
        	J.set(1, 0, vertices_[1].y - vertices_[0].y);
        	J.set(1, 1, vertices_[2].y - vertices_[0].y);
        	
        }
        return J;
    }

    double ElementMapping::jacobian( vertex x_r ) const
    {
    	/*
        std::cout << "[ElementMapping] compute jacobian determinant" << '\n';
        */
        // TODO
        DenseMatrix J = jacobian_matrix(x_r);
        double determinant = 0;
        if (border_) {
        	/*
        	std::cout << "J[0] : " << J.get(0, 0) << " , J[1] : " << J.get(1, 0) << std::endl;
        	std::cout << "calcul pow : " << pow(J.get(0, 0), 2) << " ; " << pow(J.get(1, 0), 2) << std::endl;
        	std::cout << "calcul final : " << pow(pow(J.get(0, 0), 2) + pow(J.get(1, 0), 2), 0.5) << std::endl;
        	*/
        	determinant = pow(pow(J.get(0, 0), 2) + pow(J.get(1, 0), 2), 0.5);
        	// std::cout << "det de vecteurs, on a det = " << determinant << std::endl;
        }
        
        else {
        	determinant = J.get(0, 0)*J.get(1, 1) - J.get(1, 0)*J.get(0, 1);
       		// std::cout << "On est passé par cet endoit, on a det = " << determinant << std::endl;
        }
        // return 0;
        return determinant;
    }

    /****************************************************************/
    /* Implementation of ShapeFunctions */
    /****************************************************************/
    ShapeFunctions::ShapeFunctions( int dim, int order )
        : dim_( dim ), order_( order )
    {
    	/*
        std::cout << "[ShapeFunctions] constructor in dimension " << dim << '\n';
        */
        assert(dim_ < 3);
        // TODO
    }

    int ShapeFunctions::nb_functions() const
    {
    	/*
        std::cout << "[ShapeFunctions] number of functions" << '\n';
        */
        ///*
        if (dim_ == 1) return 2;
        if (dim_ == 2) return 3;
        //*/
        // TODO
        return 0 ;
    }

    double ShapeFunctions::evaluate( int i, vertex x_r ) const
    {
    	/*
        std::cout << "[ShapeFunctions] evaluate shape function " << i << '\n';
        */
        ///*
        if (dim_ == 2) {
        	if (i==0) return 1-x_r.x - x_r.y;
        	if (i==1) return x_r.x;
        	if (i==2) return x_r.y;
        }
        if (dim_ == 1) {
        	// std::cout << "shape function de dim 1" << std::endl;
        	if (i==0) return 1-x_r.x;
        	if (i==1) return x_r.x;
        }
        //*/
        // TODO
        return 0. ; // should not be reached
    }

    vec2 ShapeFunctions::evaluate_grad( int i, vertex x_r ) const
    {
    	/*
        std::cout << "[ShapeFunctions] evaluate gradient shape function " << i << '\n';
        */
        vec2 g ;
        ///*
        if (dim_==1){
        	switch(i) {
        		case 0:
        			g.x = -1.; break;
        		case 1:
        			g.x = 1.; break;
        	}
        	g.y = 0.;
        }
        else {
        	switch(i) {
        		case 0:
        			g.x = -1.; g.y = -1.; break;
        		case 1:
        			g.x = 1.; g.y = 0.; break;
        		case 2:
        			g.x = 0.; g.y = 1.; break;
        	}
        }
        return g ;
    }

    /****************************************************************/
    /* Implementation of Finite Element functions */
    /****************************************************************/
    void assemble_elementary_matrix(
        const ElementMapping& elt_mapping,
        const ShapeFunctions& reference_functions,
        const Quadrature& quadrature,
        double (*coefficient)(vertex),
        DenseMatrix& Ke )
    {
    	/*
        std::cout << "compute elementary matrix" << '\n';
        */
        double scal;
        Ke.set_size(reference_functions.nb_functions(), reference_functions.nb_functions());
        for (int i = 0; i < reference_functions.nb_functions(); i++){
        	for(int j = 0; j < reference_functions.nb_functions(); j++){
        		Ke.set(i, j, 0);
        		scal = 0;
        		for (int q = 0; q < quadrature.nb_points(); q++){
        			vertex v = quadrature.point(q);
        			scal = dot(elt_mapping.jacobian_matrix(v).invert_2x2().transpose().mult_2x2_2(reference_functions.evaluate_grad(i, v)), elt_mapping.jacobian_matrix(v).invert_2x2().transpose().mult_2x2_2(reference_functions.evaluate_grad(j, v)));
        			Ke.add(i, j, quadrature.weight(q)*coefficient(elt_mapping.transform(v))*scal*elt_mapping.jacobian(v)); 
        				
        		}
        	}
        }
        // TODO
    }

    void local_to_global_matrix_v2(
        const Mesh& M,
        int t,
        const DenseMatrix& Ke,
        SparseMatrix& K )
    {
    	/*
        std::cout << "Ke -> K" << '\n';
        
        std::cout << "Point 0, numérotation globale : " << M.get_triangle_vertex_index(t, 0) << std::endl;
       	std::cout << "Point 1, numérotation globale : " << M.get_triangle_vertex_index(t, 1) << std::endl;
        std::cout << "Point 2, numérotation globale : " << M.get_triangle_vertex_index(t, 2) << std::endl;
        */
	///*
        int nb_lignes = Ke.height();
        int nb_colonnes = Ke.width();
        for ( int i = 0; i < nb_colonnes; i++){
        	for (int j = 0; j <= i; j++){
        		K.add(M.get_triangle_vertex_index(t, i), M.get_triangle_vertex_index(t, j), Ke.get(i, j));
        		
        	
        	}
        }
        //*/
        // TODO
    }
    
        void local_to_global_matrix(
        const Mesh& M,
        int t,
        const DenseMatrix& Ke,
        SparseMatrix& K )
    {
	for( int li = 0; li < Ke.height(); ++li ) {
		int i = M.get_triangle_vertex_index(t, li);
		for( int lj = 0; lj < Ke.width(); ++lj ) {
			int j = M.get_triangle_vertex_index(t, lj);
			K.add(i, j, Ke.get(li, lj));
		}
	}
    }

    void assemble_elementary_vector(
        const ElementMapping& elt_mapping,
        const ShapeFunctions& reference_functions,
        const Quadrature& quadrature,
        double (*source)(vertex),
        std::vector< double >& Fe )
    {
        /*
        std::cout << "compute elementary vector (source term)" << '\n';
        */
        // TODO
        double calcul;
        for (int i = 0; i < reference_functions.nb_functions(); i++) {
        	calcul = 0;
        	for (int q = 0; q < quadrature.nb_points(); q++) {
        		vertex v = quadrature.point(q);
        		
        		calcul += quadrature.weight(q)*reference_functions.evaluate(i, v)*elt_mapping.jacobian(v)*source(elt_mapping.transform(v));
        		// std::cout << "i : " << i << " , q : " << q << " calcul : " << calcul << std::endl;
        	}
        	Fe.push_back(calcul);        	
        }        		
    }

    void assemble_elementary_neumann_vector(
        const ElementMapping& elt_mapping_1D,
        const ShapeFunctions& reference_functions_1D,
        const Quadrature& quadrature_1D,
        double (*neumann)(vertex),
        std::vector< double >& Fe )
    {
    	/*
        std::cout << "compute elementary vector (neumann condition)" << '\n';
        */
        // TODO
        double calcul;
        for (int i = 0; i < reference_functions_1D.nb_functions(); i++) {
        	calcul = 0;
        	for (int q = 0; q < quadrature_1D.nb_points(); q++) {
        		vertex v = quadrature_1D.point(q);
        		calcul += quadrature_1D.weight(q)*reference_functions_1D.evaluate(i, v)*elt_mapping_1D.jacobian(v)*neumann(elt_mapping_1D.transform(v));
        		// std::cout << "i : " << i << " , q : " << q << " calcul : " << calcul << std::endl;
        	}        	
        	Fe.push_back(calcul);
    	}
    }

    void local_to_global_vector(
        const Mesh& M,
        bool border,
        int i,
        std::vector< double >& Fe,
        std::vector< double >& F )
    {
    	/*
        std::cout << "Fe -> F" << '\n';
        */
        int nb_lignes = Fe.size();
        if (border) {
        	for (int l = 0; l < nb_lignes; l++){
        		F[M.get_edge_vertex_index(i, l)] += Fe[l];
        	}	
        }
        else {
        	for (int l = 0; l < nb_lignes; l++){
        		// std::cout << "Fe[
        		F[M.get_triangle_vertex_index(i, l)] += Fe[l];
        	}
        
        // TODO
        }
    }

    void apply_dirichlet_boundary_conditions(
        const Mesh& M,
        const std::vector< bool >& attribute_is_dirichlet, /* size: nb of attributes */
        const std::vector< double >& values, /* size: nb of DOFs */
        SparseMatrix& K,
        std::vector< double >& F )
    {
        std::cout << "apply dirichlet boundary conditions" << '\n';
        // TODO
        int Pen_coeff = 10000;
        // La condition de dirichlet est imposée sur des segments, ces derniers peuvent avoir des points en commun, il faut donc gérer si on a déja pris en compte ou pas
        // Pour ne pas passer deux fois sur les points -> utilisation d'un vecteur de booléens
        std::vector <bool> processed_vertices(values.size(), false);
        assert(values.size() == M.nb_vertices());
        for (int edge = 0; edge < M.nb_edges(); edge++) {
        	int edge_attribute = M.get_edge_attribute(edge);
        	if (attribute_is_dirichlet[edge_attribute]) {
        		for (int vertice = 0; vertice < 2; vertice++){
        			int vertex_index = M.get_edge_vertex_index(edge, vertice);
        			if ( !processed_vertices[vertex_index]) {
        				processed_vertices[vertex_index] = true;
        				K.add(vertex_index, vertex_index, Pen_coeff);
        				F[vertex_index] += Pen_coeff*values[vertex_index];
        			}
        		}
        	}
        
        }
        
    }
/*
    void solve_poisson_problem(
            const Mesh& M,
            double (*diffusion_coef)(vertex),
            double (*source_term)(vertex),
            double (*dirichlet_fct)(vertex),
            double (*neumann_fct)(vertex),
            std::vector<double>& solution,
            int verbose )
    {
        std::cout << "solve poisson problem" << '\n';
        // TODO
        SparseMatrix K_glob(M.nb_vertices());
        std::vector<double> F_glob(M.nb_vertices());
        for (int triangle = 0; triangle < M.nb_triangles(); triangle++){
            ElementMapping mapping(M, false, triangle);
            ShapeFunctions shpfct(2,1);
            Quadrature quad = Quadrature::get_quadrature(2);
            DenseMatrix Ke;
            assemble_elementary_matrix(mapping, shpfct, quad, unit_fct, Ke);
            local_to_global_matrix(M, triangle, Ke, K_glob);
        }
        for (int vertice = 0; vertice < M.nb_vertices(); vertice++) {
        	ElementMapping mapping_1D(M, true, vertice);
        	ShapeFunctions shpfct_1D(1, 1);
        	Quadrature quad_1D = Quadrature::get_quadrature(2);
        	std::vector <double> Fe;
        	
        	
        
        }
        std::vector<bool> att_is_dirichlet(2, false);
        att_is_dirichlet[1] = true;
        mesh.set_attribute(unit_fct, 1, true);
        std::vector<double> imposed_values(M.nb_vertices());
        for( int i = 0; i < mesh.nb_vertices(); ++i){
            imposed_values[i] = xy_fct(M.get_vertex(i));
        }
        apply_dirichlet_boundary_conditions(M, att_is_dirichlet, imposed_values, K_glob, F_glob);
        solve(K_glob, F_glob, solution);
        std::string export_name = "pure_dirichlet";
        mesh.save(export_name+".mesh");
        save_solution(solution, export_name+".bb");
        
        
    }
*/
}
