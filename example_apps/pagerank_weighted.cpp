/**
 * @file
 * @author  Aapo Kyrola <akyrola@cs.cmu.edu>
 * @version 1.0
 *
 * @section LICENSE
 *
 * Copyright [2012] [Aapo Kyrola, Guy Blelloch, Carlos Guestrin / Carnegie Mellon University]
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 * http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 
 *
 * @section DESCRIPTION
 *
 * Simple pagerank implementation. Uses the basic vertex-based API for
 * demonstration purposes. A faster implementation uses the functional API,
 * "pagerank_functional".
 */

#include <string>
#include <fstream>
#include <cmath>

#define GRAPHCHI_DISABLE_COMPRESSION

#include "graphchi_basic_includes.hpp"
#include "util/toplist.hpp"
#include "graphchi_types.hpp"

using namespace graphchi;
 
#define THRESHOLD 1e-1    
float RANDOMRESETPROB = 0.2;

typedef float VertexDataType;
typedef struct weighted_data EdgeDataType;

struct PagerankProgram : public GraphChiProgram<VertexDataType, EdgeDataType> {
    
    /**
      * Called before an iteration starts. Not implemented.
      */
    void before_iteration(int iteration, graphchi_context &info) {
    }
    
    /**
      * Called after an iteration has finished. Not implemented.
      */
    void after_iteration(int iteration, graphchi_context &ginfo) {
    }
    
    /**
      * Called before an execution interval is started. Not implemented.
      */
    void before_exec_interval(vid_t window_st, vid_t window_en, graphchi_context &ginfo) {        
    }
    
    
    /**
      * Update the weigthed edge chivector
      * We first obtain the edge weight from the first element, sum them, then update the 
      * second item by eacg edge's weight
      */
    void update_edge_data(graphchi_vertex<VertexDataType, EdgeDataType> &v, float quota){
	    float sum = 0.0;
        for(int i=0; i < v.num_outedges(); i++) {
            graphchi_edge<EdgeDataType> * edge = v.outedge(i);
	        //We sum the weight values for each outgoing edges.
	        EdgeDataType val = edge->get_data();
	        sum += val.weight;
        }
	
        for(int i=0; i < v.num_outedges(); i++) {
            graphchi_edge<EdgeDataType> * edge = v.outedge(i);
            EdgeDataType val = edge->get_data();
            val.data= quota * val.weight / sum;
            edge->set_data(val);
#ifdef LOGOUTPUT
            if (v.id() == 1)
                std::cout <<  v.id() << " -> " << edge->vertex_id() << " with data: " << val.data<< " with weight " << val.weight << std::endl;
    #endif
             }
    }


    /**
      * Pagerank update function.
      */
    void update(graphchi_vertex<VertexDataType, EdgeDataType> &v, graphchi_context &ginfo) {
        float sum=0;
        if (ginfo.iteration == 0) {
            /* On first iteration, initialize vertex and out-edges. 
               The initialization is important,
               because on every run, GraphChi will modify the data in the edges on disk. 
             */
	        update_edge_data(v, 1.0);
            v.set_data(RANDOMRESETPROB); 
        } else {
            /* Compute the sum of neighbors' weighted pageranks by
               reading from the in-edges. */
            for(int i=0; i < v.num_inedges(); i++) {
		        EdgeDataType val= v.inedge(i)->get_data();
		        sum += val.data;
            }
            
            /* Compute my pagerank */
            float pagerank = RANDOMRESETPROB + (1 - RANDOMRESETPROB) * sum;
            
            /* Write my pagerank weighted by the weight of out-edges to
               each of my out-edges. */
	        update_edge_data(v, pagerank);
                
            /* Keep track of the progression of the computation.
               GraphChi engine writes a file filename.deltalog. */
            ginfo.log_change(std::abs(pagerank - v.get_data()));
            
            /* Set my new pagerank as the vertex value */
            v.set_data(pagerank); 
        }
    }
    
};

int main(int argc, const char ** argv) {
    graphchi_init(argc, argv);
    metrics m("pagerank");
    
    /* Parameters */
    std::string filename    = get_option_string("file"); // Base filename
    int niters              = get_option_int("niters", 4);
    bool scheduler          = false;                    // Non-dynamic version of pagerank.
    int ntop                = get_option_int("top", 20);
    RANDOMRESETPROB         = get_option_float("reset", RANDOMRESETPROB);
    
    /* Process input file - if not already preprocessed */
    int nshards             = convert_if_notexists<EdgeDataType>(filename, get_option_string("nshards", "auto"));

    /* Run */
    //graphchi_engine<float, float> engine(filename, nshards, scheduler, m); 
    graphchi_engine<VertexDataType, EdgeDataType> engine(filename, nshards, scheduler, m); 
    engine.set_modifies_inedges(false); // Improves I/O performance.
    PagerankProgram program;
    engine.run(program, niters);
        
    /* Output top ranked vertices to a file */
    std::ofstream ofile;
    ofile.open((filename + ".ntop.txt").c_str());
    std::vector< vertex_value<float> > top = get_top_vertices<float>(filename, ntop);
    std::cout << "Print top " << ntop << " vertices:" << std::endl;
    size_t num_vertices = get_num_vertices(filename);
    for(int i=0; i < (int)top.size(); i++) {
        if(i < 30)
            std::cout << (i+1) << ". " << top[i].vertex << "\t" << top[i].value / num_vertices << std::endl;
        ofile << top[i].vertex << "\t" << top[i].value /num_vertices << std::endl;
    }
    
    metrics_report(m);    
    return 0;
}

