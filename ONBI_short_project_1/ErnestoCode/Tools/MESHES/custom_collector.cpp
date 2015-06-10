
 // Begin License:
// Copyright (C) 2006-2008 Tobias Sargeant (tobias.sargeant@gmail.com).
// All rights reserved.
//
// This file is part of the Carve CSG Library (http://carve-csg.com/)
//
// This file may be used under the terms of the GNU General Public
// License version 2.0 as published by the Free Software Foundation
// and appearing in the file LICENSE.GPL2 included in the packaging of
// this file.
//
// This file is provided "AS IS" with NO WARRANTY OF ANY KIND,
// INCLUDING THE WARRANTIES OF DESIGN, MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE.
// End:


#if defined(HAVE_CONFIG_H)
#  include <carve_config.h>
#endif
#include "/usr/local/src/carve-1.4.0/common/geometry.hpp"
#include <carve/csg.hpp>
#include "myMEX.h"
#include "mex.h"
#include "MESH2carvePolyhedron.h"
#include "carvePolyhedron2MESH.h"

#include <carve/input.hpp>

#include <fstream>
#include <string>
#include <utility>
#include <set>

#include <time.h>


#define DIM 60

class Between : public carve::csg::CSG::Collector {
  Between();
  Between(const Between &);
  Between &operator=(const Between &);

public:
  std::list<carve::poly::Face<3> > faces;
  const carve::poly::Polyhedron *src_a;
  const carve::poly::Polyhedron *src_b;

  Between(const carve::poly::Polyhedron *_src_a,
          const carve::poly::Polyhedron *_src_b) : carve::csg::CSG::Collector(), src_a(_src_a), src_b(_src_b) {
  }

  virtual ~Between() {
  }

  virtual void collect(carve::csg::FaceLoopGroup *grp, carve::csg::CSG::Hooks &hooks) {
	  if( (src_a!= NULL) && (src_b!= NULL))
		  mexPrintf("poly validos .\n");
    if (grp->face_loops.head->orig_face->owner != src_a) return;
    if (grp->classificationAgainst(src_b, 1) != carve::csg::FACE_IN) return;
    if (grp->classificationAgainst(src_b, 0) != carve::csg::FACE_OUT) return;

    for (carve::csg::FaceLoop *f = grp->face_loops.head; f; f = f->next) {
    	mexPrintf("Entras aqui??\n");
      faces.push_back(carve::poly::Face<3>());
      faces.back().init(f->orig_face, f->vertices, false);
    }
  }

  virtual carve::poly::Polyhedron *done(carve::csg::CSG::Hooks &hooks) {
    return new carve::poly::Polyhedron(faces);
  }
};


void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	ALLOCATES();
	carve::poly::Polyhedron *MESH1 = NULL;
	carve::poly::Polyhedron *MESH2 = NULL;
	carve::poly::Polyhedron *MESH_OUT = NULL;

	 //Procesamiento
	  MESH1 = MESH2carvePolyhedron( prhs[0] , 1);
	  MESH2 = MESH2carvePolyhedron( prhs[1] , 0);

	  if(MESH1->faces[38].owner != MESH2)
		  mexPrintf("Owner de la Cara 38 NO es MESH2\n");

	  if(MESH1->faces[38].owner == MESH1)
	  		mexPrintf("Owner de la Cara 38 es MESH1\n");

	  plhs[0] = carvePolyhedron2MESH( MESH1 , 1);

  Between between_collector(MESH1, MESH2);
  MESH_OUT = carve::csg::CSG().compute(MESH1, MESH2, between_collector, NULL, carve::csg::CSG::CLASSIFY_NORMAL);

  //plhs[0] = carvePolyhedron2MESH( MESH_OUT );


  EXIT:
  	  delete MESH1;
  	  delete MESH2;
  	  //if(TREE != NULL) delete TREE;
  	  if(MESH_OUT != NULL) delete MESH_OUT;
  	  myFreeALLOCATES();

  return;
}



