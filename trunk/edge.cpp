#ifndef _EDGE_H_
#define _EDGE_H_

#include "edge.h"

// EDGE CONSTRUCTOR
Edge::Edge(Vec3f *vs, Vec3f *ve) {
  start_vertex = vs;
  end_vertex = ve;
  opposite = NULL;
  crease = 0;
}

// EDGE DESTRUCTOR
Edge::~Edge() { 
  // disconnect from the opposite edge
  if (opposite != NULL)
    opposite->opposite = NULL;
  // NOTE: the "prev" edge (which has a "next" pointer pointing to
  // this edge) will also be deleted as part of the triangle removal,
  // so we don't need to disconnect that
}

float Edge::Length() const 
{
  Vec3f diff = Vec3f(start_vertex->x() - end_vertex->x(), 
                     start_vertex->y() - end_vertex->y(), 
                     start_vertex->z() - end_vertex->z());
  return diff.Length();
}

#endif
