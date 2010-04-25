#ifndef EDGE_H
#define EDGE_H

#include <limits.h>
#include <stdio.h>
#include <assert.h>
#include "vectors.h"

// ===================================================================
// half-edge data structure

class Edge {

public:

  // ========================
  // CONSTRUCTORS & DESTRUCTOR
  Edge(Vec3f vs, Vec3f ve);
  ~Edge();

  // =========
  // ACCESSORS
  Vec3f getStartVertex() const { /*assert (start_vertex != NULL);*/ return start_vertex; }
  Vec3f getEndVertex() const { /*assert (end_vertex != NULL);*/ return end_vertex; }
  Edge* getOpposite() const { return opposite; }
  float getCrease() const { return crease; }
  float Length() const;

  // =========
  // MODIFIERS
  void setOpposite(Edge *e) {
    assert (opposite == NULL);
    assert (e != NULL);
    assert (e->opposite == NULL);
    opposite = e;
    e->opposite = this;
  }
  void clearOpposite() {
    if (opposite == NULL) return;
    assert (opposite->opposite == this);
    opposite->opposite = NULL;
    opposite = NULL;
  }
  void setCrease(float c) { crease = c; }

private:

  Edge(const Edge&) { assert(0); }
  Edge& operator=(const Edge&) { assert(0); exit(0); }

  // ==============
  // REPRESENTATION
  // in the half edge data adjacency data structure, the edge stores everything!
  // note: it's technically not necessary to store both vertices, but it makes
  //   dealing with non-closed meshes easier
  Vec3f start_vertex;
  Vec3f end_vertex;
  Edge *opposite;
  // the crease value is an extra field used during subdivision
  float crease;
};

// ===================================================================

#endif
