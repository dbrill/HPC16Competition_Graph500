/* -*- mode: C; mode: folding; fill-column: 70; -*- */
/* Copyright 2010-2011,  Georgia Institute of Technology, USA. */
/* See COPYING for license. */
#define _FILE_OFFSET_BITS 64
#define _THREAD_SAFE
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <assert.h>

#include "../compat.h"
#include "../graph500.h"

static int64_t maxvtx, maxdeg, nIJ;
static const struct packed_edge * restrict IJ;
static int64_t * restrict head, * restrict deg, * restrict next;

struct vertex {
  int64_t value;
  int64_t depth;
};
static int64_t count = 0;
int
create_graph_from_edgelist (struct packed_edge *IJ_in, int64_t nedge)
{

  int err = 0;

  IJ = IJ_in; //edgelist single-D array
  nIJ = nedge; // number of edges
  maxvtx = -1; //max vertex
  maxdeg = -1; //max depth
  struct vertex* ptr, *p;

  int64_t k;
  for (k = 0; k < nedge; ++k) {
        //THIS MIGHT BE WRONG. IMAGINE BOTH ARE LARGER AND V1 < V0
    if (get_v0_from_edge(&IJ[k]) > maxvtx)
      maxvtx = get_v0_from_edge(&IJ[k]);
    if (get_v1_from_edge(&IJ[k]) > maxvtx)
      maxvtx = get_v1_from_edge(&IJ[k]);
  }
    //size accounts for undirectedness meaning two "edge_entries" for each connection
    //we're mallocing with the sizeof(our struct) now for efficiency and performance and optimization

    head = malloc((sizeof(struct vertex)*(maxvtx+1)) + (2*nIJ) * sizeof(int64_t));
    //if malloc fails return error code -1
    //if (head){
    if (!head) return -1;
    ptr = (struct vertex*)head;
    //maxvtx+1 = out of bounds
    //This size allocation provides readily accessible checking of bounds
          //*deg = &head[maxvtx+1];*
    next = &head[2*maxvtx+1];

    //initialize values for head and degree
    //degree is initially zero for zeroth level
    //each head node has yet to be visited thus -1

    p = ptr;
    for (k = 0; k <= maxvtx; ++k, p++) {
      p->value = -1;
      p->depth = 0;
        //*deg[k] = 0;*
    }
//printf("\nEDGES\n");
    for (k = 0; k < nedge; ++k) {
      //get both vertices from edge
      const int64_t i = get_v0_from_edge(&IJ[k]);
      const int64_t j = get_v1_from_edge(&IJ[k]);
      int64_t t_head, t;

    //  printf("edge%d = %d<-->%d\t", k,i,j);

        // so make sure edge is existing i.e. != -1
        // also make sure we don't have self looping vertice
        //initialize deg
        //remove cycles and self edges
      if (i >= 0 && j >= 0 && i != j) {
        // translate the edges in IJ into a sequential array, next
        //  next = [V0_E0,V1_E0,V0_E1,V1_E1...V0_nIJ,V1_nIJ]

        t = 2*k;
        next[t] = -1;
        next[t++] = -1;
         /* Point at the *other* end. */


        // head[n] stores the index number of next, which contains -1
        //(end of the sequence) or another index of next.
        p = ptr+i;
        t_head = p->value;
        p->value = t;
        assert (t_head < 2*nIJ);
        next[t] = t_head;
        p->depth+=1;

        --t;
        p = ptr+j;
        t_head = p->value;
        p->value = t;
        assert (t_head < 2*nIJ);
        next[t] = t_head;
        p->depth+=1;


      }
    }
      //run through each vertex and find max degree aka depth of tree
  //  ptr = (struct vertex *)head;
    p = ptr;
    for (int64_t kg = 0; kg <= maxvtx; ++kg, p++){
      if (p->depth > maxdeg){
        maxdeg = p->depth;
      }
    }


    int64_t t;
  //  struct vertex * ptr = (struct vertex*)head;
  //  struct vertex * p;
  //  printf("\nHEAD SEQ\n");
    /*for(t=0; t < maxvtx+1; t++){
      p = &ptr[t];
  		if (p->value < 0)
  			continue;
  		else
  			printf("Vertex:%d, Value:%d, Depth:%d\n", t, p->value, p->depth);
  	} */

    return err;

}

int
make_bfs_tree (int64_t *bfs_tree_out, int64_t *max_vtx_out,
	       int64_t srcvtx)
{
  printf("******COUNT: %d*****", count);
  printf("Variables-----\nsrcvtx:%d, head[0]:%d, head[1]:%d\n", srcvtx, head[0], head[1]);
printf("Variables-----\nsrcvtx:%d, head[2]:%d, head[3]:%d\n", srcvtx, head[2], head[3]);
  count++;
  int64_t * restrict bfs_tree = bfs_tree_out;
  int err = 0;
  //# vertexes
  const int64_t nv = maxvtx+1;

  int64_t k, k1, k2;
  int64_t * restrict vlist;
  struct vertex * ptr = (struct vertex*)head;
struct vertex * p;

  *max_vtx_out = maxvtx;
  //root node
  bfs_tree[srcvtx] = srcvtx;

  vlist = malloc (nv * sizeof (*vlist));
  if (!vlist) return -1;
  //whyyy
  bfs_tree[srcvtx] = srcvtx;
  k1 = 0; k2 = 1;
  vlist[0] = srcvtx;

//set everything except for the source vertex to -1
//initialization of tree list!
          //***definitely a place for parallelization***
  for (k = 0; k < srcvtx; ++k)
    bfs_tree[k] = -1;
  for (k = srcvtx+1; k < nv; ++k)
    bfs_tree[k] = -1;

  while (k1 != k2){
    int64_t newk2 = k2;
    int64_t k;

  //  printf("\nVLIST for k%d\n", k1);

    for (k = k1; k < k2; ++k) {
      const int64_t parent = vlist[k];
      p = ptr+parent;
      //don't reassign parent until you've exhausted connections of parent
      int64_t val = p->value;
      while (val >= 0) {
        //if it's a right vertex head to the left if it's a left vertex head to the right
        // keeps search on the same degree
        const int64_t newv = ((val % 2) ? get_v1_from_edge(&IJ[val / 2]) : get_v0_from_edge(&IJ[val / 2]));
        //if there's no vertex in that spot
        if (bfs_tree[newv] < 0) {
          bfs_tree[newv] = parent;
          vlist[newk2++] = newv;
        }
        val = next[val];
      }
      //if we've run out of reachable nodes (p = -1) newk2 will not increment
      //therefore k2 will catch up, no longer being < k2 loop will break
      k1 = k2;
      k2 = newk2;
    }

  /*  int64_t t, v;
    for(t=0; t < 256; t++){
      v = vlist[t];
  		if (v < 0)
  			continue;
  		else
  			printf("Vlist:%d, Value:%d\n", t, v);
  	} */
  }
  //free malloced vlist
  free (vlist);
  //return to stderr hopefully 0
  return err;
}

void
destroy_graph (void)
{
  free (head);
}
