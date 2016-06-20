/* -*- mode: C; mode: folding; fill-column: 70; -*- */
/* Copyright 2010-2011,  Georgia Institute of Technology, USA. */
/* See COPYING for license. */
#include "compat.h"
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include <assert.h>

#include "xalloc.h"
#include "verify.h"

static int64_t count = 0;

static int
compute_levels (int64_t * level,
		int64_t nv, const int64_t * restrict bfs_tree, int64_t root)
{
  int err = 0;
/*	printf("\nBFS_TREE VERIFY\n");
	int64_t i = 0;
	for(i=0; i < nv; i++){
		//if (bfs_tree[i] < 0)
		//	continue;
		//else
			printf("Vertex:%d, Next:%d\n",i, bfs_tree[i]);
	}*/
  OMP("omp parallel shared(err)") {
    int terr;
    int64_t k;

    OMP("omp for")
				//0 - number_vertexes
      for (k = 0; k < nv; ++k)
				level[k] = (k == root? 0 : -1);

    OMP("omp for") MTA("mta assert parallel") MTA("mta use 100 streams")
      for (k = 0; k < nv; ++k) {
						//if node has been visited
					if (level[k] >= 0) continue;
						//temp error val
					terr = err;
						//if no errors thus far, and verte is not edgeless and not root
					if (!terr && bfs_tree[k] >= 0 && k != root) {
					  int64_t parent = k;
					  int64_t nhop = 0;
					  /* Run up the tree until we encounter an already-leveled vertex. */
							//if parent still has an edge and has yet to be visited
							//and traversal < num_verts
							//for loop to calculate the number of hops in the tree
					  while (parent >= 0 && level[parent] < 0 && nhop < nv) {
								//no cycle
					    assert (parent != bfs_tree[parent]);
								//move parent down to next child
					    parent = bfs_tree[parent];
								//increase number of moves
					    ++nhop;
					  }

					  if (nhop >= nv) terr = -1; /* Cycle. */
					  if (parent < 0) terr = -2; /* Ran off the end. */

					  if (!terr) {
					    /* Now assign levels until we meet an already-leveled vertex */
					    /* NOTE: This permits benign races if parallelized. */
									//USE OMP CRITICAL HERE!

							// reduce nhop by zero if parent is root
							//reduce nhop by one if parent is not root
					    nhop += level[parent];
					    parent = k;

								//while loop to assign the level to nodes in the tree
								//while level unvisited
					    while (level[parent] < 0) {
					      assert (nhop > 0);
								//modify node's level to reflect the current level in the tree
					      level[parent] = nhop--;
								//move down to the next child
					      parent = bfs_tree[parent];
					    }
							//make sure that tree was fully traversed
					    assert (nhop == level[parent]);

					    /* Internal check to catch mistakes in races... */
				#if !defined(NDEBUG)
					    nhop = 0;
					    parent = k;
					    int64_t lastlvl = level[k]+1;
					    while (level[parent] > 0) {
					      assert (lastlvl == 1 + level[parent]);
					      lastlvl = level[parent];
					      parent = bfs_tree[parent];
					      ++nhop;
					    }
				#endif
					  }
					}
					if (terr) { err = terr;	OMP("omp flush (err)"); }
			}
  }
  return err;
}

int64_t
verify_bfs_tree (int64_t *bfs_tree_in, int64_t max_bfsvtx,
		 int64_t root,
		 const struct packed_edge *IJ_in, int64_t nedge)
{
	count++;
	//printf("***RUN NUMBER %d***", count);
  int64_t * restrict bfs_tree = bfs_tree_in;
  const struct packed_edge * restrict IJ = IJ_in;

  int err;
  int64_t nedge_traversed;
  int64_t * restrict seen_edge, * restrict level;

  const int64_t nv = max_bfsvtx+1;

  /*
    This code is horrifically contorted because many compilers
    complain about continue, return, etc. in parallel sections.
  */

  if (root > max_bfsvtx || bfs_tree[root] != root)
    return -999;

  err = 0;
  nedge_traversed = 0;
  seen_edge = xmalloc_large (2 * (nv) * sizeof (*seen_edge));
  level = &seen_edge[nv];

  err = compute_levels (level, nv, bfs_tree, root);

  if (err) goto done;
	int64_t z;
	for(z = 0; z < nv; z++){
		//printf("bfs_tree[%d] = %d\n",z, bfs_tree[z]);
	}

  OMP("omp parallel shared(err)") {
    int64_t k;
    int terr = 0;
    OMP("omp for")
      for (k = 0; k < nv; ++k)
				seen_edge[k] = 0;

    OMP("omp for reduction(+:nedge_traversed)")
    MTA("mta assert parallel") MTA("mta use 100 streams")
      for (k = 0; k < nedge; ++k) {
					const int64_t i = get_v0_from_edge (&IJ[k]);
					const int64_t j = get_v1_from_edge (&IJ[k]);
				//	printf("edge %d: %d<-->%d\n", k, i, j);
					int64_t lvldiff;
					terr = err;


					if (i < 0 || j < 0) continue;
					if (i > max_bfsvtx && j <= max_bfsvtx) terr = -10;
					if (j > max_bfsvtx && i <= max_bfsvtx) terr = -11;
					if (terr) { err = terr; OMP("omp flush(err)"); }
							//error thrown OR i && j > max_bfsvtx
					if (terr || i > max_bfsvtx /* both i & j are on the same side of max_bfsvtx */)
					  continue;

					/* All neighbors must be in the tree. */


					if (bfs_tree[i] >= 0 && bfs_tree[j] < 0) terr = -12;
					if (bfs_tree[j] >= 0 && bfs_tree[i] < 0) terr = -13;
					if (terr) { err = terr; OMP("omp flush(err)"); }
					if (terr || bfs_tree[i] < 0 /* both i & j have no neighbors */)
					  continue;

					/* Both i and j are in the tree, count as a traversed edge.
					   NOTE: This counts self-edges and repeated edges.  They're
					   part of the input data.
					*/
					++nedge_traversed;
					/* Mark seen tree edges. */
					if (i != j) {
					//  printf("bfs_tree[%d](i) = %d\n", i, bfs_tree[i]);
						//printf("bfs_tree[%d](j) = %d\n", j, bfs_tree[j]);

					  if (bfs_tree[i] == j)
					    seen_edge[i] = 1;
					  if (bfs_tree[j] == i)
					    seen_edge[j] = 1;
					}
					lvldiff = level[i] - level[j];
					/* Check that the levels differ by no more than one. */
					if (lvldiff > 1 || lvldiff < -1)
					  terr = -14;
					if (terr) { err = terr; OMP("omp flush(err)"); }
      }

    if (!terr) {
      /* Check that every BFS edge was seen and that there's only one root. */
      OMP("omp for") MTA("mta assert parallel") MTA("mta use 100 streams")
			for (k = 0; k < nv; ++k) {
			  terr = err;
			  if (!terr && k != root) {
			    if (bfs_tree[k] >= 0 && !seen_edge[k]){
			      terr = -15;
					//	printf("***root = %d bfs_tree[%d] = %d and seen_edge=%d***\n", root, k, bfs_tree[k],  seen_edge[k]);
				} else 	//printf("***root = %d bfs_tree[%d] = %d and seen_edge=%d***\n", root, k, bfs_tree[k],  seen_edge[k]);
			    if (bfs_tree[k] == k)
			      terr = -16;
			    if (terr) { err = terr; OMP("omp flush(err)"); }
			  }
			}
    }
  }
 done:
	printf("YOU PASSED WITH ERROR: %d", err);
  xfree_large (seen_edge);
  if (err) return err;
  return nedge_traversed;
}
