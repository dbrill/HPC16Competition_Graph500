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

//static int64_t count = 0;

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
    int myerr = 0;
    int64_t k, i;

    OMP("omp for")
				//0 - number_vertexes
      for (i = 0; i < nv; ++i){
					//set level to 0 if it's the root otherwise set it to -1 or unvisited
					//might be able to eliminate a jump and move instruction here
				level[i] = (i == root? 0 : -1);
			}

    OMP("omp for") MTA("mta assert parallel") MTA("mta use 100 streams")
				//0 --> vertices
      for (i = 0; i < nv; ++i) {
					//if the node has no parent or has already had its level set
					//move to next iteration
				if(level[i] != -1){
					continue;
				}
				if(i != root && bfs_tree[i] >= 0 && !myerr){
					myerr = err;

					int64_t nhops = 0, temp = i;

					//while temp's depth still hasn't been set
					while(level[temp]<0 && temp > -1 && temp != bfs_tree[temp]){ //might need to add "&& bfs_tree[temp] >= 0 && nhop < nv"
						temp = bfs_tree[temp];
						nhops++;
					}
					k = i;
					while(k != temp){
						level[k] = level[temp] + nhops;
						k = bfs_tree[k];
						nhops --;
					}
					if(temp < 0)
						myerr = -1;
					if(level[root] != 0)
						myerr = -2;
					if(!myerr){
							    /* Internal check to catch mistakes in races... */
						#if !defined(NDEBUG)
							    nhops = 0;
							    temp = i;
							    int64_t lastlvl = level[i]+1;
							    while (level[temp] > 0) {
							      assert (lastlvl == 1 + level[temp]);
							      lastlvl = level[temp];
							      temp = bfs_tree[temp];
							      ++nhops;
							    }
						#endif
							  }
							}
						}
						if (myerr) { err = myerr;	OMP("omp flush (err)"); }
				}
  return err;
}

int64_t
verify_bfs_tree (int64_t *bfs_tree_in, int64_t max_bfsvtx,
		 int64_t root,
		 const struct packed_edge *IJ_in, int64_t nedge)
{

  //count++;
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
	//int64_t z;
	/*for(z = 0; z < nv; z++){
		printf("bfs_tree[%d] = %d\n",z, bfs_tree[z]);
	}*/

  OMP("omp parallel shared(err)") {
    int64_t k;
    int myerr = 0;
    OMP("omp for")
      for (k = 0; k < nv; ++k)
				//set all vertices to 0 (unseen)
				seen_edge[k] = 0;

    OMP("omp for reduction(+:nedge_traversed)")
    MTA("mta assert parallel") MTA("mta use 100 streams")
      for (k = 0; k < nedge; ++k) {
					const int64_t v0 = get_v0_from_edge (&IJ[k]);
					const int64_t v1 = get_v1_from_edge (&IJ[k]);
				//	printf("edge %d: %d<-->%d\n", k, v0, v1);
					int64_t lvldiff;
					myerr = err;


					if ((v0 < 0 || v1 < 0) || (v0 > max_bfsvtx && v1 > max_bfsvtx)) continue;
					if (v0 > max_bfsvtx || v1 > max_bfsvtx) myerr = -10;
					if (myerr) { err = myerr; OMP("omp flush(err)"); continue;}
							//error thrown OR v0 && v1 > max_bfsvtx
				//	if (myerr || v0 > max_bfsvtx /* both v0 & v1 are on the same sv0de of max_bfsvtx */)
					//  continue;

					/* All neighbors must be in the tree. */

					if(bfs_tree[v0] < 0 && bfs_tree[v1] < 0) continue;
					if (bfs_tree[v0] < 0 || bfs_tree[v1] < 0) myerr = -11;
					if (myerr) { err = myerr; OMP("omp flush(err)"); continue;}

					/* Both v0 and v1 are v0n the tree, count as a traversed edge.
					   NOTE: This counts self-edges and repeated edges.  They're
					   part of the input data.
					*/
					++nedge_traversed;
					/* Mark seen tree edges. */
					if (v0 != v1) {
					//  printf("bfs_tree[%d](v0) = %d\n", v0, bfs_tree[v0]);
						//printf("bfs_tree[%d](v1) = %d\n", v1, bfs_tree[v1]);
						if (bfs_tree[v0] == v1)
					    seen_edge[v0] = 1;
					  if (bfs_tree[v1] == v0)
					    seen_edge[v1] = 1;
					}
					lvldiff = level[v0] - level[v1];
					/* Check that the levels differ by no more than one. */
					if ((lvldiff * lvldiff) > 1)
					  myerr = -14;
					if (myerr) { err = myerr; OMP("omp flush(err)"); }
      }

    if (!myerr) {
      /* Check that every BFS edge was seen and that there's only one root. */
      OMP("omp for") MTA("mta assert parallel") MTA("mta use 100 streams")
			for (k = 0; k < nv; ++k) {
			  myerr = err;
			  if (!myerr && k != root) {
			    if (bfs_tree[k] >= 0 && !seen_edge[k]){
			      myerr = -15;
					//	printf("***root = %d bfs_tree[%d] = %d and seen_edge=%d***\n", root, k, bfs_tree[k],  seen_edge[k]);
				} else 	//printf("***root = %d bfs_tree[%d] = %d and seen_edge=%d***\n", root, k, bfs_tree[k],  seen_edge[k]);
			    if (bfs_tree[k] == k)
			      myerr = -16;
			    if (myerr) { err = myerr; OMP("omp flush(err)"); }
			  }
			}
    }
  }
 done:
	//printf("YOU PASSED WITH ERROR: %d", err);
  xfree_large (seen_edge);
  if (err) return err;
  return nedge_traversed;
}
