#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <omp.h>

#include "rbindex.h"

#define NUM_THREADS 4

/*void *rbindex_alloc(struct libavl_allocator *allocator, size_t size)
{
	struct rbindex_cache *nc;

	return malloc(size);
}

void rbindex_free(struct libavl_allocator *allocator, void *block)
{
	free(block);
}*/

void *rbindex_alloc(struct libavl_allocator *allocator, size_t size)
{
	return cache_alloc(allocator->allocator_data);
}

void rbindex_free(struct libavl_allocator *allocator, struct rb_node *nd)
{
	cache_free(allocator->allocator_data, nd);
}

void rbindex_rb_insert (struct rbindex *idx, void *item)
{
  struct rb_node *pa[RB_MAX_HEIGHT]; /* Nodes on stack. */
  unsigned char da[RB_MAX_HEIGHT];   /* Directions moved from stack nodes. */
  struct rb_table *tree;
  int k;                             /* Stack height. */

  struct rb_node *p; /* Traverses tree looking for insertion point. */
  struct rb_node *n; /* Newly inserted node. */

  tree = idx->tree;
  assert (tree != NULL && item != NULL);

  pa[0] = (struct rb_node *) &tree->rb_root;
  da[0] = 0;
  k = 1;
  for (p = tree->rb_root; p != NULL; p = p->rb_link[da[k - 1]])
    {
      int cmp = tree->rb_compare (item, p->rb_data, tree->rb_param);
      if (cmp == 0)
        return;

      pa[k] = p;
      da[k++] = cmp > 0;
    }

  n = pa[k - 1]->rb_link[da[k - 1]] =
    tree->rb_alloc->libavl_malloc (tree->rb_alloc, sizeof *n);
  if (n == NULL)
    return;

  n->rb_data = item;
  n->rb_link[0] = n->rb_link[1] = NULL;
  n->rb_color = RB_RED;
  tree->rb_count++;
  tree->rb_generation++;

  /* Determine next element. */
  if(da[k - 1] == 0){
	  list_insbefore(GET_LIST(pa[k - 1]->rb_data), item);
  }else{
	  list_insafter(GET_LIST(pa[k - 1]->rb_data), item);
  }
  idx->n++;

  while (k >= 3 && pa[k - 1]->rb_color == RB_RED)
    {
      if (da[k - 2] == 0)
        {
          struct rb_node *y = pa[k - 2]->rb_link[1];
          if (y != NULL && y->rb_color == RB_RED)
            {
              pa[k - 1]->rb_color = y->rb_color = RB_BLACK;
              pa[k - 2]->rb_color = RB_RED;
              k -= 2;
            }
          else
            {
              struct rb_node *x;

              if (da[k - 1] == 0)
                y = pa[k - 1];
              else
                {
                  x = pa[k - 1];
                  y = x->rb_link[1];
                  x->rb_link[1] = y->rb_link[0];
                  y->rb_link[0] = x;
                  pa[k - 2]->rb_link[0] = y;
                }

              x = pa[k - 2];
              x->rb_color = RB_RED;
              y->rb_color = RB_BLACK;

              x->rb_link[0] = y->rb_link[1];
              y->rb_link[1] = x;
              pa[k - 3]->rb_link[da[k - 3]] = y;
              break;
            }
        }
      else
        {
          struct rb_node *y = pa[k - 2]->rb_link[0];
          if (y != NULL && y->rb_color == RB_RED)
            {
              pa[k - 1]->rb_color = y->rb_color = RB_BLACK;
              pa[k - 2]->rb_color = RB_RED;
              k -= 2;
            }
          else
            {
              struct rb_node *x;

              if (da[k - 1] == 1)
                y = pa[k - 1];
              else
                {
                  x = pa[k - 1];
                  y = x->rb_link[0];
                  x->rb_link[0] = y->rb_link[1];
                  y->rb_link[1] = x;
                  pa[k - 2]->rb_link[1] = y;
                }

              x = pa[k - 2];
              x->rb_color = RB_RED;
              y->rb_color = RB_BLACK;

              x->rb_link[1] = y->rb_link[0];
              y->rb_link[0] = x;
              pa[k - 3]->rb_link[da[k - 3]] = y;
              break;
            }
        }
    }
  tree->rb_root->rb_color = RB_BLACK;
}

void rbindex_rb_clear(struct rbindex *eidx)
{
	struct rbindex_cache *nc;

	/* Clear free list. */
	cache_clear(eidx->nc);

//	eidx->n = 0;
	eidx->tree->rb_count = 0;
	eidx->tree->rb_generation = 0;
}

static inline int find_root(int n, int h)
{
	int half, quad, r;

	half = 1 << (h - 1);
	quad = half >> 1;
	r = n - half + 1;
	if(r < quad)
		return quad - 1 + r;
	else
		return half - 1;
}

struct cbt_info {
	int beg;
	int end;
	int root;
	int h;
};

static inline void cbt_set(struct cbt_info *t, int root, int beg, int end, int h)
{
	t->root = root;
	t->beg = beg;
	t->end = end;
	t->h = h;
}

#define QUEUE_PUSH(i) queue[(*k)++] = i
#define QUEUE_POP() queue[(*j)++]

static inline void __cbt_addchild(struct cbt_info *tree, int n, int *queue, int *j, int *k)
{
	int root, beg, end, r, h, half, quad, i, left, right;

	i = QUEUE_POP();
	beg = tree[i].beg;
	end = tree[i].end;
	root = tree[i].root;
	h = tree[i].h;

	half = 1 << (h - 1);
	quad = half >> 1;
	r = end - beg - half + 1;

	root = tree[i].root;
	// Calculate index of left child.
	if(*k < n){
		left = i * 2 + 1;
		cbt_set(&tree[left], beg + find_root(root - beg, h - 1), beg, root, h - 1);
		QUEUE_PUSH(left);
	}

	// Calculate index of left child.
	if(*k < n){
		right = left + 1;
		if(r < quad)
			cbt_set(&tree[right], root + 1 + find_root(end - root - 1, h - 2), root + 1, end, h - 2);
		else
			cbt_set(&tree[right], root + 1 + find_root(end - root - 1, h - 1), root + 1, end, h - 1);
		QUEUE_PUSH(right);
	}
}

/*void __complete_binary_tree(struct cbt_info *tree, int start)
{
	int *queue, j, k, n;

	n = tree[start].end - tree[start].beg;
	queue = malloc(sizeof(int) * n);
	queue[0] = start;
	j = 0;
	k = 1;
	do{
		__cbt_addchild(tree, n, queue, &j, &k);
	}while(k < n);
	free(queue);
}

// Calculate indices of complete binary tree from an ordered array
void complete_binary_tree_parallel(int nnodes, int *map)
{
	struct cbt_info *tree;
	int i, j, k, h0, *queue;

	tree = malloc(sizeof(struct cbt_info) * nnodes);
	queue = malloc(sizeof(int) * nnodes);

	h0 = 0;
	while((1 << h0) <= nnodes) h0++;

	cbt_set(&tree[0], find_root(nnodes, h0), 0, nnodes, h0);
	queue[0] = 0;
	i = 0;
	k = 1;
	do{
		__cbt_addchild(tree, nnodes, queue, &i, &k);
	}while(k < nnodes && k < 7);

#pragma omp parallel for schedule(static,2) num_threads(2)
	for(j = 3; j < 7; j++){
		__complete_binary_tree(tree, j);
	}

#pragma omp parallel for schedule(static,2) num_threads(2)
	for(i = 0; i < nnodes; i++)
		map[tree[i].root] = i;
//		map[i] = tree[i].root;

	free(tree);
	free(queue);
}*/

// Calculate indices of complete binary tree from an ordered array
void complete_binary_tree(int nnodes, int *map)
{
	struct cbt_info *tree;
	int i, j, k, h0, *queue;

	tree = malloc(sizeof(struct cbt_info) * nnodes);
	queue = malloc(sizeof(int) * nnodes);

	h0 = 0;
	while((1 << h0) <= nnodes) h0++;

	cbt_set(&tree[0], find_root(nnodes, h0), 0, nnodes, h0);
	queue[0] = 0;
	i = 0;
	k = 1;
	do{
		__cbt_addchild(tree, nnodes, queue, &i, &k);
	}while(k < nnodes);

	for(i = 0; i < nnodes; i++)
		map[tree[i].root] = i;
//		map[i] = tree[i].root;

	free(tree);
	free(queue);
}

// Rebuild tree index from sorted list
void rbindex_rebuild_tree(struct rbindex *eidx)
{
	int i, j, nnodes, h, half, *tree;
	struct rb_node **nodes;
	struct list_head *l;
	void **objs;

	nnodes = eidx->n;
	rbindex_rb_clear(eidx);
	if(nnodes >= eidx->nc->cache_size)
		cache_resize(eidx->nc, nnodes - eidx->nc->cache_size);
	nodes = (struct rb_node **)eidx->nc->objs;

	objs = malloc(sizeof(void *) * nnodes);
	l = eidx->ls.front;
	for(i = 0; i < nnodes; i++){
		objs[i] = GET_OBJ(l);
		l = l->next;
	}

	// Build complete binary tree
	tree = malloc(sizeof(int) * nnodes);
	complete_binary_tree(nnodes, tree);
	h = 0;
	while(1 << h <= nnodes) h++;
	half = 1 << (h - 1);

	for(j = nnodes - 1; j >= half - 1; j--){
		nodes[j]->rb_color = RB_RED;
		nodes[j]->rb_data = objs[tree[j]];
		nodes[j]->rb_link[0] = nodes[j]->rb_link[1] = NULL;
	}

	for(; j > (nnodes - 1) / 2; j--){
		nodes[j]->rb_color = RB_BLACK;
		nodes[j]->rb_data = objs[tree[j]];
		nodes[j]->rb_link[0] = nodes[j]->rb_link[1] = NULL;
	}

	if(j == nnodes / 2){
		nodes[j]->rb_color = RB_BLACK;
		nodes[j]->rb_data = objs[tree[j]];
		nodes[j]->rb_link[0] = nodes[j]->rb_link[1] = NULL;

	}else{
		nodes[j]->rb_color = RB_BLACK;
		nodes[j]->rb_data = objs[tree[j]];
		nodes[j]->rb_link[0] = nodes[(j << 1) + 1];
		nodes[j]->rb_link[1] = NULL;
	}
	j--;

	for(; j >= 0; j--){
		nodes[j]->rb_color = RB_BLACK;
		nodes[j]->rb_data = objs[tree[j]];
		nodes[j]->rb_link[0] = nodes[(j << 1) + 1];
		nodes[j]->rb_link[1] = nodes[(j << 1) + 2];
	}

	eidx->tree->rb_root = nodes[0];
	nodes[0]->rb_color = RB_BLACK;
	eidx->nc->maxnodes = nnodes;

	free(objs);
	free(tree);
}

void rbindex_destroy(struct rbindex *eidx)
{
	cache_destroy(eidx->nc);
	free(eidx->tree);
	free(eidx);
}

struct rbindex *rbindex_create(rb_comparison_func *compar, int cache_size)
{
	struct rbindex *eidx;
	int i, *pid;
	char *nd;

	eidx = malloc(sizeof(struct rbindex));

	eidx->nc = cache_create(sizeof(struct rb_node), cache_size, NULL, NULL);

	eidx->allocator.allocator_data = eidx->nc;
	eidx->allocator.libavl_malloc = rbindex_alloc;
	eidx->allocator.libavl_free = rbindex_free;

	eidx->flags = 0;
	eidx->n = 0;
	eidx->compar = compar;
	list_init(&eidx->ls);

	eidx->tree = rb_create(compar, NULL, &eidx->allocator);

	return eidx;
}

