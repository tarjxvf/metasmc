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

void rbindex_cache_extend(struct rbindex_cache *nc, int new_size)
{
	int *pid, i;
	char *nd;

	nc->nodes = realloc(nc->nodes, sizeof(struct rb_node *) * new_size);
	for(i = nc->cache_size; i < new_size; i++){
		nc->nodes[i] = nd = malloc(sizeof(struct rb_node) + sizeof(int));
		pid = (int *)(nd + sizeof(struct rb_node));
		*pid = i;

	}
	nc->cache_size = new_size;
}

void *rbindex_alloc(struct libavl_allocator *allocator, size_t size)
{
	struct rbindex_cache *nc;
	struct list *queue;
	int id;

	nc = allocator->allocator_data;
	queue = &nc->free_list;

	if(queue->front == NULL){
		/* Allocate a new node in binary indexed tree. */
		if(nc->maxnodes >= nc->cache_size){
			rbindex_cache_extend(nc, nc->cache_size + 1000);
		}
		id = nc->maxnodes++;

	}else{
		struct list_head *l;
		int *ptr;

		/* Get a existing empty node in binary indexed tree. */
		l = queue->front;
		ptr = (int *)GET_OBJ(l);
		id = *ptr;
		__list_remove(queue, l);
		free(l);
	}

#ifdef DEBUG
	fprintf(stderr, "%s: %d: id=%d, nodes[id]=%x\n", __func__, __LINE__, id, nc->nodes[id]);
#endif

	return nc->nodes[id];
}

void rbindex_free(struct libavl_allocator *allocator, struct rb_node *nd)
{
	struct rbindex_cache *nc;
	struct list_head *l;
	int id, *pidx;

	nc = allocator->allocator_data;
	id = *(int *)((char *)nd + sizeof(struct rb_node));

#ifdef DEBUG
	fprintf(stderr, "%s: %d: id=%d, nd=%x\n", __func__, __LINE__, id, nd);
#endif

	/* Add freed id to the queue. */
	l = malloc(sizeof(struct list_head) + sizeof(int));
	pidx = (int *)GET_OBJ(l);
	*pidx = id;
	__list_append(&nc->free_list, l);
}

//struct libavl_allocator rbindex_allocator = {NULL, rbindex_alloc, rbindex_free};

/*void **rbindex_rb_insert(struct rbindex *eidx, void *obj)
{
	void *next;

	next = rb_isam_find(eidx->tree, obj);
	if(next){
		list_insbefore(GET_LIST(next), obj);
		eidx->ls.n++;

	}else{
		list_append(&eidx->ls, obj);
	}

	return rb_probe(eidx->tree, obj);
}*/

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
  idx->ls.n++;

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
	struct list_head *l;
	struct list *queue;

	nc = &eidx->nc;

	/* Clear free list. */
	queue = &nc->free_list;
	l = queue->front;
	while(l){
		struct list_head *tmp;

		tmp = l->next;
		__list_remove(queue, l);
		free(l);
		l = tmp;
	}

	eidx->tree->rb_count = 0;
	eidx->tree->rb_generation = 0;
}

/*static inline int find_root(int n, int h)
{
	int half, quad, r;

	half = 1 << (h - 1);
	quad = half >> 1;
	r = n - half + 1;
	if(r < quad)
		return quad - 1 + r;
	else
		return half - 1;
}*/

// Calculate indices of complete binary tree from an ordered array
/*void complete_binary_tree(int nnodes, int *tree, int *map)
{
	int i, h, n, half, quad;
	int *beg, *end, *heights;

	beg = malloc(sizeof(int) * nnodes);
	end = malloc(sizeof(int) * nnodes);
	heights = malloc(sizeof(int) * nnodes);

	h = 0;
	while((1 << h) <= nnodes) h++;

	n = 1;
	tree[0] = find_root(nnodes, h);
	map[tree[0]] = 0;
	beg[0] = 0;
	end[0] = nnodes;
	heights[0] = h;
	i = 0;
	while(n < nnodes){
		int j, left, right, r;

		h = heights[i];
		half = 1 << (h - 1);
		quad = half >> 1;
		r = end[i] - beg[i] - half + 1;

		j = tree[i];
		// Calculate index of left child.
		if(n < nnodes){
			heights[n] = h - 1;
			left = beg[i] + find_root(j - beg[i], heights[n]);
			tree[n] = left;
			map[left] = n;
			beg[n] = beg[i];
			end[n] = j;
			n++;
		}
		// Calculate index of left child.
		if(n < nnodes){
			if(r < quad)
				heights[n] = h - 2;
			else
				heights[n] = h - 1;
			right = j + 1 + find_root(end[i] - j - 1, heights[n]);
			tree[n] = right;
			map[right] = n;
			beg[n] = j + 1;
			end[n] = end[i];
			n++;
		}
		i++;
	}

	free(heights);
	free(end);
	free(beg);
}*/

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

void __cbt_addchild(struct cbt_info *tree, int n, int *queue, int *j, int *k)
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

void __complete_binary_tree(struct cbt_info *tree, int start)
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
	memset(tree, 0, sizeof(struct cbt_info) * nnodes);
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

#pragma omp parallel for num_threads(NUM_THREADS)
	for(j = 3; j < 7; j++){
		__complete_binary_tree(tree, j);
	}

#pragma omp parallel for num_threads(NUM_THREADS)
	for(i = 0; i < nnodes; i++)
		map[i] = tree[i].root;

	free(tree);
	free(queue);
}

// Calculate indices of complete binary tree from an ordered array
void complete_binary_tree(int nnodes, int *map)
{
	struct cbt_info *tree;
	int i, j, k, h0, *queue;

	tree = malloc(sizeof(struct cbt_info) * nnodes);
	memset(tree, 0, sizeof(struct cbt_info) * nnodes);
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
		map[i] = tree[i].root;

	free(tree);
	free(queue);
}

// Rebuild tree index from sorted list
void rbindex_rebuild_tree(struct rbindex *eidx, struct rb_node **nodes)
{
	int i, nnodes, h, half, *tree;
//	struct rb_node **nodes;
	struct list_head *l;
	void **objs;

	nnodes = eidx->ls.n;
	rbindex_rb_clear(eidx);
	if(nnodes >= eidx->nc.cache_size)
		rbindex_cache_extend(&eidx->nc, nnodes);
	nodes = eidx->nc.nodes;

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

#pragma omp parallel for num_threads(NUM_THREADS)
	for(i = half - 1; i < nnodes; i++){
		nodes[i]->rb_color = RB_RED;
		nodes[i]->rb_data = objs[tree[i]];
		nodes[i]->rb_link[0] = nodes[i]->rb_link[1] = NULL;
	}

#pragma omp parallel for num_threads(NUM_THREADS)
	for(i = 0; i < half - 1; i++){
		int left, right;

		nodes[i]->rb_color = RB_BLACK;
		nodes[i]->rb_data = objs[tree[i]];

		// Link nodes
		left = i * 2 + 1;
		right = left + 1;
		nodes[i]->rb_link[0] = (left < nnodes)?nodes[left]:NULL;
		nodes[i]->rb_link[1] = (right < nnodes)?nodes[right]:NULL;
	}

//	eidx->tree = rb_create(eidx->compar, NULL, );
	eidx->tree->rb_root = nodes[0];
	nodes[0]->rb_color = RB_BLACK;
	eidx->nc.maxnodes = nnodes;

	free(objs);
	free(tree);
}

// Rebuild tree index from sorted list
/*void rbindex_rebuild_tree(struct rbindex *eidx, struct rb_node **nodes)
{
	int i, j, nnodes, *tree, *map, half, h;
//	struct rb_node **nodes;
	struct list_head *l;

	nnodes = eidx->ls.n;
//	rb_destroy(eidx->tree, NULL);
	rbindex_rb_clear(eidx);
	if(nnodes >= eidx->nc.cache_size)
		rbindex_cache_extend(&eidx->nc, nnodes);
	nodes = eidx->nc.nodes;

	tree = malloc(sizeof(int) * nnodes);
	map = malloc(sizeof(int) * nnodes);
	complete_binary_tree(nnodes, tree, map);
	h = 0;
	while(1 << h <= nnodes) h++;
	half = 1 << (h - 1);

	// Build complete binary tree
//	nodes = malloc(sizeof(struct rb_node *) * nnodes);
	l = eidx->ls.front;
	for(i = 0; i < nnodes; i++){
		j = map[i];
		nodes[j]->rb_data = GET_OBJ(l);
		l = l->next;
	}

	// Link and color nodes
	for(i = nnodes - 1; i >= 0; i--){
		// Color nodes
		if(i < half - 1)
			nodes[i]->rb_color = RB_BLACK;
		else
			nodes[i]->rb_color = RB_RED;

		// Link nodes
		if(i * 2 + 1 >= nnodes)
			nodes[i]->rb_link[0] = NULL;
		else
			nodes[i]->rb_link[0] = nodes[i * 2 + 1];
		if(i * 2 + 2 >= nnodes)
			nodes[i]->rb_link[1] = NULL;
		else
			nodes[i]->rb_link[1] = nodes[i * 2 + 2];
	}

//	eidx->tree = rb_create(eidx->compar, NUL);

	eidx->tree->rb_root = nodes[0];
	eidx->tree->rb_count = nnodes;
	nodes[0]->rb_color = RB_BLACK;
	eidx->nc.maxnodes = nnodes;

	free(tree);
	free(map);
}*/

void rbindex_destroy(struct rbindex *eidx)
{
	struct rbindex_cache *nc;
	int i;

	nc = &eidx->nc;
	rbindex_rb_clear(eidx);
	for(i = 0; i < nc->cache_size; i++)
		free(nc->nodes[i]);
	free(nc->nodes);
//	rb_destroy(eidx->tree, NULL);
	free(eidx->tree);
	free(eidx);
}

struct rbindex *rbindex_create(rb_comparison_func *compar, int cache_size)
{
	struct rbindex *eidx;
	int i, *pid;
	char *nd;

	eidx = malloc(sizeof(struct rbindex));

	eidx->nc.cache_size = cache_size;
	eidx->nc.maxnodes = 0;
	eidx->nc.nnodes = 0;
	eidx->nc.nodes = malloc(sizeof(struct rb_node *) * cache_size);
	memset(eidx->nc.nodes, 0, sizeof(struct rb_node *) * cache_size);
	for(i = 0; i < cache_size; i++){
		eidx->nc.nodes[i] = nd = malloc(sizeof(struct rb_node) + sizeof(int));
		pid = (int *)(nd + sizeof(struct rb_node));
		*pid = i;
	}

	list_init(&eidx->nc.free_list);
	eidx->allocator.allocator_data = &eidx->nc;
	eidx->allocator.libavl_malloc = rbindex_alloc;
	eidx->allocator.libavl_free = rbindex_free;

	eidx->flags = 0;
	eidx->compar = compar;
	list_init(&eidx->ls);

	eidx->tree = rb_create(compar, NULL, &eidx->allocator);

	return eidx;
}

