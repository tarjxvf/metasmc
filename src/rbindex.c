#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "rbindex.h"

void *rbindex_alloc(struct libavl_allocator *allocator, size_t size)
{
	return malloc(size);
}

void rbindex_free(struct libavl_allocator *allocator, void *block)
{
	free(block);
}

struct libavl_allocator rbindex_allocator = {rbindex_alloc, rbindex_free};

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

// Calculate indices of complete binary tree from an ordered array
void complete_binary_tree(int nnodes, int *tree, int *map)
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
}

// Rebuild tree index from sorted list
void rbindex_rebuild_tree_slow(struct rbindex *eidx, struct rb_node **nodes)
{
	int i, nnodes, h, half, *tree, *map;
//	struct rb_node **nodes;
	struct list_head *l;
	void **objs;

	nnodes = eidx->ls.n;
	rb_destroy(eidx->tree, NULL);
	objs = malloc(sizeof(void *) * nnodes);
	l = eidx->ls.front;
	for(i = 0; i < nnodes; i++){
		objs[i] = GET_OBJ(l);
		l = l->next;
	}

	// Build complete binary tree
	tree = malloc(sizeof(int) * nnodes);
	map = malloc(sizeof(int) * nnodes);
	complete_binary_tree(nnodes, tree, map);
	h = 0;
	while(1 << h <= nnodes) h++;
	half = 1 << (h - 1);

//	nodes = malloc(sizeof(struct rb_node *) * nnodes);
	for(i = nnodes - 1; i >= 0; i--){
		nodes[i] = malloc(sizeof(struct rb_node));
		nodes[i]->rb_data = objs[tree[i]];

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

	eidx->tree = rb_create(eidx->compar, NULL, &rbindex_allocator);
	eidx->tree->rb_root = nodes[0];
	nodes[0]->rb_color = RB_BLACK;

//	free(nodes);
	free(objs);
	free(tree);
}

// Rebuild tree index from sorted list
void rbindex_rebuild_tree(struct rbindex *eidx, struct rb_node **nodes)
{
	int i, j, nnodes, *tree, *map, half, h;
//	struct rb_node **nodes;
	struct list_head *l;

	nnodes = eidx->ls.n;
	rb_destroy(eidx->tree, NULL);

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
		nodes[j] =  malloc(sizeof(struct rb_node));
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

	eidx->tree = rb_create(eidx->compar, NULL, &rbindex_allocator);
	eidx->tree->rb_root = nodes[0];
	nodes[0]->rb_color = RB_BLACK;

//	free(nodes);
	free(tree);
	free(map);
}

void rbindex_destroy(struct rbindex *eidx)
{
	rb_destroy(eidx->tree, NULL);
	free(eidx);
}

struct rbindex *rbindex_create(rb_comparison_func *compar, void *param, struct libavl_allocator *allocator)
{
	struct rbindex *eidx;

	eidx = malloc(sizeof(struct rbindex));
	eidx->flags = 0;
	eidx->compar = compar;
	eidx->tree = rb_create(compar, param, allocator);
	list_init(&eidx->ls);

	return eidx;
}


