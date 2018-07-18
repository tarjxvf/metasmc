#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "evindex.h"
#include "smc.h"

#define MAXNSEC 1000000000

unsigned long long t_ev_summary = 0;
unsigned long long t_ev_tree = 0;

void
__print_event_tree(const struct evindex *evidx, struct rb_node *node, int level)
{
	struct event *ev;
  int i;

  /* You can set the maximum level as high as you like.
     Most of the time, you'll want to debug code using small trees,
     so that a large |level| indicates a ``loop'', which is a bug. */
  if (node == NULL)
    {
      fprintf (stderr, "<nil>");
      return;
    }

  ev = node->rb_data;
  fprintf (stderr, "%x(", ev);

  for (i = 0; i <= 1; i++)
    {
      if (node->rb_link[i])
        {
            __print_event_tree(evidx, node->rb_link[i], level + 1);
        }
      else
        fprintf (stderr, "NULL");

      if (i == 0)
        fputs (", ", stderr);
    }

  fprintf(stderr, ")");
}

/* Prints the entire structure of |tree| with the given |title|. */
void
print_event_tree (const struct evindex *evidx, const char *title)
{
  fprintf (stderr, "%s: ", title);
  __print_event_tree (evidx, evidx->idx->tree->rb_root, 0);
  fprintf(stderr, "\n");
}

/***** LIBAVL allocator of population-wise event index. *****/
int evindex_compar(struct event *a, struct event *b)
{
	double diff;

	diff = a->t - b->t;
	if(diff > 0){
		return 1;

	}else if(diff < 0){
		return -1;

	}else{
		return a->type - b->type;
	}
}

// Rebuild tree index from sorted list
void evindex_rebuild_tree(struct evindex *evidx)
{
	int i, j, nnodes, h, half, *tree, npop_all;
	struct rb_node **nodes;
	struct list_head *l;
	struct event *ev, *left, *right;
	struct rbindex *eidx;
	void **objs;

	eidx = evidx->idx;
	nnodes = eidx->n;
	rbindex_rb_clear(eidx);
	if(nnodes >= eidx->nc->cache_size)
		cache_resize(eidx->nc, nnodes - eidx->nc->cache_size);
	nodes = (struct rb_node **)eidx->nc->objs;

//	objs = malloc(sizeof(void *) * nnodes);
/*	l = eidx->ls.front;
	for(i = 0; i < nnodes; i++){
		objs[i] = GET_OBJ(l);
		l = l->next;
	}*/

	// Build complete binary tree
	tree = malloc(sizeof(int) * nnodes);
	complete_binary_tree(nnodes, tree);
	h = 0;
	while(1 << h <= nnodes) h++;
	half = 1 << (h - 1);

	*(eidx->ls.rear) = NULL;
	l = eidx->ls.front;
	i = 0;
	while(l){
		nodes[tree[i++]]->rb_data = GET_OBJ(l);
		l = l->next;
	}

	npop_all = evidx->npop_all;
	for(j = nnodes - 1; j >= half - 1; j--){
		nodes[j]->rb_color = RB_RED;
		ev = nodes[j]->rb_data;
		nodes[j]->rb_link[0] = nodes[j]->rb_link[1] = NULL;
		dn_set(npop_all, GET_SUMDN(ev), GET_DN(ev));
	}

	for(; j > (nnodes - 1) / 2; j--){
		nodes[j]->rb_color = RB_BLACK;
		ev = nodes[j]->rb_data;
		nodes[j]->rb_link[0] = nodes[j]->rb_link[1] = NULL;
		dn_set(npop_all, GET_SUMDN(ev), GET_DN(ev));
	}

	nodes[j]->rb_color = RB_BLACK;
	ev = nodes[j]->rb_data;
	if(j == nnodes / 2){
		nodes[j]->rb_link[0] = nodes[j]->rb_link[1] = NULL;
		dn_set(npop_all, GET_SUMDN(ev), GET_DN(ev));

	}else{
		nodes[j]->rb_link[0] = nodes[(j << 1) + 1];
		nodes[j]->rb_link[1] = NULL;
		left = nodes[(j << 1) + 1]->rb_data;
		dn_add2(npop_all, GET_SUMDN(ev), GET_DN(ev), GET_SUMDN(left));
	}
	j--;

	for(; j >= 0; j--){
		nodes[j]->rb_color = RB_BLACK;
		ev = nodes[j]->rb_data;
		nodes[j]->rb_link[0] = nodes[(j << 1) + 1];
		nodes[j]->rb_link[1] = nodes[(j << 1) + 2];
		left = nodes[(j << 1) + 1]->rb_data;
		right = nodes[(j << 1) + 2]->rb_data;
		dn_add3(npop_all, GET_SUMDN(ev), GET_DN(ev), GET_SUMDN(left), GET_SUMDN(right));
	}

	eidx->tree->rb_root = nodes[0];
	nodes[0]->rb_color = RB_BLACK;
	eidx->nc->maxnodes = nnodes;

//	free(objs);
	free(tree);
}

void evindex_seq_off(struct evindex *evidx)
{
//	struct timespec beg, end;
//	int nsec;
//	clock_gettime(CLOCK_MONOTONIC, &beg);

//	nnodes = evidx->idx->n;
	rbindex_clearflag(evidx->idx, RBINDEX_SEQUENTIAL);
	evindex_rebuild_tree(evidx);

//	clock_gettime(CLOCK_MONOTONIC, &end);
//	nsec = (end.tv_sec - beg.tv_sec) * MAXNSEC + (end.tv_nsec - beg.tv_nsec);
//	t_ev_tree += nsec;
}

/* Calculate number of lineages at time t. */
struct event *evindex_query(struct evindex *evidx, double t, int *n)
{
	struct event *ev, *evleft;
	struct rb_node *p;
	int npop, j, cmp;

	npop = evidx->npop_all;
	memset(n, 0, sizeof(int) * npop);
	p = evidx->idx->tree->rb_root;
	while(p){
		ev = (struct event *)p->rb_data;
		cmp = t >= ev->t;
		if(cmp){
			for(j = 0; j < npop; j++)
				n[j] += GET_DN(ev)[j];

			if(p->rb_link[0]){
				evleft = (struct event *)p->rb_link[0]->rb_data;
				for(j = 0; j < npop; j++)
					n[j] += GET_SUMDN(evleft)[j];
			}
			p = p->rb_link[1];

		}else{
			p = p->rb_link[0];
		}
	}

	if(cmp){
		seq_traverser l;
		l = GET_LIST(ev);
		return evindex_next(&l);
	}else{
		return ev;
	}
}

void evindex_propagate_add(int height, struct rb_node **stack, int npop, int *dn)
{
	struct event *ev;
	int i;

	// Backtrack and update summary statistics
	for(i = height - 1; i >= 0; i--){
		ev = stack[i]->rb_data;
		dn_add(npop, GET_SUMDN(ev), dn);
	}
}

void evindex_propagate_sub(int height, struct rb_node **stack, int npop, int *dn)
{
	struct event *ev;
	int i;

	// Backtrack and update summary statistics
	for(i = height - 1; i >= 0; i--){
		ev = stack[i]->rb_data;
		dn_sub(npop, GET_SUMDN(ev), dn);
	}
}

/* Deletes from |tree| and returns an item matching |item|.
   Returns a null pointer if no matching item found. */
void evindex_rb_delete(struct evindex *evidx, const void *item)
{
  struct rb_node *pa[RB_MAX_HEIGHT]; /* Nodes on stack. */
  unsigned char da[RB_MAX_HEIGHT];   /* Directions moved from stack nodes. */
  struct rb_table *tree;
  int k;                             /* Stack height. */

  struct rb_node *p;    /* The node to delete, or a node part way to it. */
  int cmp;              /* Result of comparison between |item| and |p|. */
  struct event *ev;

  tree = evidx->idx->tree;
  assert (tree != NULL && item != NULL);

  k = 0;
  p = (struct rb_node *) &tree->rb_root;
  for (cmp = -1; cmp != 0;
       cmp = tree->rb_compare (item, p->rb_data, tree->rb_param))
    {
      int dir = cmp > 0;

      pa[k] = p;
      da[k++] = dir;

      p = p->rb_link[dir];
      if (p == NULL)
        return;
    }
  item = p->rb_data;
  ev = (struct event *)item;
  if(GET_LIST(ev) == evidx->idx->cur_s)
	  rbindex_forward(evidx->idx);
  list_remove(&evidx->idx->ls, ev);
  evidx->idx->n--;
  evindex_propagate_sub(k - 1, pa + 1, evidx->npop_all, GET_DN(ev));

  if (p->rb_link[1] == NULL)
    pa[k - 1]->rb_link[da[k - 1]] = p->rb_link[0];
  else
    {
      enum rb_color t;
      struct rb_node *r = p->rb_link[1];

      if (r->rb_link[0] == NULL)
        {
          r->rb_link[0] = p->rb_link[0];
	  if(p->rb_link[0])
		  dn_add(evidx->npop_all, GET_SUMDN(r->rb_data), GET_SUMDN(p->rb_link[0]->rb_data));
          t = r->rb_color;
          r->rb_color = p->rb_color;
          p->rb_color = t;
          pa[k - 1]->rb_link[da[k - 1]] = r;
          da[k] = 1;
          pa[k++] = r;
        }
      else
        {
          struct rb_node *s;
          int j = k++;

//	  pa[j] = p;
          for (;;)
            {
              da[k] = 0;
              pa[k++] = r;
              s = r->rb_link[0];
              if (s->rb_link[0] == NULL)
                break;

              r = s;
            }

	  ev = s->rb_data;
	  evindex_propagate_sub(k - j - 1, pa + j + 1, evidx->npop_all, GET_DN(ev));
	  dn_set(evidx->npop_all, GET_SUMDN(ev), GET_DN(ev));

          da[j] = 1;
          pa[j] = s;
          pa[j - 1]->rb_link[da[j - 1]] = s;

          s->rb_link[0] = p->rb_link[0];

	  if(p->rb_link[0])
		  dn_add(evidx->npop_all, GET_SUMDN(ev), GET_SUMDN(p->rb_link[0]->rb_data));

          r->rb_link[0] = s->rb_link[1];
          s->rb_link[1] = p->rb_link[1];

	  if(p->rb_link[1])
		  dn_add(evidx->npop_all, GET_SUMDN(ev), GET_SUMDN(p->rb_link[1]->rb_data));

          t = s->rb_color;
          s->rb_color = p->rb_color;
          p->rb_color = t;
        }
    }

  // Fix double black.
  if (p->rb_color == RB_BLACK)
    {
      for (;;)
        {
          struct rb_node *x = pa[k - 1]->rb_link[da[k - 1]];
          if (x != NULL && x->rb_color == RB_RED)
            {
              x->rb_color = RB_BLACK;
              break;
            }
          if (k < 2)
            break;

          if (da[k - 1] == 0)
            {
              struct rb_node *w = pa[k - 1]->rb_link[1];

              if (w->rb_color == RB_RED)
                {
                  w->rb_color = RB_BLACK;
                  pa[k - 1]->rb_color = RB_RED;

		  dn_sub(evidx->npop_all, GET_SUMDN(pa[k - 1]->rb_data), GET_SUMDN(w->rb_data));
		  if(w->rb_link[0])
			  dn_sub(evidx->npop_all, GET_SUMDN(w->rb_data), GET_SUMDN(w->rb_link[0]->rb_data));

                  pa[k - 1]->rb_link[1] = w->rb_link[0];

		  if(w->rb_link[0])
			  dn_add(evidx->npop_all, GET_SUMDN(pa[k - 1]->rb_data), GET_SUMDN(w->rb_link[0]->rb_data));

                  w->rb_link[0] = pa[k - 1];

		  dn_add(evidx->npop_all, GET_SUMDN(w->rb_data), GET_SUMDN(pa[k - 1]->rb_data));

                  pa[k - 2]->rb_link[da[k - 2]] = w;

                  pa[k] = pa[k - 1];
                  da[k] = 0;
                  pa[k - 1] = w;
                  k++;

                  w = pa[k - 1]->rb_link[1];
                }

              if ((w->rb_link[0] == NULL
                   || w->rb_link[0]->rb_color == RB_BLACK)
                  && (w->rb_link[1] == NULL
                      || w->rb_link[1]->rb_color == RB_BLACK))
                w->rb_color = RB_RED;
              else
                {
                  if (w->rb_link[1] == NULL
                      || w->rb_link[1]->rb_color == RB_BLACK)
                    {
                      struct rb_node *y = w->rb_link[0];
                      y->rb_color = RB_BLACK;
                      w->rb_color = RB_RED;

		      dn_sub(evidx->npop_all, GET_SUMDN(w->rb_data), GET_SUMDN(y->rb_data));
		      if(y->rb_link[1])
			      dn_sub(evidx->npop_all, GET_SUMDN(y->rb_data), GET_SUMDN(y->rb_link[1]->rb_data));

                      w->rb_link[0] = y->rb_link[1];

		      if(y->rb_link[1])
			      dn_add(evidx->npop_all, GET_SUMDN(w->rb_data), GET_SUMDN(y->rb_link[1]->rb_data));

                      y->rb_link[1] = w;

		      dn_add(evidx->npop_all, GET_SUMDN(y->rb_data), GET_SUMDN(w->rb_data));

                      w = pa[k - 1]->rb_link[1] = y;
                    }

                  w->rb_color = pa[k - 1]->rb_color;
                  pa[k - 1]->rb_color = RB_BLACK;
                  w->rb_link[1]->rb_color = RB_BLACK;

		  dn_sub(evidx->npop_all, GET_SUMDN(pa[k - 1]->rb_data), GET_SUMDN(w->rb_data));
		  if(w->rb_link[0])
			  dn_sub(evidx->npop_all, GET_SUMDN(w->rb_data), GET_SUMDN(w->rb_link[0]->rb_data));

                  pa[k - 1]->rb_link[1] = w->rb_link[0];

		  if(w->rb_link[0])
			  dn_add(evidx->npop_all, GET_SUMDN(pa[k - 1]->rb_data), GET_SUMDN(w->rb_link[0]->rb_data));

                  w->rb_link[0] = pa[k - 1];

		  dn_add(evidx->npop_all, GET_SUMDN(w->rb_data), GET_SUMDN(pa[k - 1]->rb_data));

                  pa[k - 2]->rb_link[da[k - 2]] = w;
                  break;
                }
            }
          else
            {
              struct rb_node *w = pa[k - 1]->rb_link[0];

              if (w->rb_color == RB_RED)
                {
                  w->rb_color = RB_BLACK;
                  pa[k - 1]->rb_color = RB_RED;

		  dn_sub(evidx->npop_all, GET_SUMDN(pa[k - 1]->rb_data), GET_SUMDN(w->rb_data));
		  if(w->rb_link[1])
			  dn_sub(evidx->npop_all, GET_SUMDN(w->rb_data), GET_SUMDN(w->rb_link[1]->rb_data));

                  pa[k - 1]->rb_link[0] = w->rb_link[1];

		  if(w->rb_link[1])
			  dn_add(evidx->npop_all, GET_SUMDN(pa[k - 1]->rb_data), GET_SUMDN(w->rb_link[1]->rb_data));

                  w->rb_link[1] = pa[k - 1];

		  dn_add(evidx->npop_all, GET_SUMDN(w->rb_data), GET_SUMDN(pa[k - 1]->rb_data));

                  pa[k - 2]->rb_link[da[k - 2]] = w;

                  pa[k] = pa[k - 1];
                  da[k] = 1;
                  pa[k - 1] = w;
                  k++;

                  w = pa[k - 1]->rb_link[0];
                }

              if ((w->rb_link[0] == NULL
                   || w->rb_link[0]->rb_color == RB_BLACK)
                  && (w->rb_link[1] == NULL
                      || w->rb_link[1]->rb_color == RB_BLACK))
                w->rb_color = RB_RED;
              else
                {
                  if (w->rb_link[0] == NULL
                      || w->rb_link[0]->rb_color == RB_BLACK)
                    {
                      struct rb_node *y = w->rb_link[1];
                      y->rb_color = RB_BLACK;
                      w->rb_color = RB_RED;

		      dn_sub(evidx->npop_all, GET_SUMDN(w->rb_data), GET_SUMDN(y->rb_data));
		      if(y->rb_link[0])
			      dn_sub(evidx->npop_all, GET_SUMDN(y->rb_data), GET_SUMDN(y->rb_link[0]->rb_data));

                      w->rb_link[1] = y->rb_link[0];

		      if(y->rb_link[0])
			      dn_add(evidx->npop_all, GET_SUMDN(w->rb_data), GET_SUMDN(y->rb_link[0]->rb_data));

                      y->rb_link[0] = w;

		      dn_add(evidx->npop_all, GET_SUMDN(y->rb_data), GET_SUMDN(w->rb_data));

                      w = pa[k - 1]->rb_link[0] = y;
                    }

                  w->rb_color = pa[k - 1]->rb_color;
                  pa[k - 1]->rb_color = RB_BLACK;
                  w->rb_link[0]->rb_color = RB_BLACK;

		  dn_sub(evidx->npop_all, GET_SUMDN(pa[k - 1]->rb_data), GET_SUMDN(w->rb_data));
		  if(w->rb_link[1])
			  dn_sub(evidx->npop_all, GET_SUMDN(w->rb_data), GET_SUMDN(w->rb_link[1]->rb_data));

                  pa[k - 1]->rb_link[0] = w->rb_link[1];

		  if(w->rb_link[1])
			  dn_add(evidx->npop_all, GET_SUMDN(pa[k - 1]->rb_data), GET_SUMDN(w->rb_link[1]->rb_data));

                  w->rb_link[1] = pa[k - 1];

		  dn_add(evidx->npop_all, GET_SUMDN(w->rb_data), GET_SUMDN(pa[k - 1]->rb_data));

                  pa[k - 2]->rb_link[da[k - 2]] = w;
                  break;
                }
            }

          k--;
        }

    }

  tree->rb_alloc->libavl_free (tree->rb_alloc, p);
  tree->rb_count--;
  tree->rb_generation++;
}

/* Insert a new event and update node summary statistics.
 * This function is based on rb_probe of GNU libavl. */
void evindex_rb_insert(struct evindex *evidx, struct event *ev)
{
  struct rb_node *pa[RB_MAX_HEIGHT]; /* Nodes on stack. */
  unsigned char da[RB_MAX_HEIGHT];   /* Directions moved from stack nodes. */
  struct rb_table *tree;
  int k;                             /* Stack height. */

  struct rb_node *p; /* Traverses tree looking for insertion point. */
  struct rb_node *n; /* Newly inserted node. */

  tree = evidx->idx->tree;
  assert (tree != NULL && ev != NULL);

  pa[0] = (struct rb_node *) &tree->rb_root;
  da[0] = 0;
  k = 1;
  for (p = tree->rb_root; p != NULL; p = p->rb_link[da[k - 1]])
    {
      int cmp = tree->rb_compare (ev, p->rb_data, tree->rb_param);
      if (cmp == 0)
        return;

      pa[k] = p;
      da[k++] = cmp > 0;
    }

  n = pa[k - 1]->rb_link[da[k - 1]] =
    tree->rb_alloc->libavl_malloc (tree->rb_alloc, sizeof *n);
  if (n == NULL)
    return;

  n->rb_data = ev;
  n->rb_link[0] = n->rb_link[1] = NULL;
  n->rb_color = RB_RED;
  tree->rb_count++;
  tree->rb_generation++;

  /* Determine next element. */
  if(da[k - 1] == 0){
	  list_insbefore(GET_LIST(pa[k - 1]->rb_data), ev);
  }else{
	  list_insafter(GET_LIST(pa[k - 1]->rb_data), ev);
  }
  evidx->idx->n++;

  // Propagating change of summary statistics toward root
  dn_set(evidx->npop_all, GET_SUMDN(ev), GET_DN(ev));
  evindex_propagate_add(k - 1, pa + 1, evidx->npop_all, GET_DN(ev));

  while (k >= 3 && pa[k - 1]->rb_color == RB_RED)
    {
      if (da[k - 2] == 0)
        {
          struct rb_node *y = pa[k - 2]->rb_link[1];
          if (y != NULL && y->rb_color == RB_RED)
            {// recoloring
              pa[k - 1]->rb_color = y->rb_color = RB_BLACK;
              pa[k - 2]->rb_color = RB_RED;
              k -= 2;
            }
          else
            {// restructuring
              struct rb_node *x;

              if (da[k - 1] == 0)
                y = pa[k - 1];
              else
                {
                  x = pa[k - 1];
                  y = x->rb_link[1];
		  dn_sub(evidx->npop_all, GET_SUMDN(x->rb_data), GET_SUMDN(y->rb_data));
		  if(y->rb_link[0])
			  dn_sub(evidx->npop_all, GET_SUMDN(y->rb_data), GET_SUMDN(y->rb_link[0]->rb_data));

                  x->rb_link[1] = y->rb_link[0];

		  if(y->rb_link[0])
			  dn_add(evidx->npop_all, GET_SUMDN(x->rb_data), GET_SUMDN(y->rb_link[0]->rb_data));

                  y->rb_link[0] = x;

		  dn_add(evidx->npop_all, GET_SUMDN(y->rb_data), GET_SUMDN(x->rb_data));

                  pa[k - 2]->rb_link[0] = y;
                }

              x = pa[k - 2];
              x->rb_color = RB_RED;
              y->rb_color = RB_BLACK;

	      dn_sub(evidx->npop_all, GET_SUMDN(x->rb_data), GET_SUMDN(y->rb_data));
	      if(y->rb_link[1])
		      dn_sub(evidx->npop_all, GET_SUMDN(y->rb_data), GET_SUMDN(y->rb_link[1]->rb_data));

              x->rb_link[0] = y->rb_link[1];

	      if(y->rb_link[1])
		      dn_add(evidx->npop_all, GET_SUMDN(x->rb_data), GET_SUMDN(y->rb_link[1]->rb_data));

              y->rb_link[1] = x;

	      dn_add(evidx->npop_all, GET_SUMDN(y->rb_data), GET_SUMDN(x->rb_data));

              pa[k - 3]->rb_link[da[k - 3]] = y;
              break;
            }
        }
      else
        {
          struct rb_node *y = pa[k - 2]->rb_link[0];
          if (y != NULL && y->rb_color == RB_RED)
            {// recoloring
              pa[k - 1]->rb_color = y->rb_color = RB_BLACK;
              pa[k - 2]->rb_color = RB_RED;
              k -= 2;
            }
          else
            {// restructuring
              struct rb_node *x;

              if (da[k - 1] == 1)
                y = pa[k - 1];
              else
                {
                  x = pa[k - 1];
                  y = x->rb_link[0];

		  dn_sub(evidx->npop_all, GET_SUMDN(x->rb_data), GET_SUMDN(y->rb_data));
		  if(y->rb_link[1])
			  dn_sub(evidx->npop_all, GET_SUMDN(y->rb_data), GET_SUMDN(y->rb_link[1]->rb_data));

                  x->rb_link[0] = y->rb_link[1];

		  if(y->rb_link[1])
			  dn_add(evidx->npop_all, GET_SUMDN(x->rb_data), GET_SUMDN(y->rb_link[1]->rb_data));

                  y->rb_link[1] = x;

		  dn_add(evidx->npop_all, GET_SUMDN(y->rb_data), GET_SUMDN(x->rb_data));

                  pa[k - 2]->rb_link[1] = y;
                }

              x = pa[k - 2];
              x->rb_color = RB_RED;
              y->rb_color = RB_BLACK;

	      dn_sub(evidx->npop_all, GET_SUMDN(x->rb_data), GET_SUMDN(y->rb_data));
	      if(y->rb_link[0])
		      dn_sub(evidx->npop_all, GET_SUMDN(y->rb_data), GET_SUMDN(y->rb_link[0]->rb_data));

              x->rb_link[1] = y->rb_link[0];

	      if(y->rb_link[0])
		      dn_add(evidx->npop_all, GET_SUMDN(x->rb_data), GET_SUMDN(y->rb_link[0]->rb_data));

              y->rb_link[0] = x;

	      dn_add(evidx->npop_all, GET_SUMDN(y->rb_data), GET_SUMDN(x->rb_data));

              pa[k - 3]->rb_link[da[k - 3]] = y;
              break;
            }
        }
    }
  tree->rb_root->rb_color = RB_BLACK;
}

void evindex_reset(struct genealogy *G, struct evindex *evidx)
{
	list_init(&evidx->idx->ls);
	evidx->idx->n = 0;
	list_append(&evidx->idx->ls, GET_OBJ(evidx->idx->lsentinel));
	evidx->idx->n++;
	list_append(&evidx->idx->ls, GET_OBJ(evidx->idx->rsentinel));
	evidx->idx->n++;
	rbindex_rb_clear(evidx->idx);
}

void __evindex_destroy(struct genealogy *G, struct evindex *evidx)
{
	struct event *ev;

	ev = (struct event *)GET_OBJ(evidx->idx->lsentinel);
	free(GET_LIST(ev));
	ev = (struct event *)GET_OBJ(evidx->idx->rsentinel);
	free(GET_LIST(ev));

	rbindex_destroy(evidx->idx);
	free(evidx->dn);
}

void evindex_destroy(struct genealogy *G, struct evindex *evidx)
{
	__evindex_destroy(G, evidx);
	free(evidx);
}

// Sequential seek for interval endpoint.
void evindex_s_seek(struct evindex *evidx, double t)
{
	if(rbindex_isseq(evidx->idx)){
		seq_traverser lfwd, lbwd;
		struct event *ebwd, *efwd;

		lbwd = evidx->idx->cur_s;
		ebwd = evindex_cur(lbwd);
		while(ebwd->t >= t) ebwd = evindex_prev(&lbwd);
		lfwd = lbwd;
		efwd = evindex_cur(lfwd);
		while(efwd->t < t) efwd = evindex_next(&lfwd);
		evidx->idx->cur_s = lfwd;
	}
}

void evindex_init(struct genealogy *G, struct config *cfg, struct evindex *evidx)
{
	struct event *ev;
	int npop;

	evidx->idx = rbindex_create(evindex_compar, cfg->maxfrag * 4);
	evidx->npop_all = npop = cfg->npop_all;
	evidx->dn = malloc(sizeof(int) * npop);

	// Set up sentinel
	ev = alloc_event(G->cfg, EVENT_SAMP, 0);
	evidx->idx->lsentinel = GET_LIST(ev);

	list_append(&evidx->idx->ls, ev);
	evidx->idx->n++;

	ev = alloc_event(G->cfg, EVENT_GSIZ, INFINITY);
	((struct gsiz_event *)ev)->size = 1;
	evidx->idx->rsentinel = GET_LIST(ev);

	list_append(&evidx->idx->ls, ev);
	evidx->idx->n++;
}

struct evindex *evindex_create(struct genealogy *G, struct config *cfg)
{
	struct evindex *evidx;

	evidx = malloc(sizeof(struct evindex));
	evindex_init(G, cfg, evidx);

	return evidx;
}

