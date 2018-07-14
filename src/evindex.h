#ifndef EVINDEX_H
#define EVINDEX_H

#include "rbindex.h"
//#include "smc.h"
#include "global.h"

extern unsigned long long t_ev_summary;
extern unsigned long long t_ev_tree;

struct genealogy;
struct evindex {
	struct rbindex *idx;
	int npop_all;
	int *dn;
};

static inline void dn_add(int npop, int *x, int *y)
{
	int j;
	for(j = 0; j < npop; j++)
		x[j] += y[j];
}

static inline void dn_add2(int npop, int *x, int *y, int *z)
{
	int j;
	for(j = 0; j < npop; j++)
		x[j] = y[j] + z[j];
}

static inline void dn_add3(int npop, int *w, int *x, int *y, int *z)
{
	int j;
	for(j = 0; j < npop; j++)
		w[j] = x[j] + y[j] + z[j];
}

static inline void dn_sub(int npop, int *x, int *y)
{
	int j;
	for(j = 0; j < npop; j++)
		x[j] -= y[j];
}

static inline void dn_set(int npop, int *x, int *y)
{
	int j;
	for(j = 0; j < npop; j++)
		x[j] = y[j];
}

static inline void dn_clear(int npop, int *x)
{
	int j;
	for(j = 0; j < npop; j++)
		x[j] = 0;
}

void
print_event_tree (const struct evindex *evidx, const char *title);

static inline struct event *evindex_cur(seq_traverser it)
{
	return (struct event *)rbindex_cur(it);
}

static inline struct event *evindex_prev(seq_traverser *it)
{
	return (struct event *)rbindex_prev(it);
}

static inline struct event *evindex_next(seq_traverser *it)
{
	return (struct event *)rbindex_next(it);
}

static inline struct event *evindex_s_get(struct evindex *evidx)
{
	return (struct event *)GET_OBJ(evidx->idx->cur_s);
}

static inline void evindex_s_forward(struct evindex *evidx)
{
	rbindex_forward(evidx->idx);
}

static inline void evindex_s_set(struct evindex *evidx, struct event *e)
{
	rbindex_s_set(evidx->idx, e);
}

void evindex_s_seek(struct evindex *evidx, double t);

static inline void evindex_s_rewind(struct evindex *evidx)
{
	evidx->idx->cur_s = evidx->idx->ls.front->next;
//	evidx->cur_s = evidx->ls.front;
}

static inline int evindex_s_final(struct evindex *evidx)
{
	evidx->idx->cur_s = evidx->idx->rsentinel;
}

static inline int evindex_s_end(struct evindex *evidx)
{
	return evidx->idx->cur_s == NULL;
}

static inline void evindex_seq_on(struct evindex *evidx)
{
	rbindex_seq_on(evidx->idx);
}

void evindex_seq_off(struct evindex *evidx);

struct event *evindex_query(struct evindex *evidx, double t, int *n);
void evindex_propagate_add(int height, struct rb_node **stack, int npop, int *dn);
void evindex_propagate_sub(int height, struct rb_node **stack, int npop, int *dn);

static inline struct event *evindex_find(seq_traverser *it, struct evindex *evidx, struct event *key)
{
	return rbindex_find(it, evidx->idx, key);
}

static inline void __evindex_s_delete(struct evindex *evidx, struct event *e)
{
	__rbindex_s_delete(evidx->idx, (void *)e);
}

static inline void evindex_s_delete(struct evindex *evidx, struct event *e)
{
	rbindex_s_delete(evidx->idx, (void *)e);
}

void evindex_rb_delete(struct evindex *evidx, const void *item);
//static inline void evindex_rb_delete(struct evindex *evidx, struct event *e)
//{
//	rbindex_rb_delete(evidx->idx, (void *)e);
//}

static inline void evindex_delete(struct evindex *evidx, struct event *ev)
{
	if(rbindex_isseq(evidx->idx)){
		rbindex_s_delete(evidx->idx, (void *)ev);

	}else{
//		if(ev->type == EVENT_JOIN || ev->type == EVENT_SPLT){
//		}else{
			evindex_rb_delete(evidx, ev);
//			list_remove(&eidx->idx->ls, ev);
//		}
	}
}

void evindex_rb_insert(struct evindex *evidx, struct event *ev);
//static inline void evindex_rb_insert(struct evindex *evidx, struct event *e)
//{
//	rbindex_rb_insert(evidx->idx, (void *)e);
//}

static inline void evindex_s_insert(struct evindex *evidx, struct event *e)
{
	rbindex_s_insert(evidx->idx, (void *)e);
}

static inline void evindex_insert(struct evindex *evidx, struct event *ev)
{
	if(rbindex_isseq(evidx->idx)){
		evindex_s_insert(evidx, ev);

	}else{
		evindex_rb_insert(evidx, ev);
	}
}

void evindex_reset(struct genealogy *G, struct evindex *evidx);
void evindex_destroy(struct genealogy *G, struct evindex *evidx);
struct evindex *evindex_create(struct genealogy *G, struct config *cfg);

#endif
