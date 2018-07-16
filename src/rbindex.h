#ifndef RBINDEX_H
#define RBINDEX_H

#include "global.h"
#include "rb.h"
#include "cache.h"

#define RBINDEX_SEQUENTIAL	0x1	/* Set this flag to turn on sequential mode. When the index is in sequential mode, the tree index is not updated during insertion and deletion. The tree index need to be rebuilt whenever sequential mode is turned off. */

//typedef struct seq_traverser seq_traverser;
typedef struct list_head * seq_traverser;

struct rbindex {
	int flags;
	int n;	// Number of nodes in the tree

	rb_comparison_func *compar;
	struct rb_table *tree;
	struct list_head *lsentinel;
	struct list_head *rsentinel;
	struct list ls;
	struct list_head *cur_s;	// Cursor for sequential mode. New objects are inserted before cursor.

	struct cache *nc;
//	struct rbindex_cache nc;
	struct libavl_allocator allocator;
};

void rbindex_rebuild_tree(struct rbindex *eidx);

static inline void rbindex_setflag(struct rbindex *eidx, int flag)
{
	eidx->flags |= flag;
}

static inline void rbindex_clearflag(struct rbindex *eidx, int flag)
{
	eidx->flags &= ~flag;
}

static inline int rbindex_isseq(struct rbindex *eidx)
{
	return eidx->flags & RBINDEX_SEQUENTIAL;
}

//static inline void **rbindex_s_insert(struct rbindex *eidx, void *obj)
static inline void rbindex_s_insert(struct rbindex *eidx, void *obj)
{
	__list_insbefore(eidx->cur_s, GET_LIST(obj));
	eidx->n++;
}

extern struct libavl_allocator rbindex_allocator;

static inline void *rbindex_prev(seq_traverser *it)
{
	*it = __list_prev(*it);
	return GET_OBJ(*it);
}

static inline void *rbindex_cur(seq_traverser it)
{
	return GET_OBJ(it);
}

static inline void *rbindex_next(seq_traverser *it)
{
	*it = (*it)->next;
	return GET_OBJ(*it);
}

static inline void rbindex_forward(struct rbindex *eidx)
{
	eidx->cur_s = eidx->cur_s->next;
}

/* Set up cursor. */
static inline void rbindex_s_set(struct rbindex *eidx, void *obj)
{
	eidx->cur_s = GET_LIST(obj);
}

static inline seq_traverser rbindex_s_get(void *obj)
{
	return GET_LIST(obj);
}

static inline void rbindex_s_reset(struct rbindex *eidx)
{
	eidx->cur_s = eidx->ls.front->next;
}
/* Turn on sequential mode. */
static inline void rbindex_seq_on(struct rbindex *eidx)
{
	rbindex_setflag(eidx, RBINDEX_SEQUENTIAL);
}

/* Turn off sequential mode. */
/*static inline void rbindex_seq_off(struct rbindex *eidx)
{
	rbindex_clearflag(eidx, RBINDEX_SEQUENTIAL);
	rbindex_rebuild_tree(eidx);
}*/

/* Find the first object after the key. */
static inline void *rbindex_find(seq_traverser *it, struct rbindex *eidx, void *key)
{
	void *obj;

	obj = rb_isam_find(eidx->tree, key);
	eidx->cur_s = *it = GET_LIST(obj);
	return obj;
}

static inline void __rbindex_s_delete(struct rbindex *eidx, void *obj)
{
	__list_remove(&eidx->ls, GET_LIST(obj));
	eidx->n--;
}

static inline void rbindex_s_delete(struct rbindex *eidx, void *obj)
{
	if(GET_LIST(obj) == eidx->cur_s)
		rbindex_next(&eidx->cur_s);
	__list_remove(&eidx->ls, GET_LIST(obj));
	eidx->n--;
}

static inline void rbindex_rb_delete(struct rbindex *eidx, void *obj)
{
	if(rb_delete(eidx->tree, obj))
		list_remove(&eidx->ls, obj);
	eidx->n--;
}

static inline void rbindex_delete(struct rbindex *eidx, void *obj)
{
	if(rbindex_isseq(eidx))
		rbindex_s_delete(eidx, obj);
	else
		rbindex_rb_delete(eidx, obj);
}

void rbindex_rb_insert (struct rbindex *idx, void *item);
//void **rbindex_rb_insert(struct rbindex *eidx, void *obj);
void rbindex_delete(struct rbindex *eidx, void *obj);
void rbindex_rb_clear(struct rbindex *eidx);
void rbindex_destroy(struct rbindex *eidx);
struct rbindex *rbindex_create(rb_comparison_func *compar, int cache_size);

#endif
