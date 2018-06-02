#ifndef RBINDEX_H
#define RBINDEX_H

#include "global.h"
#include "rb.h"

#define RBINDEX_SEQUENTIAL	0x1	/* Set this flag to turn on sequential mode. When the index is in sequential mode, the tree index is not updated during insertion and deletion. The tree index need to be rebuilt whenever sequential mode is turned off. */

//typedef struct rb_traverser rb_traverser;
typedef struct list_head * rb_traverser;

struct rbindex {
	int flags;
	rb_comparison_func *compar;
	struct rb_table *tree;
	struct list_head *lsentinel;
	struct list_head *rsentinel;
	struct list ls;
	struct list_head *cur_s;	// Cursor for sequential mode. New objects are inserted before cursor.
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

static inline void rbindex_s_insert(struct rbindex *eidx, void *obj)
{
	list_insbefore(eidx->cur_s, obj);
}

extern struct libavl_allocator rbindex_allocator;

static inline void *rbindex_prev(rb_traverser *it)
{
	*it = __list_prev(*it);
	return GET_OBJ(*it);
}

static inline void *rbindex_cur(rb_traverser it)
{
	return GET_OBJ(it);
}

static inline void *rbindex_next(rb_traverser *it)
{
	*it = (*it)->next;
	return GET_OBJ(*it);
}

static inline void rbindex_insbefore(rb_traverser it, void *obj)
{
	list_insbefore(it, obj);
}

/* Set up cursor. */
static inline void rbindex_s_set(struct rbindex *eidx, void *obj)
{
	eidx->cur_s = GET_LIST(obj);
}

static inline rb_traverser rbindex_s_get(void *obj)
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
static inline void rbindex_seq_off(struct rbindex *eidx)
{
	rbindex_clearflag(eidx, RBINDEX_SEQUENTIAL);
	rbindex_rebuild_tree(eidx);
}

/* Find the first object after the key. */
static inline void *rbindex_find(rb_traverser *it, struct rbindex *eidx, void *key)
{
	void *obj;

	obj = rb_isam_find(eidx->tree, key);
	eidx->cur_s = *it = GET_LIST(obj);
	return obj;
}

void **rbindex_insert(struct rbindex *eidx, void *obj);
void rbindex_delete(struct rbindex *eidx, void *obj);
void rbindex_destroy(struct rbindex *eidx);
struct rbindex *rbindex_create(rb_comparison_func *compar, void *param, struct libavl_allocator *allocator);

#endif
