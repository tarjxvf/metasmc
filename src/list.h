#ifndef LIST_H
#define LIST_H

struct list_head {
	struct list_head *next;
	struct list_head **prev;
};

#define GET_LIST(obj)	((struct list_head *)((char *)(obj) - sizeof(struct list_head)))
#define GET_OBJ(l)	((char *)(l) + sizeof(struct list_head))

struct list {
	struct list_head *front;
	struct list_head **rear;
//	int n;
};

void list_init(struct list *ls);
void list_concat(struct list *dst, struct list *src);

/* List operations using object pointer. */
void list_insbefore(struct list_head *ref, void *item);
void list_insafter(struct list_head *ref, void *item);
void list_remove(struct list *ls, void *item);
void list_add(struct list *ls, void *item);	// Add an item at list head
void list_append(struct list *ls, void *item);	// Add an item at list end

/* List operations using list_head pointer. */
static inline void __list_add__(struct list *ls, struct list_head *l) 
{
	if(ls->front)
		ls->front->prev = &l->next;
	else
		ls->rear = &l->next;

	l->next = ls->front;
	ls->front = l;
	l->prev = &ls->front;
//	ls->n++;
}

/* Insert an item before an item. */
static inline void __list_insbefore__(struct list_head *ref, struct list_head *l)
{
	l->next = ref;
	l->prev = ref->prev;
	ref->prev = &l->next;
	*l->prev = l;
}

/* Insert an item before an item. */
static inline void __list_insafter__(struct list_head *ref, struct list_head *l)
{
	l->next = ref->next;
	l->prev = &ref->next;
	ref->next = l;
	l->next->prev = &l->next;
}

/* Append an item after a  list */
static inline void __list_append__(struct list *ls, struct list_head *l)
{
	l->next = NULL;
	*ls->rear = l;
	l->prev = ls->rear;
	ls->rear = &l->next;
//	ls->n++;
}

/* Get and remove last element from NON-EMPTY list. Use with causion. */
static inline struct list_head *__list_pop__(struct list *ls)
{
	struct list_head *l;
	l = (struct list_head *)ls->rear;
	ls->rear = (struct list_head **)l->prev;
	*ls->rear = NULL;
	return l;
}

/* Remove an item from the list. */
static inline void __list_remove__(struct list *ls, struct list_head *l)
{
	*l->prev = l->next;
	if(l->next)
		l->next->prev = l->prev;
	else	// If the removed item is the last one
		ls->rear = l->prev;
//	ls->n--;
//	l->prev = l->next = NULL;
	l->next = NULL;
	l->prev = NULL;
}

static inline struct list_head *__list_prev__(struct list_head *l)
{
	return (struct list_head *)l->prev;
}

void __list_add(struct list *ls, struct list_head *l);
void __list_insbefore(struct list_head *ref, struct list_head *l);
void __list_insafter(struct list_head *ref, struct list_head *l);
void __list_append(struct list *ls, struct list_head *l);
struct list_head *__list_pop(struct list *ls);
void __list_remove(struct list *ls, struct list_head *l);
struct list_head *__list_prev(struct list_head *l);

#endif

