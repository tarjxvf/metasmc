#include <stdio.h>

#include "list.h"

void __list_add(struct list *ls, struct list_head *l)
{
	__list_add__(ls, l);
}

void __list_insbefore(struct list_head *ref, struct list_head *l)
{
	__list_insbefore__(ref, l);
}

void __list_insafter(struct list_head *ref, struct list_head *l)
{
	__list_insafter__(ref, l);
}

void __list_append(struct list *ls, struct list_head *l)
{
	__list_append__(ls, l);
}

struct list_head *__list_pop(struct list *ls)
{
	return __list_pop__(ls);
}

void __list_remove(struct list *ls, struct list_head *l)
{
	__list_remove__(ls, l);
}

struct list_head *__list_prev(struct list_head *l)
{
	return __list_prev__(l);
}

void *list_prev(void *obj)
{
	return GET_OBJ(__list_prev__(GET_LIST(obj)));
}

void list_init(struct list *ls)
{
	ls->front = NULL;
	ls->rear = &ls->front;
//	ls->n = 0;
}

void list_concat(struct list *dst, struct list *src)
{
	*dst->rear = src->front;
	src->front->prev = dst->rear;
	dst->rear = src->rear;
}

/* Add an item in the front of the list. */
void list_add(struct list *ls, void *item)
{
	struct list_head *l;

	l = GET_LIST(item);
	__list_add__(ls, l);
}

/* Insert an item before an item. */
void list_insbefore(struct list_head *ref, void *item)
{
	struct list_head *l;

	l = GET_LIST(item);
	__list_insbefore__(ref, l);
}

/* Insert an item before an item. */
void list_insafter(struct list_head *ref, void *item)
{
	struct list_head *l;

	l = GET_LIST(item);
	__list_insafter__(ref, l);
}

/* Append an item after a  list */
void list_append(struct list *ls, void *item)
{
	struct list_head *l;

	l = GET_LIST(item);
	__list_append__(ls, l);
}

/* Remove an item from the list. */
void list_remove(struct list *ls, void *item)
{
	struct list_head *l;
	l = GET_LIST(item);
	__list_remove__(ls, l);
}

void list_print(struct list_head **head)
{
	struct list_head **ptr;
	ptr = head;

	printf("head=%x", ptr);
	while(*ptr){
		printf("->[ptr=%x, *ptr=%x, prev=%x, next=%x, &next=%x]", ptr, *ptr, (*ptr)->prev, (*ptr)->next, &(*ptr)->next);
		ptr = &(*ptr)->next;
	}
	printf("\n");
}

