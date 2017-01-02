#ifndef __MYSYSLIB_H__
#define __MYSYSLIB_H__

extern FILE *datfile;
void doOpen(const char *fname, const char *mode);
void doClose();
int doRead(void *buf, size_t size, size_t nr);
int safeRead(void *buf, size_t size, size_t nr);
void doWrite(const void *buf, size_t size, size_t nr);
__attribute__ ((malloc)) void * xmalloc(size_t size);
void xmemcpy(void *out, void *in, size_t size);

#endif
