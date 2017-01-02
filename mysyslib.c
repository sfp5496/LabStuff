#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//Some personalized versions of common library routines. Mostly just adds error checking automatically.


FILE *datfile;

void doOpen(const char *fname, const char *mode) {
	fprintf(stderr, "Opening %s with mode %s\n", fname, mode);
	datfile = fopen(fname, mode);
	if (!datfile) {
		perror("fopen");
		exit(2);
	}
}

void doClose() {
	if (fclose(datfile))
		perror("fclose");
}

int safeRead(void *buf, size_t size, size_t nr) {
	size_t num = fread(buf, size, nr, datfile);
	if (num < nr) {
		if (feof(datfile)) {
			fprintf(stderr, "Reached EOF\n");
			return 1;
		}
		else {
			perror("fread");
			exit(-1);
		}
	}
	return 0;
}

int doRead(void *buf, size_t size, size_t nr) {
	size_t num = fread(buf, size, nr, datfile);
	if (num < nr) {
		if (feof(datfile)) {
			fprintf(stderr, "Reached EOF, exitting...\n");
			exit(0);
		}
		else {
			perror("fread");
			exit(-1);
		}
	}
	return 0;
}

void doWrite(const void *buf, size_t size, size_t nr) {
	if (!fwrite(buf, size, nr, datfile)) {
		perror("fwrite");
		exit(-1);
	}
}


__attribute__ ((malloc)) void * xmalloc(size_t size) {
	void *p;
	p=malloc(size);
	if (!p) {
		perror("xmalloc");
		exit(-1);
	}
	return p;
}

void xmemcpy(void *out, void *in, size_t size) {
        if (!memcpy(out, in, size)) {
		perror("xmemcpy");
		exit(1);
        }
}

