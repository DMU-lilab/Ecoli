#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <time.h>

#define SIZE 4746218

void InsertSort(int a[], int n)
{
    int i, j, tmp;

    for ( i=1; i < n; i++ ) {
        tmp = a[i];
        for ( j=i; j > 0 && a[j-1] > tmp; j-- )
            a[j] = a[j-1];
        a[j] = tmp; 
    }
}

uint8_t *LoadRef( char *file )
{
    uint8_t *fa, *cur, ch;

    FILE *fp = fopen(file, "r");
    if ( !fp ) fprintf(stderr, \
            "[Err::%s::%d] Failed to open %s!\n", __func__, __LINE__, file);
    fa = cur = (uint8_t *)malloc(SIZE * sizeof(uint8_t));
    if ( !fa )
        fprintf(stderr, \
            "[Err::%s::%d] Failed to allocate memory!\n", __func__, __LINE__);
    while ( fgetc(fp) != '\n' ) ; // remove the fasta lines that starts with '>'
    while ( (ch = fgetc(fp)) != 0xff ) {
        if ( ch != '\r' && ch != '\n' ) *cur++ = ch;
    }
    fclose(fp); return fa;
}


int main( int argc, char **argv )
{
	int *pos, *cur, altnum;
	uint8_t *fa, mytype[4] = {65,67,71,84}; // A, C, G, T

	srand((unsigned)time(NULL));

	if ( argc != 3 ) {
		fprintf(stderr, "Usage: ./getsnv <refence.fa> <altnum>\n"); exit(1);
	}
	altnum = atoi(argv[2]);
	pos = cur = (int *)malloc(SIZE * sizeof(int));
	for ( int i=0; i < altnum; ++i ) *cur++ = rand() % SIZE;

	InsertSort(pos, altnum); // sort the randome position
	fa = LoadRef(argv[1]); // load the reference
	printf("Chrom\tPosition\tRef\tAlt\n");

	uint8_t alt;
	for ( int i=0; i < altnum; ++i ) {
		while (1) {
			alt = mytype[rand() % 4];
			if ( alt != fa[pos[i]] ) break;
		}
		printf("Ecoli\t%d\t%c\t%c\n", pos[i]+1, fa[pos[i]], alt);
	}
}
