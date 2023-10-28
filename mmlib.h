/* basic matrix library for  comp381                            */
/* written by ian a. mason @ une  march 15  '99                 */     
/* modified may 25 to work with subblocks                       */       


#include <stdio.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/types.h>
#include <time.h>
#include <sys/stat.h>
#include <stdlib.h>

/* fd  matrix_size row col blkrow blk int[]                          */
int get_block_row(int, int, int, int, int, int, int[]);
int set_block_row(int, int, int, int, int, int, int[]);

void init_block(int*, int*, int*, int blk, int row, int col);

void block_mult(int*, int*, int*, int blk);