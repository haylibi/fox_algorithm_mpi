/* basic matrix library for  comp381                          */
/* written by ian a. mason @ une  easter  '99                 */            
/* rows & columns are numbers 0 through dimension - 1         */
/* modified may 25 to work with subblocks                     */       

#include "mmlib.h"


int get_block_row(int fd, int matrix_size, int row, int col, int blkrow, int blk, int matrix_row[]){
  int blksize = blk * sizeof(int);
  off_t offset =  ((row * matrix_size * blk) + (blkrow * matrix_size) + (col * blk)) * sizeof(int);

  eif(row >= matrix_size){
    fprintf(stderr,"index out of range");
    return -1; }
  else if(lseek(fd, offset, 0) < 0){
    fprintf(stderr,"lseek failed for fd=%d\n",fd);
    return -1; 
  }
  else if(read(fd, matrix_row, blksize) < blksize){
    perror("read failed");
    return -1; };
  return 0;
}

int set_block_row(int fd, int matrix_size, int row, int col, int blkrow,  int blk, int matrix_row[]){
  int blksize = blk * sizeof(int);
  off_t offset =  ((row * matrix_size * blk) + (blkrow * matrix_size) + (col * blk)) * sizeof(int);
  if(offset < 0){
    fprintf(stderr,"offset overflow");
    return -1; }
  else if(row >= matrix_size){
    fprintf(stderr,"indexes out of range"); }
  else if(lseek(fd, offset, 0) < 0){
    perror("lseek failed");
    return -1; }
  else if(write(fd, matrix_row, blksize) < blksize){
      perror("write failed");
      return -1; };
  return 0;

}


void block_mult(int* c, int* a, int* b, int blk){
  int i,j,k;
  
  for (i = 0; i < blk; i++)
    for (j = 0; j < blk; j ++)
      for (k = 0; k < blk; k++)
        c[i*blk+j] += (a[i*blk+k] * b[k*blk+j]);
}