# include <stdio.h>
# include <stdlib.h>

int main(int argc,char *argv[])
{
  int i;
  FILE *fp;

  if(argc!=2){
    printf("a.out filename \n");
    exit(1);
  }
  if((fp = fopen(argv[1],"r"))==NULL){
    printf("cannot open the file \n");
    exit(1);
  }

  for(;;){
    i=fgetc(fp);
    if(i==EOF) break;
    if(i!=44){
    putchar(i);
    }
  }
}
