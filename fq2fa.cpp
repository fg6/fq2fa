#include "readnwritefaq.h"



//int readseqs(void);
int calc(void);
int myread(char* file);

int main(int argc, char *argv[])
{ 

  if (argc == 1) {
   fprintf(stderr, "Usage: %s <reads.fq/fa>\n", argv[0]);
   return 1;
  }	
  if((fp = gzopen(argv[1],"r")) == NULL){ 
    printf("ERROR main:: missing input file  !! \n");
    return 1;
  }

  string myname=myrename(argv[1],"fq2fa",".fasta");
  int err=1;

  // File type  
  int isfq=fasttype(argv[1]);

  if(!isfq)
    cout << " Error! not a fastq file! "<<endl;
  else{
     myfile.open(myname);  
     err=readfastq(argv[1],0,1);
     myfile.close();
 }

  return 0;
}
