#include "readnwritefaq.h"

int main(int argc, char *argv[])
{ 

  if (argc == 1) {
   fprintf(stderr, "Usage: %s <reads.fq/fa>\n", argv[0]);
   return 1;
  }	

  std::ifstream file(argv[1]);  
  if ( ! fexists(argv[1]) ){
    printf("ERROR main:: missing input file  !! \n");
    return 1;
  } 
  int pbformat=0;
  if (argc == 3) {
    pbformat=atoi(argv[2]);
  }

  string myname=myrename(argv[1],".fasta");
  int err=1;

  // File type  
  int isfq=fasttype(argv[1]);

  if(!isfq)
    cout << " Error! not a fastq file! "<<endl;
  else{
     myfile.open(myname);  
     err=readfastq(argv[1],0,1,pbformat);
     myfile.close();
 }

  return 0;
}
