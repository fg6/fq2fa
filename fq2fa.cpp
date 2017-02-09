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


  size_t pos = 0;
  std::string myname = argv[1];
  std::string token;
  std::string delimiter = "/";
  while ((pos = myname.find(delimiter)) != std::string::npos) {
    token = myname.substr(0, pos);
    myname.erase(0, pos + delimiter.length());
  }
  myname=myname.substr(0, myname.size()-1) + 'a';
  
  myfile.open(myname);  


  std::cout << argv[1] << std::endl;  
  
  int err=1;

  // File type  
  int isfq=fasttype(argv[1]);

  if(!isfq)
    cout << " Error! not a fastq file! "<<endl;
  else
    err=readfastq(argv[1],0,1);

  myfile.close();
 

  return 0;
}
