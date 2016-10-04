#include <vector>
#include <zlib.h>
#include <stdio.h>
#include "kseq.h"
#include <iostream>
#include <algorithm>    // sort
#include <numeric> // accumulate
#include <iomanip>  //setprecision
#include <fstream>

static gzFile fp;
static  FILE * outF;

KSEQ_INIT(gzFile, gzread)
int readseqs(void);
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
  
  outF = fopen (myname.c_str(),"w");

  std::cout << argv[1] << std::endl;  
  int mlines=myread(argv[1]); //does not work with seq on multi-lines! but a lot faster!
  //int mlines=1;
  if(mlines)
  	readseqs();
  
  fclose(outF);

  return 0;
}


// ---------------------------------------- //
int readseqs(){
// ---------------------------------------- //

  kseq_t *seq;
  int l=-1;
  seq = kseq_init(fp);


  if(1)std::cout << " ...using kseq to read file" << std::endl;

  while ((l = kseq_read(seq)) >= 0){    
	fprintf(outF,">%s %s\n", seq->name.s, seq->comment.s);  
        fprintf(outF,"%s\n", seq->seq.s);  
   }
  kseq_destroy(seq);
  gzclose(fp);
  
 return 0;
}

// ---------------------------------------- //
int myread(char* file)   //FILE *namef)
// ---------------------------------------- //
{ // won't work with seq on multilines!
  char fq[5]={"@"};
  char fa[5]={">"};
  char plus[5]={"+"};
  int readevery=1;
  int pri=0;
  std::ifstream infile(file); 

  std::string line; 
  getline(infile,line);
  if(line.at(0)==fq[0]) readevery=4;  // fastq input file
  else if(line.at(0)==fa[0]) readevery=2;  // fasta input file
  else std::cout << " Error: cannot determine if input file is fasta or fastq " << std::endl;

  getline(infile,line);
  getline(infile,line);


  if(line.at(0)!=plus[0]){ 
	if(pri)std::cout<< "Sequences on multiple lines..." << std::endl;  // seq on single line 
	return(1);
  }else{
        if(pri)std::cout<< "Sequences on single line " << std::endl;
  }

  int nseq=-1;
  int readnext=0;
  std::string thisname;
  std::string thisseq;
  infile.clear();  // start over
  infile.seekg (0, std::ios::beg);
  while (getline(infile,line))
   {
     nseq++;
     thisname="";
     thisseq="";

     if(nseq==readnext){      
      thisname=line;
      fprintf(outF,">%s\n", (thisname.erase(0,1)).c_str());
     }else if(nseq==readnext+1){
      thisseq=line;
      fprintf(outF,"%s\n", thisseq.c_str());
      readnext+=readevery;
     }else{
       if(pri)if(nseq<10)std::cout << " other: " << line[0] <<std::endl;
     }
   }


 return 0;
}

