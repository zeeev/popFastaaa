/*

This program was created at:  Fri Nov 13 20:05:10 2015
This program was created by:  zev


Contact: zev.kronenberg@gmail.com

Organization: Unviersity of Utah
    School of Medicine
    Salt Lake City, Utah


The MIT License (MIT)

Copyright (c) <2015> <zev>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.


*/

#include <string>
#include <iostream>
#include <math.h>
#include <cmath>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <unistd.h>
#include <map>
#include "Fasta.h"
#include "split.h"

struct options{
  uint                    window;
  uint                      nhap;
  uint                    length;
  uint                    number;
  uint                     start;
  std::string            pattern;
  std::string               file;
}globalOpts;

static const char *optString = "hw:f:n:p:s:";

//------------------------------- XXXXXXXXXX --------------------------------
void printVersion(void){
  cerr << "Version: " << VERSION << endl;
  cerr << "Contact: zev.kronenberg [at] gmail.com " << endl;
  cerr << "Notes  : -If you find a bug, please open a report on github!" << endl;
  cerr << endl;
}
//------------------------------- XXXXXXXXXX --------------------------------

void printHelp(void){
  //------------------------------- XXXXXXXXXX --------------------------------
  std::cerr << std::endl;
  std::cerr << " Example:                                       " << std::endl;
  std::cerr << "        samtools faidx my.fasta                 " << std::endl;
  

  std::cerr << "        PERM -f my.fasta -s \"alignment-name\" \\ " << std::endl;
  std::cerr << "        -w 200 -i  > pi.out 2> pi.err           " << std::endl;
  std::cerr << std::endl;
  std::cerr << " Required:  " << std::endl;
  //------------------------------- XXXXXXXXXX --------------------------------

  std::cerr << "          -f - <STRING> - A FASTA multiple alignment format.  " << std::endl;
  std::cerr << "                          Must be indexed by samtools faidx.  " << std::endl;
  std::cerr << "          -s - <STRING> - Start Position of window to permute." << std::endl;
  std::cerr << "          -w - <INT>    - Window size.                        " << std::endl;
  std::cerr << "          -p - <STRING> - Pattern.                            " << std::endl;
  std::cerr << "          -n - <INT>    - Number of permutations [100].       " << std::endl;
  std::cerr << std::endl;
  std::cerr << " Optional:  " << std::endl; 
   std::cerr << std::endl;
  printVersion();

}
//-------------------------------   OPTIONS   --------------------------------
int parseOpts(int argc, char** argv)
{
  int opt = 0;
  opt = getopt(argc, argv, optString);
  while(opt != -1){
    switch(opt){
    case 'h':
      {	
	break;
      }
    case 's':
      {
	globalOpts.start = atoi(optarg);
	break;
      }
    case 'p':
      {
	globalOpts.pattern = optarg;
	break;
      }
    case 'w':
      {
	string win = optarg;
	globalOpts.window = atoi( win.c_str() );
	break;
      }
    case 'f':
      {
	globalOpts.file = optarg ;
	break;
      }
    case '?':
      {
	std::cerr << "Fatal: unknown option: " << optarg << std::endl;
	exit(1);
	break;
      }
    }
    opt = getopt( argc, argv, optString ); 
  }
  return 1;
}


int ranint(int max)
{
  if(max <= 0) return 0;
  return rand() % max;
}

void pattern_match(const std::string & name,
		   const std::string obs,
		   uint offset,
		   uint perm,
		 
		   const std::string & haplotype){

  double n  = 0;
  double tp = 0;
  
  for(uint i = offset; i < (offset + globalOpts.window) - globalOpts.pattern.size() ; i++){
    n++;
    if(globalOpts.pattern == haplotype.substr(i, globalOpts.pattern.size())) tp++;
  }
  std::cerr << name << "\t" << obs << "\t" << perm << "\t" << offset << n << "\t" << tp << std::endl;
}


//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : map<std::string, std::string> ; haplotypes 
 Function does   : .
 Function returns: double

*/
void permute_regions(const std::string & name,
		     const std::string & haplotype){

  uint r = 0 ;

  pattern_match(name, "observed", globalOpts.start, 0, haplotype);
  
  for(uint i = 0; i < globalOpts.number; i++){
    r = ranint(haplotype.size() - globalOpts.window);
    pattern_match(name, "perm", r, i, haplotype);
  }

}


//-------------------------------    MAIN     --------------------------------
/*
 Comments:
*/

int main( int argc, char** argv)
{

  srand (time(NULL));
  globalOpts.window = 0  ;
  globalOpts.number = 100;
  parseOpts(argc, argv);
  
  if(globalOpts.window == 0){
    std::cerr << "FATAL: must specify window size: -w." << std::endl;
    printHelp();
    exit(1);
  }
  if(globalOpts.pattern.empty()){
    std::cerr << "FATAL: must specify a pattern: -p." << std::endl;
    printHelp();
    exit(1);
  }
  if(globalOpts.file.empty()){
    std::cerr << "FATAL: problem loading fasta: -f." << std::endl;
    printHelp();
    exit(1);
  }
  
  std::map<std::string, FastaReference * > inds;
  std::map<std::string, std::string> haplotypes;
 
  // load haplotypes
 
  FastaReference  rs;

  rs.open(globalOpts.file);
    
  globalOpts.nhap = rs.index->sequenceNames.size();
  globalOpts.length = 0;

  for(std::vector<std::string>::iterator it = rs.index->sequenceNames.begin();
      it != rs.index->sequenceNames.end(); it++){

    std::cerr << "INFO: processing haplotype: " << *it << std::endl;
    
    haplotypes[*it] = rs.getSequence(*it);
  
    if(globalOpts.length == 0){
      globalOpts.length = haplotypes[*it].size();
    }
    else if(haplotypes[*it].size() !=  globalOpts.length){
      std::cerr << "FATAL: sequence: " << *it << " was not the same length." << std::endl;
      std::cerr << "INFO: Are the sequences aligned?" << std::endl;
      printHelp();
      exit(1);
	
    }
  }

  std::cerr << "INFO: N haplotypes: " << haplotypes.size() << std::endl;
  std::cerr << "INFO: Aligned length: " << globalOpts.length << std::endl;

  for(std::map<std::string, std::string>::iterator it = haplotypes.begin();
      it != haplotypes.end(); it++){
    permute_regions(it->first, it->second);
  }
  
  return 0;
}
