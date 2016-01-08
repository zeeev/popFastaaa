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
  bool                     indel;
  uint                    window;
  uint                      nhap;
  uint                    length;
  std::string              seqid;
  std::string               file;
}globalOpts;

static const char *optString = "hiw:f:s:";

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
  

  std::cerr << "        SNP -f my.fasta -s \"alignment-name\" \\ " << std::endl;
  std::cerr << std::endl;
  std::cerr << " Required:  " << std::endl;
  //------------------------------- XXXXXXXXXX --------------------------------

  std::cerr << "          -f - <STRING> - A FASTA multiple alignment format." << std::endl;
  std::cerr << "                          Must be indexed by samtools faidx." << std::endl;
  std::cerr << "          -s - <STRING> - Any name for the first column in  " << std::endl;
  std::cerr << "                          the output file.                  " << std::endl;
  std::cerr << std::endl;
  std::cerr << " Optional:  " << std::endl; 
  std::cerr << "          -h - <FLAG>   - Print help statement     [false]." << std::endl;
  std::cerr << std::endl;
  std::cerr << " Output:   " << std::endl;
  std::cerr << "           SNP outputs a tab delimited text file with 2+N columns.   "  << std::endl;
  std::cerr << "           Where N is the number of sequences.                       " << std::endl;
  std::cerr << "           Column 1: Alignment name.                       "  << std::endl;
  std::cerr << "           Column 2: Start of window (one based).          "  << std::endl;
  std::cerr << "           Column 3-N: base at variant site.               "  << std::endl;
  std::cerr << std::endl;
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
	globalOpts.seqid = optarg;
	break;
      }
    case 'i':
      {
	globalOpts.indel = true;
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

//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : haplotypes, position, output vector, i - index

 Function does   : loads vector of bases at postion

 Function returns: returns false if there is  - or N bases

*/

bool loadBases(std::map<std::string , std::string> & haplotypes,
	       std::vector<std::string>            & basesPerHap,
	       std::map<std::string, int>          & ubases,
	       int i){
  
  for(std::map<std::string , std::string>::iterator it = haplotypes.begin();
      it != haplotypes.end(); it++){
    
    std::string baseChar = it->second.substr(i, 1);
    
    if(baseChar == "-"){
      return true;
    }
    if(baseChar == "N"){
      return true;
    }
    ubases[baseChar] = 1;

    basesPerHap.push_back(baseChar);
  }
  return false;
}

//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : haplotypes

 Function does   : prints the header

 Function returns: nothing

*/


void printHeader(std::map<std::string, std::string> & haplotypes){

  std::cout << "SEQID\tPOS";

  for(std::map<std::string , std::string>::iterator it = haplotypes.begin();
      it != haplotypes.end(); it++){
  
    std::cout << "\t" << it->first;
    
  }
  std::cout << endl;
}


//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : takes the haplotypes

 Function does   : prints the SNP positions

 Function returns: nothing

*/
void SNP(std::map<std::string, std::string> & haplotypes)
{
    
  for(unsigned int i = 0; i < globalOpts.length; i++){
    
    std::vector<std::string>   basesPerHap;
    std::map<std::string, int> uBasesCount;

    // there is a char we dont want to consider 
    if(loadBases(haplotypes, basesPerHap, uBasesCount, i)){
      continue;
    }
  // is there is no snp
  
    if(uBasesCount.size() == 1){
      continue;
    }
    std::cout << globalOpts.seqid << "\t";
    std::cout << (i+1);
    for(std::vector<std::string>::iterator iz = basesPerHap.begin(); iz != basesPerHap.end(); iz++ ){
      std::cout << "\t" << *iz ;
    }
    std::cout << endl;
  }
}

//-------------------------------    MAIN     --------------------------------
/*
 Comments:
*/

int main( int argc, char** argv)
{
  globalOpts.window = 20  ;
  parseOpts(argc, argv);
  
  globalOpts.indel  = true;

  if(globalOpts.seqid.empty()){
    std::cerr << "FATAL: must specify an alignment name: -s." << std::endl;
    printHelp();
    exit(1);
  }
  if(globalOpts.file.empty()){
    std::cerr << "FATAL: problem loading fasta: -f." << std::endl;
    printHelp();
    exit(1);
  }
  
  
  std::map<std::string, FastaReference * > inds;
  
  //sequenc name then haplotype.
  
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

  printHeader(haplotypes);
  SNP(haplotypes);
  
  return 0;
}
