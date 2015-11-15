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
  

  std::cerr << "        PI -f my.fasta -s \"alignment-name\" \\ " << std::endl;
  std::cerr << "        -w 200 -i  > pi.out 2> pi.err           " << std::endl;
  std::cerr << std::endl;
  std::cerr << " Required:  " << std::endl;
  //------------------------------- XXXXXXXXXX --------------------------------

  std::cerr << "          -f - <STRING> - A FASTA multiple alignment format." << std::endl;
  std::cerr << "                          Must be indexed by samtools faidx." << std::endl;
  std::cerr << "          -s - <STRING> - Any name for the first column in  " << std::endl;
  std::cerr << "                          the output file.                  " << std::endl;
  std::cerr << std::endl;
  std::cerr << " Optional:  " << std::endl; 
  std::cerr << "          -w - <INT>    - Window size, in base pair [20bp]." << std::endl;
  std::cerr << "          -i - <FLAG>   - Ignore indels (\"-\")    [false]." << std::endl;
  std::cerr << "          -h - <FLAG>   - Print help statement     [false]." << std::endl;
  std::cerr << std::endl;
  std::cerr << " Output:   " << std::endl;
  std::cerr << "           PI outputs a four column text file to STDOUT.   "  << std::endl;
  std::cerr << "           Column 1: Alignment name.                       "  << std::endl;
  std::cerr << "           Column 2: Start of window (one based).          "  << std::endl;
  std::cerr << "           Column 3: End of window (one based).            "  << std::endl;
  std::cerr << "           Column 4. Pi (Nucleotide diversity).            "  << std::endl;
  std::cerr << std::endl;
  std::cerr << " Citation:                                                 "  << std::endl;
  std::cerr << "           Nei, M.; Masatoshi Nei; Wen-Hsiung Li (October 1, 1979)" << std::endl;
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
 Function input  : map<std::string, int> ; haplotypes and count

 Function does   : does the math calculation.

 Function returns: double

*/

double pi(std::map<std::string, int> & hapWin){

  double piSum = 0;

  map< std::string, int > done;

  int bump = 0;


  for(std::map<std::string, int>::iterator hapA = hapWin.begin();
      hapA != hapWin.end(); hapA++){

    bump += 1;

    std::map<std::string, int>::iterator hapB = hapWin.begin();
    std::advance(hapB, bump);

    for( ; hapB != hapWin.end(); hapB++){
      if(hapA->first == hapB->first){
	done[hapA->first] = 1;
	continue;
      }
      if(done.find(hapB->first) != done.end()){
	continue;
      }

      int ndiff = 0;
      int nsame = 0;

      for(uint i = 0; i < globalOpts.window; i++){	
	if(toupper(hapA->first[i]) == toupper(hapB->first[i]) || 
	   toupper(hapA->first[i]) == 'N' || 
	   ( ( hapA->first[i] == '-' || hapB->first[i] == '-' ) && globalOpts.indel)){
	  nsame += 1;
	}
	else{
	  ndiff += 1;
	}
      }

      double f1 = double(hapA->second)/double(globalOpts.nhap);
      double f2 = double(hapB->second)/double(globalOpts.nhap);
      double perBaseDiff = double(ndiff)/double(globalOpts.window);

      piSum += f1*f2*perBaseDiff;

      //      std::cerr << "f1: " << f1 << " f2: " << f2 << " ndiff: " << ndiff << " piSum: " <<  piSum << std::endl; 

    }
    done[hapA->first] = 1;
  }
  return piSum;
}


//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  :

 Function does   : calulates PI in a sliding window

 Function returns: nothing

*/
void slideWindow(std::map<std::string, std::string> & haplotypes)
{

 uint start = 0;
  
  while(( start + globalOpts.window ) <= globalOpts.length){

    std::map<std::string, int> haplotypeWindowCount;

    for(std::map<std::string, std::string>::iterator hap = haplotypes.begin();
	hap != haplotypes.end(); hap++){

      std::string subHap = hap->second.substr(start, globalOpts.window);
      
      if(haplotypeWindowCount.find(subHap) != haplotypeWindowCount.end()){
	haplotypeWindowCount[subHap] += 1;
      }
      else{
	haplotypeWindowCount[subHap]  = 1;
      }
    }
    
    double pv = pi(haplotypeWindowCount);
    
    std::cout << globalOpts.seqid << "\t" 
	      << start + 1 << "\t"
	      << start + globalOpts.window  << "\t"
	      << pv
	      << std::endl;
    start += 1;
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

  slideWindow(haplotypes);
  
  return 0;
}
