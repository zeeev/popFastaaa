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
  int                       nhap;
  std::string              seqid;
  std::vector<std::string> files;
}globalOpts;

static const char *optString = "hiw:f:s:";

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
	globalOpts.files = split(optarg, ",");
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

  for(std::map<std::string, int>::iterator hapA = hapWin.begin();
      hapA != hapWin.end(); hapA++){
    for(std::map<std::string, int>::iterator hapB = hapWin.begin();
	hapB != hapWin.end(); hapB++){
      if(hapA->first == hapB->first){
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
      
      piSum += (double(hapA->second)/globalOpts.nhap) * (double(hapB->second)/globalOpts.nhap) * ndiff;

    }
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

  int start = 0;
  
  while(( start + globalOpts.window ) <= haplotypes.begin()->first.size()){

    std::map<std::string, int> haplotypeWindowCount;

    for(std::map<std::string, std::string>::iterator hap = haplotypes.begin();
	hap != haplotypes.end(); hap++){

      std::string subHap = hap->second.substr(start, globalOpts.window);
      
      if(haplotypeWindowCount.find(subHap) != haplotypeWindowCount.end()){
	haplotypeWindowCount[subHap] += 1;
      }
      else{
	haplotypeWindowCount[subHap]  = 0;
      }
    }
    
    double pv = pi(haplotypeWindowCount);
    
    std::cerr << globalOpts.seqid << "\t" 
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
  parseOpts(argc, argv);
  
  globalOpts.indel  = true;
  globalOpts.window = 20  ;

  
  if(globalOpts.seqid.empty()){
    std::cerr << "FATAL: must specify a seqid: -s." << std::endl;
    exit(1);
  }
  if(globalOpts.files.size() == 0){
    std::cerr << "FATAL: problem loading fastas: -f." << std::endl;
    exit(1);
  }
  
  globalOpts.nhap = globalOpts.files.size();
  
  std::map<std::string, FastaReference * > inds;
  std::map<std::string, std::string> haplotypes;
 
  // load haplotypes
 
  int hSize = -1;

  for(std::vector<std::string>::iterator fas = globalOpts.files.begin();
      fas != globalOpts.files.end(); fas++ ){
    
    FastaReference * rs;
    rs = new FastaReference;
    
    rs->open(*fas);
    
    haplotypes[*fas] = rs->getSequence(globalOpts.seqid);
    
    if(hSize == -1){
      hSize = haplotypes[*fas].size();
    }
    else{
      if(haplotypes[*fas].size() != hSize){
	std::cerr << "FATAL: aligned fasta sequences must be same length: "
		  << *fas << std::endl;
	exit(1);
      }
    }

    if( haplotypes[*fas].size() < globalOpts.window ){
      std::cerr << "FATAL: sequence cannot be smaller than the window size." << std::endl;
      exit(1);
      
    }
    
  }
  slideWindow(haplotypes);
  
  return 0;
}
