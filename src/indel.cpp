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
    uint                    offset;
    std::string              seqid;
    std::string               file;
    std::string                ref;
}globalOpts;

static const char *optString = "hiw:f:s:o:r:";

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


  std::cerr << "        INDEL -f my.fasta -s \"ref-name\" \\ " << std::endl;
  std::cerr << std::endl;
  std::cerr << " Required:  " << std::endl;
  //------------------------------- XXXXXXXXXX --------------------------------

  std::cerr << "          -f - <STRING> - A FASTA multiple alignment format." << std::endl;
  std::cerr << "                          Must be indexed by samtools faidx." << std::endl;
  std::cerr << "          -s - <STRING> - Reference sequence sequence name. " << std::endl;
  std::cerr << "          -r - <STRING> - CHROM for VCF. e.g. chr1.         " << std::endl;
  std::cerr << "          -o - <INT>    - Reference sequence offset.        " << std::endl;

  std::cerr << std::endl;
  std::cerr << " Optional:  " << std::endl;
  std::cerr << "          -h - <FLAG>   - Print help statement     [false]." << std::endl;
  std::cerr << std::endl;
  std::cerr << " Output:   " << std::endl;
  std::cerr << "           VCF4.2.                         "  << std::endl;
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
    case 'o':
        {
            std::string tmp = optarg;
            globalOpts.offset = atoi(tmp.c_str());
            break;
        }
    case 'r':
        {
            globalOpts.ref = optarg;
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
    std::cout << "##fileformat=VCFv4.2\n";
    std::cout << "##INFO=<ID=TYPE,Number=1,Type=String,Description=\"insertion or deletion\">\n";
    std::cout << "##INFO=<ID=END,Number=1,Type=Integer,Description=\"end\">\n";
    std::cout << "##INFO=<ID=SAMPLE,Number=1,Type=Integer,Description=\"sample of origin\">\n";
    std::cout << "CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
}


//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : takes the haplotypes

 Function does   : prints the SNP positions

 Function returns: nothing

*/
void SNP(std::map<std::string, std::string> & haplotypes)
{

    int refPos = 0;

    std::map<std::string, unsigned int > lastIndelBase;

    for(std::map<std::string, std::string>::iterator it = haplotypes.begin();
        it != haplotypes.end(); it++){
        lastIndelBase[it->first] = 0;
    }


  for(unsigned int i = 0; i < globalOpts.length; i++){

      if(haplotypes[globalOpts.seqid][i] != '-' ){
          refPos += 1;
      }
      else{
          refPos += 0;
      }

      for(std::map<std::string, std::string>::iterator it = haplotypes.begin();
          it != haplotypes.end(); it++){
          if(haplotypes[globalOpts.seqid][i] == '-' || haplotypes[it->first][i] == '-'){

              if(globalOpts.seqid == it->first) continue;

              unsigned int j = refPos;

              if( j <= lastIndelBase[it->first] ){
                  continue;
              }

              while(j < globalOpts.length){
                  if(haplotypes[globalOpts.seqid][j] != '-' && haplotypes[it->first][j] != '-') break;
                  j++;
              }

              lastIndelBase[it->first] = j;


              if(i - 1 <= 0){
                  std::cerr << "WARNING: Indel starts before MSA. Skipping." << std::endl;
                  continue;
              }


              std::string s1 = haplotypes[globalOpts.seqid].substr(i, j-i );
              std::string s2 = haplotypes[it->first].substr(i,        j-i );

              if(s1 == s2) continue;

              std::string type = "DEL";

              std::string ref =  haplotypes[globalOpts.seqid].substr(i -1 , j-i +1);
              std::string alt =  haplotypes[globalOpts.seqid].substr(i -1 , 1   );



              if(s1[0] == '-'){
                  type = "INS";
                  ref  = haplotypes[globalOpts.seqid].substr(i -1 , 1   );
                  alt  = haplotypes[it->first].substr(i-1, j-i + 1);
              }

              int end = refPos + (j - i) + globalOpts.offset;

              if(type == "INS"){
                  end = refPos;
              }

              std::cout << globalOpts.ref << "\t" << refPos + globalOpts.offset << "\t" << ".\t"  << ref << "\t" << alt  << "\t" << "TYPE=" << type << ";" << "END=" << end << ";" << "SAMPLE=" << it->first << std::endl;

          }
      }


  }

}

//-------------------------------    MAIN     --------------------------------
/*
 Comments:
*/

int main( int argc, char** argv)
{
  globalOpts.window = 20  ;
  globalOpts.offset = 0   ;


  parseOpts(argc, argv);

  globalOpts.indel  = true;

  if(globalOpts.ref.empty()){
      std::cerr << "FATAL: no ref name provided." << std::endl;
      printHelp();
      exit(1);
  }



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
