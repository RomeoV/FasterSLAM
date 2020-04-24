#include "read_input_file.h"
#include <fstream>
#include <iostream>
#include <sstream>

using std::string;
using std::vector;

/*****************************************************************************
 * Should have highest priority to be translated to C. Just make
 * the matrices a 2D array with fixed size for now.
 * **************************************************************************/
void read_input_file(const string s, double **lm, double **wp, size_t& N_lm, size_t& N_wp) 
{
    using std::ifstream;
    using std::istringstream;
    
    if(access(s.c_str(),R_OK) == -1) {
        std::cerr << "Unable to read input file " << s << std::endl;
        exit(EXIT_FAILURE);
    }
    ifstream in(s.c_str());
    
    int lineno = 0;
    int num_lm_elements = 0;
    int num_lm_coordinates = 0;
    int num_wp_elements = 0;
    int num_wp_coordinates = 0;

    while(in) {
        lineno++;
        string str;
        getline(in,str);
        std::istringstream line(str);
        
        vector<string> tokens;
        std::copy(std::istream_iterator<string>(line), 
             std::istream_iterator<string>(), 
             std::back_inserter<vector<string> > (tokens));
        
        if(tokens.size() ==0) {
            continue;
        }    
        else if (tokens[0][0] =='#') {
            continue;
        }    
        else if (tokens[0] == "lm") {
            if(tokens.size() != 3) {
                std::cerr<<"Wrong args for lm!"<<std::endl;
                std::cerr<<"Error occuredon line"<<lineno<<std::endl;
                std::cerr<<"line:"<<str<<std::endl;
                exit(EXIT_FAILURE);
            }   
			/* in the sample file we have e.g. lm 35 2 */ 
			/* 2 is the number of coordinates i.e. x and y */
			/* 35 is the number of landmarks */  
			num_lm_elements = strtof(tokens[1].c_str(),NULL);  
            num_lm_coordinates = strtof(tokens[2].c_str(),NULL);    
            
            N_lm = num_lm_elements; 
            (*lm) = (double*)malloc(num_lm_coordinates*num_lm_elements*sizeof(double));
            for (int r=0; r<num_lm_elements; r++) {
                lineno++;
                if (!in) {
                    std::cerr<<"EOF after reading" << std::endl;
                    exit(EXIT_FAILURE);
                }    
                getline(in,str);
                std::istringstream line(str);
                vector<string> tokens;
                std::copy(std::istream_iterator<string>(line), 
                     std::istream_iterator<string>(), 
                     std::back_inserter<vector<string> > (tokens));
                if(tokens.size() < num_lm_coordinates) {
                    std::cerr<<"invalid line for lm coordinate!"<<std::endl;
                    std::cerr<<"Error occured on line "<<lineno<<std::endl;
                    std::cerr<<"line: "<<str<<std::endl;
                    exit(EXIT_FAILURE);
                }
                
                for (unsigned c=0; c < num_lm_coordinates; c++) {
                    (*lm)[num_lm_coordinates*r + c] = strtof(tokens[c].c_str(),NULL);                
                }                
            }
        }
        else if (tokens[0] == "wp") {
            if(tokens.size() != 3) {
                std::cerr<<"Wrong args for wp!"<<std::endl;
                std::cerr<<"Error occured on line"<<lineno<<std::endl;
                std::cerr<<"line:"<<str<<std::endl;
                exit(EXIT_FAILURE);
            }
			/* in the sample file we have e.g. wp 17 2 */ 
			/* 2 is the number of coordinates i.e. x and y */
			/* 17 is the number of waypoints */    
			num_wp_elements = strtof(tokens[1].c_str(),NULL);
            num_wp_coordinates = strtof(tokens[2].c_str(),NULL);    

            N_wp = num_wp_elements; // should be 17
            (*wp) = (double*)malloc(num_wp_coordinates*num_wp_elements*sizeof(double));
            for (int r=0; r<num_wp_elements; r++) {
                lineno++;
                if (!in) {
                    std::cerr<<"EOF after reading" << std::endl;
                    exit(EXIT_FAILURE);
                }    
                getline(in,str);
                std::istringstream line(str);
                std::vector<string> tokens;
                std::copy(std::istream_iterator<string>(line), 
                     std::istream_iterator<string>(), 
                     std::back_inserter<std::vector<string> > (tokens));
                if(tokens.size() < num_wp_coordinates) {
                    std::cerr<<"invalid line for wp coordinate!"<<std::endl;
                    std::cerr<<"Error occured on line "<<lineno<<std::endl;
                    std::cerr<<"line: "<<str<<std::endl;
                    exit(EXIT_FAILURE);
                }
                
                for (int c=0; c < num_wp_coordinates; c++) {
                    (*wp)[num_wp_coordinates*r + c] = strtof(tokens[c].c_str(),NULL);                
                }                
            }
        }
        else {
            std::cerr<<"Unkwown command"<<tokens[0] <<std::endl;
            std::cerr<<"Error occured on line"<<lineno<<std::endl;
            std::cerr<<"line: "<<str<<std::endl;
            exit(EXIT_FAILURE);
        }    
    }
}
