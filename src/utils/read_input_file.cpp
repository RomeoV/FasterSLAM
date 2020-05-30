#include "read_input_file.h"
#include <fstream>
#include <iostream>
#include <sstream>

using std::string;
using std::vector;

void read_input_file(const string s, double **lm, double **wp, size_t& N_lm, size_t& N_wp) 
{
    int scale = 1;
    read_input_file_and_scale(s, scale, lm, wp, N_lm, N_wp);
}

/*****************************************************************************
 * Should have highest priority to be translated to C. Just make
 * the matrices a 2D array with fixed size for now.
 * **************************************************************************/
void read_input_file_and_scale(const string s, const int scale, double **lm, double **wp, size_t& N_lm, size_t& N_wp) 
{
    using std::ifstream;
    using std::istringstream;
    
    if(access(s.c_str(),R_OK) == -1) {
        std::cerr << "Unable to read input file" << s << std::endl;
        exit(EXIT_FAILURE);
    }
    ifstream in(s.c_str());
    
    int lineno = 0;
    int lm_rows = 0;
    int lm_cols = 0;
    int wp_rows = 0;
    int wp_cols = 0;
    
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
            lm_rows = strtof(tokens[1].c_str(),NULL);    
            lm_cols = strtof(tokens[2].c_str(),NULL);
            
            N_lm = lm_rows*scale;
            (*lm) = (double*)malloc(lm_rows*lm_cols*scale*sizeof(double));
            for (int r = 0; r<lm_rows; r++) {
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
                if(tokens.size() < lm_cols) {
                    std::cerr<<"invalid line for lm coordinate!"<<std::endl;
                    std::cerr<<"Error occured on line "<<lineno<<std::endl;
                    std::cerr<<"line: "<<str<<std::endl;
                    exit(EXIT_FAILURE);
                }
                
                for (unsigned c=0; c < lm_cols; c++) {
                    for (int i=0; i<scale; i++){
                        (*lm)[r*lm_cols*scale + c + i*lm_cols] = strtof(tokens[c].c_str(),NULL);
                    }
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
            wp_rows = strtof(tokens[1].c_str(),NULL);    
            wp_cols = strtof(tokens[2].c_str(),NULL);

            N_wp = wp_rows*scale;
            (*wp) = (double*)malloc(wp_rows*wp_cols*scale*sizeof(double));
            for (int r =0; r<wp_rows; r++) {
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
                if(tokens.size() < wp_cols) {
                    std::cerr<<"invalid line for wp coordinate!"<<std::endl;
                    std::cerr<<"Error occured on line "<<lineno<<std::endl;
                    std::cerr<<"line: "<<str<<std::endl;
                    exit(EXIT_FAILURE);
                }
                
                for (int c=0; c < wp_cols; c++) {
                    for (int i=0; i<scale; i++){
                        (*wp)[r*wp_cols*scale + c + i*wp_cols] = strtof(tokens[c].c_str(),NULL);  //Wrong, waypoints shouldnt scale!       
                    }      
                }                
            }
        }
        else {
            std::cerr << "Unkwown command" << tokens[0] << std::endl;
            std::cerr << "Error occured on line" << lineno << std::endl;
            std::cerr << "line: " << str << std::endl;
            exit(EXIT_FAILURE);
        }    
    }
}

void read_sequential_input_file(const std::string s, double **lm, size_t& N_lm, size_t& N_features) {
    using std::ifstream;
    using std::istringstream;
    
    if(access(s.c_str(),R_OK) == -1) {
        std::cerr << "Unable to read input file" << s << std::endl;
        exit(EXIT_FAILURE);
    }
    ifstream in(s.c_str());
    
    int lineno = 0;
    int lm_rows = 0;
    int lm_cols = 0;
    int wp_rows = 0;
    int wp_cols = 0;
    
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
        if(tokens.size() != 3) {
            std::cerr<<"Wrong args for lm!"<<std::endl;
            std::cerr<<"Error occuredon line"<<lineno<<std::endl;
            std::cerr<<"line:"<<str<<std::endl;
            exit(EXIT_FAILURE);
        }  
        else if (tokens[0] == "lm") {
            if(tokens.size() != 3) {
                std::cerr<<"Wrong args for lm!"<<std::endl;
                std::cerr<<"Error occuredon line"<<lineno<<std::endl;
                std::cerr<<"line:"<<str<<std::endl;
                exit(EXIT_FAILURE);
            }        
            N_features=0;
            lm_rows = strtof(tokens[1].c_str(),NULL);    
            lm_cols = strtof(tokens[2].c_str(),NULL);
            
            N_lm = lm_rows;
            (*lm) = (double*)malloc(lm_rows*lm_cols*sizeof(double));
            for (int r = 0; r<lm_rows; r++) {
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
                if(tokens.size() < lm_cols) {
                    std::cerr<<"invalid line for lm coordinate!"<<std::endl;
                    std::cerr<<"Error occured on line "<<lineno<<std::endl;
                    std::cerr<<"line: "<<str<<std::endl;
                    exit(EXIT_FAILURE);
                }
                
                for (unsigned c=0; c < lm_cols; c++) {
                    if (c == 1) {
                        int f_index = strtof(tokens[c].c_str(), NULL);
                        if (f_index > -1) {
                            if (f_index > N_features) {
                                N_features = f_index;
                            }
                        }
                        
                    }
                        
                    (*lm)[r*lm_cols + c] = strtof(tokens[c].c_str(),NULL);
                }                
            }
        }
         else {
            std::cerr << "Unkwown command" << tokens[0] << std::endl;
            std::cerr << "Error occured on line" << lineno << std::endl;
            std::cerr << "line: " << str << std::endl;
            exit(EXIT_FAILURE);
        }  
    }
}
