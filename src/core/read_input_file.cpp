#include "read_input_file.h"
#include <fstream>
#include <iostream>
#include <sstream>

void read_input_file(const string s, double *lm, double *wp) 
{
	using std::ifstream;
	using std::istringstream;
	
	if(access(s.c_str(),R_OK) == -1) {
		std::cerr << "Unable to read input file" << s << std::endl;
		exit(EXIT_FAILURE);
	}
	ifstream in(s.c_str());
    
	int lineno = 0;
	//int lm_rows = 0;
	//int lm_cols = 0;
	//int wp_rows = 0;
	//int wp_cols = 0;
    
	while(in) {
		lineno++;
		string str;
		getline(in,str);
		istringstream line(str);
        
		vector<string> tokens;
		copy(istream_iterator<string>(line), 
			 istream_iterator<string>(), 
			 back_inserter<vector<string> > (tokens));
		
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
			// our matrix is num_elements (rows) * num_coordinates (columns)
			int num_lm_coordinates = strtof(tokens[1].c_str(),NULL);
			int num_lm_elements = strtof(tokens[2].c_str(),NULL);

			for (int c =0; c<num_lm_elements; c++) {
				lineno++;
				if (!in) {
					std::cerr<<"EOF after reading" << std::endl;
					exit(EXIT_FAILURE);
				}	
				getline(in,str);
				istringstream line(str);
				vector<string> tokens;
				copy(istream_iterator<string>(line), 
                     istream_iterator<string>(), 
                     back_inserter<vector<string> > (tokens));
				if(tokens.size() < num_lm_coordinates) {
					std::cerr<<"invalid line for lm coordinate!"<<std::endl;
					std::cerr<<"Error occured on line "<<lineno<<std::endl;
					std::cerr<<"line: "<<str<<std::endl;
					exit(EXIT_FAILURE);
				}
				
				for (unsigned r=0; r < num_lm_coordinates; r++) {
					lm[r*num_lm_coordinates+c] = strtof(tokens[r].c_str(),NULL);				
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
			int num_wp_coordinates = strtof(tokens[1].c_str(),NULL);	
			int num_wp_elements = strtof(tokens[2].c_str(),NULL);

			for (int c =0; c<num_wp_elements; c++) {
				lineno++;
				if (!in) {
					std::cerr<<"EOF after reading" << std::endl;
					exit(EXIT_FAILURE);
				}	
				getline(in,str);
				istringstream line(str);
				std::vector<string> tokens;
				copy(istream_iterator<string>(line), 
                     istream_iterator<string>(), 
                     back_inserter<std::vector<string> > (tokens));
				if(tokens.size() < num_wp_elements) {
					std::cerr<<"invalid line for wp coordinate!"<<std::endl;
					std::cerr<<"Error occured on line "<<lineno<<std::endl;
					std::cerr<<"line: "<<str<<std::endl;
					exit(EXIT_FAILURE);
				}
                
				for (int r=0; r < num_wp_elements; r++) {
					wp[r*num_wp_coordinates+c] = strtof(tokens[r].c_str(),NULL);				
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
