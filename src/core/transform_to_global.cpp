#include "transform_to_global.h"
#include "pi_to_pi.h"
#include <math.h>
#include "linalg.h"

void transform_to_global(Matrix2d p, int p_size, Vector3d b) 
{
	//rotation matrix
	//since it is defined as double Matrix2d[4];
	Matrix2d rot = {cos(b[2]), -sin(b[2]), sin(b[2]), cos(b[2])};  
	// cos(b[2]) -sin(b[2])
	// sin(b[2]) cos(b[2])

	Matrix2d p_resized;
	copy(p, p_size, p_resized);
	// multiply the two matrices rot*p_resized
	// rot is 2*2 and p_resized is 2*n (here we assume 2*2)
    for(int i=0; i<2; i++){ //for each row of rot
        for(int n=0; n<2; n++){ //for each column of p_resized
			//p_resized[i][n] = rot[i][0]*p[0][n] + rot[i][1]*p[1][n];
            p_resized[i*2+n] = rot[i*2+0]*p[0+n] + rot[i*2+1]*p[2+n];
        }
    }

	// translate
	int c;
	// c<p_resized.cols()
	for (c=0; c<2; c++) {
		p[0 + c] = p_resized[0*2 + c]+b[0]; 				
		p[1*2 + c] = p_resized[1*2 + c]+b[1]; 				
	}

	double input;
	// we assume p is 2*2 for now
	// if p is a pose and not a point, p.rows() == 3
	// --> transform a list of poses [x;y;phi] so that they are global wrt a base pose
	/*if (p_size == 6){
		for (int k=0; k<p_resized.cols();k++) {
			input = p(2,k) + b(2);
   			p(2,k) = pi_to_pi(input);
		}		
	}*/
}