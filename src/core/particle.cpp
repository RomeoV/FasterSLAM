#include "particle.h"

/*****************************************************************************
 * OPTIMIZATION STATUS
 * Last Worked on: 01.04.2020
 * Done: Base Implementation, unit test, weights as an outside array
 * ToDo: Check if index necessary (increases array size by 8 Bytes atm!),
 * 		 maybe create get/set for weights for safe access (instead of pointer
 * 		 juggling).
 ****************************************************************************/

/*****************************************************************************
 * PERFORMANCE STATUS
 * Work: TBD
 * Memory moved: TBD
 * Cycles: Not measured.
 * Performance: Not measured.
 * Optimal: Not measured.
 * Status: Not started.
 ****************************************************************************/

void initParticle(Particle* p, const size_t Nf) {

	p->xf = (double*) malloc (2* Nf * sizeof (double));
    if (p->xf == NULL) {
        free (p);
		return;
	}
	// Try to allocate Pf, free structure if fail.

	p->Pf = (double*) malloc (4* Nf * sizeof (double));
	if (p->Pf == NULL) {
        free (p);
		return;
	}

	p-> Nfa = 0;
	p->Nf = Nf;
	p->set_xv = set_xv;
	p->set_Pv = set_Pv;
	p->set_xfi = set_xfi;
	p->set_Pfi = set_Pfi;	
}

void initParticle(Particle* p, const size_t Nf, int particle_index) {
	initParticle(p,Nf);
	p->index = particle_index;
}

Particle* newParticle(const size_t Nf) {
	// Try to allocate particle structure.
	Particle* p = (Particle*) malloc(sizeof(Particle));

	if (p == NULL) {
        return NULL;
	}

	initParticle(p,Nf);
	
	return p;
}

Particle* newParticle(const size_t Nf, int particle_index) {
	Particle* p = newParticle(Nf);
	p->index = particle_index;
	return p;
}

void delParticle (Particle* p) {
    // Can safely assume particle is NULL or fully built.

    if (p != NULL) {
        free (p->xf);
		free (p->Pf);
        free (p);
    }
}

void delParticle (Particle p) {
    // Only free additionally allocated memory (used for particle array!)
    free (p.xf);
	free (p.Pf);
}

void set_xv(Particle* p, Vector3d xv) {
	p->xv[0] = xv[0];
	p->xv[1] = xv[1];
	p->xv[2] = xv[2];
}

void set_Pv(Particle* p, Matrix3d Pv) {
	for (int i = 0; i<9; i++) {
		p->Pv[i] = Pv[i];
	}
}

void set_xfi(Particle* p, Vector2d xf, int index) {
	if (index > p->Nfa) {
		p->Nfa = index;
	}
	p->xf[2*index+0] = xf[0];
	p->xf[2*index+1] = xf[1];
}

void set_Pfi(Particle* p, Matrix2d Pf, int index) {
	if (index > p->Nfa) {
		p->Nfa = index;
	}
	for (int i = 0; i<4; i++) {
		p->Pf[4*index + i] = Pf[i];
	}
}
