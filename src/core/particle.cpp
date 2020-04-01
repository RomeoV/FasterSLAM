#include "particle.h"

/*****************************************************************************
 * OPTIMIZATION STATUS
 * Last Worked on: 01.04.2020
 * Done: Base Implementation, unit test
 * ToDo: Start optimizing
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

particle* newParticle(const size_t Nf) {
	// Try to allocate particle structure.
	particle* p = (particle*) malloc(sizeof(struct particle));

	if (p == NULL) {
        return NULL;
	}

	// Try to allocate xf, free structure if fail.

    p->xf = (double*) malloc (2* Nf * sizeof (double));
    if (p->xf == NULL) {
        free (p);
        return NULL;
	}
	// Try to allocate Pf, free structure if fail.

	p->Pf = (double*) malloc (4* Nf * sizeof (double));
	if (p->Pf == NULL) {
        free (p);
        return NULL;
	}

	p-> Nfa = 0;
	p->Nf = Nf;
	p->set_xv = set_xv;
	p->set_Pv = set_Pv;
	p->set_xfi = set_xfi;
	p->set_Pfi = set_Pfi;

	
	return p;
}

particle* newParticle(const size_t Nf, int particle_index) {
	particle* p = newParticle(Nf);
	p->particle_index = particle_index;
}

void delParticle (particle* p) {
    // Can safely assume particle is NULL or fully built.

    if (p != NULL) {
        free (p->xf);
		free (p->Pf);
        free (p);
    }
}

void set_xv(particle* p, Vector3d xv) {
	p->xv[0] = xv[0];
	p->xv[1] = xv[1];
	p->xv[2] = xv[2];
}

void set_Pv(particle* p, Matrix3d Pv) {
	for (int i = 0; i<9; i++) {
		p->Pv[i] = Pv[i];
	}
}

void set_xfi(particle* p, Vector2d xf, int index) {
	if (index > p->Nfa) {
		p->Nfa = index;
	}
	p->xf[2*index+0] = xf[0];
	p->xf[2*index+1] = xf[1];
}

void set_Pfi(particle* p, Matrix2d Pf, int index) {
	if (index > p->Nfa) {
		p->Nfa = index;
	}
	for (int i = 0; i<4; i++) {
		p->Pf[4*index + i] = Pf[i];
	}
}