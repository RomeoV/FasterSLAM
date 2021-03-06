#include "particle.h"
#include "linalg.h"
#include "alignment_utils.h"
#include <cstdlib>

/*****************************************************************************
 * OPTIMIZATION STATUS
 * Done: Base Implementation, unit test, weights as an outside array
 * ToDo: Check if index necessary (increases array size by 8 Bytes atm!),
 * 		 maybe create get/set for weights for safe access (instead of pointer
 * 		 juggling).
 ****************************************************************************/

/*****************************************************************************
 * PERFORMANCE STATUS
 * Work: TBD
 * Memory moved: TBD
 * Cycles: TBD
 * Performance: TBD
 * Optimal: TBD
 * Status: TBD
 ****************************************************************************/

//////////////////////////////////////////////////////////////////////////////
//! Also particle stays unaffected from base, we consider this as an immutable
//! for now.
//////////////////////////////////////////////////////////////////////////////

void initParticle_prealloc(Particle* p, const int Nf, const double* xv_initial) {
	copy(xv_initial, 3, p->xv);
	for (int i = 0; i<9;i++) {
		p->Pv[i] = 0.0;
	}
	p->xf = (double *) aligned_alloc(32, aligned_alloc_size(2*Nf*sizeof(double), 32));
    if (p->xf == NULL) {
        free (p);  // \todo is this correct? What if particle is statically allocated...
		return;
	}
	// Try to allocate Pf, free structure if fail.

	p->Pf = (double *) aligned_alloc(32, aligned_alloc_size(4*Nf*sizeof(double), 32));
	if (p->Pf == NULL) {
		free(p->xf);
        free (p);
		return;
	}

	p->Nfa = 0;
	p->Nf = Nf;
}

void initParticle(Particle* p, const int Nf, const double* xv_initial) {
	p->xv = (double*) aligned_alloc (32, aligned_alloc_size(3*sizeof(double), 32));
	p->Pv = (double*) aligned_alloc (32, aligned_alloc_size(9*sizeof(double), 32));
	initParticle_prealloc(p, Nf, xv_initial);
}

/* This method assumes that all particle have allocated the same amount of memory for xf and Pf, namely 2*Nf and 4*Nf respectively. */
void copyParticle(const Particle* p_ref, Particle* p_target) {
    //*(p_target->w) = *(p_ref->w);
    copy(p_ref->xv, 3, p_target->xv);
    copy(p_ref->Pv, 9, p_target->Pv);
    p_target->Nf = p_ref->Nf;
    p_target->Nfa = p_ref->Nfa;

    copy(p_ref->xf, 2*p_ref->Nfa, p_target->xf);
    copy(p_ref->Pf, 4*p_ref->Nfa, p_target->Pf);
}

Particle* newParticle(const int Nf, const double* xv_initial) {
	// Try to allocate particle structure.
	Particle* p = (Particle*) aligned_alloc(32, aligned_alloc_size(sizeof(Particle), 32));
	if (p == NULL) {
        return NULL;
	}

	initParticle(p, Nf, xv_initial);
	
	return p;
}

void delParticleMembers_prealloc(Particle* p) {
    free (p->xf);
    free (p->Pf);
}

void delParticleMembers (Particle* p) {
	free (p->xv);
	free (p->Pv);
    free (p->xf);
    free (p->Pf);
}

void delParticleMembersAndFreePtr (Particle* p) {
    // Can safely assume particle is NULL or fully built.

    if (p != NULL) {
		free (p->xv);
		free (p->Pv);
    }
	delParticleMembersAndFreePtr_prealloc(p);
}	

void delParticleMembersAndFreePtr_prealloc (Particle* p) {
    // Can safely assume particle is NULL or fully built.

    if (p != NULL) {
        free (p->xf);
        free (p->Pf);
        free (p);
    }
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
