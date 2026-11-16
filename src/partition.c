/*
 * partition.c — METIS 5 wrapper for FVCOM domain decomposition.
 *
 * Replaces the fork's src/partition.c, which was ported from the METIS 4
 * interface but never updated for METIS 5. Two interface bugs made any run
 * on this stack die in "TOTAL NUMBER OF ELEMENTS SET UP TO RECEIVE DATA
 * NOT EQUAL TO TOTAL NUMBER OF HALO ELEMENTS" (observed with conda metis
 * 5.x, 16 ranks; the Fortran caller mod_par.F keeps its own METIS-4-era
 * option enum):
 *
 * 1. METIS 5's C API is 0-based, but the Fortran caller passes 1-based
 *    node ids (NVG) and the old code built 1-based CSR pointers
 *    (eptr[i] = i*3+1). METIS never got METIS_OPTION_NUMBERING because
 *    FVCOM sets it at METIS-4 slot 18, while METIS 5.1's enum has it at
 *    index 19 (NIPARTS was inserted at 6, shifting everything after it).
 *    The option therefore silently landed on NSEPS, and every other
 *    FVCOM-set knob was misaligned too (NITER->NIPARTS, UFACTOR->PFACTOR,
 *    ...), while the real slots kept -1 (defaults).
 * 2. With NUMBERING left at its 0 default, METIS returned 0-based element
 *    partitions (0..nparts-1) that FVCOM consumed as 1..nparts: element
 *    set 0 was orphaned and the last rank got an empty domain.
 *
 * Fix entirely on this side of the F(ortran)/C boundary, so every caller
 * (level-1, level-2, dam paths) keeps its 1-based semantics:
 *   - convert the 1-based mesh to 0-based on input,
 *   - drive METIS 5 through the real enum values from metis.h (FVCOM's
 *     PTYPE slot coincidentally lands correctly; everything else gets a
 *     METIS default),
 *   - return 1-based element partitions.
 */
#include <stdlib.h>
#include <metis.h>
#include <stdio.h>

void partition_(int* npartsin, int* nein, int* nnin, int* ncommonin,  int* nv[], int* optionsin[], int* epartout[], int* objvalout)
{
  int i, j, nparts, edgecut, nn, ne, ncommon, objval;
  idx_t *eptr, *eind, *epart, *npart, *options;

  ncommon = *ncommonin;
  nparts = *npartsin;
  nn = *nnin;
  ne = *nein;

  /* local vars */
  eind = (malloc(sizeof(idx_t)*(3*ne)));
  eptr = (malloc(sizeof(idx_t)*(ne+1)));
  epart = (malloc(sizeof(idx_t)*(ne)));
  npart = (malloc(sizeof(idx_t)*(nn)));
  options = (malloc(sizeof(idx_t)*METIS_NOPTIONS));

  /* Fortran NVG is 1-based; METIS 5 is 0-based (NUMBERING=0 below) */
  for (j=3*ne, i=0; i<j; i++) {
      eind[i] = (*nv)[i] - 1;
  }

  for (j=ne+1, i=0; i<j; i++){
      eptr[i] = (idx_t)i*3;
  }

  /* FVCOM's option slots are METIS-4 numbered and misaligned with METIS 5;
   * take PTYPE (aligned by luck, value 1 = kway as FVCOM intends) and
   * default everything else. */
  for (i=0; i<METIS_NOPTIONS; i++) options[i] = -1;
  options[METIS_OPTION_PTYPE]   = (*optionsin)[0];
  options[METIS_OPTION_NUMBERING] = 0;

  METIS_PartMeshDual(&ne, &nn, eptr, eind, NULL, NULL, &ncommon, &nparts, NULL, options, &objval, epart, npart);
  printf("edgecut %d\n",objval);

  /* set return values: back to FVCOM's 1..nparts */
  (*objvalout) = objval;

  for (j=ne, i=0; i<j; i++) {
       (*epartout)[i] = epart[i] + 1;
  }

  /* free up the memory */
  free(epart);
  free(eptr);
  free(eind);
  free(npart);
  free(options);

  return;
}
