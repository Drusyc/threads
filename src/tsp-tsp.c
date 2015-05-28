#include <assert.h>
#include <string.h>
#include <stdint.h>

#include "tsp-types.h"
#include "tsp-genmap.h"
#include "tsp-print.h"
#include "tsp-tsp.h"
#include "tsp-lp.h"
#include "tsp-hkbound.h"

/* dernier minimum trouv� */
int minimum;

/* r�solution du probl�me du voyageur de commerce */
int present (int city, int hops, tsp_path_t path, uint64_t vpres)
{
  (void) hops;
  (void) path;
  return (vpres & (1<<city)) != 0;
}



void tsp (int hops, int len, uint64_t vpres, tsp_path_t path, long long int *cuts, tsp_path_t sol, int *sol_len , pthread_mutex_t *mutex_cut, pthread_mutex_t *mutex_min)
{
    pthread_mutex_lock(mutex_min);
    int local_min = minimum;
    pthread_mutex_unlock(mutex_min);


    if (len + cutprefix[(nb_towns-hops)] >= local_min) {
      pthread_mutex_lock(mutex_cut);
      (*cuts)++ ;
      pthread_mutex_unlock(mutex_cut);
      return;
    }

    /* calcul de l'arbre couvrant comme borne inf�rieure */
    if ((nb_towns - hops) > 6 &&
	lower_bound_using_hk(path, hops, len, vpres) >= local_min) {
      pthread_mutex_lock(mutex_cut);
      (*cuts)++ ;
      pthread_mutex_unlock(mutex_cut);
      return;
    }


    /* un rayon de coupure � 15, pour ne pas lancer la programmation
       lin�aire pour les petits arbres, plus rapide � calculer sans */
    if ((nb_towns - hops) > 22
	&& lower_bound_using_lp(path, hops, len, vpres) >= local_min) {
      pthread_mutex_lock(mutex_cut);
      (*cuts)++ ;
      pthread_mutex_unlock(mutex_cut);
      return;
    }

  
    if (hops == nb_towns) {
	    int me = path [hops - 1];
	    int dist = tsp_distance[me][0]; // retourner en 0
            if ( len + dist < local_min ) {
                    pthread_mutex_lock(mutex_min);
		    minimum = len + dist;
		    *sol_len = len + dist;
		    memcpy(sol, path, nb_towns*sizeof(int));
		    if (!quiet)
		      print_solution (path, len+dist);
                    pthread_mutex_unlock(mutex_min);
	    }
    } else {
        int me = path [hops - 1];        
        for (int i = 0; i < nb_towns; i++) {
	  if (!present (i, hops, path, vpres)) {
                path[hops] = i;
		vpres |= (1<<i);
                int dist = tsp_distance[me][i];
                tsp (hops + 1, len + dist, vpres, path, cuts, sol, sol_len, mutex_cut, mutex_min);
		vpres &= (~(1<<i));
            }
        }
    }
}

