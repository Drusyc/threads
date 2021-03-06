#include <assert.h>
#include <string.h>
#include <stdint.h>

#include "tsp-types.h"
#include "tsp-genmap.h"
#include "tsp-print.h"
#include "tsp-tsp.h"
#include "tsp-lp.h"
#include "tsp-hkbound.h"

volatile int minimum;

int present (int city, int hops, tsp_path_t path, uint64_t vpres)
{
        (void) hops;
        (void) path;
        return (vpres & (1<<city)) != 0;
}



void tsp (int hops, int len, uint64_t vpres, tsp_path_t path, long long int *cuts, tsp_path_t sol, int *sol_len , pthread_mutex_t *mutex_cut, pthread_mutex_t *mutex_min, pthread_mutex_t *mutex_printf)
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

        if ((nb_towns - hops) > 6 &&
                        lower_bound_using_hk(path, hops, len, vpres) >= local_min) {
                pthread_mutex_lock(mutex_cut);
                (*cuts)++ ;
                pthread_mutex_unlock(mutex_cut);
                return;
        }


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
                        local_min = minimum;
                        if ( len + dist < minimum) {
                                minimum = len + dist;
                                pthread_mutex_unlock(mutex_min);
                                *sol_len = len + dist;
                                memcpy(sol, path, nb_towns*sizeof(int));
                                if (!quiet) {
                                        pthread_mutex_lock(mutex_printf);
                                        print_solution (path, len+dist);
                                        pthread_mutex_unlock(mutex_printf);
                                }
                        } else {
                                pthread_mutex_unlock(mutex_min);
                        }
                }
        } else {
                int me = path [hops - 1];        
                for (int i = 0; i < nb_towns; i++) {
                        if (!present (i, hops, path, vpres)) {
                                path[hops] = i;
                                vpres |= (1<<i);
                                int dist = tsp_distance[me][i];
                                // printf de sauvage
                                //printf("hops %i len %i dist  %i gt
                                tsp (hops + 1, len + dist, vpres, path, cuts, sol, sol_len, mutex_cut, mutex_min, mutex_printf);
                                vpres &= (~(1<<i));
                        }
                }
        }
}

