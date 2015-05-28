#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <time.h>
#include <assert.h>
#include <complex.h>
#include <stdbool.h>
#include <unistd.h>
#include <pthread.h>
#include <semaphore.h>

#include "tsp-types.h"
#include "tsp-job.h"
#include "tsp-genmap.h"
#include "tsp-print.h"
#include "tsp-tsp.h"
#include "tsp-lp.h"
#include "tsp-hkbound.h"

/* macro de mesure de temps, retourne une valeur en nanosecondes */
#define TIME_DIFF(t1, t2) \
        ((t2.tv_sec - t1.tv_sec) * 1000000000ll + (long long int) (t2.tv_nsec - t1.tv_nsec))


/* tableau des distances */
tsp_distance_matrix_t tsp_distance ={};

/** Paramètres **/

/* nombre de villes */
int nb_towns=10;
/* graine */
long int myseed= 0;
/* nombre de threads */
int nb_threads=1;

/* affichage SVG */
bool affiche_sol= false;
bool affiche_progress=false;
bool quiet=false;


/* mutex */
static pthread_mutex_t mutex_jobs;      /* mutex pour exclusion mutuelle pour la liste des jobs */
static pthread_mutex_t mutex_cuts;      /* mutex pour exclusion mutuelle pour le nombre de coupes */
static pthread_mutex_t mutex_min;       /* mutex pour exclusion mutuelle pour la solution minimum */ /* tsp-tsp.h: extern int minimum */
static pthread_mutex_t mutex_printf;    /* mutex pour exclusion mutuelle pour la sortie  */

static void generate_tsp_jobs (struct tsp_queue *q, int hops, int len, uint64_t vpres, tsp_path_t path, long long int *cuts, tsp_path_t sol, int *sol_len, int depth)
{
        if (len >= minimum) {
                (*cuts)++ ;
                return;
        }

        if (hops == depth) {
                /* On enregistre du travail à faire plus tard... */
                add_job (q, path, hops, len, vpres);
        } else {
                int me = path [hops - 1];        
                for (int i = 0; i < nb_towns; i++) {
                        if (!present (i, hops, path, vpres)) {
                                path[hops] = i;
                                vpres |= (1<<i);
                                int dist = tsp_distance[me][i];
                                generate_tsp_jobs (q, hops + 1, len + dist, vpres, path, cuts, sol, sol_len, depth);
                                vpres &= (~(1<<i));
                        }
                }
        }
}

static void usage(const char *name) {
        fprintf (stderr, "Usage: %s [-s] <ncities> <seed> <nthreads>\n", name);
        exit (-1);
}

struct pthread_data {
        struct tsp_queue * q;
        //uint64_t * vpres;
        long long int * cuts;
        tsp_path_t * sol;
        int * sol_len;
        sem_t *sem;
};


/* Fonction executée par les threads 
 *
 * Accès concurrent : 
 *      Variables : - tsp_queue (liste des jobs) 
 *                  - cuts (coupes)
 *                  - minimum (cout minimum de la solution)
 *
 * */

void * machin(void * data) {
    
        struct pthread_data * pdata = (struct pthread_data *) data;

        struct tsp_queue * q = pdata->q;
        
        tsp_path_t solution;
        memset(solution, -1, MAX_TOWNS * sizeof (int));
        solution[0] = 0;

        long long int * cuts = pdata->cuts;
        tsp_path_t * sol = pdata->sol;
        int * sol_len = pdata->sol_len;

        uint64_t * vpres = malloc(sizeof(uint64_t));


        pthread_mutex_lock(&mutex_jobs);
        bool is_q_empty = empty_queue(q);
        pthread_mutex_unlock(&mutex_jobs);
        while (!is_q_empty) {
                int hops = 0, len = 0;
                *vpres  = 1;

                if (!get_job (q, solution, &hops, &len, vpres, &mutex_jobs, &mutex_printf)){
                        break;
                }

                int lower_hk = lower_bound_using_hk(solution, hops, len, *vpres);
                int lower_lp = lower_bound_using_lp(solution, hops, len, *vpres);

                pthread_mutex_lock(&mutex_min);
                int local_min = minimum;
                pthread_mutex_unlock(&mutex_min);

                // le noeud est moins bon que la solution courante
                /*
                if (minimum < INT_MAX  && (nb_towns - hops) > 10 && ( (lower_bound_using_hk(*solution, hops, len, *vpres)) >= minimum
                                        || (lower_bound_using_lp(*solution, hops, len, *vpres)) >= minimum)) return NULL;
                */

                // le noeud est moins bon que la solution courante
                if (local_min < INT_MAX
                    && (nb_towns - hops) > 10
                    && ( lower_hk >= local_min || lower_lp >= local_min)) {
                        continue;
                }
                
                tsp (hops, len, *vpres, solution, cuts, *sol, sol_len, &mutex_cuts, &mutex_min, &mutex_printf);

                pthread_mutex_lock(&mutex_jobs);
                is_q_empty = empty_queue(q);
                pthread_mutex_unlock(&mutex_jobs);
        } 

        free(vpres);
        return NULL;
}

int main (int argc, char **argv)
{
    unsigned long long perf;
    tsp_path_t path;
    uint64_t vpres=0;
    tsp_path_t sol;
    int sol_len;
    long long int cuts = 0;
    struct tsp_queue q;
    struct timespec t1, t2;

    /* lire les arguments */
    int opt;
    while ((opt = getopt(argc, argv, "spq")) != -1) {
      switch (opt) {
      case 's':
	affiche_sol = true;
	break;
      case 'p':
	affiche_progress = true;
	break;
      case 'q':
	quiet = true;
	break;
      default:
	usage(argv[0]);
	break;
      }
    }

    if (optind != argc-3)
      usage(argv[0]);

    nb_towns = atoi(argv[optind]);
    myseed = atol(argv[optind+1]);
    nb_threads = atoi(argv[optind+2]);
    assert(nb_towns > 0);
    assert(nb_threads > 0);
   
    minimum = INT_MAX;
      
    /* generer la carte et la matrice de distance */
    if (! quiet)
      fprintf (stderr, "ncities = %3d\n", nb_towns);
    genmap ();

    init_queue (&q);

    clock_gettime (CLOCK_REALTIME, &t1);

    memset (path, -1, MAX_TOWNS * sizeof (int));
    path[0] = 0;
    vpres=1;

    /* mettre les travaux dans la file d'attente */
    generate_tsp_jobs (&q, 1, 0, vpres, path, &cuts, sol, & sol_len, 3);
    no_more_jobs (&q);
   

    /** Nos petits ajouts **/

    /* init mutex */
    pthread_mutex_init(&mutex_jobs, NULL);
    pthread_mutex_init(&mutex_cuts, NULL);
    pthread_mutex_init(&mutex_min, NULL);
    pthread_mutex_init(&mutex_printf, NULL);

    pthread_t * table_threads = malloc(nb_threads * sizeof(pthread_t));
    struct pthread_data * table_data = malloc(nb_threads * sizeof(struct pthread_data));

    /* Semaphore de nbthread * /
    sem_t semu;
    sem_init(&semu, 0, nb_threads);

    pdata->q = &q;
    pdata->solution = &solution;
    pdata->cuts = &cuts;
    pdata->sol = &sol;
    pdata->sol_len = &sol_len;
    pdata->sem = &semu;


    / *
    while (!empty_queue (&q)) {
        void * statusdelaliberte;

        pthread_create(&boblethread, NULL, machin, (void *) pdata);
        pthread_join(boblethread, &statusdelaliberte); 
        
    }
    */

    /*
     * Création une seule et unique fois des n threads et il s'auto gère :
     * Il regarde eux-même dans la queue s'il reste du travail.
     * La bouche while (!empty_queue(&q)) est faite en interne par les threads
     */
    for (uint8_t i = 0; i < nb_threads; i++) {
         
       table_data[i].q = &q;
       //table_data[i].solution = &solution;
       table_data[i].cuts = &cuts;
       table_data[i].sol = &sol;
       table_data[i].sol_len = &sol_len;

       pthread_create(&(table_threads[i]), NULL, machin, (void *) (&table_data[i]));
    }//for()

    /* pas bo */
    for (uint8_t i = 0; i < nb_threads; i++) {
        pthread_join(table_threads[i],NULL);
    }

    free(table_data);
    free(table_threads);

    //sem_destroy(&semu);

    /* init mutex */
    pthread_mutex_destroy(&mutex_jobs);
    pthread_mutex_destroy(&mutex_cuts);
    pthread_mutex_destroy(&mutex_min);
    pthread_mutex_destroy(&mutex_printf);

    

    clock_gettime (CLOCK_REALTIME, &t2);

    if (affiche_sol)
      print_solution_svg (sol, sol_len);

    perf = TIME_DIFF (t1,t2);
    printf("<!-- # = %d seed = %ld len = %d threads = %d time = %lld.%03lld ms ( %lld coupures ) -->\n",
	   nb_towns, myseed, sol_len, nb_threads,
	   perf/1000000ll, perf%1000000ll, cuts);

    return 0 ;
}
