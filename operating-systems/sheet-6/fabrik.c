//#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <pthread.h>
#include <semaphore.h>
#include <stdbool.h>
#include <fcntl.h>

//lager
int betten = 0;
int schraenke = 0;
sem_t sem_lager_betten_full;
sem_t sem_lager_schraenke_full;
sem_t sem_lager_free;
sem_t sem_lager_access;

//rampe
int rbetten = 0;
int rschraenke = 0;
sem_t sem_rampe_betten_full;
sem_t sem_rampe_schraenke_full;
sem_t sem_rampe_free;
sem_t sem_rampe_access;

//geschäft 1
int g1betten = 0;
sem_t sem_g1_full;
sem_t sem_g1_free;
sem_t sem_g1_access;

//geschäft 2
int g2betten = 0;
int g2schraenke = 0;
sem_t sem_g2b_full;
sem_t sem_g2b_free;
sem_t sem_g2s_full;
sem_t sem_g2s_free;
sem_t sem_g2_access;


//getrennte threads für betten und schränke. damit beim be- und entladen jedes
//lkw nicht beides zeitgleich -> mutex

//lastwagen 2 (lkw1 fährt nur betten aus)
sem_t sem_l2_mutex;

//mitarbeiter für beladung
sem_t sem_mitarbeiter_mutex;


//thread-funktionen für Schränke
void *thread_fabrik_schraenke_func(void *data);
void *thread_mitarbeiter_schraenke_func(void *data);
void *thread_lastwagen1_betten_func(void *data);
void *thread_lastwagen2_schraenke_func(void *data);
void *thread_geschaeft1_func(void *data);
void *thread_geschaeft2_schraenke_func(void *data);

//thread-funktionen für Betten
void *thread_fabrik_betten_func(void *data);
void *thread_mitarbeiter_betten_func(void *data);
void *thread_lastwagen2_betten_func(void *data);
void *thread_geschaeft2_betten_func(void *data);

int main(int argc, char **args) {
   pthread_t thread_fabrik_schraenke;
   pthread_t thread_mitarbeiter_schraenke;
   pthread_t thread_lastwagen1_betten;
   pthread_t thread_lastwagen2_schraenke;

   pthread_t thread_geschaeft2_schraenke;

   pthread_t thread_fabrik_betten;
   pthread_t thread_mitarbeiter_betten;
   pthread_t thread_lastwagen2_betten;
   pthread_t thread_geschaeft1;
   pthread_t thread_geschaeft2_betten;

   printf("Starting up...\n");

   sem_init(&sem_lager_betten_full, 0, 0);
   sem_init(&sem_lager_schraenke_full, 0, 0);
   sem_init(&sem_lager_free, 0, 500);
   sem_init(&sem_lager_access, 0, 1);
   sem_init(&sem_rampe_betten_full, 0, 0);
   sem_init(&sem_rampe_schraenke_full, 0, 0);
   sem_init(&sem_rampe_free, 0, 30);
   sem_init(&sem_rampe_access, 0, 1);
   sem_init(&sem_g1_full, 0, 0);
   sem_init(&sem_g1_free, 0, 30);
   sem_init(&sem_g1_access, 0, 1);
   sem_init(&sem_g2b_full, 0, 0);
   sem_init(&sem_g2b_free, 0, 20);
   sem_init(&sem_g2s_full, 0, 0);
   sem_init(&sem_g2s_free, 0, 20);
   sem_init(&sem_g2_access, 0, 1);
   sem_init(&sem_l2_mutex, 0, 1);
   sem_init(&sem_mitarbeiter_mutex, 0, 1);

   printf("Semaphores created.\n");

   pthread_create(&thread_fabrik_schraenke, NULL, thread_fabrik_schraenke_func, (void*)0);
   pthread_create(&thread_mitarbeiter_schraenke, NULL, thread_mitarbeiter_schraenke_func, (void*)0);
   pthread_create(&thread_lastwagen1_betten, NULL, thread_lastwagen1_betten_func, (void*)0);
   pthread_create(&thread_lastwagen2_schraenke, NULL, thread_lastwagen2_schraenke_func, (void*)0);
   pthread_create(&thread_geschaeft2_schraenke, NULL, thread_geschaeft2_schraenke_func, (void*)0);
   pthread_create(&thread_fabrik_betten, NULL, thread_fabrik_betten_func, (void*)0);
   pthread_create(&thread_mitarbeiter_betten, NULL, thread_mitarbeiter_betten_func, (void*)0);
   pthread_create(&thread_lastwagen2_betten, NULL, thread_lastwagen2_betten_func, (void*)0);
   pthread_create(&thread_geschaeft1, NULL, thread_geschaeft1_func, (void*)0);
   pthread_create(&thread_geschaeft2_betten, NULL, thread_geschaeft2_betten_func, (void*)0);

   printf("Threads created.\n");

   pthread_join(thread_fabrik_schraenke, NULL);
   pthread_join(thread_mitarbeiter_schraenke, NULL);
   pthread_join(thread_lastwagen1_betten, NULL);
   pthread_join(thread_lastwagen2_schraenke, NULL);
   pthread_join(thread_geschaeft2_schraenke, NULL);
   pthread_join(thread_fabrik_betten, NULL);
   pthread_join(thread_mitarbeiter_betten, NULL);
   pthread_join(thread_lastwagen2_betten, NULL);
   pthread_join(thread_geschaeft1, NULL);
   pthread_join(thread_geschaeft2_betten, NULL);

   return EXIT_SUCCESS;
}


//semaphor-schema:

//1. reserviere plätze am zielort (dekrementiere ziel-free)
//2. gib belegte plätze an quelle frei (dekrementiere quell-full)
//3. locke quelle, entnimm, unlocke quelle
//4. gib neue plätze an quelle frei (inkrementiere quell-free)
//4. locke ziel, entlade, unlocke ziel
//5. belege plätze an zielort (inkrementiere ziel-full)



//produziere pro zyklus ein bett und einen schrank für das lager
void *thread_fabrik_schraenke_func(void *data) {
   while(1) {
     sem_wait(&sem_lager_free);

     sem_wait(&sem_lager_access);
     sem_post(&sem_lager_schraenke_full);
     sem_post(&sem_lager_access);
     printf("Schrank ins Lager\n");
   }
}

//transportiere pro zyklus ein bett vom lager zur rampe
void *thread_fabrik_betten_func(void *data) {
   while(1) {
     sem_wait(&sem_lager_free);

     sem_wait(&sem_lager_access);
     sem_post(&sem_lager_betten_full);
     sem_post(&sem_lager_access);
     printf("Bett ins Lager\n");
   }
}

//transportiere pro zyklus einen schrank vom lager zur rampe
void *thread_mitarbeiter_schraenke_func(void *data) {
   while(1) {
     sem_wait(&sem_rampe_free);
     sem_wait(&sem_lager_schraenke_full);

     sem_wait(&sem_mitarbeiter_mutex);

     sem_wait(&sem_lager_access);
     sem_post(&sem_lager_free);
     sem_post(&sem_lager_access);
     printf("Schrank aus dem Lager\n");

     sem_wait(&sem_rampe_access);
     sem_post(&sem_rampe_schraenke_full);
     sem_post(&sem_rampe_access);
     printf("Schrank auf die Rampe\n");

     sem_post(&sem_mitarbeiter_mutex);
   }
}

//transportiere pro zyklus ein bett vom lager zur rampe
void *thread_mitarbeiter_betten_func(void *data) {
   //endlosschleife
   while(1) {
     sem_wait(&sem_rampe_free);
     sem_wait(&sem_lager_betten_full);

     sem_wait(&sem_mitarbeiter_mutex);

     sem_wait(&sem_lager_access);
     sem_post(&sem_lager_free);
     sem_post(&sem_lager_access);
     printf("Bett aus dem Lager\n");

     sem_wait(&sem_rampe_access);
     sem_post(&sem_rampe_betten_full);
     sem_post(&sem_rampe_access);
     printf("Bett auf die Rampe\n");

     sem_post(&sem_mitarbeiter_mutex);
   }
}

//fahre ein bett pro zyklus von der rampe zu geschäft 1
void *thread_lastwagen1_betten_func(void *data) {
   while(1) {
     sem_wait(&sem_g1_free);
     sem_wait(&sem_rampe_betten_full);

     sem_wait(&sem_rampe_access);
     sem_post(&sem_rampe_free);
     sem_post(&sem_rampe_access);
     printf("Bett in L1\n");

     sem_wait(&sem_g1_access);
     sem_post(&sem_g1_full);
     sem_post(&sem_g1_access);
     printf("Bett in G1\n");
   }
}


//fahre einen schrank pro zyklus von der rampe zu geschäft 2
void *thread_lastwagen2_schraenke_func(void *data) {
   while(1) {
     sem_wait(&sem_g2s_free);
     sem_wait(&sem_rampe_schraenke_full);

     sem_wait(&sem_l2_mutex);

     sem_wait(&sem_rampe_access);
     sem_post(&sem_rampe_free);
     sem_post(&sem_rampe_access);
     printf("Schrank in L2\n");

     sem_wait(&sem_g2_access);
     sem_post(&sem_g2s_full);
     sem_post(&sem_g2_access);
     printf("Schrank in G2\n");

     sem_post(&sem_l2_mutex);
   }
}

//fahre ein bett pro zyklus von der rampe zu geschäft 2
void *thread_lastwagen2_betten_func(void *data) {
   while(1) {
     sem_wait(&sem_g2b_free);
     sem_wait(&sem_rampe_betten_full);

     sem_wait(&sem_l2_mutex);

     sem_wait(&sem_rampe_access);
     sem_post(&sem_rampe_free);
     sem_post(&sem_rampe_access);
     printf("Bett in L2\n");

     sem_wait(&sem_g2_access);
     sem_post(&sem_g2b_full);
     sem_post(&sem_g2_access);
     printf("Bett in G2\n");

     sem_post(&sem_l2_mutex);
   }
}

//verkaufe pro zyklus ein bett
void *thread_geschaeft1_func(void *data) {
   while(1) {
     sem_wait(&sem_g1_full);

     sem_wait(&sem_g1_access);
     sem_post(&sem_g1_free);
     sem_post(&sem_g1_access);
     printf("Schrankverkauf in G1\n");
   }
}

//für geschäft 2 ist sequentialität nicht vorgeschrieben, betten und schränke dürfen zeitgleich verkauft
//werden. daher ist hier, im gegensatz zu mitarbeiter und lkw2, kein mutex nötig.

//verkaufe pro zyklus einen schrank
void *thread_geschaeft2_schraenke_func(void *data) {
   while(1) {
     sem_wait(&sem_g2s_full);

     sem_wait(&sem_g2_access);
     sem_post(&sem_g2s_free);
     sem_post(&sem_g2_access);
     printf("Schrankverkauf in G2\n");
   }
}

//verkaufe pro zyklus ein bett
void *thread_geschaeft2_betten_func(void *data) {
   while(1){
     sem_wait(&sem_g2b_full);

     sem_wait(&sem_g2_access);
     sem_post(&sem_g2b_free);
     sem_post(&sem_g2_access);
     printf("Bettverkauf in G2\n");
   }
}
