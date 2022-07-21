#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include <sys/types.h>
#include "rebound.h"
#include "tools.h"

#define R_J 0.000477895  //Jupyter radius in AU
#define M_J 0.0009547919 // Jupyter mass in solar masses
#define T_MAX (pow(10,6) * 2.*M_PI)

//Select running mode
#define USE_RANDOM_INITIAL_CONDITIONS      0
#define USE_INITIAL_CONDITIONS_FROM_FILE   1
#define SIMULATE_MOONS                     0

#define SAVE_SIMULATION                    1
#define SAVE_SIMULATION_MOONS              0

#define N_MOONS 1643
#define N_MOONS_PARAM 1643 * 3
#define N_SIM_PARAM 23
#define SIM_SAVE_DT       0.01 * 2 * M_PI
#define SIM_SAVE_DT_MOONS pow(10,-4) * 2 * M_PI

struct reb_simualtion* init_random_simulation(struct reb_simulation* r);
struct reb_simulation* init_simulation_from_file(struct reb_simulation* r, char* file_name, int file_line);
struct reb_simulation* add_moons_from_file(struct reb_simulation*r, char* file_name, int file_line);
void heartbeat(struct reb_simulation* r);
void heartbeat_moons(struct reb_simulation* r);

int find_outcome(struct reb_simulation* r);    //What happens at the end of the simulation
double distance_between_particles(struct reb_simulation* r, int p1, int p2);
void find_min_distance(struct reb_simulation* r);
int is_close_encounter(struct reb_simulation* r);
int is_new_close_encounter(struct reb_simulation* r);

void read_file(char* file_name, int file_line, double* array_to_write, int N_parameters);
void write_output(struct reb_simulation* r, char* output_file);
void write_output_moon(struct reb_simulation* r, char* output_file);
char* get_name_with_pid(char* initial_string);
double min(double x, double y);
double rayleigh(double sigma);

double input_planets[N_SIM_PARAM];
double input_moons[N_MOONS_PARAM];

double min_d[3] = {1000, 1000, 1000};
double tmp_d[3];
int close_encounter_happening = 0;
int close_encounters_counter = 0;
double start_last_close_encounter = 0;
double end_last_close_encounter = 0;

int main(int argc, char* argv[]){
  pid_t pid = getpid();
  srand(pid * time(0));                                         // seeds the random number generators

  struct reb_simulation* r = reb_create_simulation();
  r->dt                = 0.01*2.*M_PI;                          // initial timestep
  r->integrator        = REB_INTEGRATOR_IAS15;                  // select the integrator
  r->heartbeat         = heartbeat;                              // select the function for managing the output
  r->exit_min_distance = 5*R_J;
  if(USE_RANDOM_INITIAL_CONDITIONS) r->exit_max_distance = 100.;    // if one body will get until 100 au from the star, the simulation will stop

  if(USE_INITIAL_CONDITIONS_FROM_FILE){
    r = init_simulation_from_file(r,"planets.txt", atoi(argv[1]));      // initiate the simulation
    if(SIMULATE_MOONS){
      reb_integrate(r, (input_planets[3] - 1) * 2 * M_PI);
      reb_output_binary(r, get_name_with_pid("snapshot"));
    }else{
      reb_integrate(r, input_planets[2] * 2 * M_PI);
    }
  }else if(USE_RANDOM_INITIAL_CONDITIONS){
    r = init_random_simulation(r);
    reb_integrate(r, T_MAX); // integrate the simulation r for 10 Myr
  }


  printf("%lf %lf %lf\n", min_d[0], min_d[1], min_d[2]);
  write_output(r, "output.out");

  struct reb_orbit o_first_planet = reb_tools_particle_to_orbit(r->G, r->particles[1], r->particles[0]);
  if(r->t <= T_MAX && find_outcome(r)==1 && min_d[0] >= 5*R_J && min_d[1] >= 5*R_J && min_d[2] >= 5*R_J && o_first_planet.e > 1.){
    write_output(r, "planets.txt");
  }

  reb_free_simulation(r);


  if(SIMULATE_MOONS){
    struct reb_simulation* r = reb_create_simulation_from_binary(get_name_with_pid("snapshot"));
    r->dt                = pow(10,-3) *2 * M_PI;                          // initial timestep
    r->integrator        = REB_INTEGRATOR_IAS15;                  // select the integrator
    r->heartbeat         = heartbeat_moons;                              // select the function for managing the output
    r->exit_min_distance = pow(10,-12);
    r->ri_ias15.min_dt   = pow(10,-5) * 2 * M_PI;
    r = add_moons_from_file(r,"moons.txt", atoi(argv[2]));      // initiate the simulation
    reb_integrate(r, (input_planets[3] + 10)* 2 * M_PI);
  }

  write_output_moon(r, get_name_with_pid("moon_output_"));

  reb_free_simulation(r);
}

void heartbeat(struct reb_simulation* r){

  if(reb_output_check(r, SIM_SAVE_DT) && USE_INITIAL_CONDITIONS_FROM_FILE){
    if(SIMULATE_MOONS){
      double percent = ( r->t / 2. / M_PI * 100 ) / input_planets[3];
      printf("PLANETS: %.1lf %%\r", percent);
    }else{
      double percent = ( r->t / 2. / M_PI * 100 ) / input_planets[2];
      printf("PLANETS: %.1lf %%\r", percent);
    }
  }

  if(SAVE_SIMULATION && reb_output_check(r, SIM_SAVE_DT) && r->t/2./M_PI > input_planets[2]-1000){
    reb_output_ascii(r,get_name_with_pid("sim_"));
  }

  if(is_close_encounter(r) == 0){
    if(close_encounter_happening == 1) end_last_close_encounter = r->t / 2. / M_PI;
    close_encounter_happening = 0;
  }else if(is_new_close_encounter(r) == 1){
    close_encounters_counter += 1;
    start_last_close_encounter = r->t / 2. / M_PI;
  }

  find_min_distance(r);
  int i;
  for(i = 0; i<3; i++){
    if(min_d[i] > tmp_d[i]){
      min_d[i] = tmp_d[i];
    }
  }
}

void heartbeat_moons(struct reb_simulation* r){

  if(reb_output_check(r, SIM_SAVE_DT_MOONS)){
    double percent = ( (r->t / 2. / M_PI  - input_planets[3] + 20) * 100 ) / (input_planets[2] - input_planets[3] + 20);
    printf("MOONS: %.8lf %%\r", percent);
  }

  if(SAVE_SIMULATION_MOONS && reb_output_check(r, SIM_SAVE_DT_MOONS)){
    reb_output_ascii(r,get_name_with_pid("sim_moons_"));
  }

}

void find_min_distance(struct reb_simulation* r){
  double d12 = distance_between_particles(r,1,2);
  double d13 = distance_between_particles(r,1,3);
  double d23 = distance_between_particles(r,2,3);
  tmp_d[0] = min(d12,d13);
  tmp_d[1] = min(d12,d23);
  tmp_d[2] = min(d13,d23);
}

double distance_between_particles(struct reb_simulation* r, int p1, int p2){
  double delta_x = r->particles[p1].x - r->particles[p2].x;
  double delta_y = r->particles[p1].y - r->particles[p2].y;
  double delta_z = r->particles[p1].z - r->particles[p2].z;
  return sqrt(pow(delta_x,2) + pow(delta_y,2) + pow(delta_z,2));
}

int is_close_encounter(struct reb_simulation* r){
  if(distance_between_particles(r,1,2)<600*R_J || distance_between_particles(r,1,3)<600*R_J) return 1;
  else return 0;
}

int is_new_close_encounter(struct reb_simulation* r){
  if(is_close_encounter(r) == 1 && close_encounter_happening == 0){
    close_encounter_happening = 1;
    return 1;
  }else{
    return 0;
  }
}

// Find the outcome of a given simulation
int find_outcome(struct reb_simulation* r){
  if(distance_between_particles(r,1,0) > 99.){
    return 1;
  }else if(distance_between_particles(r,2,0) > 99.){
    return 2;
  }else if(distance_between_particles(r,3,0) > 99.){
    return 3;
  }else{
    return 0;
  }
}

// Initiate a simulation
struct reb_simualtion* init_random_simulation(struct reb_simulation* r){
  // Add star
  struct reb_particle star = {0};
  star.m = 1;
  reb_add(r, star);

  // Add planets
  double m, a, e, inc, Omega, omega, M, f, rnd;
  m = M_J;
  input_planets[N_SIM_PARAM - 18] = a = 5.;
  input_planets[N_SIM_PARAM - 17] = e = rayleigh(0.01);
  input_planets[N_SIM_PARAM - 16] = inc = rayleigh(0.01);
  input_planets[N_SIM_PARAM - 15] = Omega = 2.*M_PI* (double)(rand()) / ((double)RAND_MAX);
  input_planets[N_SIM_PARAM - 14] = omega = 2.*M_PI* (double)(rand()) / ((double)RAND_MAX);
  M = 2.*M_PI* (double)(rand()) / ((double)RAND_MAX);
  input_planets[N_SIM_PARAM - 13] = f = reb_tools_M_to_f(e,M);
  struct reb_particle planet1 = reb_tools_orbit_to_particle(r->G, star, m, a, e, inc, Omega, omega, f);
  reb_add(r,planet1);


  rnd = (double)(rand()) / ((double)RAND_MAX);
  struct reb_orbit o_first_planet = reb_tools_particle_to_orbit(r->G, r->particles[1], r->particles[0]);
  input_planets[N_SIM_PARAM - 12] = a = o_first_planet.a * pow((rnd*0.2 + 1.2) ,2./3);
  input_planets[N_SIM_PARAM - 11] = e = rayleigh(0.01);
  input_planets[N_SIM_PARAM - 10] = inc = rayleigh(0.01);
  input_planets[N_SIM_PARAM -  9] = Omega = 2.*M_PI* (double)(rand()) / ((double)RAND_MAX);
  input_planets[N_SIM_PARAM -  8] = omega = 2.*M_PI* (double)(rand()) / ((double)RAND_MAX);
  M = 2.*M_PI* (double)(rand()) / ((double)RAND_MAX);
  input_planets[N_SIM_PARAM -  7] = f = reb_tools_M_to_f(e,M);
  struct reb_particle planet2 = reb_tools_orbit_to_particle(r->G, star, m, a, e, inc, Omega, omega, f);
  reb_add(r,planet2);

  rnd = (double)(rand()) / ((double)RAND_MAX);
  struct reb_orbit o_second_planet = reb_tools_particle_to_orbit(r->G, r->particles[2], r->particles[0]);
  input_planets[N_SIM_PARAM -  6] = a = o_second_planet.a * pow((rnd*0.2 + 1.2) ,2./3);
  input_planets[N_SIM_PARAM -  5] = e = rayleigh(0.01);
  input_planets[N_SIM_PARAM -  4] = inc = rayleigh(0.01);
  input_planets[N_SIM_PARAM -  3] = Omega = 2.*M_PI* (double)(rand()) / ((double)RAND_MAX);
  input_planets[N_SIM_PARAM -  2] = omega = 2.*M_PI* (double)(rand()) / ((double)RAND_MAX);
  M = 2.*M_PI* (double)(rand()) / ((double)RAND_MAX);
  input_planets[N_SIM_PARAM -  1] = f = reb_tools_M_to_f(e,M);
  struct reb_particle planet3 = reb_tools_orbit_to_particle(r->G, star, m, a, e, inc, Omega, omega, f);
  reb_add(r,planet3);

  reb_move_to_com(r);        // This makes sure the planetary systems stays within the computational domain and doesn't drift.

  return r;
}

struct reb_simulation* init_simulation_from_file(struct reb_simulation* r, char* file_name, int file_line){
  read_file(file_name, file_line, input_planets, N_SIM_PARAM);

  // Add star
  struct reb_particle star = {0};
  star.m = 1;
  reb_add(r, star);

  // Add planets
  double m, a, e, inc, Omega, omega, f;
  m = M_J;
  a = input_planets[N_SIM_PARAM - 18];
  e = input_planets[N_SIM_PARAM - 17];
  inc = input_planets[N_SIM_PARAM - 16];
  Omega = input_planets[N_SIM_PARAM - 15];
  omega = input_planets[N_SIM_PARAM - 14];
  f = input_planets[N_SIM_PARAM - 13];
  struct reb_particle planet1 = reb_tools_orbit_to_particle(r->G, star, m, a, e, inc, Omega, omega, f);
  reb_add(r,planet1);

  a = input_planets[N_SIM_PARAM - 12];
  e = input_planets[N_SIM_PARAM - 11];
  inc = input_planets[N_SIM_PARAM - 10];
  Omega = input_planets[N_SIM_PARAM - 9];
  omega = input_planets[N_SIM_PARAM - 8];
  f = input_planets[N_SIM_PARAM - 7];
  struct reb_particle planet2 = reb_tools_orbit_to_particle(r->G, star, m, a, e, inc, Omega, omega, f);
  reb_add(r,planet2);

  a = input_planets[N_SIM_PARAM - 6];
  e = input_planets[N_SIM_PARAM - 5];
  inc = input_planets[N_SIM_PARAM - 4];
  Omega = input_planets[N_SIM_PARAM - 3];
  omega = input_planets[N_SIM_PARAM - 2];
  f = input_planets[N_SIM_PARAM - 1];
  struct reb_particle planet3 = reb_tools_orbit_to_particle(r->G, star, m, a, e, inc, Omega, omega, f);
  reb_add(r,planet3);

  reb_move_to_com(r);        // This makes sure the planetary systems stays within the computational domain and doesn't drift.
  return r;
}

struct reb_simulation* add_moons_from_file(struct reb_simulation*r, char* file_name, int file_line){
  read_file(file_name, file_line, input_moons, N_MOONS_PARAM);
  int i;
  double m, a, e, inc, Omega, omega, f;
  for(i=0;i<N_MOONS;i++){
    m = 0;
    a = input_moons[0 + 3*i] * R_J;
    e = input_moons[1 + 3*i];
    inc = input_moons[2 + 3*i];
    Omega = 0.;
    omega = 0.;
    f = 0.;
    struct reb_particle moon = reb_tools_orbit_to_particle(r->G, r->particles[1], m, a, e, inc, Omega, omega, f);
    reb_add(r,moon);
  }

  reb_move_to_com(r);
  return r;
}

void read_file(char* file_name, int file_line, double* array_to_write, int N_parameters){
  FILE* fp;
  fp = fopen(file_name,"r");
  double garbage;
  int i;
  for(i=0;i<N_parameters*(file_line-1);i++) fscanf(fp,"%lf",&garbage);  //SKIP THE FIRST LINES
  for(i = 0; i<N_parameters; i++) fscanf(fp,"%lg",&array_to_write[i]);   //READ ONLY THE CHOSEN ONE
  fclose(fp);
}

void write_output(struct reb_simulation* r, char* output_file){
  struct reb_orbit o_first_planet = reb_tools_particle_to_orbit(r->G, r->particles[1], r->particles[0]);
  FILE* fp;
  fp = fopen(output_file,"a");
  fprintf(fp, "%d %lf %lf %lf %lf ", find_outcome(r), o_first_planet.e, r->t / 2. / M_PI, start_last_close_encounter, end_last_close_encounter);
  for(int i=-18; i<0; i++) fprintf(fp, "%.17g ", input_planets[N_SIM_PARAM + i]);
  fprintf(fp, "\n");
  fclose(fp);
}

void write_output_moon(struct reb_simulation* r, char* output_file){
  int id_moons[N_MOONS] = {0};
  double e_moons[N_MOONS];
  double a_moons[N_MOONS];
  double i_moons[N_MOONS];
  struct reb_orbit o = reb_tools_particle_to_orbit(r->G, r->particles[1], r->particles[0]);
  struct reb_orbit moon_orbit;
  for(int i = 0; i<N_MOONS; i++){
    if(distance_between_particles(r,1,i+4) < fabs(o.rhill)){
      id_moons[i] = 1;
    }
  moon_orbit = reb_tools_particle_to_orbit(r->G, r->particles[i+4], r->particles[1]);
  e_moons[i] = moon_orbit.e;
  a_moons[i] = moon_orbit.a;
  i_moons[i] = moon_orbit.inc;
  }

  FILE* fp2;
  fp2 = fopen(output_file, "w");
  for(int i=0; i<N_MOONS; i++) fprintf(fp2, "%d %lf %lf %lf %lf %lf %lf\n", id_moons[i], input_moons[0 + 3*i] * R_J, input_moons[1 + 3*i], input_moons[2 + 3*i], a_moons[i], e_moons[i], i_moons[i]);
  fclose(fp2);
}

char* get_name_with_pid(char* initial_string){
  pid_t pid = getpid();
  char buffer[50];
  sprintf(buffer, "%d", pid);
  static char str[80];
  strcpy(str, initial_string);
  strcat(str, buffer);
  return str;
}

double min(double x, double y){
  if(x >= y){
    return y;
  }else{
    return x;
  }
}

// Generate a number according to a Rayleigh distribution
double rayleigh(double sigma){
  double r = (double)(rand()) / ((double)RAND_MAX);
  return sigma * sqrt(2*log(1./(1.-r)));
}
