/* This program is based on a FORTRAN90 program by Bill Magro.

Licensing:

 This code is distributed under the GNU LGPL license.

Modified:

21 November 2007

Author:

John Burkardt.*/

# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <omp.h>

int main ( int argc, char *argv[] );
void compute ( int np, int nd, double pos[], double vel[], double mass, double f[], double *pot, double *kin );
double dist ( int nd, double r1[], double r2[], double dr[] );
void initialize ( int np, int nd, double pos[], double vel[], double acc[] );
void r8mat_uniform_ab ( int m, int n, double a, double b, int *seed, double r[] );
void update ( int np, int nd, double pos[], double vel[], double f[], double acc[], double mass, double dt );

/******************************************************************************/

int main ( int argc, char *argv[] )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for MD.

  Discussion:

    MD implements a simple molecular dynamics simulation.

    The velocity Verlet time integration scheme is used. 

    The particles interact with a central pair potential.
*/
{
  double *acc;
  double dt=0.0001; // time step
  double e0;
  double *force;
  int i;
  int id;
  double kinetic;
  double mass = 1.0;
  int nd=3;      //spatial dimension
  int np=2000;  // number of particles
  double *pos;
  double potential;
  int step;
  int step_num=100; //number of time steps
  int step_print;
  int step_print_index;
  int step_print_num;
  double *vel;

    
    acc = ( double * ) malloc ( nd * np * sizeof ( double ) );
    force = ( double * ) malloc ( nd * np * sizeof ( double ) );
    pos = ( double * ) malloc ( nd * np * sizeof ( double ) );
    vel = ( double * ) malloc ( nd * np * sizeof ( double ) );
    
/*
  This is the main time stepping loop:
    Compute forces and energies,
    Update positions, velocities, accelerations.
*/
 
  for ( step = 0; step <= step_num; step++ )
  {
    if ( step == 0 )
    {
      initialize ( np, nd, pos, vel, acc );
    }
    else
    {
      update ( np, nd, pos, vel, force, acc, mass, dt );
    }

    compute ( np, nd, pos, vel, mass, force, &potential, &kinetic );

    if ( step == 0 )
    {
      e0 = potential + kinetic;
    }
      
  }
    
    printf("potential=%f, kinetic=%f, %f\n", potential, kinetic, (potential+kinetic-e0)/e0);

/*
  Free memory.
*/
  free ( acc );
  free ( force );
  free ( pos );
  free ( vel );
/*
  Terminate.
*/
  return 0;
}
/******************************************************************************/

void compute ( int np, int nd, double pos[], double vel[], double mass, 
  double f[], double *pot, double *kin )

/******************************************************************************/
/*
  Purpose:

    COMPUTE computes the forces and energies.

  Discussion:

    The potential function V(X) is a harmonic well which smoothly
    saturates to a maximum value at PI/2:

      v(x) = ( sin ( min ( x, PI/2 ) ) )^2

    The derivative of the potential is:

      dv(x) = 2.0 * sin ( min ( x, PI/2 ) ) * cos ( min ( x, PI/2 ) )
            = sin ( 2.0 * min ( x, PI/2 ) )

  Parameters:

    Input, int NP, the number of particles.

    Input, int ND, the number of spatial dimensions.

    Input, double POS[ND*NP], the positions.

    Input, double VEL[ND*NP], the velocities.

    Input, double MASS, the mass of each particle.

    Output, double F[ND*NP], the forces.

    Output, double *POT, the total potential energy.

    Output, double *KIN, the total kinetic energy.
*/
{
  double d;
  double d2;
  double pe;
  double ki;
  int i;
  int j;
  int k;
  double PI2 = 3.141592653589793 / 2.0;
  double rij[3];

  pe = 0.0;
  ki = 0.0;

  #pragma parallel shared ( f, nd, np, pos, vel,PI2,pe, ki )  private ( i, j, k, rij, d, d2 ){
    #pragma for reduction (+:pe) reduction(+:ki)
    for ( k = 0; k < np; k++ )
    {
    /*
    Compute the potential energy and forces.
    */
        for ( i = 0; i < nd; i++ )
        {
            f[i+k*nd] = 0.0;
        }

        for ( j = 0; j < np; j++ )
        {
            if ( k != j )
            {
                d = dist ( nd, pos+k*nd, pos+j*nd, rij );
    /*
    Attribute half of the potential energy to particle J.
    */
                if ( d < PI2 )
                {
                    d2 = d;
                }
                else
                {
                    d2 = PI2;
                }

                pe = pe + 0.5 * pow ( sin ( d2 ), 2 );

                for ( i = 0; i < nd; i++ )
                {
                    f[i+k*nd] = f[i+k*nd] - rij[i] * sin ( 2.0 * d2 ) / d;
                }
        }
        }
    /*
    Compute the kinetic energy.
    */
        for ( i = 0; i < nd; i++ )
        {
        ki = ki + vel[i+k*nd] * vel[i+k*nd];
        }
    }

  }

  ki = ki * 0.5 * mass;
  
  *pot=pe;
  *kin=ki;
    
  return;
}

double dist ( int nd, double r1[], double r2[], double dr[] )

/******************************************************************************/
/*
  Purpose:

    DIST computes the displacement (and its norm) between two particles.

  Parameters:

    Input, int ND, the number of spatial dimensions.

    Input, double R1[ND], R2[ND], the positions of the particles.

    Output, double DR[ND], the displacement vector.

    Output, double D, the Euclidean norm of the displacement.
*/
{
  double d;
  int i;

  d = 0.0;
  for ( i = 0; i < nd; i++ )
  {
    dr[i] = r1[i] - r2[i];
    d = d + dr[i] * dr[i];
  }
  d = sqrt ( d );

  return d;
}
/******************************************************************************/

void initialize ( int np, int nd, double pos[], double vel[], double acc[] )

/******************************************************************************/
/*
  Purpose:

    INITIALIZE initializes the positions, velocities, and accelerations.

  Parameters:

    Input, int NP, the number of particles.

    Input, int ND, the number of spatial dimensions.

    Output, double POS[ND*NP], the positions.

    Output, double VEL[ND*NP], the velocities.

    Output, double ACC[ND*NP], the accelerations.
*/
{
  int i;
  int j;
  int seed;
/*
  Set positions.
*/
  seed = 123456789;
  r8mat_uniform_ab ( nd, np, 0.0, 10.0, &seed, pos );
/*
  Set velocities.
*/
  for ( j = 0; j < np; j++ )
  {
    for ( i = 0; i < nd; i++ )
    {
      vel[i+j*nd] = 0.0;
    }
  }
/*
  Set accelerations.
*/
  for ( j = 0; j < np; j++ )
  {
    for ( i = 0; i < nd; i++ )
    {
      acc[i+j*nd] = 0.0;
    }
  }

  return;
}
/******************************************************************************/

void r8mat_uniform_ab ( int m, int n, double a, double b, int *seed, double r[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_UNIFORM_AB returns a scaled pseudorandom R8MAT.

  Discussion:

    This routine implements the recursion

      seed = 16807 * seed mod ( 2^31 - 1 )
      unif = seed / ( 2^31 - 1 )

    The integer arithmetic never requires more than 32 bits,
    including a sign bit.

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input, double A, B, the limits of the pseudorandom values.

    Input/output, int *SEED, the "seed" value.  Normally, this
    value should not be 0.  On output, SEED has 
    been updated.

    Output, double R[M*N], a matrix of pseudorandom values.
*/
{
  int i;
  const int i4_huge = 2147483647;
  int j;
  int k;

  if ( *seed == 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8MAT_UNIFORM_AB - Fatal error!\n" );
    fprintf ( stderr, "  Input value of SEED = 0.\n" );
    exit ( 1 );
  }

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      k = *seed / 127773;

      *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

      if ( *seed < 0 )
      {
        *seed = *seed + i4_huge;
      }
      r[i+j*m] = a + ( b - a ) * ( double ) ( *seed ) * 4.656612875E-10;
    }
  }

  return;
}

void update ( int np, int nd, double pos[], double vel[], double f[], 
  double acc[], double mass, double dt )

/******************************************************************************/
/*
  Purpose:

    UPDATE updates positions, velocities and accelerations.

  Discussion:

    The time integration is fully parallel.

    A velocity Verlet algorithm is used for the updating.

    x(t+dt) = x(t) + v(t) * dt + 0.5 * a(t) * dt * dt
    v(t+dt) = v(t) + 0.5 * ( a(t) + a(t+dt) ) * dt
    a(t+dt) = f(t) / m

  Parameters:

    Input, int NP, the number of particles.

    Input, int ND, the number of spatial dimensions.

    Input/output, double POS[ND*NP], the positions.

    Input/output, double VEL[ND*NP], the velocities.

    Input, double F[ND*NP], the forces.

    Input/output, double ACC[ND*NP], the accelerations.

    Input, double MASS, the mass.

    Input, double DT, the time step.
*/
{
  int i;
  int j;
  double rmass;

  rmass = 1.0 / mass;
  
  # pragma omp parallel shared ( acc, dt, f, nd, np, pos, rmass, vel ) private ( i, j )
  # pragma omp for collapse (2)
  for ( j = 0; j < np; j++ )
  {
    for ( i = 0; i < nd; i++ )
    {
      pos[i+j*nd] = pos[i+j*nd] + vel[i+j*nd] * dt + 0.5 * acc[i+j*nd] * dt * dt;
      vel[i+j*nd] = vel[i+j*nd] + 0.5 * dt * ( f[i+j*nd] * rmass + acc[i+j*nd] );
      acc[i+j*nd] = f[i+j*nd] * rmass;
    }
  }

  return;
}

