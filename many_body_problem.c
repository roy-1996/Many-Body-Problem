/* Since the size of the log file and the trajectory file are beyond the permissible limits of Moodle, the files have been uploaded on Google Drive,
   and a shareable link for the files has been provided in the Readme file */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <omp.h>
# include <time.h>
# define N 1000

const double dT = 0.005;
const double G = 6.67E-11;


struct coord
{
	double x,y,z;
};

typedef struct coord Point;
Point position[N];							/* position[] stores the coordinates of the bodies */


struct vel
{
	double vx,vy,vz;
};

typedef struct vel Velocity;
Velocity velocity[N];						/* velocity[] stores the velocity of the bodies in x,y and z directions */


struct acc
{
	double ax,ay,az;
};

typedef struct acc Acceleration;
Acceleration A[N];							/* A[] stores the acceleration of the bodies in x,y and z directions */


struct avg_acc
{
	double ax,ay,az;
};

typedef struct avg_acc Avg_Acceleration;
Avg_Acceleration A_avg[N];					/* A_avg[] stores the average acceleration of the bodies in x,y and z directions */



int main()
{

	int i,body0,body1,k = 1,number_of_iterations = 0 ;
	double distance , cos_theta_x , cos_theta_y , cos_theta_z , new_x_ac , new_y_ac , new_z_ac ,force ;
	double begin_time,end_time,diff,time_elapsed = 0,actual_simulation_time,point1,prev_x,prev_y,prev_z,gap1,gap2;
	FILE *fp,*fp1;
	char data[100],data1[100];

	
	fp = fopen("trajectory1.txt","w");
	fp1 = fopen("log.txt","a");

	char comm[20] = "sh write.sh";
	system(comm);
	
	int num_of_processors = omp_get_num_procs();
	int num_of_threads = omp_get_max_threads();
	
	sprintf(data,"Number of processors:\t %d",num_of_processors);
	fputs(data,fp1);
	
	sprintf(data,"\nNumber of threads:\t %d",num_of_threads);
	fputs(data,fp1);

	point1 = omp_get_wtime();
	
	
	# pragma omp parallel for private(i)
	for ( i = 0 ; i < N ; i++ )
	{
		velocity[i].vx = velocity[i].vy = velocity[i].vz = 0;
		A[i].ax = A[i].ay = A[i].az = 0;
		A_avg[i].ax = A_avg[i].ay = A_avg[i].az = 0;

	}


	/* The for loop generates coordinates for the bodies in a random fashion */

	srand(time(NULL));

	//# pragma omp parallel for private(i)
	for ( i  = 0 ; i < N ; i++ )
	{
		position[i].x = rand()%100+0.0001;
		position[i].y = rand()%200+0.0001;
		position[i].z = rand()%400+0.0001;
	}


	do
	{

		begin_time = omp_get_wtime();

		printf("\nIteration number : %d",number_of_iterations);

		for ( body0 = 0 ; body0 < N ; body0++ )
		{
			/*
			A[body0].ax = A_avg[body0].ax;
			A[body0].ay = A_avg[body0].ay;
			A[body0].az = A_avg[body0].az;
			*/

			/* Calculating the new position of the body in the x-direction */
			
			prev_x = position[body0].x;
			position[body0].x += velocity[body0].vx * dT + 0.5 * A[body0].ax  * dT * dT ;
			
			if ( position[body0].x > 100 || position[body0].x < 0 )   /* When the body goes beyond the boundary. */
			{
			    if ( position[body0].x > 100 )
			    {
				gap1 = abs( position[body0].x - prev_x );
				gap2 = abs( (double)100 - prev_x );
				
				gap1 = gap1 - gap2;
				position[body0].x = (double)100 - gap1;
				
				velocity[body0].vx = -velocity[body0].vx;
				
			    }
			    else if ( position[body0].x < 0 )
			    {
			       gap1 = abs(position[body0].x);
			       position[body0].x = gap1;
			       
			       velocity[body0].vx = -velocity[body0].vx;
			    }
			    
			}
			       
			
			prev_y = position[body0].y;
			position[body0].y += velocity[body0].vy * dT + 0.5 * A[body0].ay  * dT * dT ;

			/* Calculating the new position of the body in the y-direction */
		
			if ( position[body0].y > 200 || position[body0].y < 0 )
			{
			    if ( position[body0].y > 200 )
			    {
				gap1 = abs( position[body0].y - prev_y );
				gap2 = abs( (double)200 - prev_y );
				
				gap1 = gap1 - gap2;
				position[body0].y = (double)200 - gap1;
				
				velocity[body0].vy = -velocity[body0].vy;
				
			    }
			    else if ( position[body0].y < 0 )
			    {
			       gap1 = abs(position[body0].y);
			       position[body0].y = gap1;
			       
			       velocity[body0].vy = -velocity[body0].vy;
			       
			    }
			    
			}
			
			prev_z = position[body0].z;
			position[body0].z += velocity[body0].vz * dT + 0.5 * A[body0].az  * dT * dT ;
			

			/* Calculating the new position of the body in the z-direction */
			
			if ( position[body0].z > 400 || position[body0].z < 0 )
			{
			    if ( position[body0].z > 400 )
			    {
				gap1 = abs( position[body0].z - prev_z );
				gap2 = abs( (double)400 - prev_z );
				
				gap1 = gap1 - gap2;
				position[body0].z = (double)400 - gap1;
				
				velocity[body0].vz = -velocity[body0].vz;
				
			    }
			    else if ( position[body0].z < 0 )
			    {
			       gap1 = abs(position[body0].z);
			       position[body0].z = gap1;
			       
			       velocity[body0].vz = -velocity[body0].vz;
			       
			    }
			    
			}


			/* For each body the acceleration in the x,y and z direction are being found out in the for loop below.*/
			
			# pragma omp parallel for private(body1,cos_theta_x,cos_theta_y,cos_theta_z,distance,force) reduction(+:new_x_ac,new_y_ac,new_z_ac)

			for ( body1 = 0 ; body1 < N ; body1++ )
			{
				if ( body1 != body0 )
				{

					/* distance is the Euclidean distance between the position of the two bodies. */

					distance = sqrt( pow( position[body1].x - position[body0].x,2 ) + pow( position[body1].y - position[body0].y,2 ) + pow (position[body1].z - position[body0].z,2) );

					cos_theta_x = ( position[body1].x - position[body0].x )/distance; /* cos_theta_x is the cosine of the angle made by the displacement vector with the x-axis */
					cos_theta_y = ( position[body1].y - position[body0].y )/distance; /* cos_theta_y is the cosine of the angle made by the displacement vector with the y-axis */
					cos_theta_z = ( position[body1].z - position[body0].z )/distance; /* cos_theta_z is the cosine of the angle made by the displacement vector with the z-axis */


					/*force is the gravitational force acting between the two bodies body0 and body1 */

					force = G * 1 * 1 /pow(distance,2);

					/* Adding up the contribution of the forces in the x,y and the z directions on the body body0. */

					new_x_ac += force*cos_theta_x;
					new_y_ac += force*cos_theta_y;
					new_z_ac += force*cos_theta_z;


				}

			}



			/* Calculating the average acceleration of the body body0 in each of the x,y and z directions. */ 

			A_avg[body0].ax = ( A[body0].ax + new_x_ac )/2.0;
			A_avg[body0].ay = ( A[body0].ay + new_y_ac )/2.0;
			A_avg[body0].az = ( A[body0].az + new_z_ac )/2.0;
			

			/* Finally the velocity of the body body0 in the x,y and the z directions is being computed */

			velocity[body0].vx += A_avg[body0].ax * dT;
			velocity[body0].vy += A_avg[body0].ay * dT;
			velocity[body0].vz += A_avg[body0].az * dT;


		}
	
	

		if ( number_of_iterations % 1000 == 0 )
		{
			fprintf(fp,"\nIteration number :\t%d",number_of_iterations);

			for ( i = 0 ; i < N ; i++ )
			{
				fprintf(fp,"\nPosititon of body %d: (%.20lf,%.20lf,%.20lf)",i,position[i].x,position[i].y,position[i].z);
			}

			fputs("\nTERMINATE",fp);
		}
		  

		end_time = omp_get_wtime();
		diff = end_time-begin_time;


		fprintf(fp1,"\nTime required for %d step: %lf",k,diff);

		time_elapsed += 0.005;
		k += 1;
  
		number_of_iterations++;

	}while(time_elapsed <= 3600 );




fprintf(fp1,"\nTotal simulation time = %lf",omp_get_wtime() - point1 );
fclose(fp);
fclose(fp1);

}












