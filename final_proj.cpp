//this is the final project: Lennard-Jones fluid on a graphene surface
//nuber of fluid atoms: 4000; number of carbon atom in graphene 1152; 
//Nose-Hooever thermostat is not finished due to the limit of time for this project
//In our molecular dynamics simulation systems, periodic boundary condition is utilized in 
//all direction; the velocity of graphene is set to be 0;
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <string.h>
#include <math.h>
#include <time.h>
using namespace std;
#define num_fluid_atom 4000
#define num_carbon_atom 1152
#define num_tot_atom    5152
#define delta_t         1.0    //in the unit of fs
#define x_min           0.0
#define x_max           51.0
#define y_min           0.0
#define y_max           59.0 
#define z_min           0.0
#define z_max           80.0
struct   data0 {
	double coordinate[num_tot_atom][3];
	double velo[num_tot_atom][3];
	int    atom_tag[num_tot_atom];

};
struct  force{
        double force_array[num_tot_atom][3];
};
struct  L_J{
	//tag =1 or 2; epsilon: kcal/mol; sigma: angstrom
	int tag;
	double epsilon;
	double sigma;
};
int init_sys(struct data0 *sys1); // the function which initialize the system
int update_sys(struct data0 *sys1,struct force *pot_array1, struct L_J *LJ1,struct L_J *LJ2,double cut_off_r);
int search_pbc(struct data0 *sys1,double pi[3],double pj[3],struct force *force_array1,struct L_J *LJ1,struct L_J *LJ2,double cut_off_r,int ii,int jj);
int cal_force(struct data0 *sys1,struct force *force_array1,double x_ij, double y_ij,double z_ij,double sigmaij, double epsilonij,int ii,int jj);
int init_force(struct force *force_array1);
int update_velo(struct data0 *sys1, struct force *force_array1);
int update_coor(struct data0 *sys1);
int write_file(struct data0 *sys1,int time);

int main(int nbargs, char* args[]) {
	struct data0 sys1;
	struct force force_array1;
	struct L_J LJ1,LJ2;
	init_sys(&sys1);
        int max_step,ii;
	max_step=1000;//in the unit of fs
	double cut_off_r;
	cut_off_r=7.0;
        FILE *file1;
	char ch3[66];
	//read Lennar-Jones potential parameters
	file1=fopen("LJ.file","r");
	fscanf(file1,"%d %lf %lf\n",&(LJ1.tag),&(LJ1.epsilon),&(LJ1.sigma));
	fscanf(file1,"%d %lf %lf\n",&(LJ2.tag),&(LJ2.epsilon),&(LJ2.sigma));
	fclose(file1);
	init_force(&force_array1);
	for (ii=0;ii<max_step;++ii){
		//before update the system according LJ potential, set force to 0
		cout << ii << "fs" << "  " << "OK\n";
		init_force(&force_array1);
		update_sys(&sys1,&force_array1,&LJ1,&LJ2,cut_off_r);
		write_file(&sys1,ii);
	}

	return 0;
}
// output coordinate and velocity
int write_file(struct data0 *sys1,int time){

        FILE *file2;
	FILE *file3;
	file2=fopen("coords.out","w");
	file3=fopen("coords.xyz","w");
	fprintf(file2,"%s %d %s\n","output_of_coordiante_at_time",time,"fs");
	fprintf(file3,"%d\n",5152);
	fprintf(file3,"%s\n","    ");
	int ii,kk;
	for(ii=0;ii<num_tot_atom;++ii){
		kk=ii+1;
		fprintf(file2,"%d %d %lf %lf %lf\n",kk,sys1->atom_tag[ii],sys1->coordinate[ii][0],sys1->coordinate[ii][1],sys1->coordinate[ii][2]);
		fprintf(file3,"%s %lf %lf %lf\n","C",sys1->atom_tag[ii],sys1->coordinate[ii][0],sys1->coordinate[ii][1],sys1->coordinate[ii][2]);
	}
	fclose(file2);
	fclose(file3);
	return 1;
}
//update the system
int update_sys(struct data0 *sys1,struct force *force_array1,struct L_J *LJ1,struct L_J *LJ2,double cut_off_r){
	int ii,jj,kk,ll;
	double dis_ij,pi[3],pj[3],NA,const0,epsilon11,epsilon12,epsilon22;
	double sigma11,sigma12,sigma22,sigmaij,epsilonij;
	double f_x,f_y,f_z,x_ij,y_ij,z_ij;
	NA=6.022*pow(10,-22); // avagodro number
	const0=24*4.184*1000/(NA*pow(10,-10));
	sigma11=LJ1->sigma;
	sigma22=LJ2->sigma;
	sigma12=0.5*(sigma11+sigma22);
	epsilon11=LJ1->epsilon;
	epsilon22=LJ2->epsilon;
	epsilon12=sqrt(epsilon11*epsilon22);
	int temp;
	temp=num_tot_atom-1;
	for (ii=0;ii<temp;ii++){
		//cout << ii << "  ii  OK\n";
		ll=ii+1;
		for(jj=ll;jj<num_tot_atom;jj++){
                 	//cout << jj << "  jj  OK\n";
			for (kk=0;kk<3;kk++){
				pi[kk]=sys1->coordinate[ii][kk];
				pj[kk]=sys1->coordinate[jj][kk];
			}
			x_ij=pj[0]-pi[0];
			y_ij=pj[1]-pi[1];
			z_ij=pj[2]-pi[2];
		        dis_ij=sqrt(x_ij*x_ij+y_ij*y_ij+z_ij*z_ij);
	                // if dis_ij<cut_off_r then the pair is used to calculate potential energy
			if (dis_ij<cut_off_r){
				if(sys1->atom_tag[ii]==1 && sys1->atom_tag[jj]==1){
					sigmaij=sigma11;
					epsilonij=epsilon11;
				}
				else if (sys1->atom_tag[ii]==2 && sys1->atom_tag[jj]==2) {
					sigmaij=sigma22;
					epsilonij=epsilon22;
				}
				else if (sys1->atom_tag[ii] != sys1->atom_tag[jj]){
					sigmaij=sigma12;
					epsilonij=epsilon12;
				}
				//the function to calculate and update force
				//cout << ii << " " << jj << "\n";
				cal_force(sys1,force_array1,x_ij,y_ij,z_ij,sigmaij,epsilonij,ii,jj);
				//cout << "OK1\n";
				//cout << "ok" << "\n";

			}
			else{
				//cout << ii << " " << jj << "\n";
				search_pbc(sys1,pi,pj,force_array1,LJ1,LJ2,cut_off_r,ii,jj);
				//cout << "OK2\n";
			}
			//update the coordinate
			//update_velo(sys1,force_array1);//update the veloity
			//update_coor(sys1);//update the coordinate
		}
	}
	 update_velo(sys1,force_array1);//update the veloity
	 update_coor(sys1);//update the coordinate


}
// use PBC to update potential when dis_ij> cut_off_r
int search_pbc(struct data0 *sys1,double pi[3],double pj[3],struct force *force_array1,struct L_J *LJ1,struct L_J *LJ2,double cut_off_r,int ii,int jj){
	double x_range,y_range,temp_x_ij,temp_y_ij,z_ij,dis_ij,sigma11,sigma22,sigma12,epsilon11,epsilon22,epsilon12;
	double sigmaij,epsilonij;
	x_range=x_max-x_min;
	y_range=y_max-y_min;
	z_ij=pj[2]-pj[1];
        sigma11=LJ1->sigma;
	sigma22=LJ2->sigma;
	sigma12=0.5*(sigma11+sigma22);
	epsilon11=LJ1->epsilon;
	epsilon22=LJ2->epsilon;
	epsilon12=sqrt(epsilon11*epsilon22);
	if(sys1->atom_tag[ii]==1 && sys1->atom_tag[jj]==1){
		sigmaij=sigma11;
		epsilonij=epsilon11;}
        else if (sys1->atom_tag[ii]==2 && sys1->atom_tag[jj]==2) {
		sigmaij=sigma22;
		epsilonij=epsilon22;}
        else if (sys1->atom_tag[ii] != sys1->atom_tag[jj]){
		sigmaij=sigma12;
		epsilonij=epsilon12;}
        for (ii=0;ii<2;++ii){
		for(jj=0;jj<2;++jj){
			temp_x_ij=pj[0]*x_range*pow(-1,ii)-pi[0];
			temp_y_ij=pj[1]*y_range*pow(-1,jj)-pi[1];
			dis_ij=sqrt(temp_x_ij*temp_x_ij+temp_y_ij*temp_y_ij+z_ij*z_ij);
			if (dis_ij<cut_off_r){
				cal_force(sys1,force_array1,temp_x_ij,temp_y_ij,z_ij,sigmaij,epsilonij,ii,jj);
			}
		}
	}




}
//calculate force
int cal_force(struct data0 *sys1,struct force *force_array1,double x_ij, double y_ij,double z_ij,double sigmaij, double epsilonij,int ii, int jj){
	double f_x,f_y,f_z,NA,const0;
	double dis_ij;
	NA=6.022*pow(10,23); // avagodro number
	const0=24*4.184*1000/(NA*pow(10,-10));
	dis_ij=sqrt(x_ij*x_ij+y_ij*y_ij+z_ij*z_ij);
	f_x=const0*epsilonij*(pow(sigmaij,6)/pow(dis_ij,8)-2*pow(sigmaij,12)/pow(dis_ij,14))*(x_ij); // force in the unit of N
        f_y=const0*epsilonij*(pow(sigmaij,6)/pow(dis_ij,8)-2*pow(sigmaij,12)/pow(dis_ij,14))*(y_ij);
        f_z=const0*epsilonij*(pow(sigmaij,6)/pow(dis_ij,8)-2*pow(sigmaij,12)/pow(dis_ij,14))*(z_ij);
	if (sys1->atom_tag[ii] !=2){//force on graphene is always 0 !!!!!
		force_array1->force_array[ii][0]=force_array1->force_array[ii][0]+f_x;
		force_array1->force_array[ii][1]=force_array1->force_array[ii][1]+f_y;
		force_array1->force_array[ii][2]=force_array1->force_array[ii][2]+f_z;
        }
	if (sys1->atom_tag[jj] !=2){
	        force_array1->force_array[jj][0]=force_array1->force_array[jj][0]-f_x;
	        force_array1->force_array[jj][1]=force_array1->force_array[jj][1]-f_y;
	        force_array1->force_array[jj][2]=force_array1->force_array[jj][2]-f_z;
	}

}

//update the velocity
int update_velo(struct data0 *sys1, struct force *force_array1){
	int ii;
	double mass1,mass2;
	mass1;12.001/(6.022*pow(10,23)*1000);
	mass2=12.001/(6.022*pow(10,23)*1000);
	for (ii=0;ii<num_tot_atom;++ii){
		if (sys1->atom_tag[ii] != 2)  // force on carbon in graphene remains 0 all the time
		{
			//mass of atom missing??!! 
			sys1->velo[ii][0]=sys1->velo[ii][0]+force_array1->force_array[ii][0]*delta_t/mass1*pow(10,-15); //in the unit of m/s
			sys1->velo[ii][1]=sys1->velo[ii][1]+force_array1->force_array[ii][1]*delta_t/mass1*pow(10,-15);
			sys1->velo[ii][2]=sys1->velo[ii][2]+force_array1->force_array[ii][2]*delta_t/mass1*pow(10,-15);
			

		}
	}
	return 1;
}

//update the coordinate
int update_coor(struct data0 *sys1){
	int ii;
	double temp_x,temp_y,temp_z,range_x,range_y,range_z;
        range_x=x_max-x_min;
        range_y=y_max-y_min;
        range_z=z_max-z_min;	
	for (ii=0;ii<num_tot_atom;++ii){
		//cout << sys1->velo[ii][0] << "\n";
	        temp_x=sys1->coordinate[ii][0]+sys1->velo[ii][0]*delta_t*pow(10,8);
		temp_y=sys1->coordinate[ii][1]+sys1->velo[ii][1]*delta_t*pow(10,8);
		temp_z=sys1->coordinate[ii][2]+sys1->velo[ii][2]*delta_t*pow(10,8);
		// periodic boundary condition
		if (temp_x>x_max){temp_x=temp_x-range_x;}
		if (temp_x<x_min){temp_x=temp_x+range_x;}
		if (temp_y>y_max){temp_y=temp_y-range_y;}
		if (temp_y<y_min){temp_y=temp_y+range_y;}
		if (temp_z>z_max){temp_z=temp_z-range_z;}
		if (temp_z<z_min){temp_z=temp_z+range_z;}
		sys1->coordinate[ii][0]=temp_x;
		sys1->coordinate[ii][1]=temp_y;
		sys1->coordinate[ii][2]=temp_z;
	}
	return 1;
}

//initialize the system
int init_force(struct force *force_array1){
	int ii;
	for (ii=0;ii<num_tot_atom;++ii){
		force_array1->force_array[ii][0]=0.0;
		force_array1->force_array[ii][1]=0.0;
		force_array1->force_array[ii][2]=0.0;
	}
	return 1;
}
int init_sys(struct data0 *sys1){
	//this is the code to initialize the system
        FILE * file0;
	file0=fopen("coords.inp","r");
	char ch1[31];
	fscanf(file0,"%s\n",ch1);
	double d1,d2,d3;
	int ii,jj,kk;
	fscanf(file0,"%d %lf\n",&ii,&d1);
	fscanf(file0,"%d %lf\n",&ii,&d2);
	char ch2[5];
	fscanf(file0,"%s\n",ch2);
	for (ii=0;ii<num_tot_atom;++ii){
		fscanf(file0,"%d %d %lf %lf %lf\n",&jj,&(sys1->atom_tag[ii]),&(sys1->coordinate[ii][0]),&(sys1->coordinate[ii][1]),&(sys1->coordinate[ii][2]));
		
	}
	//fscanf(file0,"%d %d %lf %lf %lf\n",&ii,&jj,&d1,&d2,&d3);
	fclose(file0);
	//give initialize velocity: 0 
	for (ii=0;ii<num_tot_atom;++ii){
		sys1->velo[ii][0]=0.0;
		sys1->velo[ii][1]=0.0;
		sys1->velo[ii][2]=0.0;
	}
	return 1;
}
