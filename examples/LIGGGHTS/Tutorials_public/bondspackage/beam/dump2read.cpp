// compile with `g++ -o dump2read dump2read.cpp`
// ussage: dump2read <liggghts-dumpfile> <outputfilename> <density of particle> <ID-Offset>
// e.g. `./dump2read post/dump1000.liggghts beam.pour 2500 0`

#include <stdio.h>
#include <sstream>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <string.h>
//#define M_PI    3.14159265358979323846f
using namespace std;
ifstream inFile;

  int main(int argc, char *argv[])
  {
    printf("Starting\n");
	
    

    double dens=atof(argv[3]);
    int offset=atoi(argv[4]);
    printf("infile:%s\noutfile:%s\n\nDensity:%f\nofffset:%d\n",argv[1],argv[2],dens,offset);

    inFile.open(argv[1]);
    
    if (!inFile) {
        printf("Unable to open file");
        return 0; // terminate with error
    }
    else {
	printf("file loaded successfull\r\n");
    }
    
    ofstream outfile (argv[2]);
    //ofstream velfile ("vel.out");
    stringstream part2 (stringstream::in | stringstream::out);
   
    if (!outfile.is_open()) {
	printf("could not open outfile");
	return 0;
    }

   char line[512]; //increase to get more chars !

   double pos[3];
   pos[0]=pos[1]=pos[2]=0;
   double velx,vely,velz,omx,omy,omz=0;

   double D=0;
   double mass=0;
   float xlo,xhi,ylo,yhi,zlo,zhi;
   float ts=0;
   int atoms=0;
   int l = sizeof(line);

   inFile.getline(line,sizeof(line)); //z1=ITEM: TIMESTEP
   inFile.getline(line,sizeof(line)); //z2=1000
	ts=atof(line);
	cout<<"timestep:"<<ts<<endl;
   //inFile.getline(line,sizeof(line)); //z3=ITEM: LABEL
   //inFile.getline(line,sizeof(line)); //z4=rem_0.6_0.2_0_1e-06
   inFile.getline(line,sizeof(line)); //z5=ITEM: NUMBER OF ATOMS
   inFile.getline(line,sizeof(line)); //z6={N}
	atoms=atoi(line);
	cout<<"Atoms:"<<atoms<<endl;
   inFile.getline(line,sizeof(line)); //z7
   inFile.getline(line,sizeof(line)); //z8

	for (int i = 0; i < l; i++)
    		{
        		if (line[i] == ' ')
        		{
            		line[i] = '\0';
            	}
    	}    
        char* startp = line;
	xlo=atof(startp);
	startp += strlen(startp) + 1;
	xhi=atof(startp);
   
	inFile.getline(line,sizeof(line)); //z7
	for (int i = 0; i < l; i++)
    		{
        		if (line[i] == ' ')
        		{
            		line[i] = '\0';
            	}
    	}    
        startp = line;
	ylo=atof(startp);
	startp += strlen(startp) + 1;
	yhi=atof(startp);
   
	inFile.getline(line,sizeof(line)); //z8
	for (int i = 0; i < l; i++)
    		{
        		if (line[i] == ' ')
        		{
            		line[i] = '\0';
            	}
    	}    
        startp = line;
	zlo=atof(startp);
	startp += strlen(startp) + 1;
	zhi=atof(startp);
	

   inFile.getline(line,sizeof(line)); //z9
	cout<<"Line:"<<line<<endl; //header
   double d=0;

   int lc=0;
   int ac=offset;

   outfile <<"LIGGGHTS data file from restart file: timestep = "<<ts<<", procs = 12\r\n\r\n";
   outfile << atoms <<" atoms\r\n\r\n";
   outfile << "2 atom types\r\n\r\n";
   outfile << "20 extra bond per atom\r\n\r\n";

   outfile << xlo << " " << xhi << " xlo xhi\r\n";
   outfile << ylo << " " << yhi << " ylo yhi\r\n";
   outfile << zlo << " " << zhi << " zlo zhi\r\n";
   outfile<<"\r\n";
   outfile<<"Atoms\r\n";
   outfile<<"\r\n";

   int type=0;
   while ( inFile.getline(line,sizeof(line)) ) {
	lc++;
	//cout<<"Line:"<<lc<<endl;
	for (int i = 0; i < l; i++)
    		{
        		if (line[i] == ' ')
        		{
            		line[i] = '\0';
            	}
    	}
    
    startp = line;

//id type x y z vx vy vz fx fy fz omegax omegay omegaz radius 

//id    
	int id = atoi(startp);
	id+=offset;
	startp += strlen(startp) + 1;
//type
	type = atof(startp);
	startp += strlen(startp) + 1;
	//startp += strlen(startp) + 1; //type is 2x in dumpfile (don't know why)
//position    
    for (int i = 0; i <= 2; i++)
    {
        d = atof(startp);
	pos[i]=d;
	startp += strlen(startp) + 1;
    }

//velocity
        d = atof(startp);
	velx=d;
        startp += strlen(startp) + 1;
        d = atof(startp);
	vely=d;
        startp += strlen(startp) + 1;
        d = atof(startp);
	velz=d;
        startp += strlen(startp) + 1;

//forces
	for (int i = 0; i < 3; i++)
    {
        startp += strlen(startp) + 1;
    }
//omega	
	d = atof(startp);
	omx=d;
        startp += strlen(startp) + 1;
	d = atof(startp);
	omy=d;
        startp += strlen(startp) + 1;
	d = atof(startp);
	omz=d;
        startp += strlen(startp) + 1;
//radius	
	D = 2*atof(startp);
	startp += strlen(startp) + 1;
//calc mass 
	mass+=(1./6.)*M_PI*D*D*D*dens; 
	//printf("D=%f,m=%f\n",D,mass);   
	ac++;
	velx=vely=velz=omx=omy=omz=0;
	//for hybrid pair styles: id type x y z dens rad
	outfile << ac << " " << type << " " << pos[0] << " " << pos[1] << " " << pos[2] << " " << D << " " << dens << " " << "0" << endl;
	//part2 << ac <<" "<<velx<<" "<<vely<<" "<<velz<<" "<<omx<<" "<<omy<<" "<<omz<< endl;
     
   } //endwhile

//outfile << endl << "Velocities" << endl<<endl;
//outfile << part2.str();


    inFile.close();
    outfile.close();
    //velfile.close();
    printf("particle mass=%f\n",mass);
    printf("done (%d)\n",ac);
    return 0;
  }

