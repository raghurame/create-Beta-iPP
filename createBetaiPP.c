#include<stdio.h>
#include<math.h>
#include<string.h>
#include<stdlib.h>
#include<unistd.h>

/*
This program will generate a datafile of beta crystal structure of isotactic polypropylene.

Structure: Trigonal
Symmetry: P3(1),
Symmetry operations:
	Chain1 : (x, y, z)
	Chain2 : (-y, x-y, z+(1/3))
	Chain3 : (-x+y, -x, z+(2/3))

angle(gamma) = 120 degrees

Source for symmetry operations: http://homepage.univie.ac.at/nikos.pinotsis/spacegroup.html

Approximations: United atom model is used
*/

int main(int argc, char const *argv[])
{
	if(argc > 1)
	{
		// Processing argc and argv
		if(strstr(argv[1], "--help"))
		{
            printf("\n[*] Input:\n    ~~~~~~\n\n\tNONE\n\n[*] Output:\n    ~~~~~~~\n\n\t\"betaipp.data\"\n\n[*] Use this program to create beta-phase iPP crystal, which is in LAMMPS data file format. Output from this program can be used as input in LAMMPS script with minimal modification.\n\n");
			exit(1);
		}
		else
		{
			printf("\nError:\n~~~~~~\n\n\tUnknown command passed. Type '--help' for more information or run the executable directly\n\n");
			exit(1);
		}
	}


	int natoms, nbonds, nangles, ndihedrals, nimpropers, natomtypes, nbondtypes, nangletypes, ndihedraltypes, nimpropertypes;
	int atomid=1, moltype;
	int xnumber, ynumber, znumber;
	int fr1, fr2, fr3;

	float xmult=0, ymult=0, zmult=0;
	float a=11.03, b=11.03, c=6.5;
	float xlo, xhi, ylo, yhi, zlo, zhi;

	FILE *output;
	output=fopen("betaipp.data", "w");

	// Real coordinates for the crystal structure (before applying symmetry operations)
	float c1x=0.2311*a, c2x=0.0823*a, c3x=0.0696*a, c4x=0.5426*a, c5x=0.4199*a, c6x=0.4169*a, c7x=0.8910*a, c8x=0.7410*a, c9x=0.7444*a;
	float c1y=0.1785*b, c2y=0.0772*b, c3y=0.0692*b, c4y=0.6813*b, c5y=0.6977*b, c6y=0.6991*b, c7y=0.4606*b, c8y=0.4004*b, c9y=0.4040*b;
	float c1z=0.5951*c, c2z=0.6740*c, c3z=0.9104*c, c4z=0.4169*c, c5z=0.4958*c, c6z=0.7328*c, c7z=0.6334*c, c8z=0.7168*c, c9z=0.9538*c;

	// Applying symmetry operations
	// Molecule 1
	float c1x_c1=c1x, c2x_c1=c2x, c3x_c1=c3x, c4x_c1=c4x, c5x_c1=c5x, c6x_c1=c6x, c7x_c1=c7x, c8x_c1=c8x, c9x_c1=c9x;
	float c1y_c1=c1y, c2y_c1=c2y, c3y_c1=c3y, c4y_c1=c4y, c5y_c1=c5y, c6y_c1=c6y, c7y_c1=c7y, c8y_c1=c8y, c9y_c1=c9y;
	float c1z_c1=c1z, c2z_c1=c2z, c3z_c1=c3z, c4z_c1=c4z, c5z_c1=c5z, c6z_c1=c6z, c7z_c1=c7z, c8z_c1=c8z, c9z_c1=c9z;

	// Molecule 2
	float c1x_c2=-c1y, c2x_c2=-c2y, c3x_c2=-c3y, c4x_c2=-c4y, c5x_c2=-c5y, c6x_c2=-c6y, c7x_c2=-c7y, c8x_c2=-c8y, c9x_c2=-c9y;
	float c1y_c2=c1x-c1y, c2y_c2=c2x-c2y, c3y_c2=c3x-c3y, c4y_c2=c4x-c4y, c5y_c2=c5x-c5y, c6y_c2=c6x-c6y, c7y_c2=c7x-c7y, c8y_c2=c8x-c8y, c9y_c2=c9x-c9y;
	float c1z_c2=c1z+(c/3), c2z_c2=c2z+(c/3), c3z_c2=c3z+(c/3), c4z_c2=c4z+(c/3), c5z_c2=c5z+(c/3), c6z_c2=c6z+(c/3), c7z_c2=c7z+(c/3), c8z_c2=c8z+(c/3), c9z_c2=c9z+(c/3);

	// Molecule 3
	float c1x_c3=-c1x+c1y, c2x_c3=-c2x+c2y, c3x_c3=-c3x+c3y, c4x_c3=-c4x+c4y, c5x_c3=-c5x+c5y, c6x_c3=-c6x+c6y, c7x_c3=-c7x+c7y, c8x_c3=-c8x+c8y, c9x_c3=-c9x+c9y;
	float c1y_c3=-c1x, c2y_c3=-c2x, c3y_c3=-c3x, c4y_c3=-c4x, c5y_c3=-c5x, c6y_c3=-c6x, c7y_c3=-c7x, c8y_c3=-c8x, c9y_c3=-c9x;
	float c1z_c3=c1z+(2*c/3), c2z_c3=c2z+(2*c/3), c3z_c3=c3z+(2*c/3), c4z_c3=c4z+(2*c/3), c5z_c3=c5z+(2*c/3), c6z_c3=c6z+(2*c/3), c7z_c3=c7z+(2*c/3), c8z_c3=c8z+(2*c/3), c9z_c3=c9z+(2*c/3);

	// Tilt angle for c axis (in degrees)
	float anglebeta=30; 
	// Converting degrees to radians
	anglebeta=anglebeta/57.2958;
	float sine, cosine;
	sine=sinf(anglebeta);
	cosine=cosf(anglebeta);

	// Gathering data from user
	printf("\nEnter the size of required crystal :-");
	printf("\n\txlo: "); scanf("%f", &xlo);
	printf("\txhi: "); scanf("%f", &xhi);
	printf("\n\tylo: "); scanf("%f", &ylo);
	printf("\tyhi: "); scanf("%f", &yhi);
	printf("\n\tzlo: "); scanf("%f", &zlo);
	printf("\tzhi: "); scanf("%f", &zhi);

	printf("\nComputing crystal size...");
	xnumber=(xhi-xlo)/a; /*ynumber=2;*/ ynumber=(yhi-ylo)/b; znumber=(zhi-zlo)/c;

	printf("\nCrystal size:\n\tX-Length: %.3f\n\tY-Length: %.3f\n\tZ-Length: %.3f", (xnumber*a), (ynumber*b), (znumber*c));
	printf("\nNumber of unit cells along:\n\tX-axis: %d\n\tY-axis: %d\n\tZ-axis: %d", xnumber, ynumber, znumber);

	// Computing natom, nbonds, nangles, ndihedrals, nimpropers.
	printf("\nCreating atoms...");
	natoms=xnumber*ynumber*znumber*27; natomtypes=3;
	printf("\nIgnoring bonds, angles, dihedrals and impropers...");
	nbonds=0; nangles=0; ndihedrals=0; nimpropers=0;
	nbondtypes=0; nangletypes=0; ndihedraltypes=0; nimpropertypes=0;
	printf("\nEnter molecule ID: "); scanf("%d", &moltype);
	printf("\nEnter the starting atom index: "); scanf("%d", &atomid);

	// Printing data file
	// Header
	fprintf(output, "Created by you v1.8.1 on today, this month, this year, current time.");

	// Data file information
	fprintf(output, "\n\n\t%d\tatoms\n\t%d\tbonds\n\t%d\tangles\n\t%d\tdihedrals\n\t%d\timpropers\n\n\t%d atom types\n\t%d bond types\n\t%d angle types\n\t%d dihedral types\n\t%d improper types\n\n\t%.2f\t%.2f\txlo xhi\n\t%.2f\t%.2f\tylo yhi\n\t%.2f\t%.2f\tzlo zhi\n\nMasses\n\n\t1\t13.0907	#CG311 CH\n\t2\t14.1707	#CG321 CH2\n\t3\t15.2507 #CG331 CH3\n\nAtoms\n\n", natoms, nbonds, nangles, ndihedrals, nimpropers, natomtypes, nbondtypes, nangletypes, ndihedraltypes, nimpropertypes, xlo, xhi, ylo, yhi, zlo, zhi);

	// Printing coordinates
	for(fr1=0; fr1<xnumber; fr1++)
	{
		for(fr2=0; fr2<ynumber; fr2++)
		{
			for(fr3=0; fr3<znumber; fr3++)
			{

				// Molecule 1

				fprintf(output, "\t%d\t%d\t1\t%.3f\t%.3f\t%.3f\n", atomid, moltype, (c1x_c1+(xmult*a))-((c1y_c1+(ymult*b))*sine), (c1y_c1+(ymult*b))*cosine, (c1z_c1+(zmult*c))); atomid++;
				fprintf(output, "\t%d\t%d\t1\t%.3f\t%.3f\t%.3f\n", atomid, moltype, (c2x_c1+(xmult*a))-((c2y_c1+(ymult*b))*sine), (c2y_c1+(ymult*b))*cosine, (c2z_c1+(zmult*c))); atomid++;
				fprintf(output, "\t%d\t%d\t1\t%.3f\t%.3f\t%.3f\n", atomid, moltype, (c3x_c1+(xmult*a))-((c3y_c1+(ymult*b))*sine), (c3y_c1+(ymult*b))*cosine, (c3z_c1+(zmult*c))); atomid++;
				fprintf(output, "\t%d\t%d\t1\t%.3f\t%.3f\t%.3f\n", atomid, moltype, (c4x_c1+(xmult*a))-((c4y_c1+(ymult*b))*sine), (c4y_c1+(ymult*b))*cosine, (c4z_c1+(zmult*c))); atomid++;
				fprintf(output, "\t%d\t%d\t1\t%.3f\t%.3f\t%.3f\n", atomid, moltype, (c5x_c1+(xmult*a))-((c5y_c1+(ymult*b))*sine), (c5y_c1+(ymult*b))*cosine, (c5z_c1+(zmult*c))); atomid++;
				fprintf(output, "\t%d\t%d\t1\t%.3f\t%.3f\t%.3f\n", atomid, moltype, (c6x_c1+(xmult*a))-((c6y_c1+(ymult*b))*sine), (c6y_c1+(ymult*b))*cosine, (c6z_c1+(zmult*c))); atomid++;
				fprintf(output, "\t%d\t%d\t1\t%.3f\t%.3f\t%.3f\n", atomid, moltype, (c7x_c1+(xmult*a))-((c7y_c1+(ymult*b))*sine), (c7y_c1+(ymult*b))*cosine, (c7z_c1+(zmult*c))); atomid++;
				fprintf(output, "\t%d\t%d\t1\t%.3f\t%.3f\t%.3f\n", atomid, moltype, (c8x_c1+(xmult*a))-((c8y_c1+(ymult*b))*sine), (c8y_c1+(ymult*b))*cosine, (c8z_c1+(zmult*c))); atomid++;
				fprintf(output, "\t%d\t%d\t1\t%.3f\t%.3f\t%.3f\n", atomid, moltype, (c9x_c1+(xmult*a))-((c9y_c1+(ymult*b))*sine), (c9y_c1+(ymult*b))*cosine, (c9z_c1+(zmult*c))); atomid++;

				// Molecule 2
				fprintf(output, "\t%d\t%d\t1\t%.3f\t%.3f\t%.3f\n", atomid, moltype, (c1x_c2+(xmult*a))-((c1y_c2+(ymult*b))*sine), (c1y_c2+(ymult*b))*cosine, (c1z_c2+(zmult*c))); atomid++;
				fprintf(output, "\t%d\t%d\t1\t%.3f\t%.3f\t%.3f\n", atomid, moltype, (c2x_c2+(xmult*a))-((c2y_c2+(ymult*b))*sine), (c2y_c2+(ymult*b))*cosine, (c2z_c2+(zmult*c))); atomid++;
				fprintf(output, "\t%d\t%d\t1\t%.3f\t%.3f\t%.3f\n", atomid, moltype, (c3x_c2+(xmult*a))-((c3y_c2+(ymult*b))*sine), (c3y_c2+(ymult*b))*cosine, (c3z_c2+(zmult*c))); atomid++;
				fprintf(output, "\t%d\t%d\t1\t%.3f\t%.3f\t%.3f\n", atomid, moltype, (c4x_c2+(xmult*a))-((c4y_c2+(ymult*b))*sine), (c4y_c2+(ymult*b))*cosine, (c4z_c2+(zmult*c))); atomid++;
				fprintf(output, "\t%d\t%d\t1\t%.3f\t%.3f\t%.3f\n", atomid, moltype, (c5x_c2+(xmult*a))-((c5y_c2+(ymult*b))*sine), (c5y_c2+(ymult*b))*cosine, (c5z_c2+(zmult*c))); atomid++;
				fprintf(output, "\t%d\t%d\t1\t%.3f\t%.3f\t%.3f\n", atomid, moltype, (c6x_c2+(xmult*a))-((c6y_c2+(ymult*b))*sine), (c6y_c2+(ymult*b))*cosine, (c6z_c2+(zmult*c))); atomid++;
				fprintf(output, "\t%d\t%d\t1\t%.3f\t%.3f\t%.3f\n", atomid, moltype, (c7x_c2+(xmult*a))-((c7y_c2+(ymult*b))*sine), (c7y_c2+(ymult*b))*cosine, (c7z_c2+(zmult*c))); atomid++;
				fprintf(output, "\t%d\t%d\t1\t%.3f\t%.3f\t%.3f\n", atomid, moltype, (c8x_c2+(xmult*a))-((c8y_c2+(ymult*b))*sine), (c8y_c2+(ymult*b))*cosine, (c8z_c2+(zmult*c))); atomid++;
				fprintf(output, "\t%d\t%d\t1\t%.3f\t%.3f\t%.3f\n", atomid, moltype, (c9x_c2+(xmult*a))-((c9y_c2+(ymult*b))*sine), (c9y_c2+(ymult*b))*cosine, (c9z_c2+(zmult*c))); atomid++;

				// Molecule 3
				fprintf(output, "\t%d\t%d\t1\t%.3f\t%.3f\t%.3f\n", atomid, moltype, (c1x_c3+(xmult*a))-((c1y_c3+(ymult*b))*sine), (c1y_c3+(ymult*b))*cosine, (c1z_c3+(zmult*c))); atomid++;
				fprintf(output, "\t%d\t%d\t1\t%.3f\t%.3f\t%.3f\n", atomid, moltype, (c2x_c3+(xmult*a))-((c2y_c3+(ymult*b))*sine), (c2y_c3+(ymult*b))*cosine, (c2z_c3+(zmult*c))); atomid++;
				fprintf(output, "\t%d\t%d\t1\t%.3f\t%.3f\t%.3f\n", atomid, moltype, (c3x_c3+(xmult*a))-((c3y_c3+(ymult*b))*sine), (c3y_c3+(ymult*b))*cosine, (c3z_c3+(zmult*c))); atomid++;
				fprintf(output, "\t%d\t%d\t1\t%.3f\t%.3f\t%.3f\n", atomid, moltype, (c4x_c3+(xmult*a))-((c4y_c3+(ymult*b))*sine), (c4y_c3+(ymult*b))*cosine, (c4z_c3+(zmult*c))); atomid++;
				fprintf(output, "\t%d\t%d\t1\t%.3f\t%.3f\t%.3f\n", atomid, moltype, (c5x_c3+(xmult*a))-((c5y_c3+(ymult*b))*sine), (c5y_c3+(ymult*b))*cosine, (c5z_c3+(zmult*c))); atomid++;
				fprintf(output, "\t%d\t%d\t1\t%.3f\t%.3f\t%.3f\n", atomid, moltype, (c6x_c3+(xmult*a))-((c6y_c3+(ymult*b))*sine), (c6y_c3+(ymult*b))*cosine, (c6z_c3+(zmult*c))); atomid++;
				fprintf(output, "\t%d\t%d\t1\t%.3f\t%.3f\t%.3f\n", atomid, moltype, (c7x_c3+(xmult*a))-((c7y_c3+(ymult*b))*sine), (c7y_c3+(ymult*b))*cosine, (c7z_c3+(zmult*c))); atomid++;
				fprintf(output, "\t%d\t%d\t1\t%.3f\t%.3f\t%.3f\n", atomid, moltype, (c8x_c3+(xmult*a))-((c8y_c3+(ymult*b))*sine), (c8y_c3+(ymult*b))*cosine, (c8z_c3+(zmult*c))); atomid++;
				fprintf(output, "\t%d\t%d\t1\t%.3f\t%.3f\t%.3f\n", atomid, moltype, (c9x_c3+(xmult*a))-((c9y_c3+(ymult*b))*sine), (c9y_c3+(ymult*b))*cosine, (c9z_c3+(zmult*c))); atomid++;

				// Incrementing Z-coordinates
				zmult++;

			}
		// reset Z-coordinates; increment Y-coordinates
		zmult=0; ymult++;
		}								
	// reset Z- and Y-coordinates; increment X-coordinates
	zmult=0; ymult=0; xmult++;
	}

	printf("\n");
	fclose(output);

	return(0);

}

