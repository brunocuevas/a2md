#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

/*
density.c

	The purpose of this program is to read wfn files in an efficient way,
	and to provide the electronic density at a given set of points.


	FAQ : brunocuevaszuviria@gmail.com

*/

/*
	symetry_index : this matrix sets the types of symetries that
	can be found until d level.
	angstrom2au : angstrom to atomic units factor
	au2angstrom : atomic units to angstroms factor
*/

const double OVERLAP = 3.0;
const int MAX_ITERS= 10;
const double LOG_TOLERANCE = 1e-3;
const int symetry_index[20][3] = {
	{ 0 , 0 , 0 }, // s
	{ 1 , 0 , 0 }, // px
	{ 0 , 1 , 0 }, // py
	{ 0 , 0 , 1 }, // pz
	{ 2 , 0 , 0 }, // dxx
	{ 0 , 2 , 0 }, // dyy
	{ 0 , 0 , 2 }, // dzz
	{ 1 , 1 , 0 }, // dxy
	{ 1 , 0 , 1 }, // dxz
	{ 0 , 1 , 1 }, // dyz
	{ 3 , 0 , 0 },
    { 0 , 3 , 0 },
    { 0 , 0 , 3 },
    { 2 , 1 , 0 },
    { 2 , 0 , 1 },
    { 0 , 2 , 1 },
    { 1 , 2 , 0 },
    { 1 , 0 , 2 },
    { 0 , 1 , 2 },
    { 1 , 1 , 1 }
};

const char atom_symbols[16][2] = {
	{'H', ' '},
	{'H', 'e'},
	{'L', 'i'},
	{'B', 'e'},
	{'B', ' '},
	{'C', ' '},
	{'N', ' '},
	{'O', ' '},
	{'F', ' '},
	{'N', 'a'},
	{'M', 'g'},
	{'A', 'l'},
	{'S', 'i'},
	{'P', ' '},
	{'S', ' '}
};

const double au2angstrom = 0.529177249;
const double angstrom2au = 1.889725989;
const double pi = 3.14159265359;

int help(){
	printf("pyDensity, by bruno cuevas zuviria \n");
	printf("\tanalysis of electronic density\n");
	printf("usage : 1) density -i WFNFILE -c COORDINATEFILE\n" );
	printf("        2) density -i WFNFILE -s SPACING -r RESOLUTION\n" );
	printf("        3) density -i WFNFILE -s SPACING -r RESOLUTION -d DXFILE\n" );
	printf("        4) density -i WFNFILE -s SPACING -r RESOLUTION -d DXFILE -x XYZFILE\n" );
	printf("        5) density -i WFNFILE -u FILENAME -p DENSITY_VALUE -r RESOLUTION\n" );
	printf("        6) density -i WFNFILE -R FILENAME -s SPACING\n" );
	printf("example : 1) density -i h2o.wfn -c coordinates > data\n");
	printf("          2) density -i h2o.wfn -r 0.01 -s 1.0 > data\n");
	printf("          3) density -i h2o.wfn -r 0.01 -s 1.0 -l -d h2o.dx\n" );
	printf("          4) density -i h2o.wfn -r 0.01 -s 1.0 -l -d h2o.dx -x h2o.xyz\n");
	printf("          5) density -i h2o.wfn -u water_01.csv -p 0.01 -r 50\n" );
	printf("          6) density -i h2o.wfn -R h2o_random.csv -s 2.0\n" );
	printf("FAQ : bruno.czuviria@upm.es\n");
}

int lookForAtom(char *atomName) {
	int i ;
	int finish_flag = 0 ;
	char checkName[2];
	if (strlen(atomName) == 1) {
		strcat(atomName, " ");
	}

	for ( i = 0; i <= 16;  i ++) {

		checkName[0] = atom_symbols[i][0];
		checkName[1] = atom_symbols[i][1];

		if (atomName[0] == checkName[0] && atomName[1] == checkName[1]) {
			finish_flag = 1;
			break;
		}
	}
	if (finish_flag == 1) {
		return i;
	} else {
		return 0;
	}
}

int parsehead(char *line, int *mo, double *oc, double *energy){

	int read_terms = 0;
	int o = 0;
	double ocup = 0;
	double ene = 0;
	read_terms = sscanf(
		line,
		"MO  %d                     OCC NO =   %lf ORB. ENERGY = %lf",
		&o, &ocup, &ene
	);
	if (read_terms != 3 ){
		read_terms = sscanf(
			line,
			"MO    %d     MO 0.0        OCC NO =    %lf  ORB. ENERGY = %lf",
			&o, &ocup, &ene
		);
	}
	if (read_terms == 3 ) {
		*mo = o;
		*oc = ocup;
		*energy = ene;
		return 0;
	} else {
		return 1;
	}

}
int fortran2c(char *number) {
	/*
	The main purpose of this function is to avoid the problem of fortran
	double precission notation. It substitutes the "D" by "E" in those
	numbers so they can be read as floats or doubles.
	*/
	int iter = 0;
	for (iter = 0; iter < strlen(number); iter ++) {
		if (number[iter] == 'D'){
			number[iter] = 'E';
		}
	}
	return 0;
}
int write_xyz(char *filename, int *ans, double *ccx, double *ccy, double *ccz, int nuclei) {
	FILE *xyz;
	xyz = fopen(filename, "w");
	int i ;
	char name[2];
	if (xyz == NULL){
		printf("@ ERROR could not write to file %s\n", filename);
	}
	fprintf(xyz, "%d\n", nuclei);
	fprintf(xyz, "%s\n", filename);
	for (i = 0; i < nuclei ; i ++){
		name[0] = atom_symbols[ans[i]][0];
		name[1] = atom_symbols[ans[i]][1];
		fprintf(
			xyz, "%s        %16.5lf %16.5lf %16.5lf\n",
			name, ccx[i]*au2angstrom, ccy[i]*au2angstrom, ccz[i]*au2angstrom
		);
	}
	return fclose(xyz) ;

}
int write_dx(char *filename, int px, int py, int pz, double resolution, double *density,
	double minx, double miny, double minz
	){
	//
	//
	FILE *dx;
	int i = 0;
	dx = fopen(filename, "w");
	if (dx == NULL) {
		printf("@ ERROR could not write to file %s\n", filename);
	}
	fprintf(dx, "# Created by pyDensity - BrunoCuevas\n" );
	printf("@ dx file open\n" );
	fprintf(dx, "object 1 class gridpositions counts %d %d %d\n", px,py,pz);
	fprintf(
		dx, "origin %3.02f %3.02f %3.02f\n",
		minx* au2angstrom, miny*au2angstrom, minz*au2angstrom
	);
	fprintf(dx, "delta %4.03f 0.000 0.000\n", resolution * au2angstrom);
	fprintf(dx, "delta 0.000 %4.03f 0.000\n", resolution * au2angstrom);
	fprintf(dx, "delta 0.000 0.000 %4.03f\n", resolution * au2angstrom);
	fprintf(dx, "object 2 class gridconnections counts %d %d %d\n", px,py,pz);
	fprintf(dx, "object 3 class array type double rank 0 items %d data follows\n", px * py * pz);
	for (i = 1 ; i <= (px*py*pz) ; i ++){
		fprintf(dx, "%e ", density[i - 1]);
		if (i % 3 == 0) {
			fprintf(dx, "\n");
		}
	}
	fprintf(dx, "\nobject \"Dataset name\" class field\n");
	return fclose(dx);
}
float set_min(double *cx, double spacing, int nc){
	/*
	Gets the minimum value of an array and substracts a given
	quantity (spacing)
	*/
	double min = 10000.0 ;
	int i = 0;
	for (i = 0 ; i < nc; i++){
		if (cx[i] < min) min = cx[i];
	}
	return min - spacing;
}
float set_max(double *cx, double spacing, int nc){
	/*
	Gets the maximum value of an array and adds a given
	quantity (spacing)
	*/
	double max = -10000.0 ;
	int i = 0;
	for (i = 0 ; i < nc; i++){
		if (cx[i] > max) max = cx[i];
	}
	return max + spacing;
}
int set_stepsize(double min_x, double max_x, double resolution){
	/*
	Sets the number of steps that are necessary to move from minx to maxx
	making steps of size resolution.
	*/
	double range = max_x - min_x ;
	return (int) (range / resolution) + 1 ;
}
int set_box(double *cx, double *cy, double *cz, double **coordinates, int nc, double spacing, double resolution){
	/*
	Sets the calculation box to create the dx file.
	*/
	printf("@ setting box\n" );

	double min_x, max_x, min_y, max_y, min_z, max_z;
	double range_x, range_y, range_z;
	int points_x, points_y, points_z;
	double ac_x, ac_y , ac_z;
	min_x = set_min(cx, spacing, nc) ; max_x = set_max(cx, spacing, nc) ;
	min_y = set_min(cy, spacing, nc) ; max_y = set_max(cy, spacing, nc) ;
	min_z = set_min(cz, spacing, nc) ; max_z = set_max(cz, spacing, nc) ;
	int i, j, k, l;
	printf("@COORD\tMAX\tMIN\n");
	printf("@x\t%lf\t%lf\n", max_x, min_x);
	printf("@y\t%lf\t%lf\n", max_y, min_y);
	printf("@z\t%lf\t%lf\n", max_z, min_z);
	printf("@\n");
	points_x = set_stepsize(min_x, max_x, resolution);
	points_y = set_stepsize(min_y, max_y, resolution);
	points_z = set_stepsize(min_z, max_z, resolution);
	ac_x = min_x;
	ac_y = min_y;
	ac_z = min_z;
	l = 0;
	printf("@ size %d\n", points_x * points_y * points_z);
	for (i = 0; i < points_x; i ++) {
		for (j = 0; j < points_y; j ++) {
			for (k = 0; k < points_z; k ++) {
				l = (i * points_z * points_y) + (j * points_z) + k;
				//printf("i = %d j = %d k = %d l = %d\n", i,j,k,l);
				coordinates[l][0] = ac_x;
				coordinates[l][1] = ac_y;
				coordinates[l][2] = ac_z;
				ac_z += resolution ;
			}
			ac_z = min_z;
			ac_y += resolution ;
		}
		ac_x += resolution ;
		ac_y = min_y;
	}
	printf("@ set grid\n");
	return points_x * points_y * points_z ;

}
double get_density_matrix(double **coefficients, double *oc, double **density, int no, int np) {
	/*
	The density matrix can be calculated only once, saving time and resources.
	If the wfn contains n orbitals and p primitives, the density matrix would be
	a p*p matrix containing the terms obtained as:
		D_ij = \sum_{k}^{n} o_k * C_ik * C_jk
	where o_k is the occupancy of the orbital k and C_ik is the coefficient for
	primitive i at molecular k.
	*/
	printf("@ orbitals %d primitives %d\n", no, np);
	int oI; // orbital iterator
	int pI1, pI2; // primtive iterators (2 levels)
	double checksum = 0.0f;
	int line_iter = 0;
	int word_iter = 0;
	for (oI = 0 ;  oI < no ; oI ++) {

		for (pI1 = 0 ; pI1 < np ; pI1 ++) for (pI2 = 0; pI2 < np ; pI2 ++ ) {
			density[pI1][pI2] += oc[oI] * coefficients[oI][pI1] * coefficients[oI][pI2];
			checksum += density[pI1][pI2];
		}
	}

	/*for (pI1 = 0 ; pI1 < np ; pI1 ++) for (pI2 = 0; pI2 < np ; pI2 ++ ) {
		if (pI2 % 5 == 0){
			printf("\n@ %d %e ", pI1, density[pI1][pI2]);
		} else {
			printf("%e ", density[pI1][pI2]);
		}
	}*/
	return checksum;
}
double gaussian(double x, double y, double z, int sym, double exponent){
	/*
	Gets the value of a gaussian basis funcion at a given point, regardless
	its normalization constant.
	*/
	double r2 = (x*x) + (y*y) + (z*z);
	double L  = pow(x, (double )symetry_index[sym][0]);
	L *= pow(y, (double )symetry_index[sym][1]);
	L *= pow(z, (double )symetry_index[sym][2]);
	return L*exp(-r2  * exponent);
}
double calculate_density(double **density_matrix, double *cx, double *cy, double *cz, int *sym, int *css, double *exponents, double *oc,
	double x, double y, double z, int np) {
	/*
	It applies the definition of electronic density, as \rho = \sum D\chi_{i}(r)\chi_{j}(r)
	*/
	int pI1;
	int pI2; // primitive iterators
	double density_output = 0.0; // density output
	int center; // symetry indexes, centr index
	double x_dif, y_dif, z_dif; // distance variables

	double *basis_functions = malloc(np * sizeof(double));


	/* ITERATION OVER PRIMITIVES */
	for ( pI1 = 0 ; pI1 < np ; pI1 ++) {
		center = css[pI1] ;
		x_dif  = x - cx[center] ;
		y_dif  = y - cy[center] ;
		z_dif  = z - cz[center] ;
		basis_functions[pI1] = gaussian(x_dif, y_dif, z_dif, sym[pI1], exponents[pI1]);
	}

	/* ITERATION OVER PRIMITIVES */
	for (pI1 = 0; pI1 < np ; pI1 ++ ){
		for (pI2 = 0; pI2 < np; pI2 ++){
			density_output += density_matrix[pI1][pI2]*basis_functions[pI1]*basis_functions[pI2];
		}
	}
	return density_output;
}

int get_bins(int **angular_bins, int theta_res, int phi_res) {
	int theta_it = 0;
	int phi_it = 0;
	int total_bins = 0;
	int bins;
	double theta_angle ;
	double phi_angle ;
	*angular_bins = (int *) malloc((theta_res) * sizeof(int));
	(*angular_bins)[0] = 1;
	(*angular_bins)[theta_res-1] = 1;
	for (theta_it = 1 ; theta_it < theta_res -1 ; theta_it ++) {
		theta_angle = pi * (double ) theta_it / (double ) theta_res;
		bins = (int ) (sin(theta_angle) * phi_res) + 1 ;
		(*angular_bins)[theta_it] = bins;
		total_bins += bins;
	}
	return total_bins + 2;
}


int get_radial_vectors(double ***v, int theta_res, int phi_res) {
	int j; int k;
	int *angular_bins ;
	int total_bins ;
	double theta_angle ;
	double phi_angle ;
	int theta_it = 0;
	int phi_it = 0;
	total_bins = get_bins(&angular_bins, theta_res, phi_res);

	*v = malloc((total_bins) * 3 * sizeof(double));
	for (j = 0; j < total_bins; j ++) {
		(*v)[j] = malloc(sizeof(double)*3);
	}
	theta_it = 0 ;
	k = 0;
	for (theta_it = 0 ; theta_it < theta_res ; theta_it ++ ){
		theta_angle = pi * (double ) theta_it / (double ) theta_res;
		if (angular_bins[theta_it] == 1) {
			(*v)[k][0] = 0.0;
			(*v)[k][1] = 0.0;
			(*v)[k][2] = cos(theta_angle);
			k ++;
		} else {
			for (phi_it = 0; phi_it < angular_bins[theta_it]; phi_it ++) {
				phi_angle = 2 * pi * (double) phi_it / (double) angular_bins[theta_it];
				(*v)[k][0] = cos(phi_angle) * sin(theta_angle);
				(*v)[k][1] = sin(phi_angle) * sin(theta_angle);
				(*v)[k][2] = cos(theta_angle);
				k ++;
			}
		}
	}
	return k;
}


int intermediate_center(
	double *buffer, double *ccx, double *ccy, double *ccz, int *ans,
	int iter1, int iter2
) {
	double x[3];
	double y[3];
	double d;
	if ((ans[iter1] != 0) && (ans[iter2] != 0)) {
		x[0] = ccx[iter1]; x[1] = ccy[iter1]; x[2] = ccz[iter1];
		y[0] = ccx[iter2]; y[1] = ccy[iter2]; y[2] = ccz[iter2];
		d = sqrt(pow(x[0]-y[0], 2.0)+pow(x[0]-y[0], 2.0)+pow(x[0]-y[0], 2.0));
		if (d < 2.0) {
			buffer[0] = (x[0] + y[0])/2.0;
			buffer[1] = (x[1] + y[1])/2.0;
			buffer[2] = (x[2] + y[2])/2.0;
			return 1;
		} else {
			return 0;
		}
	} else {
		return 0;
	}

}


int find_isosurface(
	double target_dens, char *surf, int resolution, double **density_matrix,
	double *ccx, double *ccy, double *ccz, int *sym, int *css,
	double *exponents, double *occ, int *ans, int primitives, int nuclei
) {

	/*
	FIND ISOSURFACE ALGORITHM

	The purpose of this function is to obtain a list of positions that lie at
	a isodensity surface of the value of *target_dens*. For that purpose, the
	algorithm goes as follows:

		v <- []
		for theta in 0 to pi
			for phi in 0 to 2pi
				v <- [cartesians(1.0, theta, phi)]
		for each coordinates in nuclei
			r <- 0
			criterion <- 1000
			for each vector in v
				while criterion > tolerance
					v' <- v*r
					p <- obtain density(v')
					pl <- log(p)
					dp <- obtain derivatives of pl(v')
					r <- r - ((p-log(target_dens))/dp) (Newton Raphson)
					criterion <- abs(p - target_dens)


	*/

	// OBTAINING RADIAL VECTORS
	double **v ;
	int theta_res = (int) resolution ;
	int phi_res   = 2*((int) resolution) ;
	int total_terms    ;
	total_terms = get_radial_vectors(&v, theta_res, phi_res);
	// DECLARING ITERATION VARIABLES
	int nuclei_it;
	int iv ;  int iw ;
	int it; int flag;
	// DECLARING BUFFERS FOR NEWTON RAPHSON
	double x0[3]; // stores the center position
	double x[3][3]; // stores the coordinates upon which density and derivatives
		// will be calculated
	double p[3]; // stores the log(density)
	// DECLARING HEURISTIC PARAMETERS
	double criterion = 100.0; // defined as 100, will always be too large
		// so it will alow to start the iterations
	double r = 1e-3; // radius. Solution variable of the NR problem
	double h = 1e-3; // numerical diferentation step
	double der; // stores the numerical variable
	int max_iter; // stores the number of iterations
	double min_distance = 0.0f; // voronoi partition buffer
	double dist; // voronoi partition buffer
	double xx; double yy; double zz; // voronoi partition buffer


	// OPENING FILE
	FILE *surf_file;
	surf_file = fopen(surf, "w+");
	fprintf(surf_file, "x,y,z,p\n");

	// STARTING ITERATIONS THROUGH THE CENTERS
	nuclei_it = 0;
	it = 0;
	while (nuclei_it < nuclei) {
		if (ans[nuclei_it] == 0) {
			// in this case, the atom is an hydrogen, so it can be skipped.
			// hydrogen sampling relies in the heavy atoms to which is
			// bound
			it = 0;
			nuclei_it ++;
			continue;
		}

		if (it <= nuclei_it){
			flag = intermediate_center(
				x0, ccx, ccy, ccz, ans,
				nuclei_it, it
			);
			if (flag == 0) {
				it ++;
				continue;
			}

		}


		for (iv = 0 ; iv < total_terms ; iv ++) {
			max_iter = 0; criterion = 100.0;

			r = 1e-3;
			while (max_iter < MAX_ITERS && criterion > LOG_TOLERANCE) {
				if (max_iter > 0){
					r -= (p[1]-log(target_dens))/der;
				}
				// .........................................................
				// using the finite diferences diferentiation with a central
				// element:
				//		f'(x) =  (f(x+h) - f(x-h)) / (2h)
				// .........................................................
				// center element
				x[1][0] = x0[0] + r * v[iv][0];
				x[1][1] = x0[1] + r * v[iv][1];
				x[1][2] = x0[2] + r * v[iv][2];
				// upper element
				x[0][0] = x0[0] + ((r-h) * v[iv][0]);
				x[0][1] = x0[1] + ((r-h) * v[iv][1]);
				x[0][2] = x0[2] + ((r-h) * v[iv][2]);
				// lower element
				x[2][0] = x0[0] + ((r+h) * v[iv][0]);
				x[2][1] = x0[1] + ((r+h) * v[iv][1]);
				x[2][2] = x0[2] + ((r+h) * v[iv][2]);
				// log-density calculations
				p[0] = log(calculate_density(
					density_matrix, ccx, ccy, ccz, sym, css, exponents, occ,
					x[0][0], x[0][1], x[0][2],
					primitives
				));
				p[1] = log(calculate_density(
					density_matrix, ccx, ccy, ccz, sym, css, exponents, occ,
					x[1][0], x[1][1], x[1][2],
					primitives
				));
				p[2] = log(calculate_density(
					density_matrix, ccx, ccy, ccz, sym, css, exponents, occ,
					x[2][0], x[2][1], x[2][2],
					primitives
				));
				// criterion update. it should approach zero.
				criterion = fabs(p[1] - log(target_dens));
				// numerical diferentiation
				der = (p[2] - p[0])/(2*h);
				max_iter ++ ;
			}
			// VORONOI CRITERION. sampled points should not locate at the
			// voronoi partition of other centers
			if (max_iter < MAX_ITERS) {
				min_distance = 1e3;
				for (iw = 0 ; iw < nuclei ; iw ++ ){
					if (iw != nuclei_it && ans[iw] != 0.0) {
						xx = x[1][0] - ccx[iw];
						yy = x[1][1] - ccy[iw];
						zz = x[1][2] - ccz[iw];
						dist =  sqrt(xx*xx + yy*yy + zz*zz);
						if (dist < min_distance) {
							min_distance = dist;
						}
					}
				}
				// OVERLAP allows to relax the vornoi condition
				if (r < (OVERLAP * min_distance)) {
					fprintf( surf_file,
						"%8.4e,%8.4e,%8.4e,%8.7e\n",
						x[1][0], x[1][1], x[1][2], exp(p[1])
					);
				}
			}

		}
		it ++;
		if (it > nuclei_it) {
			it = 0;
			printf("@ iteration nuclei %d, finished \n", nuclei_it);
			nuclei_it ++;

		}

	}

	fclose(surf_file);
}


int read_wfn(char *filename, double **ccx, double **ccy, double **ccz, double **ccc,
			int **sym, int **css, double **exponents, double ***coef, double ***density_matrix,
			int *primitives, int *nuclei, int *orbitals, double **occ, int **ans, char *buffer,
			char *buffer_pointer
		){
	/*
	It reads wfn files and provides the data corresponding to coordinates,
	symetries, charges, centers, exponents, coefficients, and sets the
	memory location of the density matrix
	*/
	FILE *fp;
	//char buffer[1000];
	char *term;
	//char *buffer_pointer;
	printf("@reading\n" );
	int gauss_number = 0;
	int iter = 0; int coords_iter = 0; int count = 0; // Iteration through file
	char element[2]; // Needed to parse the atoms from coordiantes. Not used later.
	int atom_number, atom_centre, centre_lines, exponent_lines;
	int i = 0;
	int j = 0;
	int o = 0; // Used to iterate in orbitals
	int atomNumber;
	char *c;
	double x, y, z;// Used to asign coordiantes
	double ocup;   // Used to asign occupancies of the orbitals
	double ene ;   // Used to asign energy to each orbital
	double charge; // Used to asign charge of each atom
	int  read_control = 0; // Used to stop reading
	//nuclei = 100;
	fp = fopen(filename, "r");
	while (fgets(buffer, 255, (FILE*)fp) != NULL && read_control == 0) {
		buffer_pointer = strdup(buffer);
		if (iter == 0) {
		}
		else if (iter == 1)
		{
			/*
			This line allows  to declare all the following variables
			*/
			sscanf(buffer_pointer,
				"GAUSSIAN             %d MOL ORBITALS    %d PRIMITIVES       %d NUCLEI ",
				orbitals, primitives, nuclei
			);
			// Simple Arrays
			*ccx = (double *) malloc((*nuclei) * sizeof(double));
			*ccy = (double *) malloc((*nuclei) * sizeof(double));
			*ccz = (double *) malloc((*nuclei) * sizeof(double));
			*ccc = (double *) malloc((*nuclei) * sizeof(double));
			*ans = (int   *) malloc((*nuclei)* sizeof(int));
			*sym = (int   *) malloc((*primitives) * sizeof(int));
			*css = (int   *) malloc((*primitives) * sizeof(int));
			*exponents = (double *) malloc((*primitives) * sizeof(double));
			*occ       = (double *) malloc((*orbitals) * sizeof(double));
			// Double arrays
			// coef does not need to set values to 0.0
			*coef =  malloc((*orbitals) * (*primitives) * sizeof(double));
			for (j = 0; j < (*orbitals); j ++) (*coef)[j] = malloc(sizeof(double)*(*primitives));
			// density matrix does
			*density_matrix = malloc((*primitives) * (*primitives) * sizeof(double));

			for (i = 0 ; i < (*primitives) ; i ++)
				(*density_matrix)[i] = malloc((*primitives) * sizeof(double));

			for (i = 0 ; i < (*primitives) ; i ++) for (j = 0; j < (*primitives) ; j ++)
				(*density_matrix)[i][j] = 0.0;

			centre_lines = ((*primitives) / 20) ;
			exponent_lines = ((*primitives) / 5) ;
			if (*primitives % 20 > 0)  centre_lines += 1 ;
			if (*primitives % 5 > 0) exponent_lines += 1 ;
			printf("@ wfn file : %d primitives %d orbitals\n", *primitives, *orbitals);
		}
		else if (iter >= 2 && iter < (*nuclei + 2))
		{
			sscanf(buffer_pointer,
				"  %s    %d    (CENTRE  %d)  %lf %lf %lf  CHARGE =  %f",
				element, &atom_number, &atom_centre, &x, &y, &z, &charge
			);
			atomNumber = lookForAtom(element);
			(*ans)[coords_iter] = atomNumber ;
			(*ccx)[coords_iter] =  x;
			(*ccy)[coords_iter] =  y;
			(*ccz)[coords_iter] =  z;
			(*ccc)[coords_iter]  = charge;
			coords_iter += 1;


		}
		else if ( (iter >= ((*nuclei) +2)) && (iter < ((*nuclei) + 2 + centre_lines)) ){
			if (iter == ((*nuclei) + 2)) i = 0;
			while (c = strsep(&buffer_pointer, " ")) {
				if (( strcmp(c, "CENTRE") != 0 && strcmp(c, "ASSIGNMENTS") != 0)) {
					if (atoi(c) != 0) {
						(*css)[i] = atoi(c) - 1;
						i += 1;
					}
				}
			}
		}
		else if ((iter >= ((*nuclei) + 2 + centre_lines)) && (iter < ((*nuclei) + 2 + 2*centre_lines))){
			if (iter == (*nuclei) + 2 + centre_lines) {
				i = 0;
			}
			while (c = strsep(&buffer_pointer, " ")) {
				if ((strcmp(c, "TYPE") != 0) && (strcmp(c, "ASSIGNMENTS") != 0)) {
					if (atoi(c) != 0) {
						(*sym)[i] = atoi(c) - 1;
						i += 1;
					}
				}
			}
		}
		else if ((iter >= (*nuclei + 2 + 2*centre_lines)) &&  (iter < exponent_lines + *nuclei + 2 + 2*centre_lines)) {
			if (iter == *nuclei + 2 + 2*centre_lines) i = 0;
			while (c = strsep(&buffer_pointer, " ")) {
				if (strcmp(c , "EXPONENTS") != 0) {
					if (atof(c) != 0.0) {
						fortran2c(c);
						(*exponents)[i] = atof(c);
						i += 1;
					}
				}
			}
		} else {
			count = parsehead(buffer_pointer, &o, &ocup, &ene);
			(*occ)[o - 1] = ocup;
			if (count == 1) {
				while (c = strsep(&buffer_pointer, " ")) {
					if (strcmp(c, "END") == 0) {
						read_control = 1;
						break;
					} else {
						if (c[0] != '\0') {
						    if (o - 1 >= *orbitals) {
						        printf("@ WARNING maximum number of orbitals reached\n");
						        break;
						    }
							fortran2c(c);
							(*coef)[o -1][i] = atof(c);
							i += 1;
						}
					}
				}
			} else {
				i = 0;
			}
		}
		iter += 1;
	}
	fclose(fp);
	return 1;
}
int main(int argc, char *argv[]) {
	// Setting options
	if (argc == 1) {
		help();
		exit(1);
	}
	char *input_filename = NULL;
	char *input_coordinates = NULL;
	char *dxname = NULL;
	char *xyzname = NULL;
	char *surf = NULL;
	int random_sample = 1;
	double target_density = 0.0;
	double spacing = 0.0;
	double resolution = 0.0;
	int logscale = 0;
	int calculate_moments = 0;
	char dd ;
	while ((dd =getopt(argc, argv, "d:i:c:s:r:k:x:u:p:Rmhl")) != -1) {
		switch (dd) {
			case 'i':
				input_filename = optarg;
				printf("@ reading wfn from %s\n", input_filename );
				break;
			case 'c':
				input_coordinates = optarg;
				printf("@ reading coordinates from %s\n", input_coordinates);
				break;
			case 's':
				spacing = atof(optarg);
				printf("@ spacing set to %4.3f\n", spacing);
				break;
			case 'r':
				resolution = atof(optarg);
				printf("@ resolution set to %4.3f\n", resolution);
				break;
			case 'd':
				dxname = optarg;
				printf("@ writting to dx file %s \n", dxname);
				break;
			case 'x':
				xyzname = optarg;
				printf("@ writting coordinates to xyz\n", xyzname);
				break;
			case 'R':
				random_sample = 1;
				printf("@ random sampling\n");
				break;
			case 'h':
				help();
				exit(1);
			case 'l':
				logscale = 1;
				break;
			case 'u':
				surf = optarg;
				break;
			case 'p':
				target_density = atof(optarg);
				break;
			case 'm':
				calculate_moments = 1;
				break;
		}
	}
	if (!input_filename){
		printf("input is not defined\n");
		help();
		exit(1);
	}
	if (input_coordinates){
		printf("@ using coordinates file\n");
	} else {
		if (spacing || resolution){
			printf("@ using grid\n");
		} else {
			if (spacing && !resolution) {
				double resolution = 0.01;
				printf("@using default resolution, 0.01 au\n");
			}
			if (resolution && !spacing){
				double spacing = 1.0;
				printf("@using default spacing, 1.0 au\n");
			}
		}
	}
	/* density is just a test */
	FILE *gp;
	char buffer[1000];
	char *buffer_pointer;
	int orbitals, primitives, nuclei;
	double x, y, z;
	// I took the terrible decission of making the second wfn to have
	// variables with the same name than the originals but with a "d" in the
	// second position. May god forgive me?
	double *ccx, *ccy, *ccz, *ccc, *occ;
	int *sym, *css ;
	int *ans ; // AtomNames
	double **density_matrix, **coef;
	double *exponents ;

	int i = 0;
	int j = 0;

	read_wfn(
		input_filename, &ccx, &ccy, &ccz, &ccc, &sym, &css, &exponents, &coef,
		&density_matrix, &primitives, &nuclei, &orbitals, &occ, &ans, buffer, buffer_pointer
	);
	/*
	Reading is finished.
	Starting calculations
	*/
	double **coordinates_buffer; // Where coordinates are stored
	double min_x, min_y, min_z, max_x, max_y, max_z; // Necessary if a box is used
	double delta_x, delta_y, delta_z;
	double pointdens; // Where density value is soted
	double *dxdensity ;
	double integral;
	double moment[3];
	moment[0] = 0.0; moment[1] = 0.0; moment[2] = 0.0;
	int points_x, points_y, points_z; // Necesary if a box is used
	// Calculation of the density matrix takes place before starting
	printf("@ calculting density matrix\n");
	double checksum = get_density_matrix(coef, occ, density_matrix, orbitals, primitives);
	printf("@ done density matrix\n");
	printf("@ calculting density matrix\n");
	int target_coordinates ;
	printf("@ running \n" );

	double txch;

	if (input_coordinates){
		/*

		Reading coordinates from file, storing them in a buffer matrix,
		and calculate density for each of the points.

		*/
		printf("@ running \n" );
		coordinates_buffer = malloc(sizeof(double) * 30000 * 3);

		for (i = 0; i < 30000; i++) {
		    coordinates_buffer[i] =malloc(sizeof(double) * 3);
		}

		for (i = 0; i < 30000; i++) {

			coordinates_buffer[i][0] = 0.0 ;
			coordinates_buffer[i][1] = 0.0 ;
			coordinates_buffer[i][2] = 0.0 ;
		}
		i = 0;
		gp = fopen(input_coordinates, "r");
		while (fgets(buffer, 255, (FILE*)gp) != NULL) {
	        buffer_pointer = strdup(buffer);
			sscanf(buffer_pointer, "%lf %lf %lf", &x, &y, &z);
			coordinates_buffer[i][0] = x;
			coordinates_buffer[i][1] = y;
			coordinates_buffer[i][2] = z;
			i += 1;
			if (i > 30000) {
				printf("@ WARNING. Buffer is complete\n");
				break;
			}
		}
		fclose(gp);
		target_coordinates = i ;
		printf("@ reading done. %d points to calculate\n", target_coordinates);
		printf("@ checksum %f \n", checksum);

//	} else if (random_sample) {
//
//		printf("@ running \n" );
//		coordinates_buffer = malloc(sizeof(double) * 1000 * 3);
//
//		min_x = set_min(ccx, spacing, nuclei);
//		max_x = set_max(ccx, spacing, nuclei);
//		min_y = set_min(ccy, spacing, nuclei);
//		max_y = set_max(ccy, spacing, nuclei);
//		min_z = set_min(ccz, spacing, nuclei);
//		max_z = set_max(ccz, spacing, nuclei);
//
//		delta_x = max_x - min_x;
//		delta_y = max_y - min_y;
//		delta_z = max_z - min_z;
//
//		for (i = 0; i < 1000; i++) {
//
//			coordinates_buffer[i] =malloc(sizeof(double) * 3);
//
//			coordinates_buffer[i][0] = ((double)rand() / (double)RAND_MAX)*delta_x + min_x;
//			coordinates_buffer[i][1] = ((double)rand() / (double)RAND_MAX)*delta_y + min_y;
//			coordinates_buffer[i][2] = ((double)rand() / (double)RAND_MAX)*delta_z + min_z;
//
//			// printf(
//			// 	"%f %f %f\n",
//			// 	coordinates_buffer[i][0],
//			// 	coordinates_buffer[i][1],
//			// 	coordinates_buffer[i][2]
//			// );
//
//		}
//
//		target_coordinates = 1000;


	} else if (dxname) {
		/*

		Setting a coordinates box according to the size of the system

		*/
        printf("hola!");
		min_x = set_min(ccx, spacing, nuclei) ; max_x = set_max(ccx, spacing, nuclei) ;
		min_y = set_min(ccy, spacing, nuclei) ; max_y = set_max(ccy, spacing, nuclei) ;
		min_z = set_min(ccz, spacing, nuclei) ; max_z = set_max(ccz, spacing, nuclei) ;

		points_x = set_stepsize(min_x, max_x, resolution);
		points_y = set_stepsize(min_y, max_y, resolution);
		points_z = set_stepsize(min_z, max_z, resolution);
		target_coordinates = points_x * points_y * points_z;
		coordinates_buffer = malloc(points_x * points_y * points_z * 3 * sizeof(double));
		for (i = 0 ; i < target_coordinates; i++ ){
			coordinates_buffer[i] = malloc(3 * sizeof(double));
			coordinates_buffer[i][0] = 0.0;
			coordinates_buffer[i][1] = 0.0;
			coordinates_buffer[i][2] = 0.0;
		}
		set_box(ccx, ccy, ccz, coordinates_buffer, nuclei, spacing, resolution);
	}
	/*
		Writting output
	*/
	if (!dxname && !surf){

		for (i = 0 ; i < target_coordinates; i ++ ){
			x = coordinates_buffer[i][0];
			y = coordinates_buffer[i][1];
			z = coordinates_buffer[i][2];

            pointdens = calculate_density(
                density_matrix, ccx, ccy, ccz, sym, css, exponents, occ, x, y, z, primitives
            ) ;
            if (logscale) {
                pointdens = log10(pointdens + 1e-14);
            }
            printf("%18.6e %18.6e %18.6e %18.6e\n", x, y, z, pointdens);
		}
	} else if (dxname) {
		dxdensity = malloc(target_coordinates * sizeof(double));
		for (i = 0 ; i < target_coordinates; i ++ ){
			x = coordinates_buffer[i][0];
			y = coordinates_buffer[i][1];
			z = coordinates_buffer[i][2];
            dxdensity[i] = calculate_density(
                density_matrix, ccx, ccy, ccz, sym, css, exponents, occ, x, y, z, primitives
            );
			if (calculate_moments) {

				integral += resolution*resolution*resolution*dxdensity[i];
				moment[0] -= resolution*resolution*resolution*dxdensity[i]*x;
				moment[1] -= resolution*resolution*resolution*dxdensity[i]*y;
				moment[2] -= resolution*resolution*resolution*dxdensity[i]*z;

			}
            if (logscale) {
                dxdensity[i] = log10(dxdensity[i] + 1e-10);
			}
		}
		printf("@ calculations finished\n" );
		write_dx(dxname, points_x, points_y, points_z, resolution, dxdensity, min_x, min_y, min_z);
	}

	if (calculate_moments){
		for (i = 0 ; i < nuclei ; i ++){
			moment[0] += ccx[i] * (double) ans[i];
			moment[1] += ccy[i] * (double) ans[i];
			moment[2] += ccz[i] * (double) ans[i];
		}
		printf("@ integral is %8.4f\n", integral);
		printf(
			"@ dipole is %8.4f %8.4f %8.4f\n", moment[0], moment[1], moment[2]
		);
	}

	if (xyzname) {
		write_xyz(xyzname, ans, ccx, ccy, ccz, nuclei);
	}
	if (surf) {
		find_isosurface(target_density, surf, resolution, density_matrix,
		ccx, ccy, ccz, sym, css, exponents, occ, ans, primitives, nuclei);
	}

	/*free(dxdensity);
	free(coordinates_buffer);
	free(density_matrix);
	free(exponents)
	free(occ);
	free(ccx); free(ccy); free(ccz);
	free(primitives);*/
	return 0;

}
