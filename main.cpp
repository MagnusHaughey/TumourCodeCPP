

# include <iostream>
# include <fstream>
# include <sstream>
# include <stdlib.h>
# include <cmath>
# include <math.h>
# include <random>
# include <cstdlib>
#include <dirent.h>
//# include <thread>
//# include <algorithm> 
//# include <vector>
//# include <sys/types.h>
//# include <sys/stat.h>
//# include <stdio.h>
//# include <unistd.h>

using namespace std;



/*******************************************************************************/



// Define global variables
const double _maxsize = 2.5e7;
const int _seed = 120;
const double _s = 0.5;
const double _ut = 2.0;
const double _ud = 0.2;
const double _ur = 1e-4;
const int div_model = 3;
const int adv_model = 3;

const double bdratio = 1.0/12.0;
const double r_surv = bdratio * log(2.0);		// Clonal survival rate
const double r_death = 0.5 * log(2.0);		// Intrinsic cell death rate
const int NEIGHBOURHOOD = 26;

double radius_double, t, max_birth, ran, r_birth;
int radius, Ntot, iter, x, y, z, cell_x, cell_y, cell_z, empty_neighbours, dir, queue, range, xrange, yrange, zrange, x_b, y_b, z_b;



/*******************************************************************************/



// Define poisson distributions
default_random_engine generator;
poisson_distribution<int> poisson_d(_ud);
poisson_distribution<int> poisson_r(_ur);
poisson_distribution<int> poisson_t(_ut);



/*******************************************************************************/



// Define a cell
class Cell
{

	public:
		int dvr;
		int res;
		int pgr;

	// Constructor for Cell object
	Cell(){}

	// Set() and get() methods
	void setDVR(int n)
	{
		this->dvr = n;
	}

	void setRES(int n)
	{
		this->res = n;
	}

	void setPGR(int n)
	{
		this->pgr = n;
	}

	void setGAs(int d , int r , int p)
	{
		this->dvr = d;
		this->res = r;
		this->pgr = p;
	}

	void newGAs()
	{
		this->dvr += poisson_d(generator);
		this->res += poisson_r(generator);
		this->pgr += poisson_t(generator);
	}


};



/*******************************************************************************/




// Define cell division in volumetric growth model with "straight line" cell displacement 
void MODEL1_divide(Cell *** tumour , int cell_x , int cell_y , int cell_z , int *Ntot , int *x_b , int *y_b , int *z_b , int radius)
{

	// Randomly select the direction in which to divide
	do
	{
		x = (int)(drand48()*3.0) - 1;
		y = (int)(drand48()*3.0) - 1;
		z = (int)(drand48()*3.0) - 1;
	}
	while ((x == 0) && (y == 0) && (z == 0));

	// Count how many cells need to be pushed in the specified direction (quantified by queue variable)
	queue = -1;

	do ++queue;
	while (tumour[cell_x + x*(queue+1)][cell_y + y*(queue+1)][cell_z + z*(queue+1)].dvr != -1 );

	// Shove cells outwards
	for (int j = 0; j < queue; j++)
	{
		tumour[cell_x + x*(queue - j + 1)][cell_y + y*(queue - j + 1)][cell_z + z*(queue - j + 1)] = tumour[cell_x + x*(queue - j)][cell_y + y*(queue - j)][cell_z + z*(queue - j)];
	}

	// Create daughter cell
	tumour[cell_x + x][cell_y + y][cell_z + z].setGAs(tumour[cell_x][cell_y][cell_z].dvr , tumour[cell_x][cell_y][cell_z].res , tumour[cell_x][cell_y][cell_z].pgr);

	*Ntot += 1;

	// Update bounds on tumour size
	if (fabs(cell_x + x*(queue + 1) - radius) > *x_b) *x_b = fabs(cell_x + x*(queue + 1) - radius);
	if (fabs(cell_y + y*(queue + 1) - radius) > *y_b) *y_b = fabs(cell_y + y*(queue + 1) - radius);
	if (fabs(cell_z + z*(queue + 1) - radius) > *z_b) *z_b = fabs(cell_z + z*(queue + 1) - radius);


	// Add new GAs to daughter cells
	tumour[cell_x][cell_y][cell_z].newGAs();
	tumour[cell_x + x][cell_y + y][cell_z + z].newGAs();

}


// Define cell division in volumetric growth model, but where cells completely surrounded by other tumour cells cannot replicate (p_divide NOT proportional to number of empty neighbours)
void MODEL2_divide(Cell *** tumour , int x_nn[] , int y_nn[] , int z_nn[] , int cell_x , int cell_y , int cell_z  , int *Ntot , int empty_neighbours , int *x_b , int *y_b , int *z_b , int radius)
{

	dir = (int)(drand48()*((double)empty_neighbours));

	// Create daughter cell
	tumour[cell_x + x_nn[dir]][cell_y + y_nn[dir]][cell_z + z_nn[dir]].setGAs(tumour[cell_x][cell_y][cell_z].dvr , tumour[cell_x][cell_y][cell_z].res , tumour[cell_x][cell_y][cell_z].pgr);
	
	*Ntot += 1;

	// Update bounds on tumour size
	if (fabs(cell_x + x_nn[dir] - radius) > *x_b) *x_b = fabs(cell_x + x_nn[dir] - radius);
	if (fabs(cell_x + y_nn[dir] - radius) > *y_b) *y_b = fabs(cell_y + y_nn[dir] - radius);
	if (fabs(cell_x + z_nn[dir] - radius) > *z_b) *z_b = fabs(cell_z + z_nn[dir] - radius);

	// Add new GAs to daughter cells
	tumour[cell_x][cell_y][cell_z].newGAs();
	tumour[cell_x + x_nn[dir]][cell_y + y_nn[dir]][cell_z + z_nn[dir]].newGAs();

}


// Define cell division in volumetric growth model, but where cells completely surrounded by other tumour cells cannot replicate (p_divide proportional to number of empty neighbours)
void MODEL3_divide(Cell *** tumour , int cell_x , int cell_y , int cell_z  , int *Ntot , int *x_b , int *y_b , int *z_b , int radius)
{


	// Randomly select the direction in which to divide
	do
	{
		x = (int)(drand48()*3.0) - 1;
		y = (int)(drand48()*3.0) - 1;
		z = (int)(drand48()*3.0) - 1;
	}
	while ((x == 0) && (y == 0) && (z == 0));

	if (tumour[cell_x + x][cell_y + y][cell_z + z].dvr == -1)		// Check if neighbour is empty
	{

		// Create daughter cell
		tumour[cell_x + x][cell_y + y][cell_z + z].setGAs(tumour[cell_x][cell_y][cell_z].dvr , tumour[cell_x][cell_y][cell_z].res , tumour[cell_x][cell_y][cell_z].pgr);

		*Ntot += 1;

		// Add new GAs to daughter cells
		tumour[cell_x][cell_y][cell_z].newGAs();
		tumour[cell_x + x][cell_y + y][cell_z + z].newGAs();

		// Update bounds on tumour size
		if (fabs(cell_x + x - radius) > *x_b) *x_b = fabs(cell_x + x - radius);
		if (fabs(cell_y + y - radius) > *y_b) *y_b = fabs(cell_y + y - radius);
		if (fabs(cell_z + z - radius) > *z_b) *z_b = fabs(cell_z + z - radius);
	}

}


void get_range( Cell *** tumour , int *range , int radius , int axis[3] )
{

	*range = 0;
	do
	{
		if ((tumour[radius + axis[0]*(*range)][radius + axis[1]*(*range)][radius + axis[2]*(*range)].dvr != -1)
			|| (tumour[radius - axis[0]*(*range)][radius - axis[1]*(*range)][radius - axis[2]*(*range)].dvr != -1)) *range += 1;
	}
	while ( (tumour[radius + axis[0]*(*range)][radius + axis[1]*(*range)][radius + axis[2]*(*range)].dvr != -1) 
		|| (tumour[radius - axis[0]*(*range)][radius - axis[1]*(*range)][radius - axis[2]*(*range)].dvr != -1) );

	if ((radius - *range - 50) > 0) *range += 50;
	else *range = radius;

}	

void compute_birthrate(Cell cell , int model_number , double *r_birth , double sel_adv)
{

	*r_birth = 0.0;

	if (model_number == 1) *r_birth = pow( (1.0 + sel_adv) , cell.dvr );

	else if (model_number == 2) *r_birth = pow( (1.0 + sel_adv) , (cell.dvr - cell.res) );

	else if (model_number == 3)
	{
		*r_birth = 1.0;
		for (int k = 1; k < cell.dvr + 1; k++)
		{
			(*r_birth) *= 1.0 + (sel_adv/(double)k);
		}
	}

	else if (model_number == 4)
	{
		if (cell.dvr <= 3) *r_birth = pow( (1.0 + sel_adv) , cell.dvr );
		else *r_birth = pow( (1.0 + sel_adv) , 3 ) * (pow( (1.0 + (sel_adv)/50.0) , (cell.dvr - 3) ));
	}

}	



/*******************************************************************************/




int main(int argc, char const *argv[])
{

	// Query number of available cores
	//unsigned concurentThreadsSupported = std::thread::hardware_concurrency();


	// Reset time and tumour size variables
	t = 0.0;
	Ntot = 0;


	// Seed random number generator
	srand48(_seed);



	//================== Initialise tumour ====================//

	// Estimate radius of resulting tumour using fitted parameters from previous simulations (slightly over-estimate)
	radius_double = pow ( (3.0*_maxsize/4.0*M_PI) , (1.0/3.0) );
	radius = (int)(1.5*radius_double);

	cout << " " << endl;

	Cell *** tumour = new Cell**[2*radius];
	for (int j = 0; j < (2*radius); j++)
	{
		tumour[j] = new Cell*[2*radius];

		for (int i = 0; i < (2*radius); i++)
		{
			// Declare cell pointers 
			tumour[j][i] = new Cell[2*radius];

			// Define cells as elements of tumour matrix
			for (int k = 0; k < (2*radius); k++)
			{
				tumour[j][i][k].setDVR(-1);
				tumour[j][i][k].setRES(-1);
				tumour[j][i][k].setPGR(-1);
			}
		}
		printf("Initialising tumour... %i%%\r", (int)((j+1)*100.0/(2*radius)));
		fflush(stdout);
	}

	cout << " " << endl;
		


	// Seed first tumour cell/gland at (x,y) = (0,0)
	tumour[radius][radius][radius].setDVR(0);
	tumour[radius][radius][radius].setRES(0);
	tumour[radius][radius][radius].setPGR(0);

	Ntot += 1;




	//================== Open data files ==================//

	DIR *dir = opendir("./DATA");
	if(dir)
	{
		system("mkdir ./DATA");
	}

	ofstream NversusT_file;
	stringstream f;
	f << "./DATA/maxsize="<< _maxsize << "_seed=" << _seed << "_s=" << _s << "_ut=" 
		<< _ut << "_ud=" << _ud << "_ur=" << _ur << "_divmodel=" << div_model << "_advmodel=" << adv_model << "_N(t).dat";
	NversusT_file.open(f.str().c_str());

	ofstream tumour_file;
	f.str("");
	f << "./DATA/maxsize="<< _maxsize << "_seed=" << _seed << "_s=" << _s << "_ut=" 
		<< _ut << "_ud=" << _ud << "_ur=" << _ur << "_divmodel=" << div_model << "_advmodel=" << adv_model << ".csv";
	tumour_file.open(f.str().c_str());

	cout << " " << endl;
	cout << "Created output files..." << endl;




	//================== Simulate tumour growth ==================//

	// Define arrays which will contain relative coordinates of empty neighbours for a chosen cell (for model 2)
	int x_nn[NEIGHBOURHOOD];
	int y_nn[NEIGHBOURHOOD];
	int z_nn[NEIGHBOURHOOD];

	iter = 0;
	max_birth = 0.0;
	x = 0;
	y = 0;
	z = 0;
	x_b = 0;
	y_b = 0;
	z_b = 0;

	do
	{
		
		++iter;

		// Randomly select one cell to divide
		cell_x = 0;
		cell_y = 0;
		cell_z = 0;

		do
		{

			cell_x = (int)((2*(x_b+1))*drand48()) + radius - x_b;
			cell_y = (int)((2*(y_b+1))*drand48()) + radius - y_b;
			cell_z = (int)((2*(z_b+1))*drand48()) + radius - z_b;

			//cout << x_b << " " << y_b << " " << z_b << " " << radius << " " << cell_x << " " << cell_y << " " << cell_z << endl;

		}
		while (tumour[cell_x][cell_y][cell_z].dvr == -1);
		

		// Compute birth rate of cell
		compute_birthrate( tumour[cell_x][cell_y][cell_z] , adv_model , &r_birth , _s );		


		// Update maximal birth and death rate of all cells 
		if (r_birth > max_birth) max_birth = r_birth;


		// Query the cell's position in the cell cycle (i.e. if it is ready to divide)
		if (drand48() < (r_birth/max_birth))
		{

			// First the cell is given the option to die
			if (drand48() < (r_surv/max_birth))
			{
				// Delete cell from tumour
				tumour[cell_x][cell_y][cell_z].setDVR(-1);
				tumour[cell_x][cell_y][cell_z].setRES(-1);
				tumour[cell_x][cell_y][cell_z].setPGR(-1);

				// Size of tumour is reduced by 1
				Ntot -= 1;
			}

			// Otherwise cell successfully divides
			else if (div_model == 1) MODEL1_divide(tumour , cell_x , cell_y , cell_z , &Ntot , &x_b , &y_b , &z_b , radius);
			else if (div_model == 2) 
			{
				empty_neighbours = 0;

				// Check for any neighbouring empty lattice sites
				for (int i = -1; i < 2; i++)
				{
					for (int j = -1; j < 2; j++)
					{
						for (int k = -1; k < 2; k++)
						{
							if (tumour[cell_x + i][cell_y + j][cell_z + k].dvr == -1) 	// if not occupied
							{
								x_nn[empty_neighbours] = i; 	// Store coordinates of empty neighbour
								y_nn[empty_neighbours] = j;
								z_nn[empty_neighbours] = k;
								empty_neighbours += 1;
							}
						}
						
					}
				}

				if (empty_neighbours != 0)
				{
					MODEL2_divide(tumour , x_nn , y_nn , z_nn , cell_x , cell_y , cell_z , &Ntot , empty_neighbours , &x_b , &y_b , &z_b , radius);
				}

			}
			else if (div_model == 3) MODEL3_divide(tumour , cell_x , cell_y , cell_z , &Ntot , &x_b , &y_b , &z_b , radius);
		}


		// Ask if the cell dies (intrinsic death rate)
		if (drand48() < (r_death/max_birth))
		{

			// Delete cell from tumour
			tumour[cell_x][cell_y][cell_z].setDVR(-1);
			tumour[cell_x][cell_y][cell_z].setRES(-1);
			tumour[cell_x][cell_y][cell_z].setPGR(-1);

			// Size of tumour is reduced by 1
			Ntot -= 1;
		}

	
		// Progress time variable
		t += 1.0/(max_birth * Ntot);


		// Write total number of cells after regular number of iterations
		if (iter%25000 == 0)
		{

			NversusT_file << t << " " << Ntot << endl;
			cout << "Iter=" << iter << ", N=" << Ntot << endl;

			// Update max_birth variable if needs be
			//max_dvr = findmax(getfield.(tumour , :dvr))[1]			# getfield.(tumour, :dvr) returns array dvr value of all cells in tumour
			//global max_birth = log(2.0) * ((1.0 + params[3])^(max_dvr))

		}


	} while (Ntot < _maxsize);

	cout << " " << endl;

/*
	// Write tumour data to file
	for (int i = 0; i < (2*radius); ++i)
	{
		for (int j = 0; j < (2*radius); ++j)
		{
			for (int k = 0; k < (2*radius); ++k)
			{
				if (tumour[i][j][k].dvr != -1)
				{
					tumour_file << i << " " << j << " " << k << " " << tumour[i][j][k].dvr << " " << tumour[i][j][k].res << " " << tumour[i][j][k].pgr << endl;
				}
			}	
		}
		printf("Writing data... %i%%\r", (int)((i+1)*100.0/(2*radius)));
		fflush(stdout);
	}
*/

	// Write tumour data to file
	tumour_file << "x coord, y coord, z coord, drivers, resistant, passengers" << endl;
	for (int i = (radius - x_b - 1); i < (radius + x_b + 2); ++i)
	{
		for (int j = (radius - y_b - 1); j < (radius + y_b + 2); ++j)
		{
			for (int k = (radius - z_b - 1); k < (radius + z_b + 2); ++k)
			{
				tumour_file << i << "," << j << "," << k << "," << tumour[i][j][k].dvr << "," << tumour[i][j][k].res << "," << tumour[i][j][k].pgr << endl;
			}	
		}
		printf("Writing data... %i%%\r", (int)((i+1)*100.0/(2*radius)));
		fflush(stdout);
	}

	NversusT_file.close();
	tumour_file.close();

	cout << "" << endl;
	cout << "Wrote ./" << f.str().c_str() << endl;
	cout << "" << endl;

	return 0;
}


















