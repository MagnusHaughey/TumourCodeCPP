

# include <iostream>
# include <fstream>
# include <sstream>
# include <stdlib.h>
# include <cmath>
# include <math.h>
# include <random>
# include <cstdlib>
//# include <algorithm> 
//# include <vector>
//# include <sys/types.h>
//# include <sys/stat.h>
//# include <stdio.h>
//# include <unistd.h>

using namespace std;



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


};


// Define global variables
const double bdratio = 0.5;
const double _maxsize = 1e7;
const int _seed = 114;
const double _s = 0.1;
const double _ut = 2.0;
const double _ud = 0.2;
const double _ur = 1e-4;
const double r_death = bdratio * log(2.0);		// Death model 1
const int model = 3;
const int NEIGHBOURHOOD = 26;

double radius_double, t, max_birth, ran;
int radius, Ntot, iter, x, y, z, cell_index_x, cell_index_y, cell_index_z, r_birth, empty_neighbours, chosen_direction, queue, range, xrange, yrange, zrange;


// Define poisson distributions
default_random_engine generator;
poisson_distribution<int> poisson_d(_ud);
poisson_distribution<int> poisson_r(_ur);
poisson_distribution<int> poisson_t(_ut);


// Define cell division in volumetric growth model with "straight line" cell displacement 
void MODEL1_divide(Cell *** tumour , int cell_index_x , int cell_index_y , int cell_index_z , int *Ntot)
{

	// Randomly select the direction in which to divide
	do{
		ran = drand48();
		if (ran < (1.0/3.0)) x = -1;
		else if (ran < (2.0/3.0)) x = 0;
		else x = 1;

		ran = drand48();
		if (ran < (1.0/3.0)) y = -1;
		else if (ran < (2.0/3.0)) y = 0;
		else y = 1;

		ran = drand48();
		if (ran < (1.0/3.0)) z = -1;
		else if (ran < (2.0/3.0)) z = 0;
		else z = 1;
	}
	while ((x == 0) && (y == 0) && (z == 0));

	// Count how many cells need to be pushed in the specified direction (quantified by queue variable)
	queue = -1;

	do ++queue;
	while (tumour[cell_index_x + x*(queue+1)][cell_index_y + y*(queue+1)][cell_index_z + z*(queue+1)].dvr != -1 );

	// Shove cells outwards
	for (int j = 0; j < queue; j++)
	{
		tumour[cell_index_x + x*(queue - j + 1)][cell_index_y + y*(queue - j + 1)][cell_index_z + z*(queue - j + 1)] = tumour[cell_index_x + x*(queue - j)][cell_index_y + y*(queue - j)][cell_index_z + z*(queue - j)];
	}

	// Create daughter cell
	tumour[cell_index_x + x][cell_index_y + y][cell_index_z + z].setDVR(tumour[cell_index_x][cell_index_y][cell_index_z].dvr);
	tumour[cell_index_x + x][cell_index_y + y][cell_index_z + z].setRES(tumour[cell_index_x][cell_index_y][cell_index_z].res);
	tumour[cell_index_x + x][cell_index_y + y][cell_index_z + z].setPGR(tumour[cell_index_x][cell_index_y][cell_index_z].pgr);

	*Ntot += 1;

	// Add new GAs to daughter cells
	tumour[cell_index_x][cell_index_y][cell_index_z].dvr += poisson_d(generator);
	tumour[cell_index_x][cell_index_y][cell_index_z].res += poisson_r(generator);
	tumour[cell_index_x][cell_index_y][cell_index_z].pgr += poisson_t(generator);

	tumour[cell_index_x + x][cell_index_y + y][cell_index_z + z].dvr += poisson_d(generator);
	tumour[cell_index_x + x][cell_index_y + y][cell_index_z + z].res += poisson_r(generator);
	tumour[cell_index_x + x][cell_index_y + y][cell_index_z + z].pgr += poisson_t(generator);

}


// Define cell division in volumetric growth model, but where cells completely surrounded by other tumour cells cannot replicate (p_divide NOT proportional to number of empty neighbours)
void MODEL2_divide(Cell *** tumour , int x_nn[] , int y_nn[] , int z_nn[] , int cell_index_x , int cell_index_y , int cell_index_z  , int *Ntot , int empty_neighbours)
{

	chosen_direction = (int)(drand48()*((double)empty_neighbours));

	// Create daughter cell
	tumour[cell_index_x + x_nn[chosen_direction]][cell_index_y + y_nn[chosen_direction]][cell_index_z + z_nn[chosen_direction]].setDVR(tumour[cell_index_x][cell_index_y][cell_index_z].dvr);
	tumour[cell_index_x + x_nn[chosen_direction]][cell_index_y + y_nn[chosen_direction]][cell_index_z + z_nn[chosen_direction]].setRES(tumour[cell_index_x][cell_index_y][cell_index_z].res);
	tumour[cell_index_x + x_nn[chosen_direction]][cell_index_y + y_nn[chosen_direction]][cell_index_z + z_nn[chosen_direction]].setPGR(tumour[cell_index_x][cell_index_y][cell_index_z].pgr);
	
	*Ntot += 1;

	// Add new GAs to daughter cells
	tumour[cell_index_x][cell_index_y][cell_index_z].dvr += poisson_d(generator);
	tumour[cell_index_x][cell_index_y][cell_index_z].res += poisson_r(generator);
	tumour[cell_index_x][cell_index_y][cell_index_z].pgr += poisson_t(generator);

	tumour[cell_index_x + x_nn[chosen_direction]][cell_index_y + y_nn[chosen_direction]][cell_index_z + z_nn[chosen_direction]].dvr += poisson_d(generator);
	tumour[cell_index_x + x_nn[chosen_direction]][cell_index_y + y_nn[chosen_direction]][cell_index_z + z_nn[chosen_direction]].res += poisson_r(generator);
	tumour[cell_index_x + x_nn[chosen_direction]][cell_index_y + y_nn[chosen_direction]][cell_index_z + z_nn[chosen_direction]].pgr += poisson_t(generator);

}


// Define cell division in volumetric growth model, but where cells completely surrounded by other tumour cells cannot replicate (p_divide proportional to number of empty neighbours)
void MODEL3_divide(Cell *** tumour , int cell_index_x , int cell_index_y , int cell_index_z  , int *Ntot)
{


	// Randomly select the direction in which to divide
	do{
		ran = drand48();
		if (ran < (1.0/3.0)) x = -1;
		else if (ran < (2.0/3.0)) x = 0;
		else x = 1;

		ran = drand48();
		if (ran < (1.0/3.0)) y = -1;
		else if (ran < (2.0/3.0)) y = 0;
		else y = 1;

		ran = drand48();
		if (ran < (1.0/3.0)) z = -1;
		else if (ran < (2.0/3.0)) z = 0;
		else z = 1;
	}
	while ((x == 0) && (y == 0) && (z == 0));

	if (tumour[cell_index_x + x][cell_index_y + y][cell_index_z + z].dvr == -1)		// Check if neighbour is empty
	{

		// Create daughter cell
		tumour[cell_index_x + x][cell_index_y + y][cell_index_z + z].setDVR(tumour[cell_index_x][cell_index_y][cell_index_z].dvr);
		tumour[cell_index_x + x][cell_index_y + y][cell_index_z + z].setRES(tumour[cell_index_x][cell_index_y][cell_index_z].res);
		tumour[cell_index_x + x][cell_index_y + y][cell_index_z + z].setPGR(tumour[cell_index_x][cell_index_y][cell_index_z].pgr);

		*Ntot += 1;

		// Add new GAs to daughter cells
		tumour[cell_index_x][cell_index_y][cell_index_z].dvr += poisson_d(generator);
		tumour[cell_index_x][cell_index_y][cell_index_z].res += poisson_r(generator);
		tumour[cell_index_x][cell_index_y][cell_index_z].pgr += poisson_t(generator);

		tumour[cell_index_x + x][cell_index_y + y][cell_index_z + z].dvr += poisson_d(generator);
		tumour[cell_index_x + x][cell_index_y + y][cell_index_z + z].res += poisson_r(generator);
		tumour[cell_index_x + x][cell_index_y + y][cell_index_z + z].pgr += poisson_t(generator);

	}

}


void get_range( Cell *** tumour , int *range , int radius , int axis[3] )
{

	*range = 0;
	do
	{
		if ((tumour[radius + axis[0]*(*range)][radius + axis[1]*(*range)][radius + axis[2]*(*range)].dvr != -1) || (tumour[radius - axis[0]*(*range)][radius - axis[1]*(*range)][radius - axis[2]*(*range)].dvr != -1)) *range += 1;
	}
	while ( (tumour[radius + axis[0]*(*range)][radius + axis[1]*(*range)][radius + axis[2]*(*range)].dvr != -1) || (tumour[radius - axis[0]*(*range)][radius - axis[1]*(*range)][radius - axis[2]*(*range)].dvr != -1) );

	if ((radius - *range - 50) > 0) *range += 50;
	else *range = radius;

}	





int main(int argc, char const *argv[])
{
	

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
	ofstream NversusT_file;
	stringstream f;
	f << "./DATA/maxsize="<< _maxsize << "_seed=" << _seed << "_s=" << _s << "_ut=" << _ut << "_ud=" << _ud << "_ur=" << _ur << "_model=" << model << "_N(t).dat";
	NversusT_file.open(f.str().c_str());

	ofstream tumour_file;
	f.str("");
	f << "./DATA/maxsize="<< _maxsize << "_seed=" << _seed << "_s=" << _s << "_ut=" << _ut << "_ud=" << _ud << "_ur=" << _ur << "_model=" << model << ".dat";
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
	int axis[3] = {0,0,0};

	do
	{
		
		++iter;

		// Randomly select one cell to divide
		cell_index_x = 0;
		cell_index_y = 0;
		cell_index_z = 0;


		// Estimate current dimensions of tumour
		axis[0] = 1;
		axis[1] = 0;
		axis[2] = 0;
		get_range ( tumour , &xrange , radius , axis );

		axis[0] = 0;
		axis[1] = 1;
		get_range ( tumour , &yrange , radius , axis );

		axis[1] = 0;
		axis[2] = 1;
		get_range ( tumour , &zrange , radius , axis );

		do
		{

			cell_index_x = (int)((2*xrange)*drand48()) + radius - xrange;
			cell_index_y = (int)((2*yrange)*drand48()) + radius - yrange;
			cell_index_z = (int)((2*zrange)*drand48()) + radius - zrange;

			//cout << range << " " << radius << " " << cell_index_x << " " << cell_index_y << " " << cell_index_z << endl;

		}
		while (tumour[cell_index_x][cell_index_y][cell_index_z].dvr == -1);
		

		// Compute birth and death rate of cell (params[3] is the selective advantage of a single driver mutation)
		//r_birth = pow( (log(2.0)*(1.0 + _s)) , tumour[cell_index_x][cell_index_y].dvr );		// Birth model 1
		r_birth = pow( (log(2.0)*(1.0 + _s)) , (tumour[cell_index_x][cell_index_y][cell_index_z].dvr - tumour[cell_index_x][cell_index_y][cell_index_z].res) );		// Birth model 2	


		// Update maximal birth and death rate of all cells 
		if (r_birth > max_birth) max_birth = r_birth;


		// Cell divides with proability r_birth/max_birth
		if (drand48() < (r_birth/max_birth))
		{
			if (model == 1) MODEL1_divide(tumour , cell_index_x , cell_index_y , cell_index_z , &Ntot);
			else if (model == 2) 
			{
				empty_neighbours = 0;

				// Check for any neighbouring empty lattice sites
				for (int i = -1; i < 2; i++)
				{
					for (int j = -1; j < 2; j++)
					{
						for (int k = -1; k < 2; k++)
						{
							if (tumour[cell_index_x + i][cell_index_y + j][cell_index_z + k].dvr == -1) 	// if not occupied
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
					MODEL2_divide(tumour , x_nn , y_nn , z_nn , cell_index_x , cell_index_y , cell_index_z , &Ntot , empty_neighbours);
				}

			}
			else if (model == 3) MODEL3_divide(tumour , cell_index_x , cell_index_y , cell_index_z , &Ntot);
		}

		else if (drand48() < (r_death/max_birth))
		{
			// Delete cell from tumour
			tumour[cell_index_x][cell_index_y][cell_index_z].setDVR(-1);
			tumour[cell_index_x][cell_index_y][cell_index_z].setRES(-1);
			tumour[cell_index_x][cell_index_y][cell_index_z].setPGR(-1);

			// Size of tumour is reduced by 1
			Ntot -= 1;
		}
	
		// Progress time variable
		t += 1.0/(max_birth * Ntot);


		// Write total number of cells after regular number of iterations
		if (iter%5000 == 0)
		{

			NversusT_file << t << " " << Ntot << endl;
			cout << "Iter=" << iter << ", N=" << Ntot << endl;

			// Update max_birth variable if needs be
			//max_dvr = findmax(getfield.(tumour , :dvr))[1]			# getfield.(tumour, :dvr) returns array dvr value of all cells in tumour
			//global max_birth = log(2.0) * ((1.0 + params[3])^(max_dvr))

		}


	} while (Ntot < _maxsize);

	cout << " " << endl;

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


	NversusT_file.close();
	tumour_file.close();

	return 0;
}


















