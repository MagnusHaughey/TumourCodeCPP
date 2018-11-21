

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
const double _maxsize = 1e6;
const int _seed = 109;
const double _s = 0.05;
const double _ut = 1.0;
const double _ud = 0.2;
const double _ur = 1e-4;
const double r_death = bdratio * log(2.0);		// Death model 1
const int model = 3;
const int NEIGHBOURHOOD = 8;

double radius_double, t, max_birth, ran;
int radius, Ntot, iter, x, y, cell_index_x, cell_index_y, r_birth, empty_neighbours, chosen_direction, queue;


// Define poisson distributions
default_random_engine generator;
poisson_distribution<int> poisson_d(_ud);
poisson_distribution<int> poisson_r(_ur);
poisson_distribution<int> poisson_t(_ut);


// Define cell division in volumetric growth model with "straight line" cell displacement 
void MODEL1_divide(Cell ** tumour , int cell_index_x , int cell_index_y , int *Ntot)
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
	}
	while ((x == 0) && (y == 0));

	// Count how many cells need to be pushed in the specified direction (quantified by queue variable)
	queue = -1;

	do ++queue;
	while (tumour[cell_index_x + x*(queue+1)][cell_index_y + y*(queue+1)].dvr != -1 );

	// Shove cells outwards
	for (int j = 0; j < queue; j++)
	{
		tumour[cell_index_x + x*(queue - j + 1)][cell_index_y + y*(queue - j + 1)] = tumour[cell_index_x + x*(queue - j)][cell_index_y + y*(queue - j)];
	}

	// Create daughter cell
	tumour[cell_index_x + x][cell_index_y + y].setDVR(tumour[cell_index_x][cell_index_y].dvr);
	tumour[cell_index_x + x][cell_index_y + y].setRES(tumour[cell_index_x][cell_index_y].res);
	tumour[cell_index_x + x][cell_index_y + y].setPGR(tumour[cell_index_x][cell_index_y].pgr);

	*Ntot += 1;

	// Add new GAs to daughter cells
	tumour[cell_index_x][cell_index_y].dvr += poisson_d(generator);
	tumour[cell_index_x][cell_index_y].res += poisson_r(generator);
	tumour[cell_index_x][cell_index_y].pgr += poisson_t(generator);

	tumour[cell_index_x + x][cell_index_y + y].dvr += poisson_d(generator);
	tumour[cell_index_x + x][cell_index_y + y].res += poisson_r(generator);
	tumour[cell_index_x + x][cell_index_y + y].pgr += poisson_t(generator);

}


// Define cell division in volumetric growth model, but where cells completely surrounded by other tumour cells cannot replicate (p_divide NOT proportional to number of empty neighbours)
void MODEL2_divide(Cell ** tumour , int x_nn[] , int y_nn[] , int cell_index_x , int cell_index_y , int *Ntot , int empty_neighbours)
{

	chosen_direction = (int)(drand48()*((double)empty_neighbours));

	// Create daughter cell
	tumour[cell_index_x + x_nn[chosen_direction]][cell_index_y + y_nn[chosen_direction]].setDVR(tumour[cell_index_x][cell_index_y].dvr);
	tumour[cell_index_x + x_nn[chosen_direction]][cell_index_y + y_nn[chosen_direction]].setRES(tumour[cell_index_x][cell_index_y].res);
	tumour[cell_index_x + x_nn[chosen_direction]][cell_index_y + y_nn[chosen_direction]].setPGR(tumour[cell_index_x][cell_index_y].pgr);
	
	*Ntot += 1;

	// Add new GAs to daughter cells
	tumour[cell_index_x][cell_index_y].dvr += poisson_d(generator);
	tumour[cell_index_x][cell_index_y].res += poisson_r(generator);
	tumour[cell_index_x][cell_index_y].pgr += poisson_t(generator);

	tumour[cell_index_x + x_nn[chosen_direction]][cell_index_y + y_nn[chosen_direction]].dvr += poisson_d(generator);
	tumour[cell_index_x + x_nn[chosen_direction]][cell_index_y + y_nn[chosen_direction]].res += poisson_r(generator);
	tumour[cell_index_x + x_nn[chosen_direction]][cell_index_y + y_nn[chosen_direction]].pgr += poisson_t(generator);

}


// Define cell division in volumetric growth model, but where cells completely surrounded by other tumour cells cannot replicate (p_divide proportional to number of empty neighbours)
void MODEL3_divide(Cell ** tumour , int cell_index_x , int cell_index_y , int *Ntot)
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
	}
	while ((x == 0) && (y == 0));

	if (tumour[cell_index_x + x][cell_index_y + y].dvr == -1)		// Check if neighbour is empty
	{

		// Create daughter cell
		tumour[cell_index_x + x][cell_index_y + y].setDVR(tumour[cell_index_x][cell_index_y].dvr);
		tumour[cell_index_x + x][cell_index_y + y].setRES(tumour[cell_index_x][cell_index_y].res);
		tumour[cell_index_x + x][cell_index_y + y].setPGR(tumour[cell_index_x][cell_index_y].pgr);

		*Ntot += 1;

		// Add new GAs to daughter cells
		tumour[cell_index_x][cell_index_y].dvr += poisson_d(generator);
		tumour[cell_index_x][cell_index_y].res += poisson_r(generator);
		tumour[cell_index_x][cell_index_y].pgr += poisson_t(generator);

		tumour[cell_index_x + x][cell_index_y + y].dvr += poisson_d(generator);
		tumour[cell_index_x + x][cell_index_y + y].res += poisson_r(generator);
		tumour[cell_index_x + x][cell_index_y + y].pgr += poisson_t(generator);

	}

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
	radius_double = pow ( (_maxsize/M_PI) , 0.5 );
	radius = (int)(5.0*radius_double);

	Cell ** tumour = new Cell*[2*radius];
	for (int i = 0; i < (2*radius); i++)
	{
		// Declare cell pointers 
		tumour[i] = new Cell[2*radius];

		// Define cells as elements of tumour matrix
		for (int j = 0; j < (2*radius); j++)
		{
			tumour[i][j].setDVR(-1);
			tumour[i][j].setRES(-1);
			tumour[i][j].setPGR(-1);
		}
	}


	// Seed first tumour cell/gland at (x,y) = (0,0)
	tumour[radius][radius].setDVR(0);
	tumour[radius][radius].setRES(0);
	tumour[radius][radius].setPGR(0);

	Ntot += 1;

	cout << "Initialised tumour..." << endl;


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

	iter = 0;
	max_birth = 0.0;
	x = 0;
	y = 0;

	do
	{
		
		++iter;

		// Randomly select one cell to divide
		cell_index_x = 0;
		cell_index_y = 0;
		do
		{
			cell_index_x = (int)((2*radius)*drand48());
			cell_index_y = (int)((2*radius)*drand48());
		}
		while (tumour[cell_index_x][cell_index_y].dvr == -1);
		

		// Compute birth and death rate of cell (params[3] is the selective advantage of a single driver mutation)
		//r_birth = pow( (log(2.0)*(1.0 + _s)) , tumour[cell_index_x][cell_index_y].dvr );		// Birth model 1
		r_birth = pow( (log(2.0)*(1.0 + _s)) , (tumour[cell_index_x][cell_index_y].dvr - tumour[cell_index_x][cell_index_y].res) );		// Birth model 2	


		// Update maximal birth and death rate of all cells 
		if (r_birth > max_birth) max_birth = r_birth;


		// Cell divides with proability r_birth/max_birth
		if (drand48() < (r_birth/max_birth))
		{
			if (model == 1) MODEL1_divide(tumour , cell_index_x , cell_index_y , &Ntot);
			else if (model == 2) 
			{
				empty_neighbours = 0;

				// Check for any neighbouring empty lattice sites
				for (int i = -1; i < 2; i++)
				{
					for (int j = -1; j < 2; j++)
					{
						if (tumour[cell_index_x + i][cell_index_y + j].dvr == -1) 	// if not occupied
						{
							x_nn[empty_neighbours] = i; 	// Store coordinates of empty neighbour
							y_nn[empty_neighbours] = j;
							empty_neighbours += 1;
						}
					}
				}

				if (empty_neighbours != 0)
				{
					MODEL2_divide(tumour , x_nn , y_nn , cell_index_x , cell_index_y , &Ntot , empty_neighbours);
				}

			}
			else if (model == 3) MODEL3_divide(tumour , cell_index_x , cell_index_y , &Ntot);
		}

		else if (drand48() < (r_death/max_birth))
		{
			// Delete cell from tumour
			tumour[cell_index_x][cell_index_y].dvr = -1;
			tumour[cell_index_x][cell_index_y].res = -1;
			tumour[cell_index_x][cell_index_y].pgr = -1;

			// Size of tumour is reduced by 1
			Ntot -= 1;
		}
	
		// Progress time variable
		t += 1.0/(max_birth * Ntot);


		// Write total number of cells after regular number of iterations
		if (iter%200 == 0)
		{

			NversusT_file << t << " " << Ntot << endl;
			cout << "Iter=" << iter << ", N=" << Ntot << endl;

			// Update max_birth variable if needs be
			//max_dvr = findmax(getfield.(tumour , :dvr))[1]			# getfield.(tumour, :dvr) returns array dvr value of all cells in tumour
			//global max_birth = log(2.0) * ((1.0 + params[3])^(max_dvr))

		}


	} while (Ntot < _maxsize);


	// Write tumour data to file
	for (int i = 0; i < (2*radius); ++i)
	{
		for (int j = 0; j < (2*radius); ++j)
		{
			if (tumour[i][j].dvr != -1)
			{
				tumour_file << i << " " << j << " " << tumour[i][j].dvr << " " << tumour[i][j].res << " " << tumour[i][j].pgr << endl;
			}
		}
	}


	NversusT_file.close();
	tumour_file.close();

	return 0;
}


















