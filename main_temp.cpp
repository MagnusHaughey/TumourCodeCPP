

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
const int _seed = 104;
const double _s = 0.025;
const double _ut = 1.0;
const double _ud = 0.1;
const double _ur = 1e-4;
const double r_death = bdratio * log(2.0);		// Death model 1

double radius_double, t, max_birth, ran;
int radius, Ntot, iter, x, y, cell_index_x, cell_index_y, r_birth, range;


// Define poisson distributions
default_random_engine generator;
poisson_distribution<int> poisson_d(_ud);
poisson_distribution<int> poisson_r(_ur);
poisson_distribution<int> poisson_t(_ut);


// Method which returns random integer within given range i.e. [lower,upper] inclusive
int uniform_range(int lower, int upper)
{

	x = 0;
	ran = drand48();

	for (int i = lower; i < (upper+1); i++)
	{
		if (ran < ((double)(i - lower + 1)/(double)(upper - lower + 1)))
		{
			x = i;
			break;
		}
	}

	return x;

}


// Define cell division in volumetric growth model, but where cells completely surrounded by other tumour cells cannot replicate (p_divide NOT proportional to number of empty neighbours)
void MODEL3_divide(Cell ** tumour, int cell_index_x , int cell_index_y , int *Ntot){


	// Randomly select the direction in which to divide
	do{
		x = uniform_range(-1 , 1);
		y = uniform_range(-1 , 1);
	}
	while ((x == 0) && (y == 0));

	if (tumour[cell_index_x + x][cell_index_y + y].dvr == -1)		// Check if neighbour is empty
	{

		// Create daughter cell
		tumour[cell_index_x + x][cell_index_y + y].dvr = tumour[cell_index_x][cell_index_y].dvr;
		tumour[cell_index_x + x][cell_index_y + y].res = tumour[cell_index_x][cell_index_y].res;
		tumour[cell_index_x + x][cell_index_y + y].pgr = tumour[cell_index_x][cell_index_y].pgr;

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
			tumour[i][j].dvr = -1;
			tumour[i][j].res = -1;
			tumour[i][j].pgr = -1;
		}
	}


	// Seed first tumour cell/gland at (x,y) = (0,0)
	tumour[radius][radius].dvr = 0;
	tumour[radius][radius].res = 0;
	tumour[radius][radius].pgr = 0;

	Ntot += 1;

	cout << "Initialised tumour..." << endl;


	//================== Open data files ==================//
	ofstream NversusT_file;
	stringstream f;
	f << "./DATA/maxsize="<< _maxsize << "_seed=" << _seed << "_s=" << _s << "_ut=" << _ut << "_ud=" << _ud << "_ur=" << _ur << "_N(t).dat";
	NversusT_file.open(f.str().c_str());

	ofstream tumour_file;
	f.str("");
	f << "./DATA/maxsize="<< _maxsize << "_seed=" << _seed << "_s=" << _s << "_ut=" << _ut << "_ud=" << _ud << "_ur=" << _ur << ".dat";
	tumour_file.open(f.str().c_str());

	cout << " " << endl;
	cout << "Created output files..." << endl;



	//================== Simulate tumour growth ==================//

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

		range = 0;
		if ((radius - Ntot - 1) < 0) range = radius; // Compute range of search
		else range = Ntot + 1;

		do
		{ 
			cell_index_x = uniform_range( radius - range , radius + range );
			cell_index_y = uniform_range( radius - range , radius + range );
		}
		while (tumour[cell_index_x][cell_index_y].dvr == -1);
		

		// Compute birth and death rate of cell (params[3] is the selective advantage of a single driver mutation)
		r_birth = pow( (log(2.0)*(1.0 + _s)) , tumour[cell_index_x][cell_index_y].dvr );		// Birth model 1	


		// Update maximal birth and death rate of all cells 
		if (r_birth > max_birth) max_birth = r_birth;


		// Cell divides with proability r_birth/max_birth
		if (drand48() < (r_birth/max_birth))
		{
			MODEL3_divide(tumour , cell_index_x , cell_index_y , &Ntot);
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


















