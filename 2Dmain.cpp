

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
const double _maxsize = 2e6;
const int _seed = 75;
const double _s = 0.25;
const double _ut = 2.0;
const double _ud = 0.1;
//const double _ur = 1e-4;
const int arising_time = 2e5;
const int div_model = 3;
const int adv_model = 6;

const bool death = false;

const double bdratio = 1.0/12.0;
const double r_surv = bdratio * log(2.0);		// Clonal survival rate
const double r_death = 0.5 * log(2.0);		// Intrinsic cell death rate
const int NEIGHBOURHOOD = 8;

double radius_double, t, max_birth, ran, r_birth;
int radius, Ntot, iter, x, y, z, cell_x, cell_y, cell_z, empty_neighbours, dir, queue, range, xrange, yrange, zrange, x_b, y_b, z_b, res_mut, Nres;

bool mutated = false;



/*******************************************************************************/



// Define poisson distributions
default_random_engine generator;
poisson_distribution<int> poisson_d(_ud);
//poisson_distribution<int> poisson_r(_ur);
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
		//res_mut = poisson_r(generator);

		this->dvr += poisson_d(generator);
		//this->res += res_mut;
		this->pgr += poisson_t(generator);

		// If resistant mutation occurs, disallow any further resistant mutations by setting average of Poisson distribution to zero
		//if (res_mut != 0) poisson_distribution<int> poisson_r(0.0);
	}


};



/*******************************************************************************/




// Define cell division in volumetric growth model with "straight line" cell displacement 
void MODEL1_divide(Cell ** tumour , int cell_x , int cell_y , int *Ntot , int *x_b , int *y_b, int radius)
{

	// Randomly select the direction in which to divide
	do
	{
		x = (int)(drand48()*3.0) - 1;
		y = (int)(drand48()*3.0) - 1;
	}
	while ((x == 0) && (y == 0));

	// Count how many cells need to be pushed in the specified direction (quantified by queue variable)
	queue = -1;

	do ++queue;
	while (tumour[cell_x + x*(queue+1)][cell_y + y*(queue+1)].dvr != -1 );

	// Shove cells outwards
	for (int j = 0; j < queue; j++)
	{
		tumour[cell_x + x*(queue - j + 1)][cell_y + y*(queue - j + 1)] = tumour[cell_x + x*(queue - j)][cell_y + y*(queue - j)];
	}

	// Create daughter cell
	tumour[cell_x + x][cell_y + y].setGAs(tumour[cell_x][cell_y].dvr , tumour[cell_x][cell_y].res , tumour[cell_x][cell_y].pgr);

	*Ntot += 1;

	// Update bounds on tumour size
	if (fabs(cell_x + x*(queue + 1) - radius) > *x_b) *x_b = fabs(cell_x + x*(queue + 1) - radius);
	if (fabs(cell_y + y*(queue + 1) - radius) > *y_b) *y_b = fabs(cell_y + y*(queue + 1) - radius);


	// Add new GAs to daughter cells
	tumour[cell_x][cell_y].newGAs();
	tumour[cell_x + x][cell_y + y].newGAs();

}


// Define cell division in volumetric growth model, but where cells completely surrounded by other tumour cells cannot replicate (p_divide NOT proportional to number of empty neighbours)
void MODEL2_divide(Cell ** tumour , int x_nn[] , int y_nn[] , int cell_x , int cell_y , int *Ntot , int empty_neighbours , int *x_b , int *y_b , int radius)
{

	dir = (int)(drand48()*((double)empty_neighbours));

	// Create daughter cell
	tumour[cell_x + x_nn[dir]][cell_y + y_nn[dir]].setGAs(tumour[cell_x][cell_y].dvr , tumour[cell_x][cell_y].res , tumour[cell_x][cell_y].pgr);
	
	*Ntot += 1;

	// Update bounds on tumour size
	if (fabs(cell_x + x_nn[dir] - radius) > *x_b) *x_b = fabs(cell_x + x_nn[dir] - radius);
	if (fabs(cell_x + y_nn[dir] - radius) > *y_b) *y_b = fabs(cell_y + y_nn[dir] - radius);

	// Add new GAs to daughter cells
	tumour[cell_x][cell_y].newGAs();
	tumour[cell_x + x_nn[dir]][cell_y + y_nn[dir]].newGAs();

}


// Define cell division in volumetric growth model, but where cells completely surrounded by other tumour cells cannot replicate (p_divide proportional to number of empty neighbours)
void MODEL3_divide(Cell ** tumour , int cell_x , int cell_y , int *Ntot , int *Nres , int *x_b , int *y_b , int radius)
{


	// Randomly select the direction in which to divide
	do
	{
		x = (int)(drand48()*3.0) - 1;
		y = (int)(drand48()*3.0) - 1;
	}
	while ((x == 0) && (y == 0));

	if (tumour[cell_x + x][cell_y + y].dvr == -1)		// Check if neighbour is empty
	{

		if (tumour[cell_x][cell_y].res == 1) *Nres += 1;

		// Create daughter cell
		tumour[cell_x + x][cell_y + y].setGAs(tumour[cell_x][cell_y].dvr , tumour[cell_x][cell_y].res , tumour[cell_x][cell_y].pgr);

		*Ntot += 1;

		// Insert resistant cell at specified time (i.e. once tumour is comprised of a certain number of cells)
		if ((mutated == false) && (*Ntot == arising_time)) {
			cout << "Added mutation at N=" << arising_time << endl;
			tumour[cell_x + x][cell_y + y].res = 1; 
			*Nres += 1; 
			mutated = true;
		}

		// Add new GAs to daughter cells
		tumour[cell_x][cell_y].newGAs();
		tumour[cell_x + x][cell_y + y].newGAs();

		// Update bounds on tumour size
		if (fabs(cell_x + x - radius) > *x_b) *x_b = fabs(cell_x + x - radius);
		if (fabs(cell_y + y - radius) > *y_b) *y_b = fabs(cell_y + y - radius);
	}

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

	else if (model_number == 5)
	{
		*r_birth = 1.0;
		for (int k = 1; k < cell.res + 1; k++)
		{
			(*r_birth) *= 1.0 + (sel_adv/(double)k);
		}
	}
	else if (model_number == 6) *r_birth = pow( (1.0 + sel_adv) , cell.res );

}	



/*******************************************************************************/




int main(int argc, char const *argv[])
{

	// Query number of available cores
	//unsigned concurentThreadsSupported = std::thread::hardware_concurrency();


	// Reset time and tumour size variables
	t = 0.0;
	Ntot = 0;
	Nres = 0;


	// Seed random number generator
	srand48(_seed);



	//================== Initialise tumour ====================//

	// Estimate radius of resulting tumour using fitted parameters from previous simulations (slightly over-estimate)
	radius_double = pow ( (_maxsize/M_PI) , (1.0/2.0) );
	radius = (int)(1.5*radius_double);

	cout << " " << endl;

	Cell ** tumour = new Cell*[2*radius];
	for (int i = 0; i < (2*radius); i++)
	{
		tumour[i] = new Cell[2*radius];

		for (int j = 0; j < (2*radius); j++)
		{

			tumour[i][j].setDVR(-1);
			tumour[i][j].setRES(-1);
			tumour[i][j].setPGR(-1);
		
		}
		printf("Initialising tumour... %i%%\r", (int)((i+1)*100.0/(2*radius)));
		fflush(stdout);
	}

	cout << " " << endl;
		


	// Seed first tumour cell/gland at (x,y) = (0,0)
	tumour[radius][radius].setDVR(0);
	tumour[radius][radius].setRES(0);
	tumour[radius][radius].setPGR(0);

	Ntot += 1;




	//================== Open data files ==================//

	DIR *dir1 = opendir("./2D_DATA");
	if(!dir1)
	{
		system("mkdir ./2D_DATA");
	}

	stringstream f;
	f << "./2D_DATA/maxsize="<< _maxsize << "_seed=" << _seed << "_s=" << _s << "_ut=" 
		<< _ut << "_ud=" << _ud << "_arisingtime=" << arising_time << "_death=" << death << "_divmodel=" << div_model << "_advmodel=" << adv_model;
	DIR *dir2 = opendir(f.str().c_str());
	if(!dir2)
	{
		f.str("");
		f << "mkdir ./2D_DATA/maxsize="<< _maxsize << "_seed=" << _seed << "_s=" << _s << "_ut=" 
		<< _ut << "_ud=" << _ud << "_arisingtime=" << arising_time << "_death=" << death << "_divmodel=" << div_model << "_advmodel=" << adv_model;
		system(f.str().c_str());
	}

	ofstream NversusT_file;
	f.str("");
	f << "./2D_DATA/maxsize="<< _maxsize << "_seed=" << _seed << "_s=" << _s << "_ut=" 
		<< _ut << "_ud=" << _ud << "_arisingtime=" << arising_time << "_death=" << death << "_divmodel=" << div_model << "_advmodel=" << adv_model << "/N(t).dat";
	NversusT_file.open(f.str().c_str());

	ofstream tumour_file;
	f.str("");
	f << "./2D_DATA/maxsize="<< _maxsize << "_seed=" << _seed << "_s=" << _s << "_ut=" 
		<< _ut << "_ud=" << _ud << "_arisingtime=" << arising_time << "_death=" << death << "_divmodel=" << div_model << "_advmodel=" << adv_model << "/tumour.csv";
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
	x_b = 0;
	y_b = 0;

	do
	{
		
		++iter;


		// Randomly select one cell to divide
		cell_x = 0;
		cell_y = 0;

		do
		{

			cell_x = (int)((2*(x_b+1))*drand48()) + radius - x_b;
			cell_y = (int)((2*(y_b+1))*drand48()) + radius - y_b;

			//cout << x_b << " " << y_b << " " << z_b << " " << radius << " " << cell_x << " " << cell_y << " " << cell_z << endl;

		}
		while (tumour[cell_x][cell_y].dvr == -1);


		// Compute birth rate of cell
		compute_birthrate( tumour[cell_x][cell_y] , adv_model , &r_birth , _s );		


		// Update maximal birth and death rate of all cells 
		if (r_birth > max_birth) max_birth = r_birth;


		// Query the cell's position in the cell cycle (i.e. if it is ready to divide)
		if (drand48() < (r_birth/max_birth))
		{


			// First the cell is given the option to die
			if ( (death) && (drand48() < (r_surv/max_birth)) )
			{
				if (tumour[cell_x][cell_y].res == 1) Nres -= 1;

				// Delete cell from tumour
				tumour[cell_x][cell_y].setDVR(-1);
				tumour[cell_x][cell_y].setRES(-1);
				tumour[cell_x][cell_y].setPGR(-1);

				// Size of tumour is reduced by 1
				Ntot -= 1;

				// Progress time variable
				t += 1.0/(max_birth * Ntot);
			}

			// Otherwise cell successfully divides
			else if (div_model == 1) MODEL1_divide(tumour , cell_x , cell_y , &Ntot , &x_b , &y_b , radius);
			else if (div_model == 2) 
			{
				empty_neighbours = 0;

				// Check for any neighbouring empty lattice sites
				for (int i = -1; i < 2; i++)
				{
					for (int j = -1; j < 2; j++)
					{

						if (tumour[cell_x + i][cell_y + j].dvr == -1) 	// if not occupied
						{
							x_nn[empty_neighbours] = i; 	// Store coordinates of empty neighbour
							y_nn[empty_neighbours] = j;
							empty_neighbours += 1;
						}						
					}
				}

				if (empty_neighbours != 0)
				{
					MODEL2_divide(tumour , x_nn , y_nn , cell_x , cell_y , &Ntot , empty_neighbours , &x_b , &y_b , radius);
				}

			}
			else if (div_model == 3)
			{

				MODEL3_divide(tumour , cell_x , cell_y , &Ntot , &Nres , &x_b , &y_b , radius);
			}


			// Progress time variable
			t += 1.0/(max_birth * Ntot);
		}


		// Ask if the cell dies (intrinsic death rate)
		if ( (tumour[cell_x][cell_y].dvr != -1) && (death) && (drand48() < (r_death/max_birth)) )
		{
			if (tumour[cell_x][cell_y].res == 1) Nres -= 1;

			// Delete cell from tumour
			tumour[cell_x][cell_y].setDVR(-1);
			tumour[cell_x][cell_y].setRES(-1);
			tumour[cell_x][cell_y].setPGR(-1);

			// Size of tumour is reduced by 1
			Ntot -= 1;

			// Progress time variable
			t += 1.0/(max_birth * Ntot);
		}



		// Write total number of cells after regular number of iterations
		if (iter%25000 == 0)
		{

			NversusT_file << t << " " << Ntot << endl;
			cout << "Iter=" << iter << ", N=" << Ntot << ", ResGA=" << Nres << endl;

			// Update max_birth variable if needs be
			//max_dvr = findmax(getfield.(tumour , :dvr))[1]			# getfield.(tumour, :dvr) returns array dvr value of all cells in tumour
			//global max_birth = log(2.0) * ((1.0 + params[3])^(max_dvr))

		}


	} while (Ntot < _maxsize);

	cout << " " << endl;





	// Write tumour data to file
	int count = 0;
	tumour_file << "x coord, y coord, z coord, drivers, resistant, passengers" << endl;
	for (int i = (radius - x_b - 1); i < (radius + x_b + 2); ++i)
	{
		for (int j = (radius - y_b - 1); j < (radius + y_b + 2); ++j)
		{

			if ( tumour[i][j].dvr != -1 ) tumour_file << i << "," << j << ",0," << tumour[i][j].dvr << "," << tumour[i][j].res << "," << tumour[i][j].pgr << endl;
		
		}
		printf("Writing data... %i%%\r", (int)((i+1)*100.0/(2*radius)));
		fflush(stdout);
	}

	NversusT_file.close();
	tumour_file.close();

	cout << "" << endl;
	cout << "Wrote " << f.str().c_str() << endl;
	cout << "" << endl;

	return 0;
}


















