

# include <iostream>
# include <fstream>
# include <sstream>
# include <stdlib.h>
# include <cmath>
# include <math.h>
# include <random>
# include <cstdlib>
# include <dirent.h>
# include <string>
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
const int _maxsize = 1e6;
const double _ut = 2.0;
const double _ud = 0.1;
const int div_model = 3;
const int adv_model = 6;

const bool death = true;

// Uncomment one of the following two lines to choose from triangular lattice (3), nearest-neighbour (4) and next-nearest-beighbour neighbourhoods (8) respectively
//const int NEIGHBOURHOOD = 3;
//const int NEIGHBOURHOOD = 4;
const int NEIGHBOURHOOD = 8;

double radius_double, t, _s, max_birth, ran, r_birth, current_diameter, chord;
int radius, Ntot, Ntagged, iter, x, y, z, cell_x, cell_y, cell_z, empty_neighbours, dir, queue, ind, length, coordX, coordY, boundary, attempts, counter;
int range, arising_time, xrange, yrange, zrange, x_b, y_b, z_b, res_mut, Nres, Nwt, direction, chosen_direction, min_length, num_mins, opposite_x, opposite_y, marker_x , marker_y;


const double bdratio = 1.0/12.0;
double r_surv = r_birth * bdratio;		// Clonal survival rate
//const double r_death = 0.5 * log(2.0);			// Intrinsic cell death rate
const double r_death = 0.0;							// Intrinsic cell death rate


// Define arrays which will contain relative coordinates of empty neighbours for a chosen cell (for model 2)
int x_nn[NEIGHBOURHOOD];
int y_nn[NEIGHBOURHOOD];
int chainX[(int)_maxsize];		// this is potentially lazy use of memory
int chainY[(int)_maxsize];
int directions[NEIGHBOURHOOD];
int previous_Nres[(int)(_maxsize)];
int previous_Nwt[(int)(_maxsize)];



bool mutated = false;
bool extended_chain = false;
bool divided = false;
bool quiet = true;
bool success = false;

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



/**********************   Only model 3 division used for initial 2D surface growth fractal dimension studies   ***********************/



/*
// Define cell division in volumetric growth model with "straight line" cell displacement 
void MODEL1_divide(Cell ** tumour , int cell_x , int cell_y , int *Ntot , int *Nres  , int *x_b , int *y_b, int radius)
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

	// Insert resistant cell at specified time (i.e. once tumour is comprised of a certain number of cells)
	if ((mutated == false) && (*Ntot == arising_time)) {
		//cout << "Added mutation at N=" << arising_time << endl;
		tumour[cell_x + x][cell_y + y].res = 1; 
		*Nres += 1; 
		mutated = true;
	}

	// Update bounds on tumour size
	if (fabs(cell_x + x*(queue + 1) - radius) > *x_b) *x_b = fabs(cell_x + x*(queue + 1) - radius);
	if (fabs(cell_y + y*(queue + 1) - radius) > *y_b) *y_b = fabs(cell_y + y*(queue + 1) - radius);


	// Add new GAs to daughter cells
	tumour[cell_x][cell_y].newGAs();
	tumour[cell_x + x][cell_y + y].newGAs();

	if (tumour[cell_x][cell_y].res == 1) *Nres += 1;

}


// Define cell division in surface growth model, where cell divides as long as >0 empty neighbours (p_divide NOT proportional to number of empty neighbours)
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
*/

// Define cell division in surface growth model, but where cells completely surrounded by other tumour cells cannot replicate (p_divide proportional to number of empty neighbours)
void MODEL3_divide(Cell ** tumour , int cell_x , int cell_y , int *Ntot , int *Nres, int *Nwt, int *Ntagged , int *x_b , int *y_b , int radius, int arising_time , int *marker_x , int *marker_y)
{


	// Randomly select the direction in which to divide
	if (NEIGHBOURHOOD == 3)  // transform to triangular lattice
	{
		do
		{
			x = (int)(drand48()*3.0) - 1;
			y = (int)(drand48()*3.0) - 1;
		}
		while ( ((x == 0) && (y == 0)) || (x*y != 0) || (((cell_x + cell_y)%2 == 0) && (y == 1)) || (((cell_x + cell_y)%2 != 0) && (y == -1)) );
	}

	else if (NEIGHBOURHOOD == 4) // only nearest neighbours
	{
		do
		{
			x = (int)(drand48()*3.0) - 1;
			y = (int)(drand48()*3.0) - 1;
		}
		while ( ((x == 0) && (y == 0)) || (x*y != 0) );
	}

	else if (NEIGHBOURHOOD == 8)		// include next-nearest neighbours
	{
		do
		{
			x = (int)(drand48()*3.0) - 1;
			y = (int)(drand48()*3.0) - 1;
		}
		while ((x == 0) && (y == 0));
	}

	//cout << x << " " << y << endl;
	

	//cout << *Ntot << " -> " << tumour[cell_x + x][cell_y + y].dvr << ": Checking " << x << " " << y << endl;

	//if ((tumour[cell_x + x][cell_y + y].dvr == -1) || (tumour[cell_x + x][cell_y + y].dvr == -2))		// Check if neighbour is empty
	if (tumour[cell_x + x][cell_y + y].dvr == -1)
	{

		if ((tumour[cell_x][cell_y].res == 1) || (tumour[cell_x][cell_y].res == 1000)) *Nres += 1;
		if ((tumour[cell_x][cell_y].res == 0) || (tumour[cell_x][cell_y].res == 1000)) *Nwt += 1;

		// Create daughter cell
		if ( tumour[cell_x][cell_y].res == 1000 )
		{
			tumour[cell_x + x][cell_y + y].setGAs(tumour[cell_x][cell_y].dvr , 1 , tumour[cell_x][cell_y].pgr);
		}
		else if ( tumour[cell_x][cell_y].res == 2000 )
		{
			tumour[cell_x + x][cell_y + y].setGAs(tumour[cell_x][cell_y].dvr , 2001 , tumour[cell_x][cell_y].pgr);
		}
		else tumour[cell_x + x][cell_y + y].setGAs(tumour[cell_x][cell_y].dvr , tumour[cell_x][cell_y].res , tumour[cell_x][cell_y].pgr);

		*Ntot += 1;

		//cout << *Ntot << " -> Division occured. " << arising_time << " " << mutated << endl;

		// Insert resistant cell at specified time (i.e. once tumour is comprised of a certain number of cells)
		if ((mutated == false) && (*Ntot == arising_time)) {
			//cout << "Added mutation at N=" << arising_time << endl;
			//if (!death) tumour[cell_x + x][cell_y + y].res = 1000; 
			//else tumour[cell_x + x][cell_y + y].res = 1;
			tumour[cell_x + x][cell_y + y].res = 1000; 
			*marker_x = cell_x + x;
			*marker_y = cell_y + y;
			*Nres += 1; 
			mutated = true;
		}


	/*
		// Insert resistant cell at specified time (i.e. once tumour is comprised of a certain number of cells)
		if ((mutated == false) && (*Ntot == arising_time)) {
			//cout << "Added mutation at N=" << arising_time << endl;
			tumour[cell_x + x][cell_y + y].res = 1000; 
			*Nres += 1; 
			mutated = true;

			// Also "tag" a neutral cell
			// Find position on opposite side of tumour to tag
			opposite_x = 0;
			opposite_y = 0;
			boundary = 0;
			attempts = 0;

			current_diameter = 1.8*pow( (arising_time/M_PI) , (1.0/2.0) );
			do
			{

				if (attempts > *Ntot) 
				{
					current_diameter *= 0.9;
					attempts = 0;
				}
				boundary = 0;

				opposite_x = (int)((2*(*x_b+1))*drand48()) + radius - *x_b;
				opposite_y = (int)((2*(*y_b+1))*drand48()) + radius - *y_b;

				chord = pow( (pow((opposite_x - cell_x),2)+pow((opposite_y - cell_y),2)) , 0.5 );

				// Check if on boundary (by checking if it has empty neighbour)
				for (int i = -1; i < 2; i++)
				{
					for (int j = -1; j < 2; j++)
					{

						if ((NEIGHBOURHOOD == 3) && ((x*y != 0) || (((cell_x + cell_y)%2 == 0) && (y == 1)) || (((cell_x + cell_y)%2 != 0) && (y == -1)))   ) continue;

						if ((NEIGHBOURHOOD == 4) && (x*y != 0)) continue;

						//cout << "Checking (" << opposite_x + i << " , "  << opposite_y + j << ") --> " << tumour[opposite_x + i][opposite_y + j].dvr << endl;
						if (tumour[opposite_x + i][opposite_y + j].dvr == -1)
						{
							boundary += 1;
						}
					}
				}
				//cout << boundary << endl;
				attempts += 1;
			}
			while ((tumour[opposite_x][opposite_y].dvr == -1) || ( (NEIGHBOURHOOD == 8) && (boundary < 3)) || (chord < current_diameter));


			// "Tag" the chosen cell 
			tumour[opposite_x][opposite_y].res = 2000; 
			*Ntagged += 1; 


			//********** Now need to make special cases for res GA = 2000 in other functions (such as selective advantage computation)

		}

	*/
		// Add new GAs to daughter cells
		tumour[cell_x][cell_y].newGAs();
		tumour[cell_x + x][cell_y + y].newGAs();

		// Update bounds on tumour size
		if (fabs(cell_x + x - radius) > *x_b) *x_b = fabs(cell_x + x - radius) + 1;
		if (fabs(cell_y + y - radius) > *y_b) *y_b = fabs(cell_y + y - radius) + 1;
	}

}

/*
// Define cell division in volumetric growth model with "minimal drag" cell displacement as described in Bartek's Nature paper
void MODEL4_divide(Cell ** tumour , int cell_x , int cell_y , int *Ntot , int *Nres  , int *x_b , int *y_b, int radius)
{

	//cout << "\n\n****************************" << endl;
	//cout << "Chosen cell (x,y) = (" << cell_x << " , " << cell_y  <<  ")" << endl;

	// Randomly select the direction in which to divide
	do
	{
		x = (int)(drand48()*3.0) - 1;
		y = (int)(drand48()*3.0) - 1;
	}
	while ((x == 0) && (y == 0));

	//cout << "Direction of division (dx,dy) = (" << x << " , " << y  <<  ")" << endl;


	// Find shortest path to an empty lattice point in the tumour
	for (int i = 0; i < NEIGHBOURHOOD; ++i)
	{
		directions[i] = 0;
	}

	queue = 0;
	chainX[queue] = cell_x + x;
	chainY[queue] = cell_y + y;

	while(1)
	{

		if (tumour[chainX[queue]][chainY[queue]].dvr == -1) break;

		//cout << "Queue = " << queue << ": chain(x,y) = " << chainX[queue] << " , " << chainY[queue] << "   ->    " << tumour[chainX[queue]][chainY[queue]].dvr << endl;

		// Search all directions for shortest path to an empty cell
		direction = 0;
		for (int i = -1; i < 2; i++)
		{
			for (int j = -1; j < 2; j++)
			{

				//cout << "Searching direction (dx,dy) = (" << i << " , " << j  <<  ")" << endl;
				if ( (i == 0) && (j == 0) )
				{
					//++direction;
					continue;
				}

				// Do not check all directions: choose 4 of them at random
				ran = drand48();
				if (ran < 0.5)
				{
					directions[direction] = _maxsize;
					++direction;
					continue;
				}

				//cout << "Searching direction (dx,dy) = (" << i << " , " << j  <<  ")" << endl;

				length = 1;
				while(1)
				{
					coordX = chainX[queue] + length*i;
					coordY = chainY[queue] + length*j;

					//cout << "The status of (dx,dy) = (" << coordX << " , " << coordY  << "   is:  " << tumour[coordX][coordY].dvr << endl;

					if (tumour[coordX][coordY].dvr == -1)
					{
						directions[direction] = length;
						//cout << "********* (dx,dy) = (" << i << " , " << j  <<  ")  ->  length = " << length << endl;
						break;
					}
					else ++length;
				}

				++direction;
			}
		}

		//cout << "Distance to empty cell in each direction: ";
		//for (int i = 0; i < NEIGHBOURHOOD; ++i)
		//{
		//	cout << directions[i] << " ";
		//}
		//cout << endl;

		extended_chain = false;
		while(extended_chain == false)
		{


			// Find which entry in directions list is smallest
			min_length = *Ntot;
			num_mins = 0;
			for (int i = 0; i < NEIGHBOURHOOD; ++i)
			{
				if (directions[i] < min_length) min_length = directions[i];
			}

			// Then count number of directions which are minimum
			for (int i = 0; i < NEIGHBOURHOOD; ++i)
			{
				if (directions[i] == min_length) ++num_mins;
			}

			//cout << "Number of minimum directions is: " << num_mins << " which have the value: " << min_length << endl;

			// If more than one direction is a minimum distance, then choose one with equal probability
			if (num_mins > 1)
			{
				ind = 0;
				while(1)
				{
					ran = drand48();
					//cout << "Checking directions[" << ind%NEIGHBOURHOOD << "]" << endl;
					//cout << "Checking if " << ran << "<" << 1.0/(float)num_mins << endl;
					if ( (directions[ind%NEIGHBOURHOOD] == min_length) && (ran < (1.0/(float)num_mins)))
					{
						chosen_direction = ind%NEIGHBOURHOOD;
						break;
					}

					++ind;
				}
			}

			else 	// Otherwise select the only minimal direction
			{
				for (int i = 0; i < NEIGHBOURHOOD; ++i)
				{
					if (directions[i] == min_length) chosen_direction = i;
				}
			}


			//cout << "Chosen direction is: " << chosen_direction << endl;




			// Before adding the next cell in the chosen direction to the chain, make sure chain is self-avoiding
			// First, find the coordinates of potential new cell 
			direction = 0;
			for (int i = -1; i < 2; i++)
			{
				for (int j = -1; j < 2; j++)
				{

					if ( (i == 0) && (j == 0) ) 
					{
						continue;
					}

					if (direction == chosen_direction)
					{
						coordX = chainX[queue] + i;
						coordY = chainY[queue] + j;

						//cout << "Previous chain coord (x,y) = " << chainX[queue] << " , " << chainY[queue] << ")" << endl;
						//cout << "Next link at coord (x,y) = " << coordX << " , " << coordY << ")" << endl;

						//chainX[queue + 1] = coordX;
						//chainY[queue + 1] = coordY;
					}

					++direction;
				}
			}

			// Second, check these coordinates are not already in the chain
			for (int i = 0; i < (queue+1); ++i)
			{
				if ( ((chainX[i] == coordX) && (chainY[i] == coordY)) && ((chainX[i] == cell_x) && (chainY[i] == cell_y)) )
				{
					// Add a large value to the length of the chosen direction so that it will be essentially eliminated
					//cout << "Potential new link in chain is already taken! ******" << endl;
					directions[chosen_direction] += (int)_maxsize;
				}

				else 
				{
					//cout << "Potential new link in chain is not taken! ******" << endl;
					extended_chain = true;
				}
			}


		}

		// Add next cell in chosen direction to the chain
		chainX[queue + 1] = coordX;
		chainY[queue + 1] = coordY; 
	

		++queue;
	}

	//cout << queue << endl;

	// Once the chain has been constructed, move all cells along one place
	divided = false;
	if (queue > 0)
	{

		// Cell divides and pushes with probability P=1/q where q=queue is the number of cells needed to push
		ran = drand48();
		//cout << 1.0/(float)queue << endl;
		if (ran < (1-(float)(queue/radius)))
		{
			divided = true;

			for (int i = 0; i < queue; ++i)
			{
				//cout << "chain(x" << i << ",y" << i << ") = (" << chainX[i] << " , " << chainY[i] << ") -> " << tumour[chainX[i]][chainY[i]].dvr << " ||| Pushing cell at (x,y)=(" << chainX[queue-i-1] << " , " << chainY[queue-i-1] << ") to (x,y)=(" << chainX[queue-i] << " , " << chainY[queue-i] << ")" << endl;
				tumour[chainX[queue-i]][chainY[queue-i]].setGAs(tumour[chainX[queue-i-1]][chainY[queue-i-1]].dvr , tumour[chainX[queue-i-1]][chainY[queue-i-1]].res , tumour[chainX[queue-i-1]][chainY[queue-i-1]].pgr);
				
				// Update bounds on tumour size
				if (fabs(chainX[queue] + 1 - radius) > *x_b) *x_b = fabs(chainX[queue] + 1 - radius);
				if (fabs(chainY[queue] + 1 - radius) > *y_b) *y_b = fabs(chainY[queue] + 1 - radius);
			}
		}
	}

	//cout << "Bounds -> x=" << *x_b << " | y=" << *y_b << endl;

	if ((divided) || (queue == 0))
	{

		// Create daughter cell
		tumour[cell_x + x][cell_y + y].setGAs(tumour[cell_x][cell_y].dvr , tumour[cell_x][cell_y].res , tumour[cell_x][cell_y].pgr);
		//cout << "Divided into (x,y) = (" << cell_x+x << " , " << cell_y+y  <<  ")" << endl;

		*Ntot += 1;

		// Insert resistant cell at specified time (i.e. once tumour is comprised of a certain number of cells)
		if ((mutated == false) && (*Ntot == arising_time)) {
			cout << "Added mutation at N=" << arising_time << endl;
			tumour[cell_x + x][cell_y + y].res = 1; 
			*Nres += 1; 
			mutated = true;
		}

		// Update bounds on tumour size
		//if (fabs(cell_x + x*(queue + 1) - radius) > *x_b) *x_b = fabs(cell_x + x*(queue + 1) - radius);
		//if (fabs(cell_y + y*(queue + 1) - radius) > *y_b) *y_b = fabs(cell_y + y*(queue + 1) - radius);


		// Add new GAs to daughter cells
		tumour[cell_x][cell_y].newGAs();
		tumour[cell_x + x][cell_y + y].newGAs();

		if (tumour[cell_x][cell_y].res == 1) *Nres += 1;
	}

}
*/


/**********************   Only model 6 birthrates (resistant .= driver mutation) used for initial 2D surface growth fractal dimension studies   ***********************/

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

	else if ((model_number == 6) && !(cell.res == 1000)) *r_birth = pow( (1.0 + sel_adv) , cell.res );
	else if ((model_number == 6) && (cell.res == 1000)) *r_birth = pow( (1.0 + sel_adv) , 1 );

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
	Nwt = 0;
	marker_x = 0;
	marker_y = 0;
	//Ntagged = 0;


	// Seed random number generator
	const int _seed = atoi(argv[2]);
	srand48(_seed);

	// Set selective advantage and arising time from command line args and level of output to stdout
	const double _s = atof(argv[3]);
	const int arising_time = atoi(argv[4]);

	if (std::string(argv[1]) == "-q") quiet = true;
	else quiet = false;
	//NEIGHBOURHOOD = atoi(argv[4]);



	//================== Initialise tumour ====================//

	// Estimate radius of resulting tumour using fitted parameters from previous simulations (slightly over-estimate)
	radius_double = pow ( (_maxsize/M_PI) , (1.0/2.0) );
	radius = (int)(1.5*radius_double);

	if (!quiet) cout << " " << endl;

	Cell ** tumour = new Cell*[2*radius];
	for (int i = 0; i < (2*radius); i++)
	{
		tumour[i] = new Cell[2*radius];

		for (int j = 0; j < (2*radius); j++)
		{

			tumour[i][j].setDVR(-1);		// use -2 as marker for cell which is empty and has *never been occupied*. if empty but once occupied, then the mutations at this position will be -1
			tumour[i][j].setRES(-1);
			tumour[i][j].setPGR(-1);
		
		}
		if (!quiet) printf(" Initialising tumour... %i%%\r", (int)((i+1)*100.0/(2*radius)));
		if (!quiet) fflush(stdout);
	}
	if (!quiet) printf(" Initialising tumour... Done.\r");
	if (!quiet) cout << " " << endl;
		


	// Seed first tumour cell/gland at (x,y) = (0,0)
	tumour[radius][radius].setDVR(0);
	tumour[radius][radius].setRES(0);
	tumour[radius][radius].setPGR(0);

	Ntot += 1;
	Nwt += 1;







	//================== Simulate tumour growth ==================//

	iter = 0;
	max_birth = 0.0;
	x = 0;
	y = 0;
	x_b = 10;
	y_b = 10;
	counter = 0;


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
		}
		while (tumour[cell_x][cell_y].dvr == -1);
		//while ((tumour[cell_x][cell_y].dvr == -1) || (tumour[cell_x][cell_y].dvr == -2));

		//cout << Ntot << " " << radius << " " << cell_x << " " << cell_y << endl;

		// Compute birth rate of cell
		compute_birthrate( tumour[cell_x][cell_y] , adv_model , &r_birth , _s );		
		
		// Update survival rate (which is a function of birth rate)
		r_surv = r_birth * bdratio;
	

		// Update maximal birth and death rate of all cells 
		if (r_birth > max_birth) max_birth = r_birth;


		// Query the cell's position in the cell cycle (i.e. if it is ready to divide)
		if (drand48() < (r_birth/max_birth))
		{


			// First the cell is given the option to die
			if ( (Ntot > 1) && (death) && (drand48() < (r_surv/max_birth)) )
			{
				if (tumour[cell_x][cell_y].res == 1) Nres -= 1;
				if (tumour[cell_x][cell_y].res == 0) Nwt -= 1;

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
/*
			else if (div_model == 1) MODEL1_divide(tumour , cell_x , cell_y , &Ntot , &Nres , &x_b , &y_b , radius);
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
*/
			else if (div_model == 3) MODEL3_divide(tumour , cell_x , cell_y , &Ntot , &Nres, &Nwt, &Ntagged , &x_b , &y_b , radius , arising_time , &marker_x , &marker_y);

			//else if (div_model == 4) MODEL4_divide(tumour , cell_x , cell_y , &Ntot , &Nres , &x_b , &y_b , radius);



			// Progress time variable
			t += 1.0/(max_birth * Ntot);
		}


		// Ask if the cell dies (intrinsic death rate)
		if ( (tumour[cell_x][cell_y].dvr != -1) && (death) && (drand48() < (r_death/max_birth)) )
		{
			if (tumour[cell_x][cell_y].res == 1) Nres -= 1;
			if (tumour[cell_x][cell_y].res == 0) Nwt -= 1;

			// Delete cell from tumour
			tumour[cell_x][cell_y].setDVR(-1);
			tumour[cell_x][cell_y].setRES(-1);
			tumour[cell_x][cell_y].setPGR(-1);

			// Size of tumour is reduced by 1
			Ntot -= 1;

			// Progress time variable
			t += 1.0/(max_birth * Ntot);
		}

		if (iter%10000 == 0)
		{

			/*
			if (Ntot >= arising_time)
			{
				previous_Nres[counter] = Nres;
				previous_Nwt[counter] = Nwt;		

				// If number of resistant cells is not growing, abort simulation
				if ((counter >= 10) && ((previous_Nres[counter]-previous_Nres[counter-10]) == 0))
				{
					if (!quiet) cout << "0 " << Nres << " " << Nwt << endl;
					else cout << "0" << endl;
					exit(0);
				}

				// If number of wild-type cells is not growing, can end the simulation and process data
				if ((counter >= (10*_s)) && ((previous_Nwt[counter]-previous_Nwt[(int)(counter-(10*_s))]) == 0))
				{
					break;
				}

				counter += 1;
			}
			*/


			if (!quiet) cout << Ntot << " " << Nres << " " << Nwt << endl;

		}


	//if (iter == 1000000) break;


	} while (Ntot < _maxsize);

	if (!quiet) cout << " " << endl;


/*
	// Check for resistant cells at outer edge -> use this as a measure of successful simulation
	success = false;
	for (int i = (radius - x_b - 1); i < (radius + x_b + 2); ++i)
	{
		for (int j = (radius - y_b - 1); j < (radius + y_b + 2); ++j)
		{

			if ( tumour[i][j].res > 0 )
			{
				// Check neighbours
				for (int dx = -1; dx < 2; dx++)
				{
					for (int dy = -1; dy < 2; dy++)
					{

						if (tumour[i+dx][j+dy].res == -1) success = true;

					}
				}
			}
		}
	}

	// Exit early if simulation deemed "unsuccessful"
	if (!success)
	{
		if (!quiet) cout << "0 " << Nres << endl;
		else cout << "0" << endl;
		exit(0);
	}
*/
/*
	
	// Exit early if number of resistant cells is small
	if ( (arising_time <= 2000) && (Nres < 2500) )
	{
		if (!quiet) cout << "0 " << Nres << endl;
		else cout << "0" << endl;
		exit(0);
	}
	else if ( (arising_time < 10000) && (Nres < 1500) )
	{
		if (!quiet) cout << "0 " << Nres << endl;
		else cout << "0" << endl;
		exit(0);
	}
	else if (Nres < 1000)
	{
		if (!quiet) cout << "0 " << Nres << endl;
		else cout << "0" << endl;
		exit(0);
	}

*/


/*

	// Check that red and blue regions are never adjascent 
	for (int i = (radius - x_b - 1); i < (radius + x_b + 2); ++i)
		{
			for (int j = (radius - y_b - 1); j < (radius + y_b + 2); ++j)
				{

			if (tumour[i][j].res == 1000) continue;	

			if ( tumour[i][j].res > 0 )
			{		

				// Check for any neighbouring empty lattice sites
				for (int x = -1; x < 2; x++)
				{
					for (int y = -1; y < 2; y++)
					{

						if ((NEIGHBOURHOOD == 3) && ((x*y != 0) || (((i + j)%2 == 0) && (y == 1)) || (((i + j)%2 != 0) && (y == -1)))   ) continue;

						if ((NEIGHBOURHOOD == 4) && (x*y != 0)) continue;

						if ( (tumour[i+x][j+y].res == 1000) || (tumour[i+x][j+y].res == 2000) ) continue;

						if ( ((i+x) < 0) || ((j+y) < 0) || ((i+x) > 2*radius ) || ((j+y) > 2*radius) ) continue;

						if (( tumour[i+x][j+y].res > 0 ) && ( tumour[i+x][j+y].res != tumour[i][j].res ))
						{
							cout << "0" << endl;
							//cout << "100" << endl;
							exit(0);
						}					
  
					}
				}	
			}
		}
	}
*/


	//================== Open data files ==================//

	stringstream f;
	f.str("");
	f << "./DATA_surfaceGrowth2D/maxSize=" << _maxsize << "_s=" << _s << "_tmut=" << arising_time << "_death=" << death << "_bdRatio=" << bdratio << "/seed=" << _seed;
	DIR *dir = opendir(f.str().c_str());
	if(!dir)
	{
		f.str("");
		f << "mkdir -p ./DATA_surfaceGrowth2D/maxSize=" << _maxsize << "_s=" << _s << "_tmut=" << arising_time << "_death=" << death << "_bdRatio=" << bdratio << "/seed=" << _seed;
		system(f.str().c_str());
	}

	//ofstream NversusT_file;
	//f.str("");
	//f << "./DATA_surfaceGrowth2D/maxSize=" << _maxsize << "_s=" << _s << "_tmut=" << arising_time << "_death=" << death << "_bdRatio=" << bdratio << "/seed=" << _seed << "/N(t).dat";
	//NversusT_file.open(f.str().c_str());

	ofstream tumour_file;
	f.str("");
	f << "./DATA_surfaceGrowth2D/maxSize=" << _maxsize << "_s=" << _s << "_tmut=" << arising_time << "_death=" << death << "_bdRatio=" << bdratio << "/seed=" << _seed << "/tumour.csv";
	tumour_file.open(f.str().c_str());
/*
	ofstream mathematica_file;
	f.str("");
	f << "./DATA_surfaceGrowth2D/_NEIGHBOURHOOD=" << NEIGHBOURHOOD << "_s=" << _s << "_arisingtime=" << arising_time << "/maxsize="<< _maxsize << "_seed=" << _seed << "_NEIGHBOURHOOD=" << NEIGHBOURHOOD << "_s=" << _s << "_ut=" 
		<< _ut << "_ud=" << _ud << "_arisingtime=" << arising_time << "_death=" << death << "_divmodel=" << div_model << "_advmodel=" << adv_model << "/plot_mathematica.csv";
	mathematica_file.open(f.str().c_str());
*/
	if (!quiet) cout << " " << endl;
	if (!quiet) cout << "Created output files..." << endl;





	// Write tumour data to file
	tumour_file << "x coord, y coord, z coord, drivers, resistant, passengers" << endl;
	/*
	for (int i = (radius - x_b - 1); i < (radius + x_b + 2); ++i)
	{
		for (int j = (radius - y_b - 1); j < (radius + y_b + 2); ++j)
		{
			if ( tumour[i][j].dvr != -1 ) tumour_file << i << "," << j << ",0," << tumour[i][j].dvr << "," << tumour[i][j].res << "," << tumour[i][j].pgr << endl;
		}
	}
	*/


	// Manually write cell at marker coordinates (for image processing program)
	tumour_file << marker_x << "," << marker_y << ",0,0,1000,0" << endl;

	for (int i = 0; i < (2*radius); i++)
	{
		for (int j = 0; j < (2*radius); j++)
		{
	
			if ( tumour[i][j].res != -1 ) tumour_file << i << "," << j << ",0," << tumour[i][j].dvr << "," << tumour[i][j].res << "," << tumour[i][j].pgr << endl;


			/*
			// Only write data to file if it is a mutated cell at the outer edge (i.e. has a '-2' neighbour), or a WT cell neighbouring a mutated cell
			// Check neighbours
			for (int di = -1; di < 2; di++)
			{
				for (int dj = -1; dj < 2; dj++)
				{
					if ( (i >= 1) && (j >= 1) && (i <= ((2*radius) - 2)) && (j <= ((2*radius) - 2)) )
					{



						if (tumour[i][j].res == 1) && (tumour[i+di][j+dj].res == 0) // mutant next to WT
						{

						}
						else if (tumour[i][j].res == 1) && (tumour[i+di][j+dj].res == -2) // mutant at outer edge
						{

						}
						else if (tumour[i][j].res == -2) && (tumour[i+di][j+dj].res == 1)
						{

						}
						else if (tumour[i][j].res == 0) && (tumour[i+di][j+dj].res == 1)

						if ( tumour[i+di][j+dj].res == 1 ) tumour_file << i << "," << j << ",0," << tumour[i][j].dvr << "," << tumour[i][j].res << "," << tumour[i][j].pgr << endl;
					}
				}
			}
			*/

			//if ( j == ((2*radius)-1) )
			//{
			//	mathematica_file << tumour[i][j].res << endl;
			//}
			//else mathematica_file << tumour[i][j].res << ",";

		}
	}


	// Manually write cell at marker coordinates (for image processing program)
	tumour_file << marker_x << "," << marker_y << ",0,0,1000,0" << endl;



	//NversusT_file.close();
	tumour_file.close();


	if (!quiet) cout << "1 " << Nres << " " << marker_x << " " << marker_y << endl;
	else cout << "1" << endl;

	return 0;
}

















