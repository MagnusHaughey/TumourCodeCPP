


/****************************************************************************************************

Worked on volumetric growth function

Now need to add in the 'sweep' functionality and add animation output 

****************************************************************************************************/




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
# include <unistd.h>
//# include <thread>
//# include <algorithm> 
//# include <vector>
//# include <sys/types.h>
//# include <sys/stat.h>
# include <stdio.h>
#include <getopt.h>



using namespace std;



/*******************************************************************************/



// Define global variables
const int _maxsize = 5e4;
const double _ut = 2.0;
const double _ud = 0.1;
const int div_model = 4;
const int adv_model = 6;
const bool death = true;



// Uncomment one of the following two lines to choose from triangular lattice (3), nearest-neighbour (4) and next-nearest-beighbour neighbourhoods (8) respectively
//const int NEIGHBOURHOOD = 3;
//const int NEIGHBOURHOOD = 4;
const int NEIGHBOURHOOD = 8;



double radius_double, t, _s, max_birth, ran, r_birth, current_diameter, chord, bdratio, optimal_direction_i, optimal_direction_j, optimal_vector_norm, scalar_prod, vector_norm, rescaled_min_length , rand_double , r_birth_WT_normalised, r_birth_mut_normalised, r_tot;
int radius, Ntot, Ntagged, iter, x, y, z, cell_x, cell_y, cell_z, empty_neighbours, dir, queue, ind, length, coordX, coordY, boundary, attempts, counter, sweeps, previous_link_direction;
int range, arising_time, xrange, yrange, zrange, x_b, y_b, z_b, res_mut, Nres, Nwt, direction, chosen_direction, min_length, num_mins, opposite_x, opposite_y, marker_x , marker_y, sweep;

char padded_sweep[6];



//const double bdratio = 1.0/6.0;
double r_surv = r_birth * bdratio;		// Clonal survival rate
//const double r_death = 0.5 * log(2.0);	// Intrinsic cell death rate
const double r_death = 0.0;			// Intrinsic cell death rate




// Define arrays which will contain relative coordinates of empty neighbours for a chosen cell (for model 2)
int x_nn[NEIGHBOURHOOD];
int y_nn[NEIGHBOURHOOD];
int chainX[(int)_maxsize];		// this is potentially lazy use of memory
int chainY[(int)_maxsize];
int directions[NEIGHBOURHOOD];
int previous_Nres[(int)(_maxsize)];
int previous_Nwt[(int)(_maxsize)];



bool chain_stuck, WT_DIVISION, MUT_DIVISION;
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





void surface_division(Cell ** tumour , int cell_x , int cell_y , int *Ntot , int *x_b , int *y_b , int radius)
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
	if (tumour[cell_x + x][cell_y + y].dvr < 0)
	{


		// Create daughter cell
		tumour[cell_x + x][cell_y + y].setGAs(tumour[cell_x][cell_y].dvr , tumour[cell_x][cell_y].res , tumour[cell_x][cell_y].pgr);

		*Ntot += 1;

		//cout << *Ntot << " -> Division occured. " << arising_time << " " << mutated << endl;

	
		// Add new GAs to daughter cells
		//tumour[cell_x][cell_y].newGAs();
		//tumour[cell_x + x][cell_y + y].newGAs();

		// Update bounds on tumour size
		if (fabs(cell_x + x - radius) > *x_b) *x_b = fabs(cell_x + x - radius) + 1;
		if (fabs(cell_y + y - radius) > *y_b) *y_b = fabs(cell_y + y - radius) + 1;
	}

}









// Volumetric growth with "minimal drag" cell displacement, following "B. Waclaw et al., Nature 525, 7568 (September 10, 2015): 261-264"
void volumetric_division(Cell ** tumour , int cell_x , int cell_y , int *Ntot , int *Nres , int *Nwt , int *x_b , int *y_b, int radius)
{

	//cout << "\n\n****************************" << endl;
	//cout << "Chosen cell (x,y) = (" << cell_x << " , " << cell_y << ")" << endl;

	//if (tumour[cell_x][cell_y].res == 1)
	//{
	//	cout << "Mutated cell" << endl;
	//}

	// If dividing cell is WT and is fully surrounded by resistant cells, exit out of function
	chain_stuck = true;
	if (tumour[cell_x][cell_y].res == 0)
	{
		for (int i = -1; i < 2; i++)
		{
			for (int j = -1; j < 2; j++)
			{

				if (tumour[cell_x + x][cell_y + y].res != 1) chain_stuck = false;

			}
		}

		// If no empty or WT neighbour found, exit function
		if (chain_stuck)
		{
			//cout << "No empty or WT neighbour found. Returning..." << endl;
			return;
		}
	}


	// Randomly select the direction in which to divide
	do
	{
		x = (int)(drand48()*3.0) - 1;
		y = (int)(drand48()*3.0) - 1;
	}
	while (((x == 0) && (y == 0)) || ((tumour[cell_x][cell_y].res == 0) && (tumour[cell_x + x][cell_y + y].res == 1)));




	//cout << "Direction of division (dx,dy) = (" << x << " , " << y << ")" << endl;


	

	// Find shortest path to an empty lattice point in the tumour
	for (int i = 0; i < NEIGHBOURHOOD; ++i)
	{
		directions[i] = 0;
	}

	queue = 0;

	//if (tumour[cell_x + x][cell_y + y][cell_z + z].stem_cell == 1) return;

	chainX[queue] = cell_x;
	chainY[queue] = cell_y;

	queue += 1;

	chainX[queue] = cell_x + x;
	chainY[queue] = cell_y + y;

	while(1)
	{

		if (tumour[chainX[queue]][chainY[queue]].res < 0) break;

		//cout << "Queue = " << queue << ": chain(x,y,z) = " << chainX[queue] << " , " << chainY[queue] << "   ->    " << tumour[chainX[queue]][chainY[queue]].cellID << endl;



		// Search all directions for shortest path to an empty cell
		direction = 0;
		for (int i = -1; i < 2; i++)
		{
			for (int j = -1; j < 2; j++)
			{

					if ( (i == 0) && (j == 0) )
					{
						//++direction;
						continue;
					}

					//cout << "Searching direction (dx,dy) = (" << i << " , " << j  <<  ")" << endl;

					// Do not check all directions: choose some of them at random
					//ran = drand48();
					//if (ran < 0.5)
					//{
					//	directions[direction] = _maxsize;
					//	++direction;
					//	continue;
					//}

					if (directions[direction] == _maxsize)
					{
						//cout << "Skipping (dx,dy,dz) = (" << i << " , " << j << " , " << k <<  ") --> direction=" << directions[direction] << endl;
						++direction;
						continue;
					}

					//cout << "Searching direction (dx,dy) = (" << i << " , " << j <<  ")" << endl;

					length = 1;
					while(1)
					{
						coordX = chainX[queue] + length*i;
						coordY = chainY[queue] + length*j;

						//cout << "The mutations of (dx,dy) = (" << coordX << " , " << coordY << "))   is:  " << tumour[coordX][coordY].res << endl;

						if (tumour[coordX][coordY].res < 0)
						{
							directions[direction] = length;
							//cout << "************************************ (dx,dy) = (" << i << " , " << j << " , " << k <<  ")  ->  length = " << length << endl;
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
			chain_stuck = true;
			for (int i = 0; i < NEIGHBOURHOOD; ++i)
			{
				if (directions[i] < min_length) min_length = directions[i];
				if (directions[i] <= _maxsize) chain_stuck = false;
			}

			if (chain_stuck == true) return;

			// Calculate probability of pushing according to minimum length
			//ran = drand48();
			//if (ran > (1.0-pow(((float)min_length/(float)(radius_double*0.25)) , 1)))
			//{
			//	return;
			//}

			// Then count number of directions which are minimum
			for (int i = 0; i < NEIGHBOURHOOD; ++i)
			{
				if (directions[i] == min_length) ++num_mins;
			}

			//cout << "Number of minimum directions is: " << num_mins << " which have the value: " << min_length << endl;



			//cout << "++++++++++++++++ chain[0]=(" << chainX[0] << "," << chainY[0] << ") \t\t chain[1]=(" << chainX[1] << "," << chainY[1] << endl;


			if (queue >= 1)
			{
				// Check in which direction the previous link in the chain is at
				direction = 0;
				previous_link_direction = 0;
				for (int i = -1; i < 2; i++)
				{
					for (int j = -1; j < 2; j++)
					{

						if ( (i == 0) && (j == 0) )
						{
							continue;
						}

						//cout << chainX[queue]+i << " " << chainX[queue-1] << "\t\t" << chainY[queue]+j << " " << chainY[queue-1] << endl;
						if ((chainX[queue]+i == chainX[queue-1]) && (chainY[queue]+j == chainY[queue-1]))
						//if (chainX[queue]+i][chainY[queue]+j].res == tumour[chainX[queue-1]][chainY[queue-1]].cellID)
						{
							previous_link_direction = direction;
 

							//cout << "Direction of previous link in chain: (" << chainX[queue-1] << "," << chainY[queue-1] << ") -> (" << chainX[queue] << "," << chainY[queue] << ")  ====  " << previous_link_direction << ", opposite direction -> " << 7-previous_link_direction << endl;
							break;
						}

						++direction;

					}
				}

				direction = 0;
				for (int i = -1; i < 2; i++)
				{
					for (int j = -1; j < 2; j++)
					{

						if ( (i == 0) && (j == 0) )
						{
							continue;
						}

						if (direction == 7-previous_link_direction)
						{
							optimal_direction_i = i;
							optimal_direction_j = j;
							optimal_vector_norm = pow(((i*i) + (j*j)) , 0.5);

							optimal_direction_i /= optimal_vector_norm;
							optimal_direction_j /= optimal_vector_norm;

							//cout << "Optimal direction = " << optimal_direction_i << " " << optimal_direction_j << endl;
						}

						++direction;

					}
				}


				// Re-scale distances vector according to relative direction to 'forward' chain direction
				direction = 0;
				for (int i = -1; i < 2; i++)
				{
					for (int j = -1; j < 2; j++)
					{

						if ( (i == 0) && (j == 0) )
						{
							continue;
						}

						vector_norm = pow(((i*i) + (j*j)) , 0.5);

						// Scalar product with unit vector pointing in 'optimal direction'
						scalar_prod = (optimal_direction_i*i/vector_norm) + (optimal_direction_j*j/vector_norm);

						//cout << "Candidate direction = " << i/vector_norm << " " << j/vector_norm << endl;

						// Rescale to within range [0,1]
						scalar_prod = (scalar_prod + 1.0)/2.0;
						directions[direction] *= 1.0 - scalar_prod;
						//directions[direction] /= scalar_prod;

						++direction;

					}
				}



				// Find new minimum after rescaling 
				rescaled_min_length = (double)*Ntot;
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

				//cout << "After rescaling: distance to empty cell in each direction: ";
				//for (int i = 0; i < NEIGHBOURHOOD; ++i)
				//{
				//	cout << directions[i] << " ";
				//}
				//cout << endl;

				//cout << "After rescaling, new number of minimum directions is: " << num_mins << " which have the value: " << min_length << endl;


				//exit(0);

			}

			// If more than one direction is a minimum distance, then choose one with equal probability
			if (num_mins > 1)
			{




				ind = 0;
				while(1)
				{
					ran = drand48();
					//if (directions[ind%NEIGHBOURHOOD] == min_length) 
					//{
						//cout << "Checking directions[" << ind%NEIGHBOURHOOD << "]" << endl;
						//cout << "Checking if " << ran << "<" << 1.0/(float)num_mins << endl;
					//}
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
						//cout << direction << endl;

						if (direction == chosen_direction)
						{
							coordX = chainX[queue] + i;
							coordY = chainY[queue] + j;
	
							//cout << "Previous chain coord (x,y,z) = " << chainX[queue] << " , " << chainY[queue] << ")" << endl;
							//cout << "Next link at coord (x,y,z) = " << coordX << " , " << coordY << " , " << ")" << endl;

							//chainX[queue + 1] = coordX;
							//chainY[queue + 1] = coordY;
						}

						++direction;

					//++direction;
				}
			}

			//for (int i = 0; i < (queue+1); ++i)
			//{
			//	cout << "++++++++++++++++ chain[" << i << "]=(" << chainX[i] << "," << chainY[i] << ")\t [" << tumour[chainX[i]][chainY[i]].res << "]" << endl;
			//}


			// Second, check these coordinates are not already in the chain, and they are consistent with the rules of pushing (resistant cells can push all cells, WT cells can only push WT cells)
			extended_chain = true;
			for (int i = 0; i < (queue+1); ++i)
			{

				if ((tumour[cell_x][cell_y].res == 0) && (tumour[coordX][coordY].res == 1))
				{
					//cout << "WT cell trying to push a mutated cell! ****** " << endl;
					extended_chain = false;
					chain_stuck = true;
					directions[chosen_direction] += (int)_maxsize;
				}

				//cout << "Checking queue entry " << i << "/" << queue << endl;
				else if ((chainX[i] == coordX) && (chainY[i] == coordY))
				{
					// Add a large value to the length of the chosen direction so that it will be essentially eliminated
					//cout << "Potential new link in chain is already taken! ****** " << endl;
					directions[chosen_direction] += (int)_maxsize;
					extended_chain = false;
				}
			}

			if (extended_chain == true)
			{
				//cout << "Potential new link in chain is not taken! ******" << endl;
			}

			if (chain_stuck == true)
			{
				return;
			}


		}

		// Add next cell in chosen direction to the chain
		chainX[queue + 1] = coordX;
		chainY[queue + 1] = coordY;
		

		++queue;
	}

	//cout << queue << endl;



	//cout << " ~~~~~~~~~~~~~~~~~~ " << endl;
	//cout << "Cell at (" << cell_x << " , " << cell_y  << ") [" << tumour[cell_x][cell_y].res << "] wants to divide into (" << cell_x+x << " , " << cell_y+y << ") [" << tumour[cell_x+x][cell_y+y].res << "]" << endl;
	//for (int i = 0; i < queue; ++i)
	//{
	//	cout << chainX[queue] << "," << chainY[queue] << endl;
	//}
	//cout << " ~~~~~~~~~~~~~~~~~~ " << endl;

	// Once the chain has been constructed, move all cells along one place
	//divided = false;
	if (queue > 0)
	{

		// Cell divides and pushes with probability P=1/q where q=queue is the number of cells needed to push
		//ran = drand48();
		//cout << "queue=" << queue << ": (float)(queue/radius)=" << ((float)queue/(float)radius) << " ---> pushing probability=" << (1.0-pow(((float)queue/(float)(radius*0.1)) , 1)) << endl;
		//if (ran < (1.0-pow(((float)queue/(float)(radius*1.0)) , 1)))
		//{
			//divided = true;


			//cout << chainX[queue]-cell_x << " " << chainY[queue]-cell_y << endl;

			for (int i = 0; i < queue-1; ++i)
			{
				//cout << "chain(x" << i << ",y" << i << ") = (" << chainX[i] << " , " << chainY[i] << ") -> " << tumour[chainX[i]][chainY[i]].res << " ||| Pushing cell at (x,y)=(" << chainX[queue-i-1] << " , " << chainY[queue-i-1] << ") to (x,y)=(" << chainX[queue-i] << " , " << chainY[queue-i] << ")" << endl;
				tumour[chainX[queue-i]][chainY[queue-i]].setGAs(tumour[chainX[queue-i-1]][chainY[queue-i-1]].dvr , tumour[chainX[queue-i-1]][chainY[queue-i-1]].res , tumour[chainX[queue-i-1]][chainY[queue-i-1]].pgr);


				// Update bounds on tumour size
				if (fabs(chainX[i] + 1 - radius) > *x_b) *x_b = fabs(chainX[i] + 1 - radius);
				if (fabs(chainY[i] + 1 - radius) > *y_b) *y_b = fabs(chainY[i] + 1 - radius);

			}
		//}
	}

	else 			// Even if queue=0, check that newly created cell increases any bounds
	{

		// Update bounds on tumour size
		if (fabs(chainX[0] + 1 - radius) > *x_b) *x_b = fabs(chainX[0] + 1 - radius);
		if (fabs(chainY[0] + 1 - radius) > *y_b) *y_b = fabs(chainY[0] + 1 - radius);


	}



	//cout << "Bounds -> x=" << *x_b << " | y=" << *y_b << endl;

	//if ((divided) || (queue == 0))
	//{

		// Create daughter cell
		//tumour[cell_x + x][cell_y + y].setID(*next_cellID);
		//tumour[cell_x + x][cell_y + y].setCellType(0);
		//inheritMutations( tumour[cell_x + x][cell_y + y] , tumour[cell_x][cell_y] , mutations );		// 2nd daughter cell inherits mutations of mother cell
		tumour[cell_x + x][cell_y + y].setGAs(tumour[cell_x][cell_y].dvr , tumour[cell_x][cell_y].res , tumour[cell_x][cell_y].pgr);
		//cout << "Divided into (x,y) = (" << cell_x+x << " , " << cell_y+y <<  ")" << endl;

		*Ntot += 1;



		if (tumour[cell_x][cell_y].res == 1) *Nres += 1;
		if (tumour[cell_x][cell_y].res == 0) *Nwt += 1;
		 


		// Add new GAs to daughter cells
		//tumour[cell_x][cell_y].newGAs();
		//tumour[cell_x + x][cell_y + y].newGAs();

			
	//}

}









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


	//=================== Parse command line arguments

	// Set quiet flag
	if (std::string(argv[1]) == "-q") quiet = true;
	else quiet = false;

	// Seed random number generator
	const int _seed = atoi(argv[2]);
	srand48(_seed);

	// Set selective advantage and arising time from command line args and level of output to stdout
	const double _s = atof(argv[3]);

	// Birth/death ratio
	//bdratio = atof(argv[4]);

	// Number of sweeps
	sweeps = atoi(argv[4]);




	//=================== Create animation directory
	stringstream f;
	stringstream g;


	f.str("");
	f << "./2D_DATA/maxSize=" << _maxsize << "_sweeps=" << sweeps << "_ut=" << _ut 
			<< "_BDratio=" << bdratio << "_s=" << _s << "/" << _seed << "/animation";
	DIR *dir = opendir(f.str().c_str());
	if(!dir)
	{
		f.str("");
		f << "mkdir -p ./2D_DATA/maxSize=" << _maxsize << "_sweeps=" << sweeps << "_ut=" << _ut 
			<< "_BDratio=" << bdratio << "_s=" << _s << "/" << _seed << "/animation";
		system(f.str().c_str());
	}





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
		if (!quiet) printf("\tInitialising tumour... %i%%\r", (int)((i+1)*100.0/(2*radius)));
		if (!quiet) fflush(stdout);
	}
	if (!quiet) printf("\tInitialising tumour... Done.\r");
	if (!quiet) cout << " " << endl;
		


	// Seed first tumour cell/gland at (x,y) = (0,0)
	tumour[radius][radius].setDVR(0);
	tumour[radius][radius].setRES(0);
	tumour[radius][radius].setPGR(0);

	Ntot += 1;
	Nwt += 1;


	//==============================================================================================//
	//===================================    Initialise cells    ===================================//
	//==============================================================================================//

	/*

	Initialise the cells on the lattice by placing the first cell at the centre and letting it divide 
	continuously (no cell death, no mutations). This will produce a ball of cells of size _maxsize.
	
	*/



	iter = 0;
	max_birth = 0.0;
	x = 0;
	y = 0;
	x_b = 0;
	y_b = 0;

	if (!quiet) printf("\r\tInitialising cells... ");
	if (!quiet) fflush(stdout);

	do
	{
		
		++iter;

		//cout << "\n" << "ITERATION #" << iter << endl;

		// Randomly select one cell to divide
		cell_x = 0;
		cell_y = 0;

		do
		{

			cell_x = (int)((2*(x_b+1))*drand48()) + radius - x_b;
			cell_y = (int)((2*(y_b+1))*drand48()) + radius - y_b;

			//cout << x_b << " " << y_b << " " << z_b << " " << 2*radius << " " << cell_x << " " << cell_y << " " << cell_z << endl;

		}
		while (tumour[cell_x][cell_y].dvr == -1);
		
		//cout << x_b << " " << y_b << " " << z_b << " " << 2*radius << " " << cell_x << " " << cell_y << " " << cell_z << " " << Ntot << endl;


		// Cell divides with probability 1
		surface_division(tumour , cell_x , cell_y , &Ntot , &x_b , &y_b , radius );

		//cout << Ntot << endl;


	} while (Ntot < _maxsize);

	Nwt = Ntot;

	if (!quiet) printf("\r\tInitialising cells... Done.");
	if (!quiet) cout << " " << endl;




	//====================== Insert mutated cells in the centre 

	for (int i = -1; i < 2; i++)
	{
		for (int j = -1; j < 2; j++)
		{
			tumour[radius + i][radius + j].res = 1;
			Nres += 1;
			Nwt -= 1;
		}
	}










	//==============================================================================================//
	//===================================    Turn on dynamics    ===================================//
	//==============================================================================================//

	/*

	Run the simulation for several sweeps with cell death and mutation turned on, using basic Gillespie algorithm.
	
	*/


	ofstream animation_file;
	f.str("");
	f << "./2D_DATA/maxSize=" << _maxsize << "_sweeps=" << sweeps << "_ut=" << _ut 
			<< "_BDratio=" << bdratio << "_s=" << _s << "/" << _seed << "/animation/000.csv";
	animation_file.open(f.str().c_str());


	for (int i = 0; i < (2*radius); i++)
	{
		for (int j = 0; j < (2*radius); j++)
		{
			if (tumour[i][j].res > -1)
			{

				// Print cell coordinates to file
				animation_file << i << "," << j << "," << tumour[i][j].res << endl;

			}

		}
	}

	animation_file.close();


	iter = 0;
	max_birth = 0.0;
	x = 0;
	y = 0;
	//x_b = 10;
	//y_b = 10;
	counter = 0;





	for (int sweep = 0; sweep < sweeps; ++sweep)
	//sweep = 0;
	//do
	{

		//++sweep;

		//============================ Kill 10% of cells
		do
		{

			// Randomly select one cell
			cell_x = 0;
			cell_y = 0;

			do
			{

				cell_x = (int)((2*(x_b+1))*drand48()) + radius - x_b;
				cell_y = (int)((2*(y_b+1))*drand48()) + radius - y_b;

				//cout << x_b << " " << y_b << " " << z_b << " " << 2*radius << " " << cell_x << " " << cell_y << " " << cell_z << endl;

			}
			while (tumour[cell_x][cell_y].res < 0);


			// Ensure mutated cells do not all die
			if (( Nres < 2) && (tumour[cell_x][cell_y].res == 1)) continue;


			// Delete cell from lattice
			if (tumour[cell_x][cell_y].res == 1)
			{

				Nres -= 1;
	
				// Delete cell from tumour
				tumour[cell_x][cell_y].setDVR(-3);
				tumour[cell_x][cell_y].setRES(-3);
				tumour[cell_x][cell_y].setPGR(-3);
			}
			if (tumour[cell_x][cell_y].res == 0)
			{
				Nwt -= 1;

				// Delete cell from tumour
				tumour[cell_x][cell_y].setDVR(-2);
				tumour[cell_x][cell_y].setRES(-2);
				tumour[cell_x][cell_y].setPGR(-2);
			}


			// Size of tumour is reduced by 1
			Ntot -= 1;


			//cout << Ntot << endl;

		} while (Ntot > (int)(0.9*_maxsize));




		//============================ Re-populate tumour
		do
		{

		
			++iter;

			/*
			if (iter == 63219)
			{
				cout << "cell_x=" << cell_x << "    ,     cell_y=" << cell_y << endl;

				ofstream animation_file;
				f.str("");
				f << "./2D_DATA/maxSize=" << _maxsize << "_sweeps=" << sweeps << "_ut=" << _ut 
						<< "_BDratio=" << bdratio << "_s=" << _s << "/" << _seed << "/animation/iter=" << iter << ".csv";
				animation_file.open(f.str().c_str());


				for (int i = 0; i < (2*radius); i++)
				{
					for (int j = 0; j < (2*radius); j++)
					{
						if (tumour[i][j].res > -1)
						{

							// Print cell coordinates to file
							animation_file << i << "," << j << "," << tumour[i][j].res << endl;

						}

					}
				}

				animation_file.close();

				cout << "Wrote spatial data file" << endl;

				cin >> ran;
			}
			*/


			// Compute normalised division rates based on number of WT and mutated cells
			r_birth_WT_normalised = (double)Nwt;
			r_birth_mut_normalised = (double)Nres * (1.0 + _s);

			//cout << "Ntot=" <<Ntot << "\tNwt=" << Nwt << "\tNres=" << Nres << "    " << r_birth_WT_normalised << " " << r_birth_mut_normalised << endl;

			r_tot = r_birth_WT_normalised + r_birth_mut_normalised;
			r_birth_WT_normalised /= r_tot;
			r_birth_mut_normalised /= r_tot;

			//cout << "Ntot=" <<Ntot << "\tNwt=" << Nwt << "\tNres=" << Nres << "    " << r_birth_WT_normalised << " " << r_birth_mut_normalised << endl;

			// Choose WT or mutated cell division
			WT_DIVISION = false;
			MUT_DIVISION = false;
			rand_double = drand48();

			if (rand_double < r_birth_WT_normalised)
			{
				WT_DIVISION = true;
			}
			else if (rand_double < r_birth_WT_normalised + r_birth_mut_normalised)
			{
				MUT_DIVISION = true;
			}

			//cout << "\nIter=" << iter << "\tNtot=" << Ntot << "\t WT_DIVISION=" << WT_DIVISION << "\tMUT_DIVISION=" << MUT_DIVISION << endl;


			//cout << iter << " " << Ntot << endl;


			// Randomly select one cell to divide
			cell_x = 0;
			cell_y = 0;

			do
			{

				cell_x = (int)((2*(x_b+1))*drand48()) + radius - x_b;
				cell_y = (int)((2*(y_b+1))*drand48()) + radius - y_b;

				//cout <<  WT_DIVISION << " " << MUT_DIVISION << " " << cell_x << " " << cell_y << " " << tumour[cell_x][cell_y].res << endl;
			}
			//while (tumour[cell_x][cell_y].dvr < 0);
			while ( !((WT_DIVISION && tumour[cell_x][cell_y].res == 0) || (MUT_DIVISION && tumour[cell_x][cell_y].res == 1)) );



			//cout << "Ntot=" <<Ntot << "\tNwt=" << Nwt << "\tNres=" << Nres << "\tWT_DIVISION=" << WT_DIVISION << "\tMUT_DIVISION=" << MUT_DIVISION << "\t cell mutations ->" << tumour[cell_x][cell_y].res << endl;

			// Compute birth rate of cell
			//compute_birthrate( tumour[cell_x][cell_y] , adv_model , &r_birth , _s );		
			
			// Update survival rate (which is a function of birth rate)
			//r_surv = r_birth * bdratio;
		

			// Update maximal birth and death rate of all cells 
			//if (r_birth > max_birth) max_birth = r_birth;



			// Cell divides
			volumetric_division(tumour , cell_x , cell_y , &Ntot , &Nres , &Nwt , &x_b , &y_b , radius);



			//cout << "\t\t (Nwt , Nres) -> (" << Nwt << "," << Nres << ")" << endl;


			/*
			// Ask if the cell dies (intrinsic death rate)
			if ( (tumour[cell_x][cell_y].dvr != -1) && (Ntot > 2) && (death) && (drand48() < (r_death/max_birth)) )
			{
				if (tumour[cell_x][cell_y].res == 1)
				{
					Nres -= 1;

					// Delete cell from tumour
					tumour[cell_x][cell_y].setDVR(-3);
					tumour[cell_x][cell_y].setRES(-3);
					tumour[cell_x][cell_y].setPGR(-3);
				}
				if (tumour[cell_x][cell_y].res == 0)
				{
					Nwt -= 1;

					// Delete cell from tumour
					tumour[cell_x][cell_y].setDVR(-2);
					tumour[cell_x][cell_y].setRES(-2);
					tumour[cell_x][cell_y].setPGR(-2);
				}

				// Size of tumour is reduced by 1
				Ntot -= 1;

				// Progress time variable
				t += 1.0/(max_birth * Ntot);
			}
			*/


		//if (iter == 1000000) break;



			/*
			if (iter%100 == 0)
			{
				sprintf(padded_iter, "%05d", iter);

				f.str("");
				f  << "./2D_DATA/maxSize=" << _maxsize << "_sweeps=" << sweeps << "_ut=" << _ut 
					<< "_BDratio=" << bdratio << "_s=" << _s << "/" << _seed << "/animation/" << padded_iter << ".csv";

				cout << f.str().c_str() << endl;


				animation_file.open(f.str().c_str());


				for (int i = 0; i < (2*radius); i++)
				{
					for (int j = 0; j < (2*radius); j++)
					{
						if (tumour[i][j].res > -1)
						{

							// Print cell coordinates to file
							animation_file << i << "," << j << "," << tumour[i][j].res << endl;

						}

					}
				}


				animation_file.close();
			}
			*/




		} while ( Ntot < _maxsize );		// Stop once we have the desired number of cells





		sprintf(padded_sweep, "%03d", sweep + 1);

		f.str("");
		f  << "./2D_DATA/maxSize=" << _maxsize << "_sweeps=" << sweeps << "_ut=" << _ut 
			<< "_BDratio=" << bdratio << "_s=" << _s << "/" << _seed << "/animation/" << padded_sweep << ".csv";


		animation_file.open(f.str().c_str());


		for (int i = 0; i < (2*radius); i++)
		{
			for (int j = 0; j < (2*radius); j++)
			{
				if (tumour[i][j].res > -1)
				{

					// Print cell coordinates to file
					animation_file << i << "," << j << "," << tumour[i][j].res << endl;

				}

			}
		}


		animation_file.close();


		if (!quiet)
		{
			cout << "\tSweep #" << sweep << "\t Ntot =" << Ntot << "\t Nwt =  " << Nwt << " \t Nres = " << Nres << endl;

		}

	}
	//while (Nwt > 0.01*_maxsize);		// Run simulation until WT cells make up only <1% of all cells

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









	//================== Open data files ==================//
/*
	ofstream tumour_file;
	f.str("");
	f << "./2D_DATA/maxSize=" << _maxsize << "_sweeps=" << sweeps << "_ut=" << _ut 
		<< "_BDratio=" << bdratio << "_s=" << _s << "/" << _seed << "/tumour.csv";
	tumour_file.open(f.str().c_str());


	ofstream mathematica_file;
	f.str("");
	f << "./DATA_surfaceGrowth2D/_NEIGHBOURHOOD=" << NEIGHBOURHOOD << "_s=" << _s << "_arisingtime=" << arising_time << "/maxsize="<< _maxsize << "_seed=" << _seed << "_NEIGHBOURHOOD=" << NEIGHBOURHOOD << "_s=" << _s << "_ut=" 
		<< _ut << "_ud=" << _ud << "_arisingtime=" << arising_time << "_death=" << death << "_divmodel=" << div_model << "_advmodel=" << adv_model << "/plot_mathematica.csv";
	mathematica_file.open(f.str().c_str());


	if (!quiet) cout << " " << endl;
	if (!quiet) cout << "\tCreated output files..." << endl;
*/




	// Cell memory -> fill in dead cells according to their previous (i.e. before dying) state
	for (int i = 0; i < (2*radius); i++)
	{
		for (int j = 0; j < (2*radius); j++)
		{

			// If cell has -2 drivers, then it used to be a WT cell, and if -3 drivers then it used to be a mutated cell
			if ( tumour[i][j].res == -2 ) tumour[i][j].res = 0;
			if ( tumour[i][j].res == -3 ) tumour[i][j].res = 1;

		}
	}


/*
	// Write tumour data to file
	tumour_file << "x coord, y coord, z coord, drivers, resistant, passengers" << endl;




	for (int i = 0; i < (2*radius); i++)
	{
		for (int j = 0; j < (2*radius); j++)
		{
	
			if ( tumour[i][j].res != -1 ) tumour_file << i << "," << j << ",0," << tumour[i][j].dvr << "," << tumour[i][j].res << "," << tumour[i][j].pgr << endl;

		}
	}




	//NversusT_file.close();
	tumour_file.close();


	if (!quiet) cout << "1 " << Nres << " " << endl;
	else cout << "1" << endl;
*/
	return 0;
}

















