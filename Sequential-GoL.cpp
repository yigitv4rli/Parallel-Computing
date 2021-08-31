#define _CRT_SECURE_NO_DEPRECATE
#include<iostream>
#include<time.h>

class Grid
{
public:

	static int procs_per_col, procs_per_row; // procs_per_col = the number of processes at every col, procs_per_row = similarly..
											 //(I know it is not it's exact meaning)
	static int map_width, map_height, grid_width, grid_height;


	Grid(int grid_num, char** _map, int ** _num_map);
	~Grid();

	void update_map();
	void update_num_map();
	void zero_num_map();

private:

	int neighbours[9];
	int grid_num;
	char **map;
	int **num_map;
	int state;

	void update_state();

};

Grid::Grid(int _grid_num, char** _map, int** _num_map) : map(_map), grid_num(_grid_num), num_map(_num_map)
{
	int grid_cord_y = grid_num / procs_per_col;
	int grid_cord_x = grid_num % procs_per_col;
	int num = 0;

	neighbours[0] = ((grid_cord_x + procs_per_col - 1) % procs_per_col) + procs_per_col * ((grid_cord_y + procs_per_row - 1) % procs_per_row);
	neighbours[1] = grid_cord_x + procs_per_col * ((grid_cord_y + procs_per_row - 1) % procs_per_row);
	neighbours[2] = ((grid_cord_x + 1)%procs_per_col) + procs_per_col * ((grid_cord_y + procs_per_row - 1) % procs_per_row);
	neighbours[3] = ((grid_cord_x + procs_per_col - 1) % procs_per_col) + procs_per_col * grid_cord_y;
	neighbours[4] = grid_num;
	neighbours[5] = ((grid_cord_x + 1) % procs_per_col) + procs_per_col * grid_cord_y;
	neighbours[6] = ((grid_cord_x + procs_per_col - 1) % procs_per_col) + procs_per_col * ((grid_cord_y + 1) % procs_per_row);
	neighbours[7] = grid_cord_x + procs_per_col * ((grid_cord_y + 1) % procs_per_row);
	neighbours[8] = ((grid_cord_x + 1) % procs_per_col) + procs_per_col * ((grid_cord_y + 1) % procs_per_row);

	update_state();
}

void Grid::update_state()
{
	int grid_cord_y = grid_num / procs_per_col;
	int grid_cord_x = grid_num % procs_per_col;
	int num = 0;

	for (int _i = 0, i = grid_cord_y * grid_height; _i < grid_height; i++, _i++)
		for (int _j = 0, j = grid_cord_x * grid_width; _j < grid_width; _j++, j++)
			if (map[i][j] == '*')
				num++;

	if ((num / ((double)grid_width * grid_height)) >= ((double)0.75))
		state = 1; //overpopulated

	else if ((num / ((double)grid_width * grid_height)) <= ((double)0.15))
		state = -1; //underpopulated

	else
		state = 0; //neither
}

Grid::~Grid()
{
}

void Grid::update_num_map()
{
	int curr_grid = 0;
	int curr_grid_x = 0;
	int curr_grid_y = 0;
	int own_grid = neighbours[4];
	int own_grid_x = (own_grid % procs_per_col) * grid_width;
	int own_grid_y = (own_grid / procs_per_col) * grid_height;

	for (int k = 0; k < 9; k++)
	{
		if (k == 4)
			continue;


		curr_grid = neighbours[k];
		curr_grid_x = (curr_grid % procs_per_col) * grid_width;
		curr_grid_y = (curr_grid / procs_per_col) * grid_height;

		for (int i = 0; i < grid_height; i++)
			for (int j = 0; j < grid_width; j++)
				if (map[i + own_grid_y][j + own_grid_x] == '*')
					num_map[i + curr_grid_y][j + curr_grid_x] += 1;
	}
}

void Grid::update_map()
{
	int num_of_neighbours = 0;
	int own_grid = neighbours[4];
	int own_grid_x = (own_grid % procs_per_col) * grid_width;
	int own_grid_y = (own_grid / procs_per_col) * grid_height;

	for (int _i = 0, i = own_grid_y; _i < Grid::grid_height; _i++, i++)
	{
		for (int _j = 0, j = own_grid_x; _j < Grid::grid_width; _j++, j++)
		{
			num_of_neighbours = num_map[i][j];

			//  Any live cell with two or three live neighbours survives.
			if (map[i][j] == '*' && (num_of_neighbours == 2 || num_of_neighbours == 3))
			{
				//No change is necessary.
			}

			//Any dead cell with three live neighbours becomes a live cell.
			else if (map[i][j] == '-' && num_of_neighbours == 3)
			{

				if (state != 1) //If overpopulated no birth will occur.
				{
					map[i][j] = '*';
				}
			}

			//All other live cells die in the next generation. Similarly, all other dead cells stay dead.
			else if (map[i][j] == '*')
			{
				if (state == -1)//If underpopulated no death will occur.
				{
					//No change is necessary.
				}

				else
				{
					map[i][j] = '-';
				}
			}

			//The dead will stay dead so no change is necessary.
		}
	}

	update_state();
}

void Grid::zero_num_map()
{
	for (int i = 0; i < map_height; i++)
		for (int j = 0; j < map_width; j++)
			num_map[i][j] = 0;
}

int Grid::procs_per_col = 0, Grid::procs_per_row = 0, Grid::map_width = 0, Grid::map_height = 0;
int Grid::grid_width = 0, Grid::grid_height = 0;

int main()
{

	int loop;
	
	std::cin >> Grid::map_width >> Grid::map_height >> Grid::procs_per_col >> Grid::procs_per_row >> loop;
	Grid::grid_height = Grid::map_height / Grid::procs_per_row; 
	Grid::grid_width = Grid::map_width / Grid::procs_per_col;

	clock_t t1, t2;

	t1 = clock();
	FILE* incoming = fopen("input.txt", "r");

	if (incoming == NULL) {
		perror("Unable to open file");
		exit(1);
	}

	//Creation of grid
	char** grid = new char* [Grid::map_height];
	for (int i = 0; i < Grid::map_height; i++) {
		grid[i] = new char[Grid::map_width];
	}

	int row = 0;
	while (fgets(grid[row], Grid::map_width + 2, incoming)) {
		row++;
	}

	fclose(incoming);

	Grid **grid_arr = new Grid*[Grid::procs_per_col * Grid::procs_per_row];
	
	int** neighbour_count_map = new int*[Grid::map_height];

	for (int i = 0; i < Grid::map_height; i++)
		neighbour_count_map[i] = new int[Grid::map_width];
	
	for (int i = 0; i < Grid::procs_per_col * Grid::procs_per_row; i++)
		grid_arr[i] = new Grid(i, grid, neighbour_count_map);
	

	while (loop--)
	{
		grid_arr[0]->zero_num_map();

		for (int i = 0; i < Grid::procs_per_col * Grid::procs_per_row; i++)
		{
			grid_arr[i]->update_num_map();
		}

		for (int i = 0; i < Grid::procs_per_col * Grid::procs_per_row; i++)
		{
			grid_arr[i]->update_map();
		}
	}

	FILE* outgoing = fopen("output.txt", "w");

	if (outgoing == NULL) {
		std::cout << "Unable to create file" << '\n';
		exit(EXIT_FAILURE);
	}

	for (int i = 0; i < Grid::map_height; i++) {
		fprintf(outgoing, "%s", grid[i]);
	}

	fclose(outgoing);

	t2 = clock();

	std::cout << "Time taken:" << ((float)(t2 - t1)) / CLOCKS_PER_SEC << std::endl;
	return 0;
}
