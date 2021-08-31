#include <iostream>
#include <mpi.h>

class NeighbourGrid 
{
public:

    char*** get_grid() { return grids; }

    NeighbourGrid(char **grids_[9], const int& num_of_rows_, const int& num_of_cols_)
    :num_of_rows(num_of_rows_), num_of_cols(num_of_cols_), grids(grids_)
    {
        for (int i = 0; i < num_of_rows; i++)
            for (int j = 0; j < num_of_cols; j++)
                if (grids[4][i][j] == '*')
                    alive_cells++;

        update_state();
    }

    ~NeighbourGrid()
    {
        for (int k = 0; k < 9; k++)
        {
            for (int i = 0; i < num_of_rows; i++)
            {
                delete[] grids[k][i];
            }

            delete[] grids[k];

        }

        grids = NULL;
    }

    void take_one_step()
    {
        int num_of_neighbours = 0;
        alive_cells = 0;

        for (int i = 0; i < num_of_rows; i++)
        {
            for (int j = 0; j < num_of_cols; j++)
            {
                num_of_neighbours = get_neighbours(i, j);

                //  Any live cell with two or three live neighbours survives.
                if (grids[4][i][j] == '*' && (num_of_neighbours == 2 || num_of_neighbours == 3))
                {
                    //No change is necessary.
                    alive_cells++;
                }

                //Any dead cell with three live neighbours becomes a live cell.
                else if (grids[4][i][j] == '-' && num_of_neighbours == 3)
                {

                    if (grid_state != 1) //If overpopulated no birth will occur.
                    {
                        grids[4][i][j] = '*';
                        alive_cells++;
                    }
                }

                //All other live cells die in the next generation. Similarly, all other dead cells stay dead.
                else if (grids[4][i][j] == '*')
                {
                    if (grid_state == -1)//If underpopulated no death will occur.
                    {
                        //No change is necessary.
                        alive_cells++;
                    }

                    else
                    {
                        grids[4][i][j] = '-';
                    }
                }

                //The dead will stay dead so no change is necessary.

            }
        }

        update_state();
    }

private:

    int get_neighbours(const int& i, const int& j)
    {
        int neighbours_num = 0;

        for (int k = 0; k < 9; k++)
            if (grids[k][i][j] == '*' && k != 4)
                neighbours_num++;

        return neighbours_num;
    }

    void update_state()
    {
        if (alive_cells / ((double)num_of_cols * num_of_rows) >= ((double)0.75))
            grid_state = 1; //overpopulated

        else if (alive_cells / ((double)num_of_cols * num_of_rows) <= ((double)0.15))
            grid_state = -1; //underpopulated

        else
            grid_state = 0; //neither
    }

    int num_of_rows, num_of_cols, alive_cells = 0;
    int grid_state; // -1 = underpopulated, 0 = neither, 1 = overpopulated.
    char ***grids;

};

void exchange_info(char ***Grids, int neighbour_ranks[9], int num_of_rows, int num_of_cols)
{
    //The function responsible for sending and receiving neighbour information.

    MPI_Request request0 = MPI_REQUEST_NULL, request1 = MPI_REQUEST_NULL;
    MPI_Request req_arr[2];
    int res;

    for (int k = 0; k < 9; k++) //k = 0 send to upper left, k = 1 send to up ...
    {

        if (k != 4)
        {

            for (int i = 0; i < num_of_rows; i++)
            {
                MPI_Isend(Grids[4][i], num_of_cols, MPI_CHAR, neighbour_ranks[k], neighbour_ranks[4], MPI_COMM_WORLD, &request0);
                MPI_Irecv(Grids[8 - k][i], num_of_cols, MPI_CHAR, neighbour_ranks[8 - k], neighbour_ranks[8 - k], MPI_COMM_WORLD, &request1);
                //Suppose k = 0 then it sends to upper left neighbours_rank[0] and,
                //it receives from lower right which is neighbours_rank[8-k]
                //Think neighbours_rank as 
                //{up_left,  up,  up_right,
                //    left,  mid, right,
                // low_left, low, low_right}
                req_arr[0] = request0, req_arr[1] = request1;
                res = MPI_Waitall(2, req_arr, MPI_STATUSES_IGNORE);

                if (res != MPI_SUCCESS)
                {
                    std::cout << "An error occured in func. exchange_info at rank: " << neighbour_ranks[4] << " loop: " << k << std::endl;
                }
            }

        }

        MPI_Barrier(MPI_COMM_WORLD);
    }

}

int main(int *argc, char *argv)
{
    int width = -1, height = -1;
    int commsize, rank;
    int procs_per_row = 0, procs_per_col = 0;
    int loop = 0;
    MPI_Datatype viewtype;
    MPI_File fh;

    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &commsize);

    if(rank == 0)
        std::cin >> width >> height >> procs_per_col >> procs_per_row >> loop; //Only the 0th process can take inputs.
                                                                       //rank == procs_per_row * procs_per_col must be true.

    double t1, t2;
    t1 = MPI_Wtime();

    MPI_Bcast(&width, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&height, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&procs_per_col, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&procs_per_row, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&loop, 1, MPI_INT, 0, MPI_COMM_WORLD);

    //procs_rooted = processes_rooted.

    int num_of_cols = (width / procs_per_col), num_of_rows = (height / procs_per_row); 
    //Calculate number of columns and rows in each grid.

    int global_start_col = (rank % procs_per_col)*num_of_cols, global_start_row = (rank/procs_per_col)*num_of_rows; 
    //Calculate the upper left coordinates of each grid with respect to the input file.

    int global_end_col = (global_start_col + num_of_cols - 1), global_end_row = (global_start_row + num_of_rows - 1); 
    //Calculate the lower right coordinates of each grid with respect to the input file (not neccesary when calling function).

    int input_size[2] = { num_of_rows, num_of_cols }, start_indexes[2] = { global_start_row, global_start_col }, 
        global_size[2] = {height, width+2};//Neccesary values for the MPI_File_set_view function.


    MPI_File_open(MPI_COMM_WORLD, "input.txt", MPI_MODE_RDWR, MPI_INFO_NULL, &fh);
    // DOSYAYI AÇABILMESI IÇIN EXE DOSYASIYLA AYNI YERDE DOSYANIN BULUNMASI GEREKIYOR YOKSA KOD TRASH VALUE BASIYOR YA DA DEADLOCK YIYOR.
    
    char **grid = new char*[num_of_rows];

    for (int j = 0; j < num_of_rows; j++)
        grid[j] = new char[num_of_cols];

    MPI_Type_create_subarray(2, global_size, input_size, start_indexes, MPI_ORDER_C, MPI_CHAR, &viewtype);
    MPI_Type_commit(&viewtype);
    MPI_File_set_view(fh, 0, MPI_CHAR, viewtype, "native", MPI_INFO_NULL); 
    // Sets the view such that every process will read it's own grid from the file.

    for (int i = 0; i < num_of_rows; i++)
    {
        MPI_File_read(fh, grid[i], num_of_cols, MPI_CHAR, MPI_STATUS_IGNORE);
        //Reading row by row since each row is stored discontinuously.
    }
    MPI_File_close(&fh);

    //PREPEARING FOR RULE APPLICATION.

    char **Grids[9];
    
    for (int k = 0; k < 9; k++)
    {
        Grids[k] = new char*[num_of_rows];

        for (int i = 0; i < num_of_rows; i++)
            Grids[k][i] = new char[num_of_cols];

    }

    for (int i = 0; i < num_of_rows; i++)
        for (int j = 0; j < num_of_cols; j++)
            Grids[4][i][j] = grid[i][j];

    for (int i = 0; i < num_of_rows; i++)
        delete[] grid[i];

    delete[] grid;

    //Map the ranks from 1d to 2d cordinates than calculate it's neighbours that way.
    int x_cord = rank % procs_per_col, y_cord = rank / procs_per_col;
    int up_left = ((y_cord + procs_per_row - 1) % procs_per_row) * procs_per_col + ((x_cord + procs_per_col - 1) % procs_per_col),
        up = ((y_cord + procs_per_row - 1) % procs_per_row) * procs_per_col + x_cord,
        up_right = ((y_cord + procs_per_row - 1) % procs_per_row) * procs_per_col + ((x_cord + procs_per_col + 1) % procs_per_col),
        left = y_cord * procs_per_col + ((x_cord + procs_per_col - 1) % procs_per_col),
        right = y_cord * procs_per_col + ((x_cord + procs_per_col + 1) % procs_per_col),
        low_left = ((y_cord + procs_per_row + 1) % procs_per_row) * procs_per_col + ((x_cord + procs_per_col - 1) % procs_per_col),
        low = ((y_cord + procs_per_row + 1) % procs_per_row) * procs_per_col + x_cord,
        low_right = ((y_cord + procs_per_row + 1) % procs_per_row) * procs_per_col + ((x_cord + procs_per_col + 1) % procs_per_col);
        
    int neighbour_ranks[9] = { up_left, up, up_right, left, rank, right, low_left, low, low_right };

    //START RULE APPLICATION.

    NeighbourGrid all_grids(Grids, num_of_rows, num_of_cols);

    MPI_Barrier(MPI_COMM_WORLD);

    double total_time = 0;

    while(loop--)
    {
        double t1, t2;
        t1 = MPI_Wtime();
        exchange_info(all_grids.get_grid(), neighbour_ranks, num_of_rows, num_of_cols);
        MPI_Barrier(MPI_COMM_WORLD);
        t2 = MPI_Wtime();
        total_time += t2 - t1;

        all_grids.take_one_step();
        MPI_Barrier(MPI_COMM_WORLD);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    //PRINT OUTPUT ON THE FILE output.txt.

    MPI_File nh;
    MPI_Datatype new_viewtype;

    if ((rank + 1) % procs_per_col == 0)
        input_size[1]++; //For printing an extra new line.

    global_size[1]--;

    MPI_File_open(MPI_COMM_WORLD, "output.txt", MPI_MODE_RDWR, MPI_INFO_NULL, &nh);

    MPI_Type_create_subarray(2, global_size, input_size, start_indexes, MPI_ORDER_C, MPI_CHAR, &new_viewtype);
    MPI_Type_commit(&new_viewtype);
    MPI_File_set_view(nh, 0, MPI_CHAR, new_viewtype, "native", MPI_INFO_NULL);

    for (int i = 0; i < num_of_rows; i++)
    {
        MPI_File_write(nh, Grids[4][i], num_of_cols, MPI_CHAR, MPI_STATUS_IGNORE);

        if ((rank + 1) % procs_per_col == 0)
            MPI_File_write(nh, &("\n"), 1, MPI_CHAR, MPI_STATUS_IGNORE);
    }

    MPI_File_close(&nh);

    t2 = MPI_Wtime();
    if (rank == 0)
        std::cout << "Time taken to process: " << t2-t1 << std::endl;

    if (rank == 0)
        std::cout << "Time taken on exchange_info: " << total_time << std::endl;


    MPI_Finalize();

    return 0;
}
