#include <stdio.h>

#define ITER 10
#define ROWS 3
#define COLS 9

int main( int argc, char *argv[] ) {
	if (argc != 3) {
		printf("Usage: ./diff file_output_ok file_output_check \n");
		return 1;
	}
	
	FILE* file_ok = 		fopen(argv[1], "r");
	FILE* file_check = 	fopen(argv[2], "r");

	int ok_ints[ITER][ROWS][COLS];
	int check_ints[ITER][ROWS][COLS];

	int i, j, iters;
	char text[7];
	for (iters = 0; iters < ITER; ++iters) {
		for (i = 0; i < ROWS; ++i) {
			
			fscanf(file_ok,			"%s", text);
			fscanf(file_check,	"%s", text);
			for (j = 0; j < COLS; ++j) {
				fscanf(file_ok, 		"%d", &ok_ints[iters][i][j]);
				fscanf(file_check,	"%d", &check_ints[iters][i][j]);
			}
		}
	}
	//To this point, ok_ints and check_ints hold all the integers
	//for the correctness check, which is done below.
	
	int x, y, k, z;
	int found;
	for (iters = 0; iters < ITER; ++iters) {
		found = 0;
		for (x = 0; x < ROWS; ++x) {
			for (y = 0; y < ROWS; ++y) {
				//Check if coordinates and angles are the same and score 
				//doesn't differ too much
				if (check_ints[iters][x][3] == ok_ints[iters][y][3] 		&& 
						check_ints[iters][x][4] == ok_ints[iters][y][4] 		&&
						check_ints[iters][x][5] == ok_ints[iters][y][5] 		&& 
						check_ints[iters][x][6] == ok_ints[iters][y][6] 		&& 
						check_ints[iters][x][7] == ok_ints[iters][y][7] 		&& 
						check_ints[iters][x][8] == ok_ints[iters][y][8] 		&& 
						check_ints[iters][x][2] <= ok_ints[iters][y][2] + 2 &&
						check_ints[iters][x][2] >= ok_ints[iters][y][2]	- 2) {
							
							++found;
				}
			}
		}
		if (found != 3) {
			printf("\033[31;1m  WRONG! Differences found at solution: %d\033[0m\n", iters+1);
			return 1;
		}		
	}
	printf("\033[32;1m  GOOD! Files don't differ!  \033[0m\n");
	
	/*
	for (k = 0; k < ITER; ++k) {
		for (x = 0; x < ROWS; ++x) {
			for (y = 0; y < COLS; ++y) {
				printf("%d %d ", ok_ints[k][x][y], check_ints[k][x][y]);
			}
		}
		printf("\nline ended\n");
	}*/
}
