#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

int find_value(int i, int j, int filter[3][3], int start, int n, int m, int mat[n][m], int depl, int factor, int vect[2]){
	int sum = 0;
	if(i != start && j != 0 && i < n - 1 && j < m - 1){
		sum += mat[i-1][j-1] * filter[0][0];
		sum += mat[i-1][j] * filter[0][1];
		sum += mat[i-1][j+1] * filter[0][2];
		sum += mat[i][j-1] * filter[1][0];
		sum += mat[i][j] * filter[1][1];
		sum += mat[i][j+1] * filter[1][2];
		sum += mat[i+1][j-1] * filter[2][0];
		sum += mat[i+1][j] * filter[2][1];
		sum += mat[i+1][j+1] * filter[2][2];
	}else if(i == start && j > 0 && j < m - 1){
		sum += vect[j-1] * filter[0][0];
		sum += vect[j] * filter[0][1];
		sum += vect[j+1] * filter[0][2];
		sum += mat[i][j-1] * filter[1][0];
		sum += mat[i][j] * filter[1][1];
		sum += mat[i][j+1] * filter[1][2];
		sum += mat[i+1][j-1] * filter[2][0];
		sum += mat[i+1][j] * filter[2][1];
		sum += mat[i+1][j+1] * filter[2][2];
	} else if(j == 0 && i > start && i < n - 1){
		sum += 0 * filter[0][0];
		sum += mat[i-1][j] * filter[0][1];
		sum += mat[i-1][j+1] * filter[0][2];
		sum += 0 * filter[1][0];
		sum += mat[i][j] * filter[1][1];
		sum += mat[i][j+1] * filter[1][2];
		sum += 0 * filter[2][0];
		sum += mat[i+1][j] * filter[2][1];
		sum += mat[i+1][j+1] * filter[2][2];
	}else if(i == start && j == 0){
		sum += 0 * filter[0][0];
		sum += vect[j] * filter[0][1];
		sum += vect[j+1] * filter[0][2];
		sum += 0 * filter[1][0];
		sum += mat[i][j] * filter[1][1];
		sum += mat[i][j+1] * filter[1][2];
		sum += 0 * filter[2][0];
		sum += mat[i+1][j] * filter[2][1];
		sum += mat[i+1][j+1] * filter[2][2];
	}else if(i == n - 1 && j > 0 && j < m - 1){
		sum += mat[i-1][j-1] * filter[0][0];
		sum += mat[i-1][j] * filter[0][1];
		sum += mat[i-1][j+1] * filter[0][2];
		sum += mat[i][j-1] * filter[1][0];
		sum += mat[i][j] * filter[1][1];
		sum += mat[i][j+1] * filter[1][2];
		sum += vect[j + m - 1] * filter[2][0];
		sum += vect[j + m] * filter[2][1];
		sum += vect[j + m + 1] * filter[2][2];
	}else if(i > start && i < n - 1 && j == m-1){
		sum += mat[i-1][j-1] * filter[0][0];
		sum += mat[i-1][j] * filter[0][1];
		sum += 0 * filter[0][2];
		sum += mat[i][j-1] * filter[1][0];
		sum += mat[i][j] * filter[1][1];
		sum += 0 * filter[1][2];
		sum += mat[i+1][j-1] * filter[2][0];
		sum += mat[i+1][j] * filter[2][1];
		sum += 0 * filter[2][2];
	}else if(i == n - 1 && j == m - 1){
		sum += mat[i-1][j-1] * filter[0][0];
		sum += mat[i-1][j] * filter[0][1];
		sum += 0 * filter[0][2];
		sum += mat[i][j-1] * filter[1][0];
		sum += mat[i][j] * filter[1][1];
		sum += 0 * filter[1][2];
		sum += vect [j + m - 1] * filter[2][0];
		sum += vect[j + m] * filter[2][1];
		sum += 0 * filter[2][2];
	}else if(i == start && j == m - 1){
		sum += vect[j - 1] * filter[0][0];
		sum += vect[j] * filter[0][1];
		sum += 0 * filter[0][2];
		sum += mat[i][j-1] * filter[1][0];
		sum += mat[i][j] * filter[1][1];
		sum += 0 * filter[1][2];
		sum += mat[i+1][j-1] * filter[2][0];
		sum += mat[i+1][j] * filter[2][1];
		sum += 0 * filter[2][2];
	}else if(j == 0 && i == n - 1){
		sum += 0 * filter[0][0];
		sum += mat[i-1][j] * filter[0][1];
		sum += mat[i-1][j+1] * filter[0][2];
		sum += 0 * filter[1][0];
		sum += mat[i][j] * filter[1][1];
		sum += mat[i][j+1] * filter[1][2];
		sum += 0 * filter[2][0];
		sum += vect[j + m] * filter[2][1];
		sum += vect[j + m + 1] * filter[2][2];
	}

	sum = ceil((double)sum / factor);
	sum += depl;

	if(sum > 255){
		sum = 255;
	}
	if(sum < 0){
		sum = 0;
	}
	return sum;
}

int main(int argc, char * argv[]) {
	if(argc != 4){
		printf("Lista de argumente incorecte\n");
		return 0;
	}
	int rank;
	int nProcesses;
	MPI_Init(&argc, &argv);
	MPI_Status status, status2, Istatus;
	MPI_Request request;

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nProcesses);

	char* topologie = strdup(argv[1]);
	int vecini[nProcesses - 1];
	
	FILE* f_top, *f_img, *f_stat;
	f_top = fopen(topologie, "r");

	int sobel[3][3] = {{1,0,-1},{2,0,-2}, {1,0,-1}};
	int mr[3][3] = {{-1,-1,-1}, {-1,9,-1}, {-1,-1,-1}};
	
	//tag cu care parintele trimite blocurile catre copii
	int sondaj = 1;
	//tag penru a trimite width si heght
	int wh = 5;
	//tag linii
	int linii_ss = 2;
	//tag reply copil->parinte cu bucata prelucrata
	int reply = 3;
	//tag finish
	int finish = 4;
	//tag statistica reply
	int stat = 8;
	//tag pentru linii
	int line_tag = 6;
	//tag pentru a trimite filtrul
	int tag_filtru = 7;
	
	int parinte = -1;
	int numarul_de_linii = 0;
	int i, j, nr, width, height;
	int lineCount = 0;
	int line_count = 0, rest = 0;
	int k = 0;

	char* token;
	//char* filter;
	char* line = calloc(100, sizeof(char));

	//citim topologia => fiecare proces isi afla vecinii
	while(fgets(line, 100, f_top) != NULL){
		if(lineCount == rank){
			token = strtok(line, " ");
			token = strtok(NULL, " ");
			while(token != NULL){
				vecini[k] = atoi(token);
				k++;
				token = strtok(NULL, " ");
			}
			break;
		}

		lineCount++;
	}

	int nrVecini = k;

	if(rank == 0){
		char* img_in = strdup(argv[2]);
		char* stat_out = strdup(argv[3]);
		f_img = fopen(img_in, "r");
		f_stat = fopen(stat_out, "w");
		
		//citim numarul de imagini 
		fgets(line, 100, f_img);
		int n = atoi(line);
		int zero = 0;

		//citim fiecare imagine
		while(n != 0){
			fgets(line, 100, f_img);
			char s[1] = {' '};
		
			//primul token = filtrul ce va fi aplicat
			token = strtok(line, s);
			char filter[strlen(token)];
			strcpy(filter, token);

			for(i = 0; i < nrVecini; i++){
				MPI_Send(filter, strlen(filter) , MPI_CHAR, vecini[i], tag_filtru, MPI_COMM_WORLD);
			}
			
			//al doilea token = imaginea nefiltrata
	   		token = strtok(NULL, s);
		   	char *original_img = strdup(token);
		    
		    //al treila token = imaginea filtrata
		    token = strtok(NULL, s);
		    char * filtered_img = strdup(token);
		   
		    //deschidem poza originala
			FILE *picture = fopen(original_img, "r");
			char *str = calloc(strlen(filtered_img) - 1, sizeof(char));
			strncpy(str, filtered_img, strlen(filtered_img) - 1);
		    FILE * out_picture = fopen(str, "w");
		    
		    char* line2 = calloc(100, sizeof(char));
		    char* token;

		    for(i = 0; i < 4; i++){
		    	fgets(line2, 100, picture);
		    	fputs(line2, out_picture);
		    	if(i == 2){
		    		token = strtok(line2, " ");
		    		width = atoi(token);
		    		token = strtok(NULL, " ");
		    		height = atoi(token);
		       	}
		    }

		    int v[2] = {width, height};
		    for(i = 0; i < nrVecini; i++){
		    	MPI_Send(v, 2 , MPI_INT, vecini[i], wh, MPI_COMM_WORLD);
		    }

		   int mat[height][width];
		   printf("0 incepe citirea matricii %s\n", original_img);

		   	for(i = 0; i < height; i++){
		   		for(j = 0; j < width; j++){
		   			fscanf(picture, "%d", &nr);
		   			mat[i][j] = nr;
		   		}
		   	}

		   	int vect[2*width];
		   	line_count = height / nrVecini;
			rest = height - nrVecini * line_count;

			//trimitem primului copil

			printf("0 trimite lui %d primului copil matricea\n", vecini[0]);
			MPI_Send(mat, height * width , MPI_INT, vecini[0], sondaj, MPI_COMM_WORLD);

			//trimitem liniile de start si de stop
			int start, stop;
			start = 0;
			stop = line_count;
			MPI_Send(&start, 1 , MPI_INT, vecini[0], linii_ss, MPI_COMM_WORLD);
			MPI_Send(&stop, 1 , MPI_INT, vecini[0], linii_ss, MPI_COMM_WORLD);

			for(i = 0; i < width; i++){
				vect[i] = 0;
				vect[i + width] = mat[line_count][i];
			}

			MPI_Send(vect, 2*width , MPI_INT, vecini[0], line_tag, MPI_COMM_WORLD);

			if(nrVecini > 1){
				//trimitem celorlalti copii
				for(i = 1; i < nrVecini - 1; i++){
					printf("0 trimite matricea lui %d\n", vecini[i]);
					MPI_Send(mat, height * width , MPI_INT, vecini[i], sondaj, MPI_COMM_WORLD);

					//trimite liniile de start si de stop
					start = i * line_count;
					stop = (i + 1) * line_count;
					MPI_Send(&start, 1 , MPI_INT, vecini[i], linii_ss, MPI_COMM_WORLD);
					MPI_Send(&stop, 1 , MPI_INT, vecini[i], linii_ss, MPI_COMM_WORLD);

					//trimite liniile de sus si de jos
					for(j = 0; j < width; j++){
						vect[j] = mat[i * line_count - 1][j];
						vect[j + width] = mat[(1 + i) * line_count][j];
					}

					MPI_Send(vect, 2*width , MPI_INT, vecini[i], line_tag, MPI_COMM_WORLD);
				}
				//trimitem ultimului copil
				printf("0 trimiteu ultimului copil %d \n", vecini[nrVecini - 1]);
				MPI_Send(mat, height * width , MPI_INT, vecini[nrVecini - 1], sondaj, MPI_COMM_WORLD);

				start = (nrVecini - 1) * line_count;
				stop = height;
				MPI_Send(&start, 1 , MPI_INT, vecini[nrVecini - 1], linii_ss, MPI_COMM_WORLD);
				MPI_Send(&stop, 1 , MPI_INT, vecini[nrVecini - 1], linii_ss, MPI_COMM_WORLD);
				
				for(i = 0; i < width; i++){
					vect[i + width] = 0;
					vect[i] = mat[line_count * (nrVecini - 1) - 1][i];
				} 

				MPI_Send(vect, 2*width, MPI_INT, vecini[nrVecini - 1], line_tag, MPI_COMM_WORLD);
			}

			int new_mat[height][width];
			int mat_fin[height][width];
			int k, size;

			//asteapta rezultatele prelucrarii de la copii
			int aux = nrVecini;
			while(aux > 0){
				MPI_Recv(new_mat, height*width, MPI_INT,MPI_ANY_SOURCE, reply, MPI_COMM_WORLD, &status);
				
				int idx = 0;
				for(i = 0; i < nrVecini; i++){
					if(status.MPI_SOURCE == vecini[i]){
						idx = i;
						break;
					}
				}

				int inf, sup;
				if (idx == 0){
					inf = 0;
					sup = line_count;
				}else if(idx == nrVecini - 1){
					inf = line_count*(nrVecini - 1);
					sup = height;
				}else{
					inf = idx*line_count;
					sup = (idx + 1) * line_count;
				}

				for(i = inf; i < sup; i++){
					for(j = 0; j < width; j++){
						mat_fin[i][j] = new_mat[i][j];
					}
				}

				printf("0 primeste de la copilul %d\n", vecini[idx]);
				aux--;
			}

			for(i = 0; i < height; i++){
				for(j = 0; j < width; j++){
					fprintf(out_picture, "%d\n", mat_fin[i][j]);
				}
			}

		    fclose(picture);
		    fclose(out_picture);
		    n--;
		    printf("\n\n\n");
		}

		//trimitem mesajul de finish
		for(i = 0; i < nrVecini; i++){
			char str[4] = "gata";
			MPI_Send(str, 1, MPI_INT, vecini[i], finish, MPI_COMM_WORLD);
		}

		int statistica[nProcesses];
		for(i = 0; i < nProcesses; i++){
			statistica[i] = 0;
		}

		int stat_aux[nProcesses];
		for(i = 0; i < nrVecini; i++){
			MPI_Recv(stat_aux, nProcesses, MPI_INT, vecini[i], stat, MPI_COMM_WORLD, &status);
			for(j = 0; j < nProcesses; j++){
				if(statistica[j] == 0 && stat_aux[j] != 0){
					statistica[j] = stat_aux[j];
				}
			}
		}
		statistica[rank] = 0;

		for(i = 0; i < nProcesses; i++){
			fprintf(f_stat, "%d: %d\n", i, statistica[i]);
		}

		fclose(f_img);
		fclose(f_stat);
	} else{
		while(1){
			//filtrul de aplicat
			int statistica[nProcesses];
			for(i = 0; i < nProcesses; i++){
				statistica[i] = 0;
			}

			char *filter_aux =  calloc(15, sizeof(char));
			MPI_Recv(filter_aux, 15, MPI_CHAR, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			
			if(status.MPI_TAG == finish){
				for(i = 0; i < nrVecini; i++){
					if(vecini[i] != parinte){
						MPI_Send(filter_aux, 1, MPI_INT, vecini[i], finish, MPI_COMM_WORLD);
					}
				}

				if(nrVecini == 1){
					statistica[rank] = numarul_de_linii;
					MPI_Send(statistica, nProcesses, MPI_INT, parinte, stat, MPI_COMM_WORLD);
				}else{
					int stat_aux[nProcesses];
					for(i = 0; i < nrVecini; i++){
						if(vecini[i] != parinte){
							MPI_Recv(stat_aux, nProcesses, MPI_INT, vecini[i], stat, MPI_COMM_WORLD, &status);
							for(j = 0; j < nProcesses; j++){
								if(statistica[j] == 0 && stat_aux[j] != 0){
									statistica[j] = stat_aux[j];
								}
							}
						}
					}
					statistica[rank] = 0;
					MPI_Send(statistica, nProcesses, MPI_INT, parinte, stat, MPI_COMM_WORLD);
				}

				break;
			}

			parinte = status.MPI_SOURCE;

			char filter[strlen(filter_aux)];
			strcpy(filter, filter_aux);
			
			for(i = 0; i < nrVecini; i++){
				if(vecini[i] != parinte){
					MPI_Send(filter, strlen(filter) , MPI_CHAR, vecini[i], tag_filtru, MPI_COMM_WORLD);
				}
			}

			//height si width imagine
			int v[2];
			MPI_Recv(v, 2, MPI_INT, parinte, wh, MPI_COMM_WORLD, &status);
			
			for(i = 0; i < nrVecini; i++){
				if(vecini[i] != parinte){
					MPI_Send(v, 2 , MPI_INT, vecini[i], wh, MPI_COMM_WORLD);
				}
			}
			
			width = v[0];
			height = v[1];
			int start = 0;
			int stop = 0;

			//asteptam matricea, linia de start si stop si liniile de bordura
			int *mat_recv_aux = calloc(width*height, sizeof(int));
			int vect[2*width];
		 
			MPI_Recv(mat_recv_aux, width*height, MPI_INT, parinte, sondaj, MPI_COMM_WORLD, &status);
			MPI_Recv(&start, 1, MPI_INT, parinte, linii_ss, MPI_COMM_WORLD, &status2);
			MPI_Recv(&stop, 1, MPI_INT, parinte, linii_ss, MPI_COMM_WORLD, &status2);
			MPI_Recv(vect, width*2, MPI_INT, parinte, line_tag, MPI_COMM_WORLD, &status2);
		
			
			//daca procesul nu e frunza trebuie sa imparta copiilor ce a primit
			if(nrVecini > 1){
				line_count = (stop - start) / (nrVecini - 1);
			}else{
				numarul_de_linii += stop - start;
			}

			int line_start = 0;
			int line_stop = 0;
			int vector[2*width];
			
			int copy[height][width];

			for(i = 0; i < height; i++){
				for(j = 0; j < width; j++){
					copy[i][j] = mat_recv_aux[i*width + j];
				}
			}
			
			//trimitem primului copil (daca nu este parintele)
			//daca primul vecin este si parinte am doua posibilitati: frunza (nu exista alti vecini, deci nu trimit nimanui) sau vecinul 1 este primul copil
			if(vecini[0] != parinte){
				MPI_Send(copy, height * width , MPI_INT, vecini[0], sondaj, MPI_COMM_WORLD);
				line_start = start;
				line_stop = start + line_count;
				MPI_Send(&line_start, 1 , MPI_INT, vecini[0], linii_ss, MPI_COMM_WORLD);
				MPI_Send(&line_stop, 1 , MPI_INT, vecini[0], linii_ss, MPI_COMM_WORLD);
				
				for(j = 0; j < width; j++){
					vector[j] = vect[j];
					vector[j + width] = copy[start + line_count][j];
				}
				MPI_Send(vector, 2*width , MPI_INT, vecini[0], line_tag, MPI_COMM_WORLD);
				
			}else if (nrVecini > 1){
				printf("procesul %d trimite lui %d %d linii \n", rank, vecini[1], line_count);
				MPI_Send(copy, height * width , MPI_INT, vecini[1], sondaj, MPI_COMM_WORLD);
				line_start = start;
				line_stop = start + line_count;
				MPI_Send(&line_start, 1 , MPI_INT, vecini[1], linii_ss, MPI_COMM_WORLD);
				MPI_Send(&line_stop, 1 , MPI_INT, vecini[1], linii_ss, MPI_COMM_WORLD);
				
				for(j = 0; j < width; j++){
					vector[j] = vect[j];
					
					if(nrVecini == 2){
						vector[j + width] = vect[j + width];
					}else{
						vector[j + width] = copy[line_stop][j];
					}
				}
				
				MPI_Send(vector, 2*width , MPI_INT, vecini[1], line_tag, MPI_COMM_WORLD);
				
				//trimitem celorlalti copii
				for(i = 2; i < nrVecini - 1; i++){
					if(vecini[i] != parinte){
						printf("procesul %d trimite lui %d \n", rank, vecini[i]);
						MPI_Send(copy, height * width , MPI_INT, vecini[i], sondaj, MPI_COMM_WORLD);

						line_start = start + (i - 1) * line_count;
						line_stop = start + i * line_count ;
						MPI_Send(&line_start, 1 , MPI_INT, vecini[i], linii_ss, MPI_COMM_WORLD);
						MPI_Send(&line_stop, 1 , MPI_INT, vecini[i], linii_ss, MPI_COMM_WORLD);

						for(j = 0; j < width; j++){
							vector[j] = copy[line_start - 1][j];
							vector[j + width] = copy[line_stop][j];
						}
						
						MPI_Send(vector, 2*width , MPI_INT, vecini[i], line_tag, MPI_COMM_WORLD);
					}
				}

				//trimitem ultimului copil
				if(nrVecini >= 3){
					if(vecini[nrVecini - 1] != parinte){
						printf("procesul %d trimite lui %d \n", rank, vecini[nrVecini - 1]);
						MPI_Send(copy, height * width , MPI_INT, vecini[nrVecini - 1], sondaj, MPI_COMM_WORLD);

						line_start = start + line_count * (nrVecini - 2);
						line_stop = stop;
						MPI_Send(&line_start, 1 , MPI_INT, vecini[nrVecini - 1], linii_ss, MPI_COMM_WORLD);
						MPI_Send(&line_stop, 1 , MPI_INT, vecini[nrVecini - 1], linii_ss, MPI_COMM_WORLD);
						
						for(j = 0; j < width; j++){
							vector[j] = copy[line_start - 1][j];
							vector[j + width] = vect[j + width];
						} 

						MPI_Send(vector, 2*width , MPI_INT, vecini[nrVecini - 1], line_tag, MPI_COMM_WORLD);
				}
			}
		}

			//int *mat_aux_copil = calloc(height * width, sizeof(int));
			//int* new_mat = calloc(height * width, sizeof(int));
			int new_mat[height][width];
			int mat_aux_copil[height][width];

			if(nrVecini == 1){
				for(i = start; i < stop; i++){
					for(j = 0; j < width; j++){
						if(strcmp(filter, "sobel") == 0){
							mat_aux_copil[i][j] = find_value(i, j, sobel, start, stop, width, copy, 127, 1, vect);
			   			}else{
			   				mat_aux_copil[i][j] = find_value(i, j, mr, start, stop, width, copy, 0, 1, vect);
			   			}
					}
				}

				printf("frunza %d trimite parintelui %d\n", rank, parinte);
				MPI_Send(mat_aux_copil, height*width, MPI_INT, parinte, reply, MPI_COMM_WORLD);
			}

			if(nrVecini != 1){
				int aux = nrVecini - 1;
				while(aux > 0){
					MPI_Recv(mat_aux_copil, height*width, MPI_INT, MPI_ANY_SOURCE, reply, MPI_COMM_WORLD, &status);
					int idx = 0;
					for(i = 0; i < nrVecini; i++){
						if(status.MPI_SOURCE == vecini[i]){
							idx = i;
							break;
						}
					}

					int inf, sup;
					if ((idx == 0 && parinte == vecini[1]) || (idx == 1 && parinte == vecini[0])){
						inf = start;
						sup = start + line_count;
					}else if(idx == nrVecini - 1){
						inf = start + line_count*(nrVecini - 2);
						sup = stop;
					}else{
						inf = start + (idx - 1)*line_count;
						sup = start + idx * line_count;
					}

					for(k = inf; k < sup; k++){
						for(j = 0; j < width; j++){
							new_mat[k][j] = mat_aux_copil[k][j];
						}
					}
					printf("parintele %d primeste de la copilul %d\n", rank, vecini[idx]);
					aux--;
				}

				printf("copilul %d trimite parintelui %d\n", rank, parinte);
				MPI_Send(new_mat, height*width, MPI_INT, parinte, reply, MPI_COMM_WORLD);
			}

		}
	}

	fclose(f_top);
	printf("Bye from %d\n", rank);
	MPI_Finalize();
	return 0;
}