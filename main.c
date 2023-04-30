#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>


//Population
#define POPULATION 100

// The maximum genes for one ant
#define MAX_GENES 20

// The four ways the ants can go, these are used for the ant's gene creation
const char GENES[] = {'l','r','u','d'};

// Mutation genes if the given gene is mutated select from these genes instead.
const char MUT_GENES_L[] = {'r','u','d'};
const char MUT_GENES_R[] = {'l','u','d'};
const char MUT_GENES_U[] = {'r','l','d'};
const char MUT_GENES_D[] = {'r','u','l'};

// Default food location is at (10,10), this can be change by the user
int gFood_x = 10;
int gFood_y = 10;

// Helper variable to check for the perfect ant
int gPerfect_ant = 0;

// Helper varaible to check for the generations
int gGenerations = 0;

// Maximum generation to run.
int gMax_generation = 50;

// Max is 100% then all child will mutate a random gene
int gMutationprob = 0;

//Struct for an ant
struct Ant
{
    char gene[MAX_GENES];
    int fitness;
};

// Log to file

char *filename = "ants_log.txt";
FILE *gfp;


// Gives a random number between min and max 
int random_num_gen(int min, int max){
    return (rand() % (max - min + 1)) + min;
}

// Gives random genes for each ant given the maximum genes possible for one.
void rand_gene(char *gene_arr){
    // The start and end for the random generation, this is based on the genes
    int start = 0;
    int end = (sizeof(GENES) / sizeof(char)) - 1;
    int random_number;

    // Generate the genes randomly for each gene
    for(int z = 0; z < MAX_GENES; z++){
        random_number = random_num_gen(start,end);
        gene_arr[z] = GENES[random_number];
    }

    return;
}

// Get fitness
int get_fitness(char *gene_arr){
    
    // Distance from food coords, temporary needed to calculate fitness;
    int tmp_dist_x, tmp_dist_y;

    // Fitness for x and y separately
    int fitness_x, fitness_y;

    // Combined fitness
    int fitness;

    // Ant coordinants from genes.
    int tmp_x = 0;
    int tmp_y = 0;

    // Calculate coordinates for the ant given the genes
    for(int xy = 0; xy < MAX_GENES; xy++){
        switch (gene_arr[xy])
        {
        case 'l':
            tmp_x--;
            break;
        case 'r':
            tmp_x++;
            break;
        case 'u':
            tmp_y++;
            break;
        case 'd':
            tmp_y--;
            break;
        default:
            printf("Illegal character in gene!\n");
            break;
        }
    }

    // Distance from the food
    tmp_dist_x = abs(gFood_x-tmp_x);
    tmp_dist_y = abs(gFood_y-tmp_y);

    // The fitness is basically the ant's distance from the food, thus the lower the fitnes the closer the ant is to the food
    fitness_x = tmp_dist_x * 10;
    fitness_y = tmp_dist_y * 10;
    fitness = fitness_x + fitness_y;

    return fitness;
}

// "Map check" since the maximum genes limit the reachable food.
// The food should be reachable for ants
void map_check(){
    //Calculate the max reach, considering negative coords.
    int gFood_x_local = gFood_x;
    int gFood_y_local = gFood_y;
    int max_to_reach = (int) sqrt((double)(gFood_x_local*=gFood_x_local)) + sqrt((double)(gFood_y_local*=gFood_y_local));
    if(max_to_reach > MAX_GENES)
    {
        printf("With the given food coordinents the ants wont reach the food!\n");
        printf("The maximum gense to move that the ants have is: %d\n", MAX_GENES);
        printf("The food is : %d moves away from the start!\n", max_to_reach);
        if(gfp != NULL)
        {
            fprintf(gfp,"The food isn't reachable by the ants! \n");
            fprintf(gfp,"Max moves: %d | Food distance: %d \n", MAX_GENES, max_to_reach);
        }
    }
    return;
}

// Generate the start population
void generate_start_population(struct Ant pop_arr[]){
        if(gfp != NULL)
        {
            fprintf(gfp,"Generate the start population! \n");
        }
        for(int i = 0; i < POPULATION; i++)
        {
            if(gfp != NULL)
            {
                fprintf(gfp,"#%d ANT GENES: ",i);
            }
            // Generate random genes
            rand_gene(pop_arr[i].gene);
            for(int xg=0; xg < MAX_GENES; xg++)
            {
                if(gfp != NULL)
                {
                    fputc(pop_arr[i].gene[xg],gfp);
                }
            }

            // Get the starting fitness
            pop_arr[i].fitness = get_fitness(pop_arr[i].gene);
            if(gfp != NULL)
            {
                fprintf(gfp," FITNESS: %d \n",pop_arr[i].fitness);
            }
        }
}

int check_mut_prob(int tmp_mut_in){
    int tmp_mut;
    if(tmp_mut_in < 0){
        tmp_mut = 0;
    }
    else if(tmp_mut_in > 100){
        tmp_mut = 100;
    }
    else{
        tmp_mut = tmp_mut_in;
    }
    return tmp_mut;
}

// Arguments handler
void arg_handlr(int arg_count, char *arg_arr[]){
    if(arg_count == 2){
        if(strcmp(arg_arr[1],"-h") == 0){
            printf("If the program is run with the '-i' flag. The parameters can be set 'interactively'\n");
            printf("Program parameters that can be set with arguments: \n");
            printf(" <INT:MAX_GENERATIONS> <INT:MUTATION_PROBABILITY> <INT:FOOD_X> <INT:FOOD_Y> \n");
            printf("Note: mutation probability is a number between 0 and 100\n");
            printf("Note: both food coordinates has to be set!\n");
            if(gfp != NULL)
            {
                fprintf(gfp,"Program were run with the '-h' flag.\n");
                fprintf(gfp,"Program execution end!\n");
            }
            fclose(gfp);
            exit(0);
        }
        else if (strcmp(arg_arr[1],"-i") == 0)
        {
            printf("Set the maximum generation: \n>");
            scanf("%d",&gMax_generation);
            printf("Maximum generation is set to: %d\n",gMax_generation);
            printf("Set the mutation probability (0-100)%%\n>");
            int tmp_mut_prob;
            scanf("%d",&tmp_mut_prob);
            gMutationprob = check_mut_prob(tmp_mut_prob);
            printf("Mutation probability is set to: %d%%\n",gMutationprob);
            char set_food[1];
            printf("Do you want to set the food coordinates?(y/n):\n>");
            scanf("%s",&set_food[0]);
            if(strcmp(set_food,"y") == 0 || strcmp(set_food,"Y") == 0){
                printf("Set the x coordinates for the food\n>");
                scanf("%d",&gFood_x);
                printf("Set the y coordinates for the food\n>");
                scanf("%d",&gFood_y);
            }else{
                printf("Default food coordinates are used! \n");
            }
            printf("Food coordinates (%d,%d)\n",gFood_x,gFood_y);

            return;
        }
        else{
            gMax_generation = atoi(arg_arr[1]);
        }
        
    }
    if(arg_count == 3){
        gMax_generation = atoi(arg_arr[1]);
        int tmp_mut_prob = atoi(arg_arr[2]);
        gMutationprob = check_mut_prob(tmp_mut_prob);
    }
    if(arg_count == 5){
        gMax_generation = atoi(arg_arr[1]);
        int tmp_mut_prob = atoi(arg_arr[2]);
        if(tmp_mut_prob < 0){
            gMutationprob = 0;
        }
        else if(tmp_mut_prob > 100){
            gMutationprob = 100;
        }
        else{
            gMutationprob = tmp_mut_prob;
        }
        gFood_x = atoi(arg_arr[3]);
        gFood_y = atoi(arg_arr[4]);
    }

    return;
}

// For array sort
void swap_arr_items(struct Ant *a, struct Ant *b){
    struct Ant tmp = *a;
    *a = *b;
    *b = tmp;
}

// Sort the array
void arr_sort(struct Ant ant_arr[])
{
int tmp_indx = 0;
    for(int i = 0; i < POPULATION - 1; i++){
        tmp_indx = i;
        for(int v = i + 1; v < POPULATION; v++)
            {
                if(ant_arr[v].fitness < ant_arr[tmp_indx].fitness)
                {
                    tmp_indx = v;
                }
            }
            if(tmp_indx != i){
                swap_arr_items(&ant_arr[tmp_indx], &ant_arr[i]);
            }
    }
}

// Mutate the gene based on the possibility
short mutate_gene(char gene_arr[]){
    int gene_mutation = random_num_gen(1,100);
    int gene_to_mutate;
    int new_gene;
    short is_mutated = 0;
    // If the gene_mutation is smaller then the probability for the mutation
    // Thene mutate a random gene
    if(gfp != NULL)
    {
        for(int za=0; za < MAX_GENES; za++)
        {
            fputc(gene_arr[za],gfp);
        }
    }
    if(gene_mutation < gMutationprob+1)
    {
        // Get the random gene to mutate
        gene_to_mutate = random_num_gen(0,MAX_GENES-1);
        // The new genes can be selected from the mutation char arrays 
        new_gene = random_num_gen(0,2);

        switch (gene_arr[gene_to_mutate])
        {
        case 'l':
            gene_arr[gene_to_mutate] = MUT_GENES_L[new_gene];
            break;
        case 'r':
            gene_arr[gene_to_mutate] = MUT_GENES_R[new_gene];
            break;
        case 'u':
            gene_arr[gene_to_mutate] = MUT_GENES_U[new_gene];
            break;
        case 'd':
            gene_arr[gene_to_mutate] = MUT_GENES_D[new_gene];
            break;
        default:
            printf("The gene can't be mutated..\n");
            break;
        }
        
        if(gfp != NULL)
        {
            fprintf(gfp," GENE: %d is mutated! MUTATED GENES: ",gene_to_mutate + 1);
            for(int za=0; za < MAX_GENES; za++)
            {
                fputc(gene_arr[za],gfp);
            }
        }
        is_mutated = 1;
    }
    return is_mutated;
}
// Generate the next population
void next_gen_pop(struct Ant arr_1[], struct Ant arr_2[])
{
    int whl_a, whl_b, whl_i;

    // Get the point where the genes are sliced.
    int rand_slice_point = random_num_gen(1, MAX_GENES - 1);

   
    if(gfp != NULL)
    {
        fprintf(gfp,"===========================\n");
        fprintf(gfp,"Generation: %d\n",gGenerations + 1);
    }
    whl_a = 0;
    whl_b = 0;

    while((whl_a < POPULATION) && (whl_b < POPULATION/2)){

        // a1-b2
        struct Ant child_1;
        // b1-a2
        struct Ant child_2;

        // a1-a2
        arr_1[whl_a] = arr_2[whl_b];
        
        if(gfp != NULL)
        {
            fprintf(gfp,"#%d ANT GENES: ",whl_a);
        }
        if(mutate_gene(arr_1[whl_a].gene))
        {
            arr_1[whl_a].fitness = get_fitness(arr_1[whl_a].gene);
        }
        else{
            if(gfp != NULL)
            {
                fprintf(gfp," GENE DIDN'T MUTATE! ");
            }
        }
        if(gfp != NULL)
        {
            fprintf(gfp," FITNESS: %d\n",arr_1[whl_a].fitness);
        }
        
        // Increment the array index (select the next ant)
        whl_a++;

        // b1-b2
        arr_1[whl_a] = arr_2[whl_b+1];
        
        if(gfp != NULL)
        {
            fprintf(gfp,"#%d ANT GENES: ",whl_a);
        }

        if(mutate_gene(arr_1[whl_a].gene))
        {
            arr_1[whl_a].fitness = get_fitness(arr_1[whl_a].gene);
        }
        else{
            if(gfp != NULL)
            {
                fprintf(gfp," GENE DIDN'T MUTATE! ");
            }
        }
        if(gfp != NULL)
        {
            fprintf(gfp," FITNESS: %d\n",arr_1[whl_a].fitness);
        }

        // Increment the array index (select the next ant)
        whl_a++;

        // Make the first half of the genes
        for(int i = 0; i < rand_slice_point; i++)
        {   
            child_1.gene[i] = arr_2[whl_b].gene[i];
            child_2.gene[i] = arr_2[whl_b+1].gene[i];
        } 

        // Make the second half of the genes
        for(int y = rand_slice_point; y < MAX_GENES; y++)
        {   
            child_1.gene[y] = arr_2[whl_b+1].gene[y];
            child_2.gene[y] = arr_2[whl_b].gene[y];
        } 

        // Get the fitness for the new childs
        child_1.fitness = get_fitness(child_1.gene);
        child_2.fitness = get_fitness(child_2.gene);

        // Add the new childs to the array
        arr_1[whl_a] = child_1;
        
        if(gfp != NULL)
        {
            fprintf(gfp,"#%d ANT GENES: ",whl_a);
        }

        if(mutate_gene(arr_1[whl_a].gene))
        {
            arr_1[whl_a].fitness = get_fitness(arr_1[whl_a].gene);
        }
        else{
            if(gfp != NULL)
            {
                fprintf(gfp," GENE DIDN'T MUTATE! ");
            }
        }
        if(gfp != NULL)
        {
            fprintf(gfp," FITNESS: %d\n",arr_1[whl_a].fitness);
        }

        // Increment the array index (select the next ant)
        whl_a++;

        arr_1[whl_a] = child_2;

        if(gfp != NULL)
        {
            fprintf(gfp,"#%d ANT GENES: ",whl_a);
        }

        if(mutate_gene(arr_1[whl_a].gene))
        {
            arr_1[whl_a].fitness = get_fitness(arr_1[whl_a].gene);
        }
        else{
            if(gfp != NULL)
            {
                fprintf(gfp," GENE DIDN'T MUTATE! ");
            }
        }
        if(gfp != NULL)
        {
            fprintf(gfp," FITNESS: %d\n",arr_1[whl_a].fitness);
        }


        // Increment the array index (select the next ant)
        whl_a++;
        
        whl_b +=2;
    }

    // Sort the newly created array
    arr_sort(arr_1);

    // The first ant in the array has to have the best fitness because of the sort
    // if it's 0 then the perfect ant is found
    if(arr_1[0].fitness == 0){
        gPerfect_ant = 1;
    }

    // Print out the current generation, and the best fitness. 
    printf("Generation: %d | Best fitness (the lower, the better): %d |\n", gGenerations+1, arr_1[0].fitness);
    
    // Increment the generation count
    gGenerations++;

    // Check for perfect ant and if max generation is reached
    if((gPerfect_ant != 1) && (gGenerations < gMax_generation))
    {

        // Pass the first 50 child to the helper array to create the new ants
        whl_i = 0;
        while(whl_i < 50){
            arr_2[whl_i] = arr_1[whl_i];
            whl_i++;
        }

        // Call this function again.
        next_gen_pop(arr_1, arr_2);
    }
    else{
        // Print some basic information out, why the program is stopped.
        if(gPerfect_ant == 1)
        {
            printf("Perfect ant's gene: ");
            if(gfp != NULL)
            {
                fprintf(gfp,"===========================\n");
                fprintf(gfp,"Perfect ant is found in generation: %d\n",gGenerations);
                fprintf(gfp,"Prfect ant's gene: ");
                
            }
            for(int xa=0; xa < MAX_GENES;xa++)
            {
                printf("%c", arr_1[0].gene[xa]);
                if(gfp != NULL)
                {
                    fputc(arr_1[0].gene[xa],gfp);
                }
            }
            
            printf("\n");
            if(gfp != NULL)
            {
                fprintf(gfp,"\n");
            }
        }
        else if(!(gGenerations < gMax_generation))
        {
            printf("Maximum generation limit: %d is reached. Program execution ends!\n",gMax_generation);
            if(gfp != NULL)
            {
                fprintf(gfp,"Maximum generation limit: %d is reached\n",gMax_generation);
            }
        }

        return;
    }
}

int main(int argc, char *argv[]){

    // Logger filename
    gfp = fopen(filename,"w");
    if(gfp != NULL)
    {
        fprintf(gfp,"Program execution start!\n");
    }
    // For each program run different seed for rand()
    srand(time(0));
    
    int whl_i;

    // Init the arrays that are used
    struct Ant ant_pop_a[POPULATION];
    struct Ant ant_pop_b[POPULATION/2];
    
    // Argument check
    arg_handlr(argc,argv);

    // Print some basic information into the logfile.
    if(gfp != NULL)
    {
        fprintf(gfp,"---------------------------\n");
        fprintf(gfp,"Program parameters!\n");
        fprintf(gfp,"MAX GENERATIONS: %d\n",gMax_generation);
        fprintf(gfp,"MUTATION PROBABILITY: %d\n",gMutationprob);
        fprintf(gfp,"FOOD X COORDINATE: %d\n",gFood_x);
        fprintf(gfp,"FOOD Y COORDINATE: %d\n",gFood_y);
        fprintf(gfp,"---------------------------\n");
    }

    // Map check
    map_check();
    
    // Generate the start population
    generate_start_population(ant_pop_a);

    // The lover the fitness the better!
    arr_sort(ant_pop_a);

    // Pass the first 50 ant to the helper array to create the new ants later.
    whl_i = 0;
    while(whl_i < 50){
        ant_pop_b[whl_i] = ant_pop_a[whl_i];
        whl_i++;
    }

    // Generate the next generation.
    // Return if the perfect ant is found (it's gene stops on the food)
    // Or return if the max generation (MAX_GENERATION) 
    next_gen_pop(ant_pop_a, ant_pop_b);

    if(gfp != NULL)
    {

        printf("For more details check the generated file: %s\n",filename);
        fprintf(gfp,"Program execution end!\n");
    }

    fclose(gfp);
    return 0;
}
