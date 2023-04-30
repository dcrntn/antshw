# Ant colony
Generates an ant colony with the given population and uses a genetic algorithm to try to reach the "perfect ant"!
This is a homework project and useless outside its purpose! 

## Installation
1. Clone or download this repo.
2. Compile the "main.c" file with any C compiler e.g: gcc.
3. Run the application.

## Usage
In this example the compiled filename is "ants".
```sh
# Run the program itself
./ants

# Run with the help flag
# This way the program will only display what flag(s) and parameters can be used or set.
./ants -h

# Run it "interactively" - This way you can set the parameters in the console
./ants -i

# Run it with parameters <INT:MAX_GENERATIONS> <INT:MUTATION_PROBABILITY> <INT:FOOD_X> <INT:FOOD_Y>
# Both food parameter needs to be set simultaneously
./ants 200 85 -14 6
```

## Output
The program will generate a "ants_log.txt" file that logs the following:
- Program parameterss
- Start population with each ant's genes
- Next generation populations with generation number, ants' genes, mutated genes
- The reason for the program to end.
