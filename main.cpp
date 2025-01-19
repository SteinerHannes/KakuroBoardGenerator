#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <unordered_set>
#include <utility>
#include <algorithm>

using namespace std;

const int MIN_DIGIT = 1;
const int MAX_DIGIT = 9;

struct Solution {
    vector<vector<int>> grid; // 2D grid representing the solution
    int fitness;

    Solution(int gridSize) {
        grid = vector<vector<int>>(gridSize, vector<int>(gridSize, 0));
        fitness = 0;
    }

     // Overloading < operator
    bool operator<(const Solution& obj1) const
    {
        return !(fitness < obj1.fitness);
    }
};

// Generate a random solution ensuring no duplicates in rows or columns
Solution generateRandomSolution(int gridSize) {
    Solution sol(gridSize);

    for (int i = 0; i < gridSize; ++i) {
        for (int j = 0; j < gridSize; ++j) {
            if(i==0 || j==0){
                sol.grid[i][j] = -1; // Unfillable cell
            }
            else if (rand() % 5 == 0) {
                sol.grid[i][j] = -1; 
            } else {
                sol.grid[i][j] = rand() % (MAX_DIGIT - MIN_DIGIT + 1) + MIN_DIGIT;
            }
        }
    }

    return sol;
}

// Evaluate fitness of a solution
// int evaluateFitness(const Solution &sol) {
//     int gridSize = sol.grid.size();
//     int fitness = 0;

//     // Check rows and columns for uniqueness
//     for (int i = 0; i < gridSize; ++i) {
//         unordered_set<int> rowSet, colSet;
//         bool validRow = true, validCol = true;

//         for (int j = 0; j < gridSize; ++j) {
//             // Row uniqueness
//             if (sol.grid[i][j] != -1) {
//                 if (rowSet.count(sol.grid[i][j])) {
//                     validRow = false;
//                 } else {
//                     rowSet.insert(sol.grid[i][j]);
//                 }
//             }
//             else{
//                 rowSet.clear();
//             }

//             // Column uniqueness
//             if (sol.grid[j][i] != -1) {
//                 if (colSet.count(sol.grid[j][i])) {
//                     validCol = false;
//                 } else {
//                     colSet.insert(sol.grid[j][i]);
//                 }
//             }
//             else{
//                 colSet.clear();
//             }
//         }

//         if (validRow & validCol){
//             fitness += 100; // Reward valid row
//         }
//         else{
//             fitness -= 100; // Punish invalid row
//         }
//     }

//     return fitness;
// }

int evaluateFitness(const Solution &sol) {
    int gridSize = sol.grid.size();
    int fitness = 0;
    bool validBoard = true;

    // Check rows and columns
    for (int i = 0; i < gridSize; ++i) {
        // Check each row and column
        for (int isColumn = 0; isColumn <= 1; ++isColumn) {
            vector<int> sequence;
            int consecutiveBlocks = 0;
            int numbersInGroup = 0;
            
            for (int j = 0; j < gridSize; ++j) {
                int value = isColumn ? sol.grid[j][i] : sol.grid[i][j];
                
                if (value == -1) {
                    // Check for too many consecutive blocks
                    consecutiveBlocks++;
                    if (consecutiveBlocks > 4) {
                        fitness -= 50; // Punish too many consecutive blocks
                    }
                    
                    // Reset number sequence when hitting a block
                    sequence.clear();
                    numbersInGroup = 0;
                } else {
                    consecutiveBlocks = 0;
                    numbersInGroup++;
                    
                
                    
                    // Check for repeating numbers in current sequence
                    bool foundRepeat = false;
                    for (int num : sequence) {
                        if (num == value) {
                            fitness -= 600; // Punish repeated numbers
                            foundRepeat = true;
                            validBoard = false;
                            break;
                        }
                    }
                    
                    // if (!foundRepeat) {
                    //     fitness += 10; // Reward unique numbers
                    // }
                    
                    sequence.push_back(value);
                    
                    // Check if group size exceeds 9
                    if (numbersInGroup > 6) {
                        fitness -= 60; // Punish groups larger than 9
                    }
                }
            }
            
            // Reward if the line meets all criteria
            if (fitness > 0) {
                fitness += 60;
            }
        }
    }

    if(validBoard){
        fitness += 1000;
    }

    return fitness;
}

// Helper function to check if a value is unique in a row and column
bool isUnique(const Solution &sol, int row, int col, int value) {
     int gridSize = sol.grid.size();

    // Check row for duplicates within contiguous groups
    unordered_set<int> rowSet;
    for (int j = 0; j < gridSize; ++j) {
        if (sol.grid[row][j] == -1) {
            rowSet.clear(); // Reset when encountering a block
        } else if (sol.grid[row][j] == value && j != col) {
            return false; // Duplicate found in the same group
        } else {
            rowSet.insert(sol.grid[row][j]);
        }
    }

    // Check column for duplicates within contiguous groups
    unordered_set<int> colSet;
    for (int i = 0; i < gridSize; ++i) {
        if (sol.grid[i][col] == -1) {
            colSet.clear(); // Reset when encountering a block
        } else if (sol.grid[i][col] == value && i != row) {
            return false; // Duplicate found in the same group
        } else {
            colSet.insert(sol.grid[i][col]);
        }
    }

    return true;
}

Solution crossover(const Solution &parent1, const Solution &parent2) {
    int gridSize = parent1.grid.size();
    Solution child(gridSize);
    
    // First row and column are always -1
    for (int i = 0; i < gridSize; i++) {
        child.grid[0][i] = -1;
        child.grid[i][0] = -1;
    }

    // Keep track of numbers in each column between blocks
    vector<unordered_set<int>> colSets(gridSize);
    vector<bool> colHasBlock(gridSize, false);  // Track if column has -1

    for (int i = 1; i < gridSize; ++i) {
        unordered_set<int> rowSet;
        bool hasBlockInRow = false;

        for (int j = 1; j < gridSize; ++j) {
            // Randomly choose between placing a block or a number
            bool placeBlock = (rand() % 10 == 0);  // 10% chance for a block
            
            if (placeBlock) {
                child.grid[i][j] = -1;
                rowSet.clear();  // Reset row set after block
                colSets[j].clear();  // Reset column set after block
                hasBlockInRow = true;
                colHasBlock[j] = true;
            } else {
                // Try to get valid value from either parent
                const Solution &selectedParent = (rand() % 2 == 0) ? parent1 : parent2;
                vector<int> validValues;
                
                // Collect all valid values
                for (int val = 1; val <= 9; ++val) {
                    if (rowSet.count(val) == 0 && colSets[j].count(val) == 0) {
                        validValues.push_back(val);
                    }
                }
                
                if (validValues.empty()) {
                    // If no valid values, place a block
                    child.grid[i][j] = -1;
                    rowSet.clear();
                    colSets[j].clear();
                    hasBlockInRow = true;
                    colHasBlock[j] = true;
                } else {
                    // Place a random valid value
                    int valueIndex = rand() % validValues.size();
                    int value = validValues[valueIndex];
                    child.grid[i][j] = value;
                    rowSet.insert(value);
                    colSets[j].insert(value);
                }
            }
        }
    }

    // Ensure no more than 4 consecutive blocks
    for (int i = 1; i < gridSize - 4; ++i) {
        for (int j = 1; j < gridSize; ++j) {
            if (child.grid[i][j] == -1 && 
                child.grid[i+1][j] == -1 && 
                child.grid[i+2][j] == -1 && 
                child.grid[i+3][j] == -1 && 
                child.grid[i+4][j] == -1) {
                child.grid[i+4][j] = rand() % 9 + 1;  // Replace last block with number
            }
        }
    }

    return child;
}

// Mutate a solution ensuring no duplicates in the row or column
void mutate(Solution &sol) {
    int gridSize = sol.grid.size();
    int r = 0;
    int c = 0;

    while (r == 0 || c == 0 || sol.grid[r][c] == -1)
    {
        r = rand() % gridSize;
        c = rand() % gridSize;
    }
    
    int newValue;
            
    // if(rand() % 5 == 0) {
    //     newValue = -1; 
    // } else {
        newValue = rand() % (MAX_DIGIT - MIN_DIGIT + 1) + MIN_DIGIT;
    //}

    sol.grid[r][c] = newValue;    
}

// Main evolutionary algorithm
Solution evolutionaryAlgorithm(int gridSize, int populationSize, int generations) {
    vector<Solution> population;

    // Initialize random population
    for (int i = 0; i < populationSize; ++i) {
        population.push_back(generateRandomSolution(gridSize));
    }

    Solution bestSolution = population[0];
    bestSolution.fitness = evaluateFitness(bestSolution);

    for (int gen = 0; gen < generations; ++gen) {
        // Evaluate fitness
        for (auto &sol : population) {
            sol.fitness = evaluateFitness(sol);
            if (sol.fitness > bestSolution.fitness) {
                bestSolution = sol;
                
                //     for (const auto &row : bestSolution.grid) {
                //         for (int cell : row) {
                //             cout << (cell == -1 ? "X" : to_string(cell)) << " ";
                //         }
                //         cout << endl;
                //     }
                // cout << endl;
            }
        }

        // Selection
        vector<Solution> newPopulation;

        //sort population by fitness
        std::sort(population.begin(), population.end());

        for (int i = 0; i < populationSize; ++i) {
            Solution &parent1 = population[i];
            Solution &parent2 = population[i];
            
            if(i == populationSize - 1){
                parent2 = population[0];
            }
            else{
            }

            Solution child = crossover(parent1, parent2);

            if (rand() % 100 < 30) { // Mutation chance 10%
                mutate(child);
            }

            newPopulation.push_back(child);
        }

        population = newPopulation;

        cout << "Generation " << gen << " best fitness: " << bestSolution.fitness << endl;
    }

    return bestSolution;
}

int main() {
    srand(time(0));

    // Define board size
    int gridSize = 12;
;
    // Run evolutionary algorithm
    Solution best = evolutionaryAlgorithm(gridSize, 500, 19);

    // Print the best solution
    for (const auto &row : best.grid) {
        for (int cell : row) {
            cout << (cell == -1 ? "X" : to_string(cell)) << " ";
        }
        cout << endl;
    }

    return 0;
}