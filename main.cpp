#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <unordered_set>
#include <utility>
#include <algorithm>
#include <iomanip>
#include <fstream>

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


int evaluateFitness(const Solution &sol) {
    int gridSize = sol.grid.size();
    int fitness = 0;
    bool validBoard = true;

    // Check rows and columns
    for (int i = 1; i < gridSize; ++i) {
        // Check each row and column
        for (int isColumn = 0; isColumn <= 1; ++isColumn) {
            vector<int> sequence;
            int consecutiveBlocks = 0;
            int numbersInGroup = 0;
            
            for (int j = 1; j < gridSize; ++j) {
                int value = isColumn ? sol.grid[j][i] : sol.grid[i][j];
                
                if (value == -1) {
                    // Check for too many consecutive blocks
                    consecutiveBlocks++;
                    if (consecutiveBlocks > 5) {
                        fitness -= 50; // Punish too many consecutive blocks
                       
                    }
                    
                    // Reset number sequence when hitting a block
                    sequence.clear();
                    numbersInGroup = 0;
                } else {
                    consecutiveBlocks = 0;
                    numbersInGroup++;
                    //cout << "Cons blocks 0" << endl;
                
                    
                    // Check for repeating numbers in current sequence
                    bool foundRepeat = false;
                    for (int num : sequence) {
                        if (num == value) {
                            fitness -= 60; // Punish repeated numbers
                            foundRepeat = true;
                            validBoard = false;
                            //cout << "repeated numbers " << num << endl;
                            break;
                        }
                        
                    }
                    
                    // if (!foundRepeat) {
                    //     fitness += 10; // Reward unique numbers
                    // }
                    
                    sequence.push_back(value);
                    
                    // Check if group size exceeds 9
                    if (numbersInGroup > 9) {
                        fitness -= 60; // Punish groups larger than 9
                        validBoard = false;
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

bool isValidInChild(const Solution &child, int row, int col, int value) {
    if (value == -1) return true; // Always valid for block cells

    int gridSize = child.grid.size();

    // Check row for duplicates within contiguous groups
    unordered_set<int> rowSet;
    for (int j = 0; j < gridSize; ++j) {
        if (child.grid[row][j] == -1) {
            rowSet.clear(); // Reset when encountering a block
        } else if (child.grid[row][j] == value) {
            return false; // Duplicate found in the same group
        } else {
            rowSet.insert(child.grid[row][j]);
        }
    }

    // Check column for duplicates within contiguous groups
    unordered_set<int> colSet;
    for (int i = 0; i < gridSize; ++i) {
        if (child.grid[i][col] == -1) {
            colSet.clear(); // Reset when encountering a block
        } else if (child.grid[i][col] == value) {
            return false; // Duplicate found in the same group
        } else {
            colSet.insert(child.grid[i][col]);
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

    // Randomly choose sections from parents
    for (int i = 1; i < gridSize; ++i) {
        for (int j = 1; j < gridSize; ++j) {
            if (rand() % 100 < 70) {
                // 80% chance to copy directly from one of the parents
                const Solution &selectedParent = (rand() % 2 == 0) ? parent1 : parent2;
                child.grid[i][j] = selectedParent.grid[i][j];
                continue; 
            }

            // 20% chance or if the copied value was invalid, search for a valid value
            bool valueFound = false;
            for (int attempts = 0; attempts < 10; ++attempts) { // Limit retries
                const Solution &selectedParent = (rand() % 2 == 0) ? parent1 : parent2;
                int r = 1 + (rand() % (gridSize - 1));
                int c = 1 + (rand() % (gridSize - 1));
                int value = selectedParent.grid[r][c];

                if (isValidInChild(child, i, j, value)) {
                    child.grid[i][j] = value;
                    valueFound = true;
                    break;
                }
            }

            // Fallback: If no valid value found, copy directly from the same position in one of the parents
            if (!valueFound) {
                const Solution &selectedParent = (rand() % 2 == 0) ? parent1 : parent2;
                child.grid[i][j] = selectedParent.grid[i][j];
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
    int gen = 0;

    while (bestSolution.fitness < 1000) {
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
        

        for (int i = 0; i < populationSize; ++i) {
            int idx1 = rand() % populationSize;
            int idx2 = rand() % populationSize;
                    
            Solution &parent1 = population[idx1];
            Solution &parent2 = population[idx2];
            

            Solution child = crossover(parent1, parent2);

            if (rand() % 100 < 40) { // Mutation chance 10%
                mutate(child);
            }

            newPopulation.push_back(child);
        }

        population = newPopulation;

        cout << "Generation " << gen << " best fitness: " << bestSolution.fitness << endl;
        gen++;
    }

    return bestSolution;
}

vector<vector<string>> deriveEmptyBoard(const vector<vector<int>> &grid) {
    int gridSize = grid.size();
    vector<vector<string>> emptyBoard(gridSize, vector<string>(gridSize, "E"));

    for (int i = 0; i < gridSize; ++i) {
        for (int j = 0; j < gridSize; ++j) {
            if (grid[i][j] == -1) { // Potential clue cell
                int horizontalSum = 0, verticalSum = 0;

                // Compute horizontal sum
                for (int k = j + 1; k < gridSize && grid[i][k] != -1; ++k) {
                    horizontalSum += grid[i][k];
                }

                // Compute vertical sum
                for (int k = i + 1; k < gridSize && grid[k][j] != -1; ++k) {
                    verticalSum += grid[k][j];
                }

                // Construct clue cell label, ensuring "/" is always included
                string clue = to_string(verticalSum) + "/" + to_string(horizontalSum);
                if (verticalSum == 0) clue = "/" + to_string(horizontalSum); // Only horizontal
                if (horizontalSum == 0) clue = to_string(verticalSum) + "/"; // Only vertical

                if (clue == "0/0" ||  clue == "/0" ||  clue == "0/") {
                    clue = "#"; // Replace unmodifiable cells
                }

                emptyBoard[i][j] = clue;
            } else if (grid[i][j] != -1) {
                emptyBoard[i][j] = ""; // Replace modifiable cells
            }
        }
    }

    return emptyBoard;
}

// Function to print the derived empty board with aligned columns
void printEmptyBoard(const vector<vector<string>> &emptyBoard) {
    // Determine the maximum width of any cell
    size_t maxWidth = 0;
    for (const auto &row : emptyBoard) {
        for (const auto &cell : row) {
            maxWidth = max(maxWidth, cell.size());
        }
    }

    // Print each row with aligned columns
    for (const auto &row : emptyBoard) {
        for (size_t j = 0; j < row.size(); ++j) {
            cout << std::setw(maxWidth) << row[j]; // Align each cell
            if (j < row.size() - 1) {
                cout << " "; // Add space between columns
            }
        }
        cout << endl;
    }
}

void saveEmptyBoardToFile(const vector<vector<string>> &emptyBoard, const string &filename) {
    // Open the file for writing
    ofstream outFile(filename);

    if (!outFile) {
        cerr << "Error opening file for writing!" << endl;
        return;
    }

    // Determine the maximum width of any cell
    size_t maxWidth = 0;
    for (const auto &row : emptyBoard) {
        for (const auto &cell : row) {
            maxWidth = max(maxWidth, cell.size());
        }
    }

    // Write each row with aligned columns to the file
    for (const auto &row : emptyBoard) {
        for (size_t j = 0; j < row.size(); ++j) {
            outFile << std::setw(maxWidth) << row[j]; // Align each cell
            if (j < row.size() - 1) {
                outFile << " "; // Add space between columns
            }
        }
        outFile << endl;
    }

    outFile.close(); // Close the file
}

int main() {
    srand(time(0));

    // Define board size
    int gridSize = 12;
;
    // Run evolutionary algorithm
    Solution best = evolutionaryAlgorithm(gridSize, 500, 3000);

    // Print the best solution
    for (const auto &row : best.grid) {
        for (int cell : row) {
            cout << (cell == -1 ? "X" : to_string(cell)) << " ";
        }
        cout << endl;
    }

    // Derive empty board
    vector<vector<string>> emptyBoard = deriveEmptyBoard(best.grid);

    // Print empty board
    printEmptyBoard(emptyBoard);

    return 0;
}