#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <unordered_set>
#include <utility>
#include <algorithm>
#include <iomanip>
#include <fstream>
#include <set>

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

Solution generateRandomSolution(int gridSize) {
    Solution sol(gridSize);
    bool hasNumber = false;
    
    // First pass: generate grid
    for (int i = 0; i < gridSize; ++i) {
        for (int j = 0; j < gridSize; ++j) {
            if (i == 0 || j == 0) {
                sol.grid[i][j] = -1; // Unfillable cell
            }
            else if (rand() % 5 == 0) {
                sol.grid[i][j] = -1;
            } else {
                sol.grid[i][j] = rand() % (MAX_DIGIT - MIN_DIGIT + 1) + MIN_DIGIT;
                hasNumber = true;
            }
        }
    }
    
    // If no numbers were generated, force at least one
    if (!hasNumber) {
        int i = 1 + (rand() % (gridSize - 1));
        int j = 1 + (rand() % (gridSize - 1));
        sol.grid[i][j] = rand() % (MAX_DIGIT - MIN_DIGIT + 1) + MIN_DIGIT;
    }
    
    // Second pass: fix isolated numbers
    bool needsFixing;
    do {
        needsFixing = false;
        for (int i = 1; i < gridSize; ++i) {
            for (int j = 1; j < gridSize; ++j) {
                if (sol.grid[i][j] != -1) {
                    // Check if number is isolated
                    bool isIsolated = true;
                    // Check adjacent cells (not diagonal)
                    if ((i > 0 && sol.grid[i-1][j] != -1) ||
                        (i < gridSize-1 && sol.grid[i+1][j] != -1) ||
                        (j > 0 && sol.grid[i][j-1] != -1) ||
                        (j < gridSize-1 && sol.grid[i][j+1] != -1)) {
                        isIsolated = false;
                    }
                    
                    if (isIsolated) {
                        // Fix isolation by randomly making one adjacent cell a number
                        vector<pair<int,int>> adjacent;
                        if (i > 0) adjacent.push_back({i-1, j});
                        if (i < gridSize-1) adjacent.push_back({i+1, j});
                        if (j > 0) adjacent.push_back({i, j-1});
                        if (j < gridSize-1) adjacent.push_back({i, j+1});
                        
                        if (!adjacent.empty()) {
                            int idx = rand() % adjacent.size();
                            sol.grid[adjacent[idx].first][adjacent[idx].second] = 
                                rand() % (MAX_DIGIT - MIN_DIGIT + 1) + MIN_DIGIT;
                            needsFixing = true;
                        }
                    }
                }
            }
        }
    } while (needsFixing);
    return sol;
}


int evaluateFitness(const Solution &sol) {
    int gridSize = sol.grid.size();
    int fitness = 0;
    bool validBoard = true;
    /*
      4 4
    4 1 3 
    4 3 1 

    find pairs, +1 for each, scale rewards 
    
    */

    // Check rows and columns
    for (int i = 1; i < gridSize; ++i) {
        // Check each row and column
        for (int isColumn = 0; isColumn <= 1; ++isColumn) {
            vector<int> sequence;
            int numbersInGroup = 0;
            
            for (int j = 1; j < gridSize; ++j) {
                int value = isColumn ? sol.grid[j][i] : sol.grid[i][j];
                
                if (value == -1) {
                    // Reset number sequence when hitting a block
                    sequence.clear();
                    if(numbersInGroup == 1){
                        fitness -= 1;
                    }
                    numbersInGroup = 0;
                } else {
                    numbersInGroup++;
                  
                    // Check for repeating numbers in current sequence
                    for (int num : sequence) {
                        if (num == value) {
                            fitness -= 2; // Punish repeated numbers
                            validBoard = false;
                            //cout << "repeated numbers " << num << endl;
                            break;
                        }
                        
                    }
                    
                    sequence.push_back(value);

                    /*  
                    0 0 0 0 0 0 0 0 0 0 0
                    0 0 x x x x x x x x x 
                    0 x x x x x x x x x x

                    iterate through board + save indices of numbers
                    get first number in board 
                    depth first search through all indices, collect visited 
                    if visited smaller than all numbers, not reachable
                    start search again with number at smallest not visited index

                    */
                    // Check if group size exceeds 9
                    if (numbersInGroup > 9) {
                        fitness -= 1; // Punish groups larger than 9 or smaller than 2
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
    int gridSize = child.grid.size();

     // If trying to place -1, check for consecutive numbers rule
    if (value == -1) {
        // Only need to check if previous cell exists and contains a number
        if (col > 0 && child.grid[row][col-1] != -1) {
            // Count consecutive numbers before this position
            int consecutiveCount = 0;
            int checkCol = col - 1;
            
            // Keep going left until we hit a -1 or the start of the row
            while (checkCol >= 0 && child.grid[row][checkCol] != -1) {
                consecutiveCount++;
                checkCol--;
            }
            
            // If only one number is found, -1 cannot be placed
            if (consecutiveCount == 1) {
                //cout << "one consecutive" << endl;
                return false;
            }
        }
        return true;
    }

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

    // vertical or horizontal
    const int orientation = (rand() % 2 == 0) ? 0 : 1;
    const int percentageP1 = rand() % 100;
    const int partP1 = gridSize * (percentageP1/100.0);

    //cout << partP1 << endl;

    if(orientation == 0){ //horizontal
        for (int i = 0; i < gridSize; ++i) {
            for (int j = 0; j < gridSize; ++j) {
                if(i <= partP1){
                    child.grid[i][j] = parent1.grid[i][j];
                }else{
                    child.grid[i][j] = parent2.grid[i][j];
                }
            }

        }
    }else{//vertical
         for (int i = 0; i < gridSize; ++i) {
            for (int j = 0; j < gridSize; ++j) {
                if(j <= partP1){
                    child.grid[i][j] = parent1.grid[i][j];
                }else{
                    child.grid[i][j] = parent2.grid[i][j];
                }
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

    const int mutations = (rand() % 100 / 100) * (gridSize-1)*(gridSize-1); //number of cells to mutate (perc * numOfCells)
    std::set<std::pair<int, int>> coordinates;
    
    for(int i = 0; i < mutations; i++){
        while (coordinates.find({r, c}) != coordinates.end()){ // skip already used coordinates
            r = rand() % (gridSize - 1) + 1;
            c = rand() % (gridSize - 1) + 1;
        }
        coordinates.insert({r,c});
        int newValue = rand() % (MAX_DIGIT - MIN_DIGIT + 1) + MIN_DIGIT;
        sol.grid[r][c] = newValue;    
    }    
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

    for (auto &sol : population) {
        sol.fitness = evaluateFitness(sol);
        cout << "parent fitness " << sol.fitness << endl;
            
        if (sol.fitness > bestSolution.fitness) {
            bestSolution = sol; // save to file?
        }
    }

    int gen = 0;

    while (gen < generations) { // do not check for fitness, stop if valid board generated

        for (auto &sol : population) {
            sol.fitness = evaluateFitness(sol);   
            if (sol.fitness > bestSolution.fitness) {
                bestSolution = sol; // save to file?
                cout << "New best solution" << endl;
            }
        }
        if(population[0].fitness > bestSolution.fitness){
            bestSolution = population[0];
            
        }

        cout << "Generation " << gen << " best fitness: " << bestSolution.fitness << endl;
        // extra validity check for best board, stop if valid

        // Selection
        vector<Solution> children;
        // iterate through population: parent1 = i, parent2 = random
        // better fitness means more likely to be parent

        for (int i = 0; i < populationSize; ++i) {
            
            int idx2 = 0;

            if(rand() % 2 == 0){
                idx2 = rand() % populationSize;
            }
            else{
                idx2 = rand() % (populationSize/2);
            }
        
            Solution &parent1 = population[i];
            Solution &parent2 = population[idx2];
            
            Solution child = crossover(parent1, parent2);

            if (rand() % 100 < 60) { // Mutation chance 60%
                mutate(child); 
            }

            child.fitness = evaluateFitness(child);

            // check if perfect solution has been created, then return

            children.push_back(child);
        }

        vector<Solution>* totalPopulation = new vector<Solution>;

        totalPopulation->reserve(population.size() + children.size());  // Optional but improves performance
        totalPopulation->insert(totalPopulation->end(), population.begin(), population.end());
        totalPopulation->insert(totalPopulation->end(), children.begin(), children.end()); 

        sort(totalPopulation->begin(), totalPopulation->end(), [](const Solution& a, const Solution& b) {
                return a.fitness > b.fitness;
            }); 
        
        population.clear();
        population = vector<Solution>(totalPopulation->begin(), totalPopulation->begin() + 500);

        delete totalPopulation;
        gen++;

        // fitness check
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
    int gridSize = 10;
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