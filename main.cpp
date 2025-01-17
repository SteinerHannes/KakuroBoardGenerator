#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <unordered_set>
#include <utility>

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
};

// Generate a random solution ensuring no duplicates in rows or columns
Solution generateRandomSolution(int gridSize) {
    Solution sol(gridSize);
    vector<unordered_set<int>> rowUsed(gridSize);
    vector<unordered_set<int>> colUsed(gridSize);

    for (int i = 0; i < gridSize; ++i) {
        for (int j = 0; j < gridSize; ++j) {
            if (rand() % 5 == 0) {
                sol.grid[i][j] = -1; // Unfillable cell
            } else {
                int value;
                do {
                    value = rand() % (MAX_DIGIT - MIN_DIGIT + 1) + MIN_DIGIT;
                } while (rowUsed[i].count(value) || colUsed[j].count(value));

                sol.grid[i][j] = value;
                rowUsed[i].insert(value);
                colUsed[j].insert(value);
            }
        }
    }

    return sol;
}

// Evaluate fitness of a solution
int evaluateFitness(const Solution &sol) {
    int gridSize = sol.grid.size();
    int fitness = 0;

    // Check rows and columns for uniqueness
    for (int i = 0; i < gridSize; ++i) {
        unordered_set<int> rowSet, colSet;
        bool validRow = true, validCol = true;

        for (int j = 0; j < gridSize; ++j) {
            // Row uniqueness
            if (sol.grid[i][j] != -1) {
                if (rowSet.count(sol.grid[i][j])) {
                    validRow = false;
                } else {
                    rowSet.insert(sol.grid[i][j]);
                }
            }

            // Column uniqueness
            if (sol.grid[j][i] != -1) {
                if (colSet.count(sol.grid[j][i])) {
                    validCol = false;
                } else {
                    colSet.insert(sol.grid[j][i]);
                }
            }
        }

        if (validRow) fitness += 1; // Reward valid rows
        if (validCol) fitness += 1; // Reward valid columns
    }

    return fitness;
}

// Perform crossover between two parents
Solution crossover(const Solution &parent1, const Solution &parent2) {
    int gridSize = parent1.grid.size();
    Solution child(gridSize);
    vector<unordered_set<int>> rowUsed(gridSize);
    vector<unordered_set<int>> colUsed(gridSize);

    for (int i = 0; i < gridSize; ++i) {
        for (int j = 0; j < gridSize; ++j) {
            int value = (rand() % 2 == 0) ? parent1.grid[i][j] : parent2.grid[i][j];

            // Ensure no duplicates in rows and columns
            if (value != -1 && (rowUsed[i].count(value) || colUsed[j].count(value))) {
                value = -1; // Mark as unfillable if duplicate
            }

            child.grid[i][j] = value;
            if (value != -1) {
                rowUsed[i].insert(value);
                colUsed[j].insert(value);
            }
        }
    }

    return child;
}

// Mutate a solution ensuring no duplicates in the row or column
void mutate(Solution &sol) {
    int gridSize = sol.grid.size();
    int r = rand() % gridSize;
    int c = rand() % gridSize;

    if (sol.grid[r][c] != -1) {
        unordered_set<int> used;

        // Collect all used values in the row and column
        for (int i = 0; i < gridSize; ++i) {
            if (sol.grid[r][i] != -1) used.insert(sol.grid[r][i]);
            if (sol.grid[i][c] != -1) used.insert(sol.grid[i][c]);
        }

        int newValue;
        do {
            newValue = rand() % (MAX_DIGIT - MIN_DIGIT + 1) + MIN_DIGIT;
        } while (used.count(newValue));

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

    for (int gen = 0; gen < generations; ++gen) {
        // Evaluate fitness
        for (auto &sol : population) {
            sol.fitness = evaluateFitness(sol);
            if (sol.fitness > bestSolution.fitness) {
                bestSolution = sol;
            }
        }

        // Selection
        vector<Solution> newPopulation;
        for (int i = 0; i < populationSize; ++i) {
            Solution &parent1 = population[rand() % populationSize];
            Solution &parent2 = population[rand() % populationSize];
            Solution child = crossover(parent1, parent2);

            if (rand() % 100 < 10) { // Mutation chance 10%
                mutate(child);
            }

            newPopulation.push_back(child);
        }

        population = newPopulation;
    }

    return bestSolution;
}

int main() {
    srand(time(0));

    // Define board size
    int gridSize = 5;
;
    // Run evolutionary algorithm
    Solution best = evolutionaryAlgorithm(gridSize, 50, 1000);

    // Print the best solution
    for (const auto &row : best.grid) {
        for (int cell : row) {
            cout << (cell == -1 ? "X" : to_string(cell)) << " ";
        }
        cout << endl;
    }

    return 0;
}