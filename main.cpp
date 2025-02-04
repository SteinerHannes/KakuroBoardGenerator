#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <unordered_set>
#include <algorithm>
#include <iomanip>
#include <fstream>
#include <set>
#include <queue>

constexpr bool debugMutations = false;
constexpr bool debugScore = false;
constexpr bool debugOutput = true;

constexpr auto gridSizes = { 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20 };
constexpr int populationSizeFactor = 100;
constexpr int nGenerationsFactor = 1000;

constexpr int punishmentLargeGroups = 1;
constexpr int punishmentRepetition = 3;
constexpr int rewardCorrectBoard = 100;
constexpr int rewardSingleGroup = 20;
constexpr int punishmentPerGroup = 10;
constexpr int rewardNoIsolatedNumbers = 10;
constexpr int punishmentIsolatedNumbers = 10;
constexpr int punishmentThreeByThree = 10;
constexpr int punishmentPairs = 4;
constexpr int punishmentNotSymmetric = 1;
constexpr int punishmentUnusedRowOrColumn = 10;

constexpr int MIN_DIGIT = 1;
constexpr int MAX_DIGIT = 9;

struct Solution {
    std::vector<std::vector<int>> grid; // 2D grid representing the solution
    int fitness;

    explicit Solution(const int gridSize) {
        grid = std::vector(gridSize, std::vector(gridSize, 0));
        fitness = 0;
    }

    bool operator<(const Solution& obj1) const
    {
        return fitness >= obj1.fitness;
    }
};

Solution generateRandomSolution(const int gridSize) {
    Solution sol(gridSize);
    bool hasNumber = false;
    
    // First pass: generate grid
    for (int i = 0; i < gridSize; ++i) {
        for (int j = 0; j < gridSize; ++j) {
            if (i == 0 || j == 0 || rand() % 5 == 0) {
                sol.grid[i][j] = -1;
            } else {
                sol.grid[i][j] = rand() % (MAX_DIGIT - MIN_DIGIT + 1) + MIN_DIGIT;
                hasNumber = true;
            }
        }
    }
    
    // If no numbers were generated, force at least one
    if (!hasNumber) {
        const int i = 1 + (rand() % (gridSize - 1));
        const int j = 1 + (rand() % (gridSize - 1));
        sol.grid[i][j] = rand() % (MAX_DIGIT - MIN_DIGIT + 1) + MIN_DIGIT;
    }

    return sol;
}

int countConnectedGroups(const Solution &sol) {
    const int rows = static_cast<int>(sol.grid.size());
    const int cols = static_cast<int>(sol.grid[0].size());

    // Directions for vertical and horizontal movement
    const std::vector<std::pair<int, int>> directions = {{-1, 0}, {1, 0}, {0, -1}, {0, 1}};

    // Visited set to track visited cells
    std::vector visited(rows, std::vector(cols, false));

    // BFS function to explore the connected component
    auto bfs = [&](int startX, int startY) {
        std::queue<std::pair<int, int>> q;
        q.emplace(startX, startY);
        visited[startX][startY] = true;

        while (!q.empty()) {
            auto [x, y] = q.front();
            q.pop();

            // Explore all 4 possible neighbors (up, down, left, right)
            for (auto& dir : directions) {
                int newX = x + dir.first;
                int newY = y + dir.second;

                if (newX >= 0 && newX < rows && newY >= 0 && newY < cols &&
                    sol.grid[newX][newY] >= 1 && sol.grid[newX][newY] <= 9 && !visited[newX][newY]) {
                    visited[newX][newY] = true;
                    q.emplace(newX, newY);
                    }
            }
        }
    };

    int groupCount = 0;

    // Iterate through the grid and start BFS whenever we find an unvisited number
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            // Start a BFS if the current cell contains a number (1-9) and hasn't been visited
            if (sol.grid[i][j] >= 1 && sol.grid[i][j] <= 9 && !visited[i][j]) {
                bfs(i, j); // Explore the connected component starting from (i, j)
                ++groupCount; // Increment the connected group count
            }
        }
    }

    return groupCount;
}

int countNumbersWithoutNeighbours(const Solution &sol) {
    const int rows = static_cast<int>(sol.grid.size());
    if (rows == 0) return 0;
    const int cols = static_cast<int>(sol.grid[0].size());

    int count = 0;

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            if (sol.grid[i][j] > 0) { // Only consider numbers (1-9)
                bool hasHorizontalNeighbor = false;
                bool hasVerticalNeighbor = false;

                // Check left and right neighbors for horizontal
                if (j > 0 && sol.grid[i][j - 1] > 0) {
                    hasHorizontalNeighbor = true;
                }
                if (j < cols - 1 && sol.grid[i][j + 1] > 0) {
                    hasHorizontalNeighbor = true;
                }

                // Check top and bottom neighbors for vertical
                if (i > 0 && sol.grid[i - 1][j] > 0) {
                    hasVerticalNeighbor = true;
                }
                if (i < rows - 1 && sol.grid[i + 1][j] > 0) {
                    hasVerticalNeighbor = true;
                }

                // If either condition fails, increment the count
                if (!hasHorizontalNeighbor || !hasVerticalNeighbor) {
                    ++count;
                }
            }
        }
    }

    return count;
}

int evaluateFitness(const Solution &sol) {
    const int gridSize = (int)sol.grid.size();
    int fitness = 0;
    bool validBoard = true;


    // Check rows and columns
    for (int i = 1; i < gridSize; ++i) {
        bool rowUsed = false;
        bool colUsed = false;
        // Check each row and column
        for (int isColumn = 0; isColumn <= 1; ++isColumn) {
            std::vector<int> sequence;
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
                            fitness -= punishmentRepetition; // Punish repeated numbers
                            validBoard = false;
                            //cout << "repeated numbers " << num << endl;
                            break;
                        }
                        
                    }
                    
                    sequence.push_back(value);

                    auto punishmentFactor = 0.0;

                    if (numbersInGroup < 4) {
                        punishmentFactor = 0;
                    }  else if (numbersInGroup == 4 ) {
                        punishmentFactor = 2;
                    } else if ( numbersInGroup == 5) {
                        punishmentFactor = 3;
                    } else if (numbersInGroup == 6) {
                        punishmentFactor = 4;
                    } else if (numbersInGroup == 7) {
                        punishmentFactor = 6;
                    } else if (numbersInGroup == 8) {
                        punishmentFactor = 10;
                    } else if (numbersInGroup == 9) {
                        punishmentFactor = 16;
                    } else {
                        punishmentFactor = 26;
                    }
                    if (punishmentFactor > 0) {
                        fitness -= (int)(punishmentLargeGroups * punishmentFactor);
                        validBoard = false;
                    }
                }

                // Check everything below only once per grid
                if (isColumn) {
                    // Check if the board is symmetric across the main diagonal
                    if ((sol.grid[i][j] == -1 && sol.grid[j][i] != -1) || (sol.grid[i][j] != -1 && sol.grid[j][i] == -1)) {
                        fitness -= punishmentNotSymmetric;
                    }

                    // Check if there are 3x3 squares. If so punish for each
                    if (sol.grid[i][j] != -1 && sol.grid[i-1][j] != -1 && sol.grid[i][j-1] != -1 && sol.grid[i-1][j-1] != -1) {
                        fitness -= punishmentThreeByThree;
                        validBoard = false;
                    }

                    // Check for swappable pairs
                    if (sol.grid[i][j] != -1 && sol.grid[i-1][j] != -1 && sol.grid[i][j-1] != -1 && sol.grid[i-1][j-1] != -1) {
                        if(sol.grid[i][j] == sol.grid[i-1][j] || sol.grid[i][j] == sol.grid[i][j-1] || sol.grid[i][j] == sol.grid[i-1][j-1]){
                            fitness -= punishmentPairs;
                            validBoard = false;
                        }
                    }

                    if (sol.grid[i][j] != -1) {
                        rowUsed = true;
                    }
                    if (sol.grid[j][i] != -1) {
                        colUsed = true;
                    }
                }
            }
            // Check only once per grid
            if (isColumn) {
                if (!rowUsed) {
                    fitness -= punishmentUnusedRowOrColumn;
                    validBoard = false;
                }
                if (!colUsed) {
                    fitness -= punishmentUnusedRowOrColumn;
                    validBoard = false;
                }
            }
        }
    }

    // Check for connected groups
    int connectedGroups = countConnectedGroups(sol);
    if (connectedGroups == 1) {
        fitness += rewardSingleGroup;
    } else {
        fitness -= punishmentPerGroup * (connectedGroups - 1);
        validBoard = false;
    }

    // Check for numbers without horizontal or vertical neighbors
    int isolatedNumbers = countNumbersWithoutNeighbours(sol);
    if(isolatedNumbers == 0){
        fitness += rewardNoIsolatedNumbers;
    }else{
        fitness -= punishmentIsolatedNumbers * isolatedNumbers;
        validBoard = false;
    }


    if(validBoard){
        fitness += rewardCorrectBoard;
    }

    return fitness;
}

Solution crossover(const Solution &parent1, const Solution &parent2) {
    const int gridSize = static_cast<int>(parent1.grid.size());
    Solution child(gridSize);

    // vertical or horizontal
    const int orientation = (rand() % 2 == 0) ? 0 : 1;
    const int percentageP1 = rand() % 100;
    const int partP1 = gridSize * (percentageP1 / 100);

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
    } else {//vertical
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
    const int gridSize = static_cast<int>(sol.grid.size());
    int r = rand() % (gridSize - 1) + 1;
    int c = rand() % (gridSize - 1) + 1;

    //number of cells to mutate (min 2 cells and max 20% of the board)
    const int mutations = 2 + rand() % (gridSize * gridSize / 5 - 2);
    if (debugMutations) {
        std::cout << "mutations: " << mutations << std::endl;
    }

    std::set<std::pair<int, int>> coordinates;
    
    for(int i = 0; i < mutations; i++){
        while (coordinates.find({r, c}) != coordinates.end()){ // skip already used coordinates
            r = rand() % (gridSize - 1) + 1;
            c = rand() % (gridSize - 1) + 1;
        }
        coordinates.insert({r,c});
        // Create new value (-1) or random number (1-9) with chance 20%
        const int newValue = (rand() % 10 < 2) ? -1 : rand() % (MAX_DIGIT - MIN_DIGIT + 1) + MIN_DIGIT;
        if (debugMutations) {
            std::cout << "mutating " << r << " " << c << " to " << newValue << std::endl;
        }
        sol.grid[r][c] = newValue;    
    }    
}

// Main evolutionary algorithm
Solution evolutionaryAlgorithm(const int gridSize, const int populationSize, const int generations) {
    std::vector<Solution> population;
    population.reserve(populationSize);

    // Initialize random population
    for (int i = 0; i < populationSize; ++i) {
        population.push_back(generateRandomSolution(gridSize));
    }

    Solution bestSolution = population[0];
    bestSolution.fitness = evaluateFitness(bestSolution);

    for (auto &sol : population) {
        sol.fitness = evaluateFitness(sol);
        if (sol.fitness > bestSolution.fitness) {
            bestSolution = sol;
        }
    }

    // Top 10% of population retained as elites
    const int eliteCount = std::max(1, populationSize / 10);
    int gen = 0;
    int bestSolutionGeneration = 0;

    while (gen < generations && gen < bestSolutionGeneration + nGenerationsFactor * 2) {
        if (debugScore) {
            std::cout << "Generation " << gen << " best fitness: " << bestSolution.fitness
                 << " average fitness: " << population[populationSize / 2].fitness << std::endl;
        }

        // Sort population by fitness (descending order)
        sort(population.begin(), population.end(), [](const Solution &a, const Solution &b) {
            return a.fitness > b.fitness;
        });

        // Save the elites
        const std::vector elites(population.begin(), population.begin() + eliteCount);

        // Roulette wheel selection (probabilistic selection based on fitness)
        std::vector<double> cumulativeProbabilities(populationSize);
        double totalFitness = 0.0;
        for (const auto &sol : population) {
            totalFitness += std::max(1.0, static_cast<double>(sol.fitness + 1)); // Ensure all fitness values are positive
        }

        cumulativeProbabilities[0] = std::max(1.0, static_cast<double>(population[0].fitness + 1)) / totalFitness;
        for (int i = 1; i < populationSize; ++i) {
            cumulativeProbabilities[i] = cumulativeProbabilities[i - 1] +
                std::max(1.0, static_cast<double>(population[i].fitness + 1)) / totalFitness;
        }

        // Create children using selection and crossover
        std::vector<Solution> children;
        while (children.size() < populationSize - eliteCount) {
            const double rand1 = static_cast<double>(rand()) / RAND_MAX;
            const double rand2 = static_cast<double>(rand()) / RAND_MAX;

            // Select two parents based on roulette wheel probabilities
            const auto it1 = std::lower_bound(cumulativeProbabilities.begin(), cumulativeProbabilities.end(), rand1);
            const auto it2 = std::lower_bound(cumulativeProbabilities.begin(), cumulativeProbabilities.end(), rand2);

            const int parent1Idx = (int)std::distance(cumulativeProbabilities.begin(), it1);
            const int parent2Idx = (int)std::distance(cumulativeProbabilities.begin(), it2);

            Solution child = crossover(population[parent1Idx], population[parent2Idx]);

            if (rand() % 100 < 80) { // 80% chance to mutate
                mutate(child);
            }

            child.fitness = evaluateFitness(child);
            children.push_back(child);
        }

        // Combine elites and children to form the next generation
        population = elites;
        population.insert(population.end(), children.begin(), children.end());

        // Update the best solution if needed
        for (const auto &sol : population) {
            if (sol.fitness > bestSolution.fitness) {
                bestSolution = sol;
                bestSolutionGeneration = gen;
            }
        }

        gen++;
    }

    return bestSolution;
}

std::vector<std::vector<std::string>> deriveEmptyBoard(const std::vector<std::vector<int>> &grid) {
    const int gridSize = static_cast<int>(grid.size());
    std::vector emptyBoard(gridSize, std::vector<std::string>(gridSize, "E"));

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
                std::string clue = std::to_string(verticalSum) + "/" + std::to_string(horizontalSum);
                if (verticalSum == 0) clue = "/" + std::to_string(horizontalSum); // Only horizontal
                if (horizontalSum == 0) clue = std::to_string(verticalSum) + "/"; // Only vertical

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
void printEmptyBoard(const std::vector<std::vector<std::string>> &emptyBoard) {
    // Determine the maximum width of any cell
    size_t maxWidth = 0;
    for (const auto &row : emptyBoard) {
        for (const auto &cell : row) {
            maxWidth = std::max(maxWidth, cell.size());
        }
    }

    // Print each row with aligned columns
    for (const auto &row : emptyBoard) {
        for (const auto &cell : row) {
            std::cout << std::setw(static_cast<int>(maxWidth)) << cell << " "; // Align each cell and add space
        }
        std::cout << std::endl;
    }
}

void saveEmptyBoardToFile(const std::vector<std::vector<std::string>> &emptyBoard, const std::string &filename) {
    std::ofstream outFile(filename);
    if (!outFile) {
        std::cerr << "Error opening file for writing!" << std::endl;
        return;
    }

    size_t maxWidth = 0;
    for (const auto &row : emptyBoard) {
        for (const auto &cell : row) {
            maxWidth = std::max(maxWidth, cell.size());
        }
    }

    for (const auto &row : emptyBoard) {
        for (size_t j = 0; j < row.size(); ++j) {
            outFile << std::setw(static_cast<int>(maxWidth)) << row[j];
            if (j < row.size() - 1) {
                outFile << " ";
            }
        }
        outFile << std::endl;
    }

    outFile.close();
}

int main() {
    srand(time(nullptr));

    for (auto size: gridSizes) {
        auto populationSize = size * populationSizeFactor;
        auto generations = size * nGenerationsFactor;
        // Run evolutionary algorithm
        Solution best = evolutionaryAlgorithm(size, populationSize , generations);
        // Derive empty board
        std::vector<std::vector<std::string>> emptyBoard = deriveEmptyBoard(best.grid);

        if (debugOutput) {
            for (const auto &row : best.grid) {
                for (int cell : row) {
                    std::cout << (cell == -1 ? "X" : std::to_string(cell)) << " ";
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;
            for (const auto &row : best.grid) {
                for (const int cell : row) {
                    std::cout << (cell == -1 ? "X" : "O") << " ";
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;
            printEmptyBoard(emptyBoard);
        }

        saveEmptyBoardToFile(emptyBoard, "output_" + std::to_string(size) + ".txt");
    }

    return 0;
}