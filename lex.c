#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <stdint.h>
#include <stdbool.h>

// ============================================================================
// UTILITY FUNCTIONS
// ============================================================================

/**
 * Compute factorials up to n and store in fact array
 * @param n Maximum factorial to compute
 * @param fact Array to store factorials (must be size n+1)
 * @return 0 on success, -1 on overflow
 */
int compute_factorials(int n, unsigned long long *fact) {
    fact[0] = 1ULL;
    for (int i = 1; i <= n; ++i) {
        // Check for overflow before multiplication
        if (ULLONG_MAX / fact[i-1] < (unsigned long long)i) {
            return -1;
        }
        fact[i] = fact[i-1] * (unsigned long long)i;
    }
    return 0;
}

/**
 * Calculate P(a,b) = a! / (a-b)! (number of permutations)
 * @param fact Precomputed factorial array
 * @param a Total number of elements
 * @param b Number of elements to choose
 * @return Number of permutations P(a,b)
 */
unsigned long long P(unsigned long long *fact, int a, int b) {
    if (b <= 0) return 1ULL;
    return fact[a] / fact[a - b];
}

// ============================================================================
// PERMUTATION RANKING/UNRANKING
// ============================================================================

/**
 * Rank a k-permutation in lexicographic order among all k-permutations of [1..n]
 * @param perm The permutation to rank (elements 1..n)
 * @param n Total number of elements
 * @param k Length of permutation
 * @param fact Precomputed factorial array
 * @return Rank in range [0, P(n,k)-1]
 */
unsigned long long rank_k_perm(const int *perm, int n, int k, unsigned long long *fact) {
    bool used[256]; // Supports n up to 255 (practical limits)
    memset(used, 0, sizeof(bool) * (n + 1));
    unsigned long long rank = 0ULL;
    
    for (int i = 0; i < k; ++i) {
        int current_element = perm[i];
        int smaller_unused = 0;
        
        // Count how many smaller elements are still unused
        for (int v = 1; v < current_element; ++v) {
            if (!used[v]) {
                smaller_unused++;
            }
        }
        
        // Number of k-permutations with a smaller choice at this position
        unsigned long long permutations_with_smaller = P(fact, n - i - 1, k - i - 1);
        rank += (unsigned long long)smaller_unused * permutations_with_smaller;
        used[current_element] = true;
    }
    return rank;
}

/**
 * Unrank a k-permutation (inverse of rank_k_perm)
 * @param rank Rank in range [0, P(n,k)-1]
 * @param n Total number of elements
 * @param k Length of permutation
 * @param fact Precomputed factorial array
 * @param perm_out Output array to store the permutation
 */
void unrank_k_perm(unsigned long long rank, int n, int k, unsigned long long *fact, int *perm_out) {
    bool used[256];
    memset(used, 0, sizeof(bool) * (n + 1));
    
    for (int i = 0; i < k; ++i) {
        unsigned long long permutations_remaining = P(fact, n - i - 1, k - i - 1);
        unsigned long long element_index = (permutations_remaining == 0) ? 0 : (rank / permutations_remaining);
        rank = (permutations_remaining == 0) ? rank : (rank % permutations_remaining);

        // Pick the element_index-th unused number (0-based) from 1..n
        int chosen_element = -1;
        for (int v = 1; v <= n; ++v) {
            if (used[v]) continue;
            if (element_index == 0) {
                chosen_element = v;
                break;
            }
            element_index--;
        }
        perm_out[i] = chosen_element;
        used[chosen_element] = true;
    }
}

// ============================================================================
// PERMUTATION GENERATION
// ============================================================================

/**
 * Generate next permutation in lexicographic order (in-place)
 * @param a Array to permute (elements a[0..n-1])
 * @param n Length of array
 * @return 1 if next permutation produced, 0 if last permutation (no next)
 */
int next_permutation(int *a, int n) {
    // Find the largest index i such that a[i] < a[i+1]
    int i = n - 2;
    while (i >= 0 && a[i] >= a[i + 1]) {
        i--;
    }
    if (i < 0) {
        return 0; // Last permutation
    }
    
    // Find the largest index j such that a[i] < a[j]
    int j = n - 1;
    while (a[j] <= a[i]) {
        j--;
    }
    
    // Swap a[i] and a[j]
    int temp = a[i];
    a[i] = a[j];
    a[j] = temp;
    
    // Reverse the suffix a[i+1..n-1]
    int left = i + 1;
    int right = n - 1;
    while (left < right) {
        temp = a[left];
        a[left] = a[right];
        a[right] = temp;
        left++;
        right--;
    }
    return 1;
}

// ============================================================================
// MAIN ALGORITHM: LEXICOGRAPHICALLY SMALLEST UNIVERSAL CYCLE
// ============================================================================

/**
 * Generate lexicographically smallest universal cycle for (n-1)-permutations of [1..n]
 * @param n Number of elements (must be >= 1)
 * @param out_len Pointer to receive the length of the returned sequence
 * @return Pointer to int array containing the sequence (values 1..n), or NULL on error
 *         Caller must free() the returned pointer
 */
int *generate_lex_uc(int n, size_t *out_len) {
    if (n <= 0) {
        return NULL;
    }

    // Handle trivial cases
    if (n == 1) {
        *out_len = 1;
        int *result = malloc(sizeof(int));
        if (!result) {
            perror("malloc failed for n=1 case");
            return NULL;
        }
        result[0] = 1;
        return result;
    }
    
    if (n == 2) {
        *out_len = 2;
        int *result = malloc(2 * sizeof(int));
        if (!result) {
            perror("malloc failed for n=2 case");
            return NULL;
        }
        result[0] = 1;
        result[1] = 2;
        return result;
    }

    // Algorithm parameters
    int k = n - 2; // Vertex tuple length
    
    // Allocate and compute factorials
    unsigned long long *factorials = malloc((n + 1) * sizeof(unsigned long long));
    if (!factorials) {
        perror("malloc failed for factorials");
        return NULL;
    }
    
    if (compute_factorials(n, factorials) != 0) {
        fprintf(stderr, "Error: Factorial overflow or computation error for n=%d\n", n);
        free(factorials);
        return NULL;
    }

    // Calculate graph dimensions
    // Number of vertices: P(n, k) = n! / (n-k)! ; here n-k = 2 => vertices = n! / 2
    unsigned long long vertices_ull = P(factorials, n, k);
    if (vertices_ull > (unsigned long long)SIZE_MAX) {
        fprintf(stderr, "Error: Vertex count too large for memory (n=%d)\n", n);
        free(factorials);
        return NULL;
    }
    size_t num_vertices = (size_t)vertices_ull;

    // Number of edges (n-1)-permutations = P(n, n-1) = n!
    unsigned long long edges_ull = P(factorials, n, n - 1);
    if (edges_ull > (unsigned long long)SIZE_MAX) {
        fprintf(stderr, "Error: Edge count too large for memory (n=%d)\n", n);
        free(factorials);
        return NULL;
    }
    size_t num_edges = (size_t)edges_ull;

    // Allocate graph data structures
    // Each vertex has exactly n - k = 2 outgoing symbols (distinct), so store flat array of size V*2
    int *adjacency_flat = calloc(num_vertices * 2, sizeof(int));
    if (!adjacency_flat) {
        perror("calloc failed for adjacency_flat");
        free(factorials);
        return NULL;
    }
    
    unsigned char *adjacency_count = calloc(num_vertices, sizeof(unsigned char)); // Should be 2 at end
    if (!adjacency_count) {
        perror("calloc failed for adjacency_count");
        free(adjacency_flat);
        free(factorials);
        return NULL;
    }

    // Store vertex tuples for quick access (V x k ints)
    int *vertex_tuples = malloc(num_vertices * k * sizeof(int));
    if (!vertex_tuples) {
        perror("malloc failed for vertex_tuples");
        free(adjacency_flat);
        free(adjacency_count);
        free(factorials);
        return NULL;
    }
    
    unsigned char *vertices_initialized = calloc(num_vertices, sizeof(unsigned char));
    if (!vertices_initialized) {
        perror("calloc failed for vertices_initialized");
        free(vertex_tuples);
        free(adjacency_flat);
        free(adjacency_count);
        free(factorials);
        return NULL;
    }

    // Build adjacency by iterating through all permutations of 1..n (lexicographic)
    // For each perm p[0..n-1] take v = p[0..k-1] and a = p[n-2], append a to adj[v_index]
    // This enumerates all (n-1)-permutations exactly once
    int *current_permutation = malloc(n * sizeof(int));
    if (!current_permutation) {
        perror("malloc failed for current_permutation");
        free(vertices_initialized);
        free(vertex_tuples);
        free(adjacency_flat);
        free(adjacency_count);
        free(factorials);
        return NULL;
    }
    
    // Initialize with identity permutation
    for (int i = 0; i < n; ++i) {
        current_permutation[i] = i + 1;
    }

    bool first_permutation_processed = false;
    while (1) {
        // Extract vertex tuple v = p[0..k-1] and edge symbol a = p[n-2]
        unsigned long long vertex_index = rank_k_perm(current_permutation, n, k, factorials);
        
        if (!vertices_initialized[vertex_index]) {
            // Copy tuple to vertex storage
            for (int i = 0; i < k; ++i) {
                vertex_tuples[vertex_index * k + i] = current_permutation[i];
            }
            vertices_initialized[vertex_index] = 1;
        }
        
        int edge_symbol = current_permutation[n - 2];
        
        // Append edge symbol to adjacency list
        unsigned char adjacency_position = adjacency_count[vertex_index];
        if (adjacency_position >= 2) {
            fprintf(stderr, "Warning: Unexpected adjacency_count > 2 for vertex %llu\n", vertex_index);
            // Continue but avoid overflow
        } else {
            adjacency_flat[vertex_index * 2 + adjacency_position] = edge_symbol;
            adjacency_count[vertex_index] = adjacency_position + 1;
        }

        if (!first_permutation_processed) {
            first_permutation_processed = true;
        }
        
        if (!next_permutation(current_permutation, n)) {
            break;
        }
    }
    free(current_permutation);

    // Verify that every vertex has exactly 2 adjacencies and sort them for lexicographic order
    for (size_t i = 0; i < num_vertices; ++i) {
        if (adjacency_count[i] != 2) {
            fprintf(stderr, "Error: Vertex %zu has adjacency_count %u (expected 2) - algorithm error\n", 
                    i, adjacency_count[i]);
            // Continue to allow graceful failure
        } else {
            // Sort the two adjacency values so that smallest is taken first (lexicographic order)
            int adj0 = adjacency_flat[i * 2 + 0];
            int adj1 = adjacency_flat[i * 2 + 1];
            if (adj0 > adj1) {
                adjacency_flat[i * 2 + 0] = adj1;
                adjacency_flat[i * 2 + 1] = adj0;
            }
        }
    }

    // Prepare data structures for Hierholzer's algorithm
    // - path: stack of vertex indices (size <= E + 1)
    // - next_edge_idx[v]: index in [0,2) of next outgoing edge to consume
    // - edge_stack: stack of symbols we took to get from each vertex in path
    // - circuit: vector to collect popped edge labels (symbols)
    size_t path_capacity = num_edges + 2;
    size_t *vertex_path = malloc(path_capacity * sizeof(size_t));
    if (!vertex_path) {
        perror("malloc failed for vertex_path");
        goto cleanup_fail;
    }
    size_t path_size = 0;

    unsigned char *next_edge_index = calloc(num_vertices, sizeof(unsigned char));
    if (!next_edge_index) {
        perror("calloc failed for next_edge_index");
        goto cleanup_fail;
    }

    int *edge_symbol_stack = malloc((num_edges + 4) * sizeof(int));
    if (!edge_symbol_stack) {
        perror("malloc failed for edge_symbol_stack");
        goto cleanup_fail;
    }
    size_t edge_stack_size = 0;

    int *eulerian_circuit = malloc((num_edges + 4) * sizeof(int));
    if (!eulerian_circuit) {
        perror("malloc failed for eulerian_circuit");
        goto cleanup_fail;
    }
    size_t circuit_size = 0;

    // Start vertex = tuple (1,2,...,n-2)
    int *start_tuple = malloc(k * sizeof(int));
    if (!start_tuple) {
        perror("malloc failed for start_tuple");
        goto cleanup_fail;
    }
    
    for (int i = 0; i < k; ++i) {
        start_tuple[i] = i + 1;
    }
    unsigned long long start_vertex_index = rank_k_perm(start_tuple, n, k, factorials);
    free(start_tuple);

    // Push start vertex onto path
    vertex_path[path_size++] = (size_t)start_vertex_index;

    // Hierholzer's algorithm with lexicographic outgoing order
    while (path_size > 0) {
        size_t current_vertex = vertex_path[path_size - 1];
        
        if (next_edge_index[current_vertex] < adjacency_count[current_vertex]) {
            // Take smallest unused edge (lexicographic order)
            int edge_symbol = adjacency_flat[current_vertex * 2 + next_edge_index[current_vertex]];
            next_edge_index[current_vertex]++;

            // Compute next vertex tuple = current_vertex[1..k-1] + (edge_symbol)
            int next_vertex_tuple[64]; // Safe for reasonable n (k <= 64)
            if (k > 64) {
                fprintf(stderr, "Error: k=%d too large (max 64)\n", k);
                goto cleanup_fail;
            }
            
            for (int i = 0; i < k - 1; ++i) {
                next_vertex_tuple[i] = vertex_tuples[current_vertex * k + i + 1];
            }
            next_vertex_tuple[k - 1] = edge_symbol;
            unsigned long long next_vertex_index = rank_k_perm(next_vertex_tuple, n, k, factorials);

            // Push next vertex onto path
            vertex_path[path_size++] = (size_t)next_vertex_index;
            // Record the symbol taken
            edge_symbol_stack[edge_stack_size++] = edge_symbol;
        } else {
            // Backtrack: no more edges from current vertex
            path_size--;
            if (edge_stack_size > 0) {
                // Append last taken symbol to circuit
                eulerian_circuit[circuit_size++] = edge_symbol_stack[--edge_stack_size];
            }
        }
    }

    // Circuit contains edge-labels in reverse; reverse it to get correct order
    for (size_t i = 0; i < circuit_size / 2; ++i) {
        int temp = eulerian_circuit[i];
        eulerian_circuit[i] = eulerian_circuit[circuit_size - 1 - i];
        eulerian_circuit[circuit_size - 1 - i] = temp;
    }

    // Build final sequence: prefix (start tuple 1..k) + circuit
    unsigned long long factorial_n = factorials[n];
    size_t final_length = (size_t)factorial_n + (size_t)k;
    int *output_sequence = malloc(final_length * sizeof(int));
    if (!output_sequence) {
        perror("malloc failed for output_sequence");
        goto cleanup_fail;
    }

    // Add prefix (start tuple 1..k)
    for (int i = 0; i < k; ++i) {
        output_sequence[i] = i + 1;
    }
    
    // Append circuit
    for (size_t i = 0; i < circuit_size; ++i) {
        output_sequence[k + i] = eulerian_circuit[i];
    }

    *out_len = final_length;

    // Cleanup and return
    free(adjacency_flat);
    free(adjacency_count);
    free(vertex_tuples);
    free(vertices_initialized);
    free(factorials);
    free(vertex_path);
    free(next_edge_index);
    free(edge_symbol_stack);
    free(eulerian_circuit);
    return output_sequence;

cleanup_fail:
    // Free allocated memory on error
    if (adjacency_flat) free(adjacency_flat);
    if (adjacency_count) free(adjacency_count);
    if (vertex_tuples) free(vertex_tuples);
    if (vertices_initialized) free(vertices_initialized);
    if (factorials) free(factorials);
    if (vertex_path) free(vertex_path);
    if (next_edge_index) free(next_edge_index);
    if (edge_symbol_stack) free(edge_symbol_stack);
    if (eulerian_circuit) free(eulerian_circuit);
    return NULL;
}

// ============================================================================
// MAIN FUNCTION
// ============================================================================

int main(int argc, char **argv) {
    // Default value for n
    int n = 4;
    if (argc >= 2) {
        n = atoi(argv[1]);
    }
    
    if (n <= 0) {
        fprintf(stderr, "Error: n must be positive (got %d)\n", n);
        return 1;
    }

    // Generate the lexicographically smallest universal cycle
    size_t sequence_length = 0;
    int *universal_cycle = generate_lex_uc(n, &sequence_length);
    if (!universal_cycle) {
        fprintf(stderr, "Error: generate_lex_uc failed for n=%d\n", n);
        return 2;
    }
    
    // Account for the wraparound (remove the last n-2 elements that overlap with the beginning)
    sequence_length = sequence_length - (n - 2);

    // Output results
    printf("n=%d, length=%zu\n", n, sequence_length);
    for (size_t i = 0; i < sequence_length; ++i) {
        printf("%d", universal_cycle[i]);
    }
    printf("\n");

    // Cleanup
    free(universal_cycle);
    return 0;
}
