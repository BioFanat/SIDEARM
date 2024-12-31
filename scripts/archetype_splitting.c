#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <fcntl.h>
#include <unistd.h>
#include <ctype.h>
#include <stdbool.h>

#define BUFFER_SIZE 8192
#define MAX_LINE_LENGTH 1024
#define MAX_CLUSTERS 512  // Safe upper bound for TF clusters
#define MAX_PATH_LENGTH 256

typedef struct {
    char name[64];
    FILE *plus_handle;
    FILE *minus_handle;
} Cluster;

typedef struct {
    Cluster *clusters;
    int count;
    char **names;  // For sorting and printing
} ClusterManager;

// Initialize cluster manager
ClusterManager* init_cluster_manager() {
    ClusterManager *cm = malloc(sizeof(ClusterManager));
    cm->clusters = malloc(MAX_CLUSTERS * sizeof(Cluster));
    cm->names = malloc(MAX_CLUSTERS * sizeof(char*));
    cm->count = 0;
    return cm;
}

// Get or create cluster handles
Cluster* get_cluster(ClusterManager *cm, const char *cluster_name, const char *output_dir) {
    // Search for existing cluster
    for (int i = 0; i < cm->count; i++) {
        if (strcmp(cm->clusters[i].name, cluster_name) == 0) {
            return &cm->clusters[i];
        }
    }
    
    // Create new cluster
    if (cm->count >= MAX_CLUSTERS) {
        fprintf(stderr, "Maximum number of clusters exceeded\n");
        exit(1);
    }
    
    Cluster *cluster = &cm->clusters[cm->count];
    // Copy and sanitize cluster name (replace '/' with '_')
    strncpy(cluster->name, cluster_name, 63);
    cluster->name[63] = '\0';
    for (char *p = cluster->name; *p; p++) {
        if (*p == '/') *p = '_';
    }
    
    // Create cluster directory
    char dir_path[MAX_PATH_LENGTH];
    int ret = snprintf(dir_path, MAX_PATH_LENGTH, "%s/%s", output_dir, cluster_name);
    if (ret < 0 || ret >= MAX_PATH_LENGTH) {
        fprintf(stderr, "Path too long for cluster directory\n");
        exit(1);
    }
    mkdir(dir_path, 0755);
    
    // Create plus strand file
    char file_path[MAX_PATH_LENGTH];
    ret = snprintf(file_path, MAX_PATH_LENGTH, "%s/%s_plus.bed", dir_path, cluster_name);
    if (ret < 0 || ret >= MAX_PATH_LENGTH) {
        fprintf(stderr, "Path too long for plus strand file\n");
        exit(1);
    }
    cluster->plus_handle = fopen(file_path, "w");
    if (!cluster->plus_handle) {
        fprintf(stderr, "Failed to create plus strand file: %s\n", file_path);
        exit(1);
    }
    
    // Create minus strand file
    ret = snprintf(file_path, MAX_PATH_LENGTH, "%s/%s_minus.bed", dir_path, cluster_name);
    if (ret < 0 || ret >= MAX_PATH_LENGTH) {
        fprintf(stderr, "Path too long for minus strand file\n");
        exit(1);
    }
    cluster->minus_handle = fopen(file_path, "w");
    if (!cluster->minus_handle) {
        fprintf(stderr, "Failed to create minus strand file: %s\n", file_path);
        fclose(cluster->plus_handle);  // Clean up previously opened file
        exit(1);
    }
    
    // Store cluster name for sorting
    cm->names[cm->count] = strdup(cluster_name);
    cm->count++;
    
    return cluster;
}

// Process a single line
void process_line(char *line, ClusterManager *cm, const char *output_dir) {
    char contig[32], cluster[64], strand;
    int start, stop;
    float score;
    
    // Parse line
    char *token = strtok(line, "\t");
    if (!token) return;
    strncpy(contig, token, 31);
    contig[31] = '\0';
    
    token = strtok(NULL, "\t");
    if (!token) return;
    start = atoi(token);
    
    token = strtok(NULL, "\t");
    if (!token) return;
    stop = atoi(token);
    
    token = strtok(NULL, "\t");
    if (!token) return;
    // Copy and sanitize cluster name
    strncpy(cluster, token, 63);
    cluster[63] = '\0';
    for (char *p = cluster; *p; p++) {
        if (*p == '/') *p = '_';
    }
    
    token = strtok(NULL, "\t");
    if (!token) return;
    score = atof(token);
    
    token = strtok(NULL, "\t");
    if (!token) return;
    strand = token[0];
    
    // Get cluster handles and write output
    Cluster *c = get_cluster(cm, cluster, output_dir);
    FILE *out = (strand == '+') ? c->plus_handle : c->minus_handle;
    fprintf(out, "%s\t%d\t%d\t%c\t%f\n", contig, start, stop, strand, score);
}

// Compare function for qsort
int compare_strings(const void *a, const void *b) {
    return strcmp(*(const char **)a, *(const char **)b);
}

// Print sorted clusters
void print_clusters(ClusterManager *cm) {
    printf("\nProcessed clusters:\n");
    qsort(cm->names, cm->count, sizeof(char*), compare_strings);
    for (int i = 0; i < cm->count; i++) {
        printf("%s\n", cm->names[i]);
    }
}

// Clean up resources
void cleanup_cluster_manager(ClusterManager *cm) {
    for (int i = 0; i < cm->count; i++) {
        fclose(cm->clusters[i].plus_handle);
        fclose(cm->clusters[i].minus_handle);
        free(cm->names[i]);
    }
    free(cm->clusters);
    free(cm->names);
    free(cm);
}

int main(int argc, char **argv) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <bed_file> <output_dir>\n", argv[0]);
        return 1;
    }
    
    // Create output directory
    mkdir(argv[2], 0755);
    
    // Initialize cluster manager
    ClusterManager *cm = init_cluster_manager();
    
    // Process input file
    FILE *input = fopen(argv[1], "r");
    if (!input) {
        perror("Failed to open input file");
        return 1;
    }
    
    // Skip header if present
    char first_line[MAX_LINE_LENGTH];
    if (fgets(first_line, MAX_LINE_LENGTH, input)) {
        if (first_line[0] != 't' &&  // not track
            first_line[0] != '#' &&  // not comment
            first_line[0] != 'b') {  // not browser
            rewind(input);
        } else if (first_line[0] == 't' && strncmp(first_line, "track", 5) != 0) {
            rewind(input);  // starts with 't' but isn't "track"
        } else if (first_line[0] == 'b' && strncmp(first_line, "browser", 7) != 0) {
            rewind(input);  // starts with 'b' but isn't "browser"
        }
    }
    
    char buffer[BUFFER_SIZE];
    char line_buffer[MAX_LINE_LENGTH];
    int line_pos = 0;
    
    char *line = NULL;
    size_t len = 0;
    ssize_t read;

    while ((read = getline(&line, &len, input)) != -1) {
        if (read > 0) {
            // Remove newline if present
            if (line[read-1] == '\n') {
                line[read-1] = '\0';
            }
            process_line(line, cm, argv[2]);
        }
    }

    free(line);
    
    // Process last line if needed
    if (line_pos > 0) {
        line_buffer[line_pos] = '\0';
        process_line(line_buffer, cm, argv[2]);
    }
    
    fclose(input);
    print_clusters(cm);
    cleanup_cluster_manager(cm);
    
    return 0;
}