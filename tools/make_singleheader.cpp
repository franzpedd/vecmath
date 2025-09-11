#include <iostream>
#include <string.h>

typedef struct content_node_t
{
    const char* file;
    int start;
    int end;
} content_node_t;

void fprintf_content_node(FILE* output, content_node_t* node)
{
    FILE* file = fopen(node->file, "r+");
    if (!file) return;

    char buffer[256];
    int currentLine = 0;
    
    while (currentLine < node->end) {
        memset(buffer, 0, sizeof(buffer));
        fgets(buffer, sizeof(buffer), file);

        if (currentLine >= node->start && currentLine < node->end) {
            fprintf(output, "%s", buffer);
        }
        currentLine++;
    }

    fclose(file);
}

void create_header_only()
{
    // define
    char header[] = 
    {
        "#ifndef VECMATH_INCLUDED\n"
        "#define VECAMTH_INCLUDED\n\n"
    };

    char libraries[] = {
        "#include <math.h>\n"
        "#include <string.h>\n\n"
        "#define VECMATH_REQUESTING_HEADER_ONLY // will make functions static \n\n"
    };

    char separator[] = "// functions definitions\n\n";

    content_node_t types; types.start = 3; types.end = 116; types.file = "../vecmath_types.h";
    content_node_t defines; defines.start = 3; defines.end = 19; defines.file = "../vecmath_defines.h";
    content_node_t basic; basic.start = 4; basic.end = 862; basic.file = "../vecmath_basic_op.c";
    content_node_t vec; vec.start = 3; vec.end = 729; vec.file = "../vecmath_vec_op.c";
    content_node_t mat; mat.start = 5; mat.end = 1841; mat.file = "../vecmath_mat_op.c";
    content_node_t quat; quat.start = 3; quat.end = 425; quat.file = "../vecmath_quat_op.c";
    content_node_t util; util.start = 3; util.end = 130; util.file = "../vecmath_util.c";
    
    char footer[] = {
        "\n\n#endif // VECAMTH_INCLUDED\n"
    };

    // write
    FILE* outputFile = fopen("../headeronly/vecmath.h", "w");
    if (!outputFile) {
        printf("File headeronly/vecmath.h is non-existent, please created and re-run\n");
        return;
    }

    fprintf(outputFile, "%s", header);
    fprintf(outputFile, "%s", libraries);
    fprintf(outputFile, "%s", separator);

    fprintf_content_node(outputFile, &types);
    fprintf_content_node(outputFile, &defines);
    fprintf_content_node(outputFile, &basic);
    fprintf_content_node(outputFile, &vec);
    fprintf_content_node(outputFile, &mat);
    fprintf_content_node(outputFile, &quat);
    fprintf_content_node(outputFile, &util);

    fprintf(outputFile, "%s", footer);
    fclose(outputFile);
}

int main(int args, char** argv)
{
    create_header_only();
    return 0;
}