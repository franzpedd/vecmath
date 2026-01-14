#define _CRT_SECURE_NO_WARNINGS

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

void create_header_file()
{
    // define
    char header[] = 
    {
        "#ifndef VECMATH_INCLUDED\n"
        "#define VECMATH_INCLUDED\n\n"
        "#define VECMATH_REQUESTING_HEADER_ONLY // will make functions static \n\n"
    };

    char separator[] = "// functions definitions\n\n";

    // headers
    content_node_t types; types.start = 0; types.end = 128; types.file = "../vecmath_types.h";
    content_node_t defines; defines.start = 0; defines.end = 32; defines.file = "../vecmath_defines.h";
    content_node_t basic; basic.start = 6; basic.end = 100; basic.file = "../vecmath_basic_op.h";
    content_node_t vec; vec.start = 6; vec.end = 112; vec.file = "../vecmath_vec_op.h";
    content_node_t mat; mat.start = 6; mat.end = 122; mat.file = "../vecmath_mat_op.h";
    content_node_t quat; quat.start = 6; quat.end = 60; quat.file = "../vecmath_quat_op.h";
    content_node_t ray; ray.start = 6; ray.end = 20; ray.file = "../vecmath_ray_op.h";
    content_node_t util; util.start = 6; util.end = 47; util.file = "../vecmath_util.h";
    
    char footer[] = {
        "#endif // VECMATH_INCLUDED\n\n"
        "/// @brief prevents circular dependency\n"
        "#ifdef VECMATH_IMPLEMENTATION\n"
        "#undef VECMATH_IMPLEMENTATION\n"
        "#include \"vecmath.c\"\n"
        "#endif\n"
    };

    // write
    FILE* outputFile = fopen("../headeronly/vecmath.h", "w");
    if (!outputFile) {
        printf("File headeronly/vecmath.h is non-existent, please created and re-run\n");
        return;
    }

    fprintf(outputFile, "%s", header);
    fprintf(outputFile, "%s", separator);

    fprintf_content_node(outputFile, &types);
    fprintf_content_node(outputFile, &defines);
    fprintf_content_node(outputFile, &basic);
    fprintf_content_node(outputFile, &vec);
    fprintf_content_node(outputFile, &mat);
    fprintf_content_node(outputFile, &quat);
    fprintf_content_node(outputFile, &util);
    fprintf_content_node(outputFile, &ray);

    fprintf(outputFile, "%s", footer);
    fclose(outputFile);
}

void create_source_file()
{
    // define
    char header[] = 
    {
        "#include \"vecmath.h\"\n\n"
        "#include <math.h>\n"
        "#include <string.h>\n\n"
    };

    char separator[] = "// functions implementation\n\n";

    content_node_t basic; basic.start = 4; basic.end = 870; basic.file = "../vecmath_basic_op.c";
    content_node_t vec; vec.start = 3; vec.end = 754; vec.file = "../vecmath_vec_op.c";
    content_node_t mat; mat.start = 5; mat.end = 1901; mat.file = "../vecmath_mat_op.c";
    content_node_t quat; quat.start = 5; quat.end = 425; quat.file = "../vecmath_quat_op.c";
    content_node_t ray; ray.start = 4; ray.end = 64; ray.file = "../vecmath_ray_op.c";
    content_node_t util; util.start = 5; util.end = 155; util.file = "../vecmath_util.c";

    // write
    FILE* outputFile = fopen("../headeronly/vecmath.c", "w");
    if (!outputFile) {
        printf("File headeronly/vecmath.c is non-existent, please created and re-run\n");
        return;
    }

    fprintf(outputFile, "%s", header);
    fprintf(outputFile, "%s", separator);

    fprintf_content_node(outputFile, &basic);
    fprintf_content_node(outputFile, &vec);
    fprintf_content_node(outputFile, &mat);
    fprintf_content_node(outputFile, &quat);
    fprintf_content_node(outputFile, &util);
    fprintf_content_node(outputFile, &ray);

    fclose(outputFile);
}

int main(int args, char** argv)
{
    create_header_file();
    create_source_file();
    return 0;
}