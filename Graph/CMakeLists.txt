aux_source_directory(. GRAPH_SRC_LIST)
add_library(Graph ${GRAPH_SRC_LIST})
target_compile_definitions(Graph PRIVATE GENG_MAIN=geng_main) # for NAUTY utility "geng.c" which is in library
