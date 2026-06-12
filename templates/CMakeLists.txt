cmake_minimum_required(VERSION 3.20)
project(nom_projet VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Flags de compilation
# Release : cmake -B build -DCMAKE_BUILD_TYPE=Release
# Debug   : cmake -B build -DCMAKE_BUILD_TYPE=Debug
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Debug)
endif()

add_compile_options(-Wall -Wextra -Wpedantic)

# Répertoire d'inclusion partagé entre les TPs
include_directories(include)

# ---------------------------------------------------------------------------
# Bibliothèque principale
# Décommenter et compléter au fil des TPs pour exposer le code réutilisable.
# ---------------------------------------------------------------------------
# file(GLOB_RECURSE LIB_SOURCES "src/*.cpp")
# add_library(${PROJECT_NAME}_lib STATIC ${LIB_SOURCES})
# target_include_directories(${PROJECT_NAME}_lib PUBLIC include)

# ---------------------------------------------------------------------------
# Exécutable principal (optionnel)
# ---------------------------------------------------------------------------
# add_executable(${PROJECT_NAME} src/main.cpp)
# target_link_libraries(${PROJECT_NAME} PRIVATE ${PROJECT_NAME}_lib)

# ---------------------------------------------------------------------------
# Tests — ajouter un bloc par TP
# Compilation : cmake --build build
# Exécution   : ctest --test-dir build --output-on-failure
# ---------------------------------------------------------------------------
enable_testing()

# Exemple TP01 — décommenter et adapter :
# add_executable(test_tp01 subjects/TP01_nom/test_skeleton.cpp)
# target_include_directories(test_tp01 PRIVATE include)
# # target_link_libraries(test_tp01 PRIVATE ${PROJECT_NAME}_lib)
# add_test(NAME TP01 COMMAND test_tp01)

# Exemple TP02 :
# add_executable(test_tp02 subjects/TP02_nom/test_skeleton.cpp)
# target_include_directories(test_tp02 PRIVATE include)
# add_test(NAME TP02 COMMAND test_tp02)
