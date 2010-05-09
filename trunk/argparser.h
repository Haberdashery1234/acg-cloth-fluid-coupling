#ifndef __ARG_PARSER_H__
#define __ARG_PARSER_H__

#include <cstring>
#include <cstdio>
#include <cassert>
#include <cstdlib>
#include <string>

#include "vectors.h"

// ================================================================================
// ================================================================================

class ArgParser
{

public:

    ArgParser()
    {
        DefaultValues();
    }

    ArgParser(int argc, char *argv[])
    {
        DefaultValues();

        for (int i = 1; i < argc; i++)
        {
            if (!strcmp(argv[i],"-cloth"))
            {
                i++;
                assert (i < argc);
                cloth_file = argv[i];
            }
            else if (!strcmp(argv[i],"-fluid"))
            {
                i++;
                assert (i < argc);
                fluid_file = argv[i];
            }
            else if (!strcmp(argv[i],"-both"))
            {
                i++;
                assert (i < argc);
                fluid_file = argv[i];
                i++;
                assert (i < argc);
                cloth_file = argv[i];
            }
            else if (!strcmp(argv[i],"-size"))
            {
                i++;
                assert (i < argc);
                width = height = atoi(argv[i]);
            }
            else if (!strcmp(argv[i],"-timestep"))
            {
                i++;
                assert (i < argc);
                timestep = atof(argv[i]);
                assert (timestep > 0);
            }
            else if (!strcmp(argv[i],"-collision"))
            {
                i++;
                assert (i < argc);
                collision = argv[i];
                assert (strcmp(collision, "fluid") || strcmp(collision, "both"));
            }
            else
            {
                printf ("whoops error with command line argument %d: '%s'\n",i,argv[i]);
                assert(0);
            }
        }
    }

    // ===================================
    // ===================================

    void DefaultValues()
    {
        collision = "none";
        width = 500;
        height = 500;

        timestep = 0.01;
        animate = false;

        particles = true;
        velocity = true;
        force = true;

        edge_velocity = 0;
        dense_velocity = 0;

        surface = false;
        isosurface = 0.7;

        wireframe = false;
        bounding_box = true;
        cubes = false;
        pressure = false;
        updraft = false;

        gravity = Vec3f(0,-9.8,0);
    }

    // ===================================
    // ===================================
    // REPRESENTATION
    // all public! (no accessors)

    std::string cloth_file;
    std::string fluid_file;
    std::string collision;
    int width;
    int height;

    // animation control
    float timestep;
    bool animate;
    Vec3f gravity;

    // display option toggles
    // (used by both)
    bool particles;
    bool velocity;
    bool surface;
    bool bounding_box;

    // used by cloth
    bool force;
    bool wireframe;
    bool updraft;

    // used by fluid
    int edge_velocity;
    int dense_velocity;
    float isosurface;
    bool cubes;
    bool pressure;
};

#endif
