/*
 *  glfem.h
 *  Library for LEPL1110 : Finite Elements for dummies
 *
 *  Copyright (C) 2017 UCL-EPL : Vincent Legat
 *  All rights reserved.
 *
 *  Pour GLFW (version utilisée 3.1)
 *  Pour l'installation de la librairie, voir http://www.glfw.org/
 *
 */

#ifndef _GLFEM_H_
#define _GLFEM_H_

#define GLFW_INCLUDE_GLU
#include <GLFW/glfw3.h>

#include "fem.h"

/**
 * Draws a colored element.
 * 
 * @param x Array of x coordinates of the element's nodes.
 * @param y Array of y coordinates of the element's nodes.
 * @param u Array of values at the element's nodes.
 * @param n Number of nodes in the element.
 */
void        glfemDrawColorElement(float *x, float *y, double *u, int n);

/**
 * Draws an element.
 * 
 * @param x Array of x coordinates of the element's nodes.
 * @param y Array of y coordinates of the element's nodes.
 * @param n Number of nodes in the element.
 */
void 		    glfemDrawElement(float *x, float *y, int n);

/**
 * Draws nodes.
 * 
 * @param x Array of x coordinates of the nodes.
 * @param y Array of y coordinates of the nodes.
 * @param n Number of nodes.
 */
void 		    glfemDrawNodes(double* x, double* y,int n);

/**
 * Gets the current action.
 * 
 * @return The current action.
 */
int         glfemGetAction();

/**
 * Reshapes the window.
 * 
 * @param theNodes The nodes to be reshaped.
 * @param width The new width of the window.
 * @param height The new height of the window.
 */
void 		    glfemReshapeWindows(femNodes *theNodes, int width, int height);

/**
 * Plots a field on the mesh.
 * 
 * @param theMesh The mesh to plot the field on.
 * @param u Array of values at the mesh's nodes.
 */
void 		    glfemPlotField(femMesh *theMesh, double *u);

/**
 * Plots the mesh.
 * 
 * @param theMesh The mesh to be plotted.
 */
void 		    glfemPlotMesh(femMesh *theMesh);

/**
 * Plots the domain.
 * 
 * @param theDomain The domain to be plotted.
 */
void        glfemPlotDomain(femDomain *theDomain);

/**
 * Displays a message.
 * 
 * @param message The message to be displayed.
 */
void glfemMessage(char *message);

/**
 * Draws a message at a specified position.
 * 
 * @param h The horizontal position.
 * @param v The vertical position.
 * @param message The message to be drawn.
 */
void glfemDrawMessage(int h, int v, char *message);

/**
 * Sets the raster size.
 * 
 * @param width The width of the raster.
 * @param height The height of the raster.
 */
void glfemSetRasterSize(int width, int height);

/**
 * Initializes the GLFW window.
 * 
 * @param windowName The name of the window.
 * @return The initialized GLFW window.
 */
GLFWwindow* glfemInit(char *windowName);

/**
 * Key callback function for GLFW.
 * 
 * @param self The GLFW window.
 * @param key The key that was pressed.
 * @param scancode The system-specific scancode of the key.
 * @param action The action (press, release, repeat).
 * @param mods Bit field describing which modifier keys were held down.
 */
static void glfemKeyCallback(GLFWwindow* self,int key,int scancode,int action,int mods);




#endif