/*
 * Title: Sextactic and Type-9 Points on the Fermat Cubic and Associated Objects
 * Authors: Łukasz Merta and Maciej Zięba
 */

option(noredefine);
option(noloadLib);

LIB "elim.lib";

// Load necessary procedures for the code below
<"src/procedures.cpp";

/*****************************************************************************/
/*!
*  Type 9 points on the Fermat cubic
*
*****************************************************************************/

// Compute the 9-type points on the Fermat cubic
< "src/points.cpp";

// Compute polynomial of degree 24 on which lies 9-type points
< "src/ideal_of_points.cpp";

/*****************************************************************************/
/*!
*  Tangents to the Fermat cubic
*
*****************************************************************************/

// Compute tangents to the Fermat cubic at 6-type points
< "src/tangents_6_type.cpp";

// Compute tangents to the Fermat cubic at 9-type points
< "src/tangents_9_type.cpp";

/*****************************************************************************/
/*!
*  Type 9 points and a peculiar arrangement of conics
*
*****************************************************************************/

// Compute conics passing through 6-type points with multiplicity 3
< "src/conics_6_type.cpp";

// Compute conics passing through 9-type points with multiplicity 3
< "src/conics_9_type.cpp";




