#pragma once

/*
 *  LBM method with complex domains
 *  
 *  Method 1)   any set of lattice points
 *              advantage: avoid large void spaces
 *              drawback:  needs indirection lists for neighborhood information
 *
 *  Method 2)   a collection of rectangular domains
 *              drawback:  less flexible for complex geometries
 *              advantage: does not require indirection lists
 *
 *  Currently, we have two types of cell: liquid and solid.
 *  All border cells of the domain are initialized as solid.
 *  Moving populations from liquid cell end up in solid cells after propagate
 *  has been called inside the collide loop.
 *  Then, in a second loop, we have to reflect those populations
 *  (, treat wall speeds) and propage them.
 *
 *  In a setup with a collection of rectangular domains,
 *  we need 1 point (or 1 layer) of overlap.
 *  The cells of overlap describe the same point, but in both domains
 *  so we have to treat them as liquid in one but not in the other domain.
 *  
 *   Sketch: connect two wall-enclosed domains
 *  
 *  +---------------------------------+
 *  |wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww|
 *  |W Domain #0                     W|
 *  |W                               W|
 *  |W                               W+-------------------------+
 *  |W                               Wwwwwwwwwwwwwwwwwwwwwwwwwww|
 *  |W                                |              Domain #1 W|
 *  |W                                |                        W|
 *  |W                                |                        W|
 *  |W                                |                        W|
 *  |W                                |                        W|
 *  |Wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww                       W|
 *  +---------------------------------+W                       W|
 *                                    |W                       W|
 *                                    |W                       W|
 *                                    |WwwwwwwwwwwwwwwwwwwwwwwwW|
 *                                    +-------------------------+
 *  
 *  
 *  
 */
