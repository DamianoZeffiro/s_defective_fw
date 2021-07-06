# Fast cluster detection in networks by first-order optimization
Given a graph _G_ an s-defective clique is a subset of vertices with at most s edges missing to form a clique. It can be proved that the problem of finding maximal s-defective cliques admits a continuous cubic formulation of the form

         min f(x,y)
    s.t. (x,y) in \Delta_{|V| - 1} x D,
    
with \Delta_{|V| - 1} the |V| - 1 dimensional simplex, for V set of vertices in G, and D a suitable polytope. The variable x represents a set of weights for the vertices, while the variable y represents a set of edges to be added to an s-defective clique in order to get a clique. 

This repository contains two first order optimization methods to find a local solution for the problem above: the Frank Wolfe method with in face directions (FDFW), and the Frank Wolfe method for s-defective clique FWdc. Both are variants of the classic Frank Wolfe (FW) method, which at every iteration moves toward the vertex of the domain maximizing the scalar product with the current gradient. The FDFW is a generic scheme for constrained optimization on compact and convex sets which at every step considers a direction in the minimal face containing the current point as an alternative to the FW direction. The FWdc instead is a tailored variant designed to exploit the block product domain property of the problem. It updates the x variable with a FDFW step, and the y variable with a FW step having maximal feasible stepsize. 

After a finite number of iterations, both methods find a point (x, y) such that the support of x is an s-defective clique.  The FWdc showed a better performance in all the instances we tested. 

## Reference paper

I. Bomze, F. Rinaldi, D. Zeffiro (2021). _Fast cluster detection in networks by first-order optimization_. Pre-print available at <https://arxiv.org/abs/2103.15907>.

In order to speed up the methods, we applied the Short Step Chain procedure described in F. Rinaldi, D. Zeffiro (2020), _Avoiding bad steps in Frank Wolfe variants_, pre_print available at https://arxiv.org/abs/2012.12737.

## Authors

* Immanuel M. Bomze  (e-mail: [immanuel.bomze@univie.ac.at](mailto: immanuel.bomze@univie.ac.at))
* Francesco Rinaldi (e-mail: [rinaldi@math.unipd.it](mailto:rinaldi@math.unipd.it))
* Damiano Zeffiro (e-mail: [zeffiro@math.unipd.it](mailto:zeffiro@math.unipd.it))

## Licensing

s_defective_fw is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
s_defective_fw is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with s_defective_fw. If not, see <http://www.gnu.org/licenses/>.

## How to use the FDFW and the FWdc

1. This directory should contain the following files:
    * `main.m`,
    * `FDFW.m`,
    * `FWdc.m`,
    * `Qaugp.m`,
    * `buildface.m`,
    * `clique_init2.m`,
    * `embedfs.m`,
    * `graphcreator2.m`,
    * `inittable.m`,
    * `missingedgecount.m`,
    * `partials.m`,
    * `updatetable.m`,
    * `merge.m`,
    * `cliquenames.mat`,
    * `README.md`,
beside the folder “instances” with graphs from the 2nd DIMACS implementation challenge. 

2. See the file `main.m` for an example of how to call the FDFW and the FWdc to find s-defective cliques starting from random points, and print a summary of the results to a table.

3. The main algorithms are written in ‘FDFW.m’ and ‘FWdc.m’, while the remaining functions are auxiliary. ‘cliquenames.mat’ contains the list of instances tested in the article. 
