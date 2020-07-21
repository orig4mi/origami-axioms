# Axioms

`axioms` is a Python implementations of the so-called origami axioms, which are elementary fold operations on a sheet of paper that model the geometry of origami constructions [[1]](#ref1). Each axiom solves a set of incidence constraints between given points and lines in the plane with a single fold line.

## Usage

Import `axioms` as a module. It defines classes for points and lines, and each of the eight axioms is defined as a function.

The classes are: 

- `Point(x, y)`: x and y are the coordinates of the point.
- `Line(a, b, c)`: the line is the defined as ax + by + c = 0.

The axioms are:

1. `axiom_1(p, q)`: given two distinct points p and q, fold to place p onto q.
1. `axiom_2(m, n)`: given two distinct lines m and n, fold to align m and n.
1. `axiom_3(m)`: fold along a given line m.
1. `axiom_4(p, q)`: given two distinct points p and q, fold along a line passing through p and q.
1. `axiom_5(p, m)`: given a line m and a point q, fold along a line passing through p to reflect half of m onto its other half.
1. `axiom_6(p, q, m)`: given a line m, a point p not on m and a point q, fold along a line passing through q to place p onto m.
1. `axiom_7(p, q, m, n)`: given two lines m and n, a point p not on m, and a point q not on n, where m and n are distinct or pa and q are distinct, fold to place p onto m, and q onto n.
1. `axiom_8(p, m,n )`: given two line m and n, and a point not on m, fold to place p onto m, and to reflect half of n onto its other half.

If the module is run as a script, then it will compute the cubic root of 2 by applying axiom 7 [[2]](#ref2).

## References

<a name="ref1">[1]</a> J. C. Lucero, ["On the elementary single-fold operations of origami: reflections and incidence constraints on the plane"](http://forumgeom.fau.edu/FG2017volume17/FG201725.pdf), Forum Geometricorum 17, 207-221 (2017).  
<a name="ref2">[2]</a> J. C. Lucero, ["Existence of a solution for Beloch's fold"](https://www.tandfonline.com/doi/abs/10.1080/0025570X.2019.1526591?journalCode=umma20), Mathematics Magazine 29, 24-31 (2019). 

