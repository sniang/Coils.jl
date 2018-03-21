using PyPlot

push!(LOAD_PATH, pwd())
using Coils
using CoilsPlot

#=
If it doesn't work, use the command lines :
julia -e 'Pkg.update()'
julia -e 'Pkg.add("LightGraphs")'
julia -e 'Pkg.add("ProgressMeter")'

Introduction
The coil design method works on a predefined grid. The grid is defined by a set
of points in space (vertices) and a directed graph g defining which are connected
with one another. A connection between two vertices means, that a wire of a coil
can be laid in a straight line between them.
It is absolutely possible to supply vertices (an array of arrays of x, y, z
coordinates) and the graph manually. There is also fuction cuboid_system that
defines a cuboid system, given it size (x, y and z) and the number of tiles along
x, y and z directions.
=#

g, vertex_positions = cuboid_system([1, 1, 1], [3, 3, 3])

#=
In CoilsPlot there is a number of functions for plotting the objects.
=#

figure(figsize = (6, 6))
plot_vertices(vertex_positions)
plot_edges(g, vertex_positions, standalone = false)

#=
The method optimizes the field in a finite set of Points Of Interest (poi).
They can be supplied manually or a helper function cuboid_poi may be used.
It has to be given the size, the position of the centre, the number of points
(along each direction).
The points may be created either in the whole volume (filled = true) or just on
the surface.
=#

poi = cuboid_poi([0.75, 0.75, 0.75], [0.0, 0.0, 0.0], [10, 10, 10], filled = false)
figure(figsize = (6, 6))
plot_system(g, vertex_positions, poi)

#=
The goal field is defined as a function of space x (an array of size 3 - x, y
and z), returning an array of size 3 - the x, y and z components of the field.
=#

Bgoal(x) = [0, 100e-6, 0]

#=
Next, the graph is analysed to find cells (or tiles) - smallest loops to be found in the graph. Each will be treated as a tile coil.
=#
cells = find_cells(g)
figure(figsize = (8, 8))
plot_vertices(vertex_positions, labels = false)
plot_cells(cells, vertex_positions)

#=
Now the matrix of the system is calculated. It is a matrix of proportionality
constants between the current in each cell (columns) and the magnetic field in
each of the points of interest, in each direction (rows). It is calculated using
the Biot-Savart law.
=#

M = system_matrix(poi, vertex_positions, cells)

#=
We evaluate the goal magnetic field in each of the points of interest and
concatenate the values into one big vector.
=#

Bpoi = vcat(Bgoal.(poi)...)

#=
Now we solve the system. We look for the optimal currents in each of the tile
coils, that will reproduce the goal magnetic field in the points of interest
best. It is a linear least-squares problem.
=#

optI = M \ Bpoi
figure(figsize = (10, 10))
plot_vertices(vertex_positions, labels = false)
current_norm = 1000 / maximum(abs, optI)
plot_cells(cells, vertex_positions, optI * current_norm)

#=
We proceed to simplify the system. The first step is for each edge to  add the
currents of the adjacent tiles. We get the net current flowing along each edge.
=#

edgecurrents = find_edgecurrents(g, cells, optI)
figure(figsize = (8, 8))
plot_vertices(vertex_positions, labels = false)
plot_edge_currents(g, edgecurrents, vertex_positions)

#=
The result is a complicated net of current. We decompose the net into what we call simple loops.
=#
simpleloops, simpleloopscurrents = find_all_simpleloops(g, edgecurrents)
#=
Now we plot the simple loops, together with the corresponding currents.
=#
figure(figsize = (6, 6))
plot_vertices(vertex_positions, labels = false)
plot_loop(simpleloops[1], vertex_positions, color = "orange")
plot_loop(simpleloops[2], vertex_positions, color = "orange")
title("Current: $(simpleloopscurrents[1] * current_norm)")
figure(figsize = (6, 6))
plot_vertices(vertex_positions, labels = false)
plot_loop(simpleloops[3], vertex_positions, color = "green")
plot_loop(simpleloops[4], vertex_positions, color = "green")
title("Current: $(simpleloopscurrents[3] * current_norm)")
figure(figsize = (6, 6))
plot_vertices(vertex_positions, labels = false)
plot_loop(simpleloops[5], vertex_positions, color = "blue")
plot_loop(simpleloops[6], vertex_positions, color = "blue")
title("Current: $(simpleloopscurrents[5] * current_norm)")
figure(figsize = (6, 6))
plot_vertices(vertex_positions, labels = false)
plot_loop(simpleloops[7], vertex_positions, color = "indigo")
plot_loop(simpleloops[8], vertex_positions, color = "indigo")
title("Current: $(simpleloopscurrents[7] * current_norm)")
figure(figsize = (6, 6))
plot_vertices(vertex_positions, labels = false)
plot_loop(simpleloops[9], vertex_positions, color = "grey")
plot_loop(simpleloops[10], vertex_positions, color = "grey")
title("Current: $(simpleloopscurrents[9] * current_norm)")

#=
We can characterise the performance of the design. For example by plotting the
histogram of the difference of the goal field and the one generated with the design.
=#
plot_deviation_histogram(poi, simpleloops, simpleloopscurrents, vertex_positions, Bgoal)
#=
Or we can plot the map of the difference between the goal field an the one
produced by the system.
=#
plot_loops_field(cells, optI, vertex_positions, [1, 2, 3], 0;
        n = 50, spanA = [-0.5, 0.5], spanB = [-0.5, 0.5], Bref = Bgoal, levels = linspace(-10, 10, 21))
