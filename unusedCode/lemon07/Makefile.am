EXTRA_DIST += \
	lemon/Makefile \
	lemon/lemon.pc.in

pkgconfig_DATA += lemon/lemon.pc

lib_LTLIBRARIES += lemon/libemon.la

lemon_libemon_la_SOURCES = \
	lemon/arg_parser.cc \
	lemon/lp_base.cc \
	lemon/lp_utils.cc \
	lemon/lp_skeleton.cc \
	lemon/base.cc \
	lemon/color.cc \
	lemon/eps.cc \
	lemon/bits/mingw32_time.cc \
	lemon/random.cc

lemon_libemon_la_CXXFLAGS = $(GLPK_CFLAGS) $(CPLEX_CFLAGS) $(SOPLEX_CXXFLAGS)
lemon_libemon_la_LDFLAGS = $(GLPK_LIBS) $(CPLEX_LIBS) $(SOPLEX_LIBS)

if HAVE_GLPK
lemon_libemon_la_SOURCES += lemon/lp_glpk.cc
lemon_libemon_la_SOURCES += lemon/mip_glpk.cc
endif

if HAVE_CPLEX
lemon_libemon_la_SOURCES += lemon/lp_cplex.cc
lemon_libemon_la_SOURCES += lemon/mip_cplex.cc
endif

if HAVE_SOPLEX
lemon_libemon_la_SOURCES += lemon/lp_soplex.cc
endif

lemon_HEADERS += \
	lemon/arg_parser.h \
	lemon/bellman_ford.h \
	lemon/bfs.h \
	lemon/bin_heap.h \
	lemon/bipartite_matching.h \
	lemon/bpugraph_adaptor.h \
	lemon/bucket_heap.h \
	lemon/capacity_scaling.h \
	lemon/circulation.h \
	lemon/color.h \
	lemon/config.h \
	lemon/concept_check.h \
	lemon/cost_scaling.h \
	lemon/counter.h \
	lemon/csp.h \
	lemon/cycle_canceling.h \
	lemon/dag_shortest_path.h \
	lemon/dfs.h \
	lemon/dijkstra.h \
	lemon/dinitz_sleator_tarjan.h \
	lemon/dist_log.h \
	lemon/dim2.h \
	lemon/dimacs.h \
	lemon/dynamic_tree.h \
	lemon/edge_set.h \
	lemon/edmonds_karp.h \
	lemon/elevator.h \
	lemon/eps.h \
	lemon/error.h \
	lemon/euler.h \
	lemon/fib_heap.h \
	lemon/floyd_warshall.h \
	lemon/fredman_tarjan.h \
	lemon/full_graph.h \
	lemon/graph_adaptor.h \
	lemon/graph_reader.h \
	lemon/graph_to_eps.h \
	lemon/graph_utils.h \
	lemon/graph_writer.h \
	lemon/grid_ugraph.h \
	lemon/goldberg_tarjan.h \
	lemon/gomory_hu_tree.h \
	lemon/hao_orlin.h \
	lemon/hypercube_graph.h \
	lemon/iterable_maps.h \
	lemon/johnson.h \
	lemon/kruskal.h \
	lemon/lemon_reader.h \
	lemon/lemon_writer.h \
	lemon/list_graph.h \
	lemon/lp.h \
	lemon/lp_base.h \
	lemon/lp_cplex.h \
	lemon/lp_glpk.h \
	lemon/lp_skeleton.h \
	lemon/lp_soplex.h \
	lemon/lp_utils.h \
	lemon/map_iterator.h \
	lemon/maps.h \
	lemon/math.h \
	lemon/matrix_maps.h \
	lemon/max_matching.h \
	lemon/min_cost_arborescence.h \
	lemon/min_cost_flow.h \
	lemon/min_cost_max_flow.h \
	lemon/min_mean_cycle.h \
	lemon/mip_glpk.h \
	lemon/mip_cplex.h \
	lemon/nagamochi_ibaraki.h \
	lemon/nauty_reader.h \
	lemon/network_simplex.h \
	lemon/path.h \
	lemon/path_utils.h \
	lemon/planarity.h \
	lemon/polynomial.h \
	lemon/preflow.h \
	lemon/prim.h \
	lemon/pr_bipartite_matching.h \
	lemon/radix_heap.h \
	lemon/radix_sort.h \
	lemon/random.h \
	lemon/refptr.h \
	lemon/simann.h \
	lemon/smart_graph.h \
	lemon/static_graph.h \
	lemon/steiner.h \
	lemon/sub_graph.h \
	lemon/suurballe.h \
	lemon/tabu_search.h \
	lemon/time_measure.h \
	lemon/tolerance.h \
	lemon/topology.h \
	lemon/ugraph_adaptor.h \
	lemon/unionfind.h

bits_HEADERS += \
	lemon/bits/alteration_notifier.h \
	lemon/bits/array_map.h \
	lemon/bits/base_extender.h \
	lemon/bits/bezier.h \
	lemon/bits/debug_map.h \
	lemon/bits/default_map.h \
	lemon/bits/edge_set_extender.h \
	lemon/bits/graph_adaptor_extender.h \
	lemon/bits/graph_extender.h \
	lemon/bits/invalid.h \
	lemon/bits/item_reader.h \
	lemon/bits/item_writer.h \
	lemon/bits/lp_id.h \
	lemon/bits/map_extender.h \
	lemon/bits/mingw32_time.h \
	lemon/bits/path_dump.h \
	lemon/bits/traits.h \
	lemon/bits/utility.h \
	lemon/bits/variant.h \
	lemon/bits/vector_map.h

concept_HEADERS += \
	lemon/concept_check.h \
	lemon/concepts/bpugraph.h \
	lemon/concepts/graph.h \
	lemon/concepts/graph_components.h \
	lemon/concepts/heap.h \
	lemon/concepts/maps.h \
	lemon/concepts/matrix_maps.h \
	lemon/concepts/path.h \
	lemon/concepts/ugraph.h
