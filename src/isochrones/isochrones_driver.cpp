/*PGR-GNU*****************************************************************
File: isochrones_driver.cpp

Copyright (c) 2020 Vjeran Crnjak
vjeran@crnjak.xyz

------

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

 ********************************************************************PGR-GNU*/

#include "drivers/isochrones/isochrones_driver.h"

#include <deque>
#include <set>
#include <sstream>
#include <vector>

#include "cpp_common/basic_vertex.h"
#include "cpp_common/pgr_alloc.hpp"
#include "cpp_common/pgr_assert.h"

// // agg_cost at node, node Id, predecessor node id
// typedef std::tuple<double, int64_t, int64_t> pq_el;
// agg_cost at node, node Id
typedef std::tuple<double, int64_t> pq_el;

struct Pred {
  int64_t edgeId;
  int64_t predId;
  double cost;
};

std::vector<std::vector<const pgr_edge_t *>>
construct_adjacency_matrix(size_t n, const pgr_edge_t *edges,
                           size_t total_edges) {
  std::vector<std::vector<const pgr_edge_t *>> adj(n);
  for (size_t i = 0; i < total_edges; ++i) {
    if (edges[i].cost >= 0.) {
      adj[edges[i].source].push_back(&edges[i]);
    }
    if (edges[i].reverse_cost >= 0.) {
      adj[edges[i].target].push_back(&edges[i]);
    }
  }
  return adj;
}

std::pair<std::vector<Pred>, std::vector<double>>
dijkstra(size_t n, const std::vector<std::vector<const pgr_edge_t *>> &adj,
         int64_t startVertex, double drivingDistance) {
  std::set<pq_el> q; // priority queue
  std::vector<Pred> predecessors(n);
  std::vector<double> distances(n, std::numeric_limits<double>::infinity());
  predecessors[0] = {nullptr, -1, 0.};
  q.insert({0., startVertex});
  while (!q.empty()) {
    double dist;
    int64_t nodeId;
    std::tie(dist, nodeId) = *q.begin();
    if (dist >= drivingDistance) {
      break;
    }
    q.erase(q.begin());
    // if (predecessors[nodeId] != -1) {
    //   // Node was already visited and the shortest path
    //   // was found (predecessor too).
    //   continue;
    // }
    // std::cout << nodeId << " " << dist << std::endl;
    for (auto &&e : adj[nodeId]) {
      int64_t target = e->target == nodeId ? e->source : e->target;
      double cost = e->target == nodeId ? e->reverse_cost : e->cost;
      double aggDist = dist + cost;
      if (distances[target] > aggDist) {
        q.erase({distances[target], target});
        distances[target] =
            aggDist > drivingDistance ? drivingDistance : aggDist;
        predecessors[target] = {e->id, nodeId, distances[target] - dist};
        q.insert({distances[target], e->target});
      }
    }
  }
  return make_pair(predecessors, distances);
}

void do_pgr_many_to_isochrones(pgr_edge_t *data_edges, size_t total_edges,
                               int64_t *start_vertex, size_t s_len,
                               double distance, bool directedFlag,
                               bool equiCostFlag,
                               General_path_element_t **return_tuples,
                               size_t *return_count, char **log_msg,
                               char **notice_msg, char **err_msg) {
  std::ostringstream log;
  std::ostringstream err;
  std::ostringstream notice;

  try {
    pgassert(total_edges != 0);
    pgassert(!(*log_msg));
    pgassert(!(*notice_msg));
    pgassert(!(*err_msg));
    pgassert(!(*return_tuples));
    pgassert(*return_count == 0);
    pgassert((*return_tuples) == NULL);

    std::vector<int64_t> start_vertices(start_vertex, start_vertex + s_len);

    // extract_vertices sorts the vertices by id.
    auto vertices(pgrouting::extract_vertices(data_edges, total_edges));
    // TODO remap vertex ids to 0-N.

    if (directedFlag) {
      pgrouting::DirectedGraph digraph(vertices, gType);

      digraph.insert_edges(data_edges, total_edges, true);

      paths = pgr_drivingDistance(digraph, start_vertices, distance,
                                  equiCostFlag, log);
    } else {
      pgrouting::UndirectedGraph undigraph(vertices, gType);

      undigraph.insert_edges(data_edges, total_edges, true);

      paths = pgr_drivingDistance(undigraph, start_vertices, distance,
                                  equiCostFlag, log);
    }

    size_t count(count_tuples(paths));

    if (count == 0) {
      log << "\nNo return values were found";
      *notice_msg = pgr_msg(log.str().c_str());
      return;
    }
    *return_tuples = pgr_alloc(count, (*return_tuples));
    auto trueCount(collapse_paths(return_tuples, paths));
    *return_count = trueCount;

    *log_msg = log.str().empty() ? *log_msg : pgr_msg(log.str().c_str());
    *notice_msg =
        notice.str().empty() ? *notice_msg : pgr_msg(notice.str().c_str());
  } catch (AssertFailedException &except) {
    (*return_tuples) = pgr_free(*return_tuples);
    (*return_count) = 0;
    err << except.what();
    *err_msg = pgr_msg(err.str().c_str());
    *log_msg = pgr_msg(log.str().c_str());
  } catch (std::exception &except) {
    (*return_tuples) = pgr_free(*return_tuples);
    (*return_count) = 0;
    err << except.what();
    *err_msg = pgr_msg(err.str().c_str());
    *log_msg = pgr_msg(log.str().c_str());
  } catch (...) {
    (*return_tuples) = pgr_free(*return_tuples);
    (*return_count) = 0;
    err << "Caught unknown exception!";
    *err_msg = pgr_msg(err.str().c_str());
    *log_msg = pgr_msg(log.str().c_str());
  }
}
