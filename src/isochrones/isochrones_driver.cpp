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

#include <set>
#include <sstream>
#include <unordered_map>
#include <vector>

#include "drivers/isochrones/isochrones_driver.h"

#include "cpp_common/basic_vertex.h"
#include "cpp_common/pgr_alloc.hpp"
#include "cpp_common/pgr_assert.h"

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

struct Pred {
  int64_t edge_id;
  int64_t pred_id;
  double cost; // cost of the edge by going from pred_id to target, cost can be
               // partial if this is the leaf (driving_distance is respected).
               // target can be either pgr_edge_t.{source,target}, depending on
               // which cost is used pgr_edge_t.{cost,reverse_cost}.
};

void dijkstra(int64_t start_vertex, double driving_distance,
              const std::vector<std::vector<const pgr_edge_t *>> &adj,
              std::vector<Pred> *predecessors, std::vector<double> *distances) {
  size_t n = adj.size();
  predecessors->assign(n, {-1, -1, std::numeric_limits<double>::infinity()});
  distances->assign(n, std::numeric_limits<double>::infinity());
  (*predecessors)[0] = {-1, -1, 0.};

  typedef std::tuple<double, int64_t> pq_el; // <agg_cost at node, node id>
  std::set<pq_el> q;                         // priority queue
  q.insert({0., start_vertex});
  while (!q.empty()) {
    double dist;
    int64_t node_id;
    std::tie(dist, node_id) = *q.begin();
    if (dist >= driving_distance) {
      break;
    }
    q.erase(q.begin());
    for (auto &&e : adj[node_id]) {
      int64_t target = e->target == node_id ? e->source : e->target;
      double cost = e->target == node_id ? e->reverse_cost : e->cost;
      double agg_cost = dist + cost;
      if ((*distances)[target] > agg_cost) {
        q.erase({(*distances)[target], target});
        (*distances)[target] =
            agg_cost > driving_distance ? driving_distance : agg_cost;
        double edge_cost =
            agg_cost > driving_distance ? (driving_distance - dist) : cost;
        // Storing the partial travel edge cost.
        (*predecessors)[target] = {e->id, node_id, edge_cost};
        q.emplace((*distances)[target], target);
      }
    }
  }
}

std::unordered_map<int64_t, int64_t> remap_edges(pgr_edge_t *data_edges,
                                                 size_t total_edges) {
  std::unordered_map<int64_t, int64_t> mapping;
  int64_t id = 0;
  for (size_t i = 0; i < total_edges; ++i) {
    pgr_edge_t *e = data_edges + i;
    int64_t source_id, target_id;
    auto it = mapping.find(e->source);
    // better with if-initialization-statement
    if (it != mapping.end()) {
      source_id = it->second;
    } else {
      source_id = id++;
      mapping[e->source] = source_id;
    }
    it = mapping.find(e->target);
    if (it != mapping.end()) {
      target_id = it->second;
    } else {
      target_id = id++;
      mapping[e->target] = target_id;
    }
    e->source = source_id;
    e->target = target_id;
  }
  return mapping;
}

std::vector<General_path_element_t>
do_many_dijkstras(pgr_edge_t *data_edges, size_t total_edges,
                  int64_t *start_vertex, size_t s_len, double distance) {
  // Extracting vertices and mapping the ids from 0 to N-1. Remapping is done
  // so that data structures used can be simpler (arrays instead of maps).
  std::unordered_map<int64_t, int64_t> mapping =
      // modifying data_edges source/target fields.
      remap_edges(data_edges, total_edges);
  size_t nodes_count = mapping.size();
  std::vector<General_path_element_t> results;
  std::vector<int64_t> start_vertices(start_vertex, start_vertex + s_len);
  auto adj =
      construct_adjacency_matrix(mapping.size(), data_edges, total_edges);
  // Storing the result of dijkstra call and reusing the memory for each vertex.
  std::vector<Pred> predecessors(nodes_count);
  std::vector<double> distances(nodes_count);
  for (int64_t start_v : start_vertices) {
    auto it = mapping.find(start_v);
    // If start_v did not appear in edges then it has no particular mapping but
    // pgr_drivingDistance result includes one row for this node.
    if (it == mapping.end()) {
      General_path_element_t r;
      r.seq = 0;
      r.start_id = start_v;
      r.end_id = start_v;
      r.node = start_v;
      // -2 tags the unmapped starting vertex and won't use the reverse_mapping
      // because mapping does not exist. -2 is changed to -1 later.
      r.edge = -2;
      r.cost = 0.0;
      r.agg_cost = 0.0;
      results.push_back(r);
      continue;
    }
    // Calling the dijkstra algorithm and storing the results in predecessors
    // and distances.
    dijkstra(it->second,
             /* driving_distance */ distance, adj, &predecessors, &distances);
    // Appending the row results.
    for (int64_t node_id = 0; node_id < predecessors.size(); ++node_id) {
      const Pred &pred = predecessors[node_id];
      if (std::isinf(pred.cost) || pred.cost > distance) {
        // Node not reached or reached too late is skipped.
        continue;
      }
      General_path_element_t r;
      // r.seq is filled later. r.node is remapped later.
      r.start_id = start_v;
      r.end_id = start_v;
      r.edge = pred.edge_id;
      r.node = node_id;
      r.cost = pred.cost;
      r.agg_cost = distances[node_id];
      results.push_back(r);
    }
  }
  // The 0 to N-1 ids need to be mapped back to their previous values.
  std::vector<int64_t> reverse_mapping(mapping.size());
  for (auto &m : mapping) {
    reverse_mapping[m.second] = m.first;
  }
  int seq_cnt = 0;
  int64_t prev_start_id = -1;
  for (size_t i = 0; i < results.size(); ++i) {
    auto &r = results[i];
    if (r.edge == -2) {
      // This is a starting vertex that has no edges in graph.
      r.seq = 0;
      // Changing back the -2 edge indicator to -1.
      r.edge = -1;
      // Force seq_cnt reset in next iteration and avoid accidental id clash
      // with the mapping.
      prev_start_id = -1;
      continue;
    }
    if (r.start_id != prev_start_id) {
      prev_start_id = r.start_id;
      seq_cnt = 0;
    }
    r.seq = seq_cnt++;
    r.node = reverse_mapping[r.node];
  }
  return results;
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

    auto results = do_many_dijkstras(data_edges, total_edges, start_vertex,
                                     s_len, distance);

    size_t count(results.size());
    if (count == 0) {
      log << "\nNo return values were found";
      *notice_msg = pgr_msg(log.str().c_str());
      return;
    }
    *return_tuples = pgr_alloc(count, (*return_tuples));
    *return_count = count;

    // Moving results to allocated return_tuples array.
    for (size_t i = 0; i < count; ++i) {
      (*return_tuples)[i] = results[i];
    }

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
