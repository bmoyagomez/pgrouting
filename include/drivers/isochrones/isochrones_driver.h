/*PGR-GNU*****************************************************************
File: isochrones_driver.h

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

#ifndef INCLUDE_DRIVERS_ISOCHRONES_MANY_TO_ISOCHRONES_H_
#define INCLUDE_DRIVERS_ISOCHRONES_MANY_TO_ISOCHRONES_H_

/* for size_t */
#ifdef __cplusplus
#include <cstddef>
#else
#include <stddef.h>
#endif

#include "c_types/general_path_element_t.h"
#include "c_types/pgr_edge_t.h"

typedef struct {
  int64_t start_id;
  int64_t edge;
  double start_perc;
  double end_perc;
  double start_cost;
  double end_cost;
} Isochrones_path_element_t;

#ifdef __cplusplus
extern "C" {
#endif

void do_pgr_many_to_isochrones(pgr_edge_t *edges, size_t total_edges,
                               int64_t *start_vertex, size_t s_len,
                               double *distances, size_t d_len,
                               bool remove_duplicates,
                               Isochrones_path_element_t **return_tuples,
                               size_t *return_count, char **log_msg,
                               char **notice_msg, char **err_msg);

#ifdef __cplusplus
}
#endif

#endif // INCLUDE_DRIVERS_ISOCHRONES_MANY_TO_ISOCHRONES_H_
