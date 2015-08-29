/*PGR

Copyright (c) 2015 Celia Virginia Vergara Castillo
vicky_vergara@hotmail.com

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

*/


CREATE OR REPLACE FUNCTION _pgr_ksp(sql text, start_vid bigint, end_vid bigint, k integer, has_rcost boolean, directed boolean,
  OUT seq integer, OUT path_id integer, OUT path_seq bigint, OUT node bigint, OUT edge bigint, OUT cost float, OUT agg_cost float)
  RETURNS SETOF RECORD AS
    '$libdir/librouting-2.1', 'kshortest_path'
    LANGUAGE c STABLE STRICT;


-- V2 the graph is directed and there are no heap paths 
CREATE OR REPLACE FUNCTION pgr_ksp(sql text, start_vid integer, end_vid integer, k integer, has_rcost boolean)
  RETURNS SETOF pgr_costresult3 AS
  $BODY$  
  DECLARE
  has_reverse boolean;
  BEGIN
      has_reverse =_pgr_parameter_check('ksp', sql::text, false);

      if (has_reverse != has_rcost) then 
         if (has_reverse) then -- raise NOTICE 'has_reverse_cost set to false but reverse_cost column found, Ignoring';
         else raise EXCEPTION 'has_reverse_cost set to true but reverse_cost not found';
         end if;
      end if;

      return query SELECT ((row_number() over()) -1)::integer  as seq,  path_id::integer as id1, node::integer as id2, edge::integer as id3, cost 
            FROM _pgr_ksp(sql::text, start_vid, end_vid, k, has_reverse, true) where path_id < k;
  END
  $BODY$
  LANGUAGE plpgsql VOLATILE
  COST 100
  ROWS 1000;

--V3 DEFAULTS directed:=true heap_paths:=false
CREATE OR REPLACE FUNCTION pgr_ksp(sql text, start_vid bigint, end_vid bigint, k integer,
  directed boolean default true, heap_paths boolean default false,
  OUT seq integer, OUT path_id integer, OUT path_seq bigint, OUT node bigint, OUT edge bigint, OUT cost float, OUT agg_cost float)
  RETURNS SETOF RECORD AS
  $BODY$
  DECLARE
  has_rcost boolean;
  BEGIN
      has_rcost =_pgr_parameter_check('ksp', sql::text, true);
      if heap_paths = false then
         return query SELECT *
                FROM _pgr_ksp(sql::text, start_vid, end_vid, k, has_rcost, directed) a where a.path_id < k;
      else
         return query SELECT *
                FROM _pgr_ksp(sql::text, start_vid, end_vid, k, has_rcost, directed);
      end if;
  END
  $BODY$
  LANGUAGE plpgsql VOLATILE
  COST 100
  ROWS 1000;


