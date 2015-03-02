///////////////////////////////////////////////////////////////////////////////
///  GraphElements.h
///  <TODO: insert file description here>
///
///  @remarks <TODO: insert remarks here>
///
///  @author Yan Qi @date 5/28/2010
///
///  $Id: GraphElements.h 65 2010-09-08 06:48:36Z yan.qi.asu $
///////////////////////////////////////////////////////////////////////////////
#ifndef SRC_KSP_SRC_BASEPATH_H_
#define SRC_KSP_SRC_BASEPATH_H_


#include <string>
#include <cassert>
#include <deque>
#include <iostream>
#include <limits>


#include "BaseEdge.h"
#include "BaseVertex.h"



/**************************************************************************
*  BasePath
*  <TODO: insert class description here>
*  edge oriented
*  modified version from the original that is vertex oriented
*
*  @remarks <TODO: insert remarks here>
*
*  @author Yan Qi @date 6/6/2010
*  @modified Vicky Vergara @date Feb/2015
**************************************************************************/
class BasePath {
 protected:
        double m_dWeight;
        int m_start_id;
        std::deque<BaseEdge> m_vtEdgesList;

 public:
        BasePath(unsigned int  start_id, unsigned int sink_id)
            :m_dWeight(0),
            m_start_id(start_id) {
                m_vtEdgesList.clear();
        }

        BasePath()
           :m_dWeight(0),
            m_start_id(-1) {
                m_vtEdgesList.clear();
        }

        explicit BasePath(unsigned int  start_id)
            :m_dWeight(0),
             m_start_id(start_id) {
                m_vtEdgesList.clear();
        }

        BasePath(const std::deque<BaseEdge>& edges_list, double weight)
            :m_dWeight(weight),
             m_start_id(edges_list[0].getStart()) {
                m_vtEdgesList.assign(edges_list.begin(), edges_list.end());
        }

        ~BasePath(void) {}

        double Weight() const { return m_dWeight;}
        void Weight(double val) { m_dWeight = val;}

        POS size() const { return m_vtEdgesList.size();}
        bool isEmpty() const { return size() == 0;}

        int getNID(POS i) const { 
		assert(i < m_vtEdgesList.size());
                return m_vtEdgesList[i].ID();
        }
        int getOriginalID(POS i) const {return m_vtEdgesList[i].originalID();}
        int GetVertex(int i) {
                return m_vtEdgesList.at(i).getStart();
        }

        BaseEdge GetEdge(POS i) { return m_vtEdgesList.at(i); }

        BaseEdge operator[](POS i) const {
		assert(i < m_vtEdgesList.size());
                return (m_vtEdgesList[i]);
        }
/*
        BaseEdge* operator[](POS i) {
                return m_vtEdgesList.at(i);
        }
*/
        bool FromTo(POS from, POS to) const {
             if (size() == 0) return false;
             return from == m_vtEdgesList[0].getStart()
                    && to == m_vtEdgesList[ size()-1 ].getEnd();
        }

       bool isEqual(const  BasePath &largerPath) const {
            if (size() > largerPath.size()) return false;
            for (POS i = 0 ; i < size() ; i++) {
                if (!m_vtEdgesList[i].ID() == largerPath.m_vtEdgesList[i].ID()) return false;
            }
            return true;
       }


        bool EdgesLessComapre(const BasePath &p2) const {
              POS limit = (size() < p2.size()) ? size() : p2.size();
              for (POS i = 0 ; i < limit; i++) {
                   if (m_vtEdgesList[i].ID() < p2.m_vtEdgesList[i].ID())
                       return true;
              }
              return false;
         }

        void push_back(BaseEdge edge) {
            if (size() == 0) m_start_id = edge.getStart();
            m_vtEdgesList.push_back(edge);
            m_dWeight += edge.Weight();
        }

        void push_front(BaseEdge edge) {
            m_start_id = edge.getStart();
            m_vtEdgesList.push_front(edge);
            m_dWeight += edge.Weight();
        }

        void clear() {
            m_vtEdgesList.clear();
            m_dWeight = 0;
        }

        // retuns true & a subPath from the begining of the path with upTo edges
        //        (number of nodes is upTo +1
        // returns false and an empty subpath when
        //  the path is empty
        //  its requierd more elements than it has

        bool subPath(BasePath &sub_path, POS upTo) {
                if (m_vtEdgesList.size() == 0) return false;
                if (upTo >= size()) return false;
                sub_path.m_start_id = m_start_id;
                for (POS i = 0; i < upTo; i++) {
                   sub_path.push_back(m_vtEdgesList[i]);
                }
                if (sub_path.size() == upTo) return true;
                sub_path.clear();
                return false;
        }

        void append(const BasePath &trail) {
            for (POS i=0; i < trail.size(); i++)
                push_back(trail.m_vtEdgesList[i]);
        }


        // display the content
        void PrintOut(std::ostream& out_stream) const {
                out_stream << "Cost: " << m_dWeight << " Length: " << m_vtEdgesList.size() << std::endl;
                for (POS i = 0; i < m_vtEdgesList.size(); i++) {
                        m_vtEdgesList[i].PrintOut(out_stream);
                        out_stream << "->";
                }
                out_stream << std::endl
                  <<  "*********************************************" << std::endl;
        }
};

#endif  // SRC_KSP_SRC_BASEPATH_H_