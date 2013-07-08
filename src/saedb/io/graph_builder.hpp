#pragma once
#include <cstdint>
#include <cstddef>
#include <map>
#include <vector>
#include <utility>
#include <sstream>
#include <cstring>
#include <functional>

#include "../serialization/serialization_includes.hpp"

#include "type_info.hpp"

using sae::serialization::ISerializeStream;
using sae::serialization::OSerializeStream;

namespace sae {
namespace io {

    typedef std::uint64_t vid_t;
    typedef std::uint64_t eid_t;
    using std::size_t;
    using std::map;
    using std::vector;
    using std::make_pair;

    enum GraphStorageType {
        MemoryMappedGraph = 0,
        GraphStorageTypeCount
    };

    struct edge_with_data {
        vid_t source;
        vid_t target;
        eid_t global_id;
        eid_t local_id;
        std::string data;
        uint32_t type_name;
        edge_with_data(vid_t source, vid_t target, eid_t global_id, eid_t local_id, std::string data, uint32_t type_name) :
            source(source), target(target), global_id(global_id), local_id(local_id), data(data), type_name(type_name) {}
    };

    struct vertex_with_data {
        vid_t global_id;
        vid_t local_id;
        uint32_t type_name;
        std::string data;
        vertex_with_data(vid_t global_id, vid_t local_id, uint32_t type_name, std::string data) :
            global_id(global_id), local_id(local_id), type_name(type_name), data(data) {
        }
    };

    struct data_type {
        uint32_t data_size;
        DataTypeAccessor* data_type_accessor;
        uint32_t count;
        data_type(int data_size, DataTypeAccessor* data_type_accessor) :
            data_size(data_size), data_type_accessor(data_type_accessor) {
            count = 0;
        }
    };

    struct GraphWriter {
        virtual void AppendVertex(vid_t, vid_t, uint32_t, std::string) = 0;
        virtual void AppendEdge(vid_t, vid_t, eid_t, eid_t, std::string, uint32_t) = 0;
        virtual void AppendVertexDataType(uint32_t data_size, DataTypeAccessor* accessor, uint32_t count) = 0;
        virtual void AppendEdgeDataType(uint32_t data_size, DataTypeAccessor* accessor, uint32_t count) = 0;
        virtual void Close() = 0;
        virtual ~GraphWriter(){}
    };

    GraphWriter* CreateMemoryMappedGraphWriter(const char*, vid_t, eid_t, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t*, uint32_t*, uint32_t*, uint32_t*);

    template<typename vkey_t>
    struct GraphBuilder
    {
        template <typename T>
        void AddVertex(vkey_t key, const char* data_type_name, T data) {
            int data_type_rank = -1;
            for (int i=0, size = vertex_data_types.size(); i<size; i++) {
                if (strcmp(vertex_data_types[i].data_type_accessor->getTypeName(), data_type_name) == 0) {
                    data_type_rank = i;
                    break;
                }
            }
            if (data_type_rank == -1) {
                std::cout << "ERROR: type " << data_type_name << " is undefined." << std::endl;
                return;
            }
            auto id = map(key);
            auto local_id = vertex_data_types[data_type_rank].count ++;

            std::string code = sae::serialization::convert_to_string(data);

            vertices[id] = vertex_with_data(id, local_id, data_type_rank, code);
        }

        template <typename T>
        void AddEdge(vkey_t source, vkey_t target, const char* data_type_name, T data) {
            int data_type_rank = -1;
            for (int i=0, size = edge_data_types.size(); i<size; i++) {
                if (strcmp(edge_data_types[i].data_type_accessor->getTypeName(), data_type_name) == 0) {
                    data_type_rank = i;
                    break;
                }
            }
            if (data_type_rank == -1) {
                std::cout << "ERROR: type " << data_type_name << " is undefined." << std::endl;
                return;
            }
            auto sid = map(source);
            auto tid = map(target);
            
            eid_t id = edges.size();
            auto local_id = edge_data_types[data_type_rank].count ++;
            
            std::string code = sae::serialization::convert_to_string(data);

            edges.push_back(edge_with_data(sid, tid, id, local_id, code, data_type_rank));
        }

        virtual vid_t VertexCount() {
            return vertices.size();
        }

        virtual eid_t EdgeCount() {
            return edges.size();
        }

        uint32_t VertexTypeTotalSize() {
            uint32_t result = 0;
            for (auto p : vertex_data_types) {
                result += p.data_type_accessor->Size();
            }
            return result;
        }

        uint32_t EdgeTypeTotalSize() {
            uint32_t result = 0;
            for (auto p: edge_data_types) {
                result += p.data_type_accessor->Size();
            }
            return result;
        }

        virtual bool Save(const char * prefix, GraphStorageType type = MemoryMappedGraph) {
            if (type == MemoryMappedGraph) {
                uint32_t * vertex_type_count = new uint32_t[vertex_data_types.size()];
                for (int i=0, size = vertex_data_types.size(); i<size; i++) {
                    vertex_type_count[i] = vertex_data_types[i].count;
                }

                uint32_t * edge_type_count = new uint32_t[edge_data_types.size()];
                for (int i=0, size = edge_data_types.size(); i<size; i++) {
                    edge_type_count[i] = edge_data_types[i].count;
                }

                uint32_t * vertex_data_size = new uint32_t[vertex_data_types.size()];
                for (int i=0, size = vertex_data_types.size(); i<size; i++) {
                    vertex_data_size[i] = vertex_data_types[i].data_size;
                }

                uint32_t * edge_data_size = new uint32_t[edge_data_types.size()];
                for (int i=0, size = edge_data_types.size(); i<size; i++) {
                    edge_data_size[i] = edge_data_types[i].data_size;
                }

                GraphWriter* writer = CreateMemoryMappedGraphWriter(prefix, vertices.size(), edges.size(), vertex_data_types.size(), edge_data_types.size(), VertexTypeTotalSize(), EdgeTypeTotalSize(), vertex_type_count, edge_type_count, vertex_data_size, edge_data_size);

                if (!writer) return false;
                for (auto data_type: vertex_data_types) {
                    writer->AppendVertexDataType(data_type.data_size, data_type.data_type_accessor, data_type.count);
                }
                for (auto data_type: edge_data_types) {
                    writer->AppendEdgeDataType(data_type.data_size, data_type.data_type_accessor, data_type.count);
                }
                for (auto& v : vertices) {
                    writer->AppendVertex(v.global_id, v.local_id, v.type_name, v.data);
                }
                for (auto& e: edges) {
                    writer->AppendEdge(e.source, e.target, e.global_id, e.local_id, e.data, e.type_name);
                }
                writer->Close();

                delete writer;
                delete [] vertex_type_count;
                delete [] edge_type_count;
                delete [] vertex_data_size;
                delete [] edge_data_size;
            }
        }

        virtual DataTypeAccessor* CreateType(const char* tn) {
            return DataTypeAccessorFactory(tn);
        }

        virtual void SaveVertexDataType(DataTypeAccessor* dt, int data_type_size) {
            vertex_data_types.push_back(data_type(data_type_size, dt));
        }

        virtual void SaveEdgeDataType(DataTypeAccessor* dt, int data_type_size) {
            edge_data_types.push_back(data_type(data_type_size, dt));
        }

        /*
         * Associate each data type with a serialize/deserialize function. 
         *
         * After call Associate() with each data type, Graph Builder saves all the functions
         * for serialization. And in later loading/saving of Graph, builder will get related
         * serialize/deserialize function at first, then use them as a wrapper for underlying 
         * serialization.
         * 
         * Parameter: 
         * data_type_rank: index of specific data type
         * 
         */
        template <typename T>
        void Associate(int data_type_rank) {
            auto iit = serialize_func_map.find(data_type_rank);
            if (iit == serialize_func_map.end()) {
                serialize_func_map[data_type_rank] = [&](OSerializeStream& out, void* d) {
                    out << *((T*)d);
                };
            }

            auto oit = deserialize_func_map.find(data_type_rank);
            if (oit == deserialize_func_map.end()) {
                deserialize_func_map[data_type_rank] = [&](ISerializeStream& in, void* d) {
                    in >> *((T*)d);
                };
            }
        }

        // no destructor?
    private:
        std::vector<vertex_with_data> vertices;
        std::vector<edge_with_data> edges;

        std::vector<data_type> vertex_data_types;
        std::vector<data_type> edge_data_types;

        std::map<vkey_t, vid_t> vid_map;

        /* map data type rank to corresponding serializa/deserialize function */
        std::map<int, std::function<void(OSerializeStream&, void*)> > serialize_func_map;   /* serialize function map */
        std::map<int, std::function<void(ISerializeStream&, void*)> > deserialize_func_map; /* deserialize function map */

        vid_t map(vkey_t key) {
            vid_t result;
            auto it = vid_map.find(key);
            if (it == vid_map.end()) {
                result = vertices.size();
                vertices.push_back(vertex_with_data(0, 0, 0, std::string("")));
                vid_map.insert(make_pair(key, result));
            } else {
                result = it->second;
            }

            return result;
        }
    };

}
}
