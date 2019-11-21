#ifndef LAYER_CONTAINER_H_
#define LAYER_CONTAINER_H_

#include <iostream>  // DEBUGING

#include <assert.h>
#include <utility>
#include <vector>
#include <unordered_map>

#include <Kernels/precision.hpp>

#include "Initializer/binning/EncodingConstants.h"
#include "Initializer/binning/Condition.h"
#include "Initializer/binning/ConditionalTable.h"
#include "Initializer/binning/IndexTable.h"
#include "Initializer/binning/specific_types.h"

using conditional_table_t = std::unordered_map< ConditionalKey, IndexTable, ConditionalHash<ConditionalKey> >;

namespace seissol {
    namespace initializers {
        class LayerContainer {

        public:
            void set_conditional_table(conditional_table_t i_table) {
                m_cond_table = i_table;
            }

            conditional_table_t& get_conditional_table() {
                //assert(!m_cond_table.empty() && "conditional table hasn't been initialized");
                return m_cond_table;
            }

            void free_conditional_table() {
                for (auto& index_table: m_cond_table) {
                    for (auto ptr_to_basic_indices_obj: index_table.second.variable_indices) {
                        if (ptr_to_basic_indices_obj != nullptr) {
                            delete ptr_to_basic_indices_obj;
                        }
                    }
                }
            }

            conditional_table_t& get_table_reference_to_init() {
                return m_cond_table;
            }

        private:
            conditional_table_t m_cond_table{};
        };
    }
}

#endif  //LAYER_CONTAINER_H_
