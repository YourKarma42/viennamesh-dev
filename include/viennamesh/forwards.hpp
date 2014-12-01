#ifndef _VIENNAMESH_FORWARDS_HPP_
#define _VIENNAMESH_FORWARDS_HPP_

#include "viennamesh/backend/api.h"
#include <string>

namespace viennamesh
{

  class context_handle;

  class abstract_data_handle;
  template<typename DataT>
  class data_handle;

  class algorithm_handle;




  namespace result_of
  {
    template<typename T>
    struct data_information;
  }

  static std::string local_binary_format()
  {
    return __VERSION__;
  }
}

#endif