/*=============================================================================
    Copyright (c) 2001-2006 Joel de Guzman

    Use, modification and distribution is subject to the Boost Software
    License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
    http://www.boost.org/LICENSE_1_0.txt)
==============================================================================*/
#include <boost/fusion/sequence/container/list/list.hpp>

#define FUSION_SEQUENCE list
#define FUSION_NO_PRIOR
#define FUSION_TRAVERSAL_TAG forward_traversal_tag
#include "./iterator.hpp"

int
main()
{
    test();
    return boost::report_errors();
}



